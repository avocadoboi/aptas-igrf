#include "aptas-igrf.h"

#include "math_constants.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef CPP_ASSOC_LEGENDRE_TESTING
#include <array>
#include <cmath>
#endif

#define MAX_HARMONIC_DEGREE 13

#define G_COEFFICIENT_COUNT MAX_HARMONIC_DEGREE*(MAX_HARMONIC_DEGREE + 3)/2
#define H_COEFFICIENT_COUNT MAX_HARMONIC_DEGREE*(MAX_HARMONIC_DEGREE + 1)/2

static double g_coefficients[G_COEFFICIENT_COUNT];
static double g_coefficients_rate[G_COEFFICIENT_COUNT];

static double h_coefficients[H_COEFFICIENT_COUNT];
static double h_coefficients_rate[H_COEFFICIENT_COUNT];

IGRFError load_IGRF_coefficients(char const* const file_name) {
  FILE* const file = fopen(file_name, "rt");
  if (file == NULL) {
    return IGRFError_FileNotFound;
  }
  unsigned g_index = 0;
  unsigned h_index = 0;
  for (unsigned n = 1; n <= MAX_HARMONIC_DEGREE; ++n) {
    for (unsigned m = 0; m <= n; ++m) {
      if (fscanf(file, "%lf%lf", g_coefficients + g_index, g_coefficients_rate + g_index) != 2) {
        fclose(file);
        return IGRFError_UnexpectedFileFormat;
      }
      // printf("g_%i^%i: %lf %lf\n", n, m, g_coefficients[g_index], g_coefficients_rate[g_index]);
      ++g_index;
      if (m != 0) {
        if (fscanf(file, "%lf%lf", h_coefficients + h_index, h_coefficients_rate + h_index) != 2) {
          fclose(file);
          return IGRFError_UnexpectedFileFormat;
        }
        // printf("h_%i^%i: %lf %lf\n", n, m, h_coefficients[h_index], h_coefficients_rate[h_index]);
        ++h_index;
      }
    }
  }
  return IGRFError_Success;
}

typedef struct {
  // Starting values for the recursive formulas for the associated legendre functions.
  double P_nm_initial[1+2+3];
  double P_nm_derivative_initial[1+2+3];
  /*
  {
    { P_{n-1}^0, P_{n-2}^0 },
    { P_{n-1}^1, P_{n-2}^1 },
  }
  */
  double P_up[2][2];
  // Derivative of the legendre functions above, with respect to mu = cos(theta).
  double P_up_derivative[2][2];
  /*
  { P_n^{m-1}, P_n^{m-2} }
  */
  double P_left[2];
  double P_left_derivative[2];
} associated_legendre_state_t;

static void init_associated_legendre_state(associated_legendre_state_t* const legendre_state, double const cos_theta, double const sin_theta) {
  // Initialize the associated legendre state struct, containing the initial values that the recurrence relations
  // start from as well as arrays containing P_n^{m-1} and P_n^{m-2} (P_left) that are updated for each recurrence step
  // and similarly for the downward propagation in the two leftmost columns (m = 0 and m = 1).
  memcpy(legendre_state->P_nm_initial, (double[1+2+3]){
    // 0
    1,
    // 1                              2
    cos_theta,                        sin_theta,
    // 3                              4                                  5
    0.5f*(3*cos_theta*cos_theta - 1), (double)sqrt_3*cos_theta*sin_theta, 0.5f*(double)sqrt_3*sin_theta*sin_theta
  }, sizeof(legendre_state->P_nm_initial));

  memcpy(legendre_state->P_nm_derivative_initial, (double[1+2+3]){
    // 0
    0,
    // 1                    2
    -sin_theta,             cos_theta,
    // 3                    4                                                          5
    -3*cos_theta*sin_theta, (double)sqrt_3*(cos_theta*cos_theta - sin_theta*sin_theta), (double)sqrt_3*sin_theta*cos_theta
  }, sizeof(legendre_state->P_nm_derivative_initial));

  // P_left and P_left_derivative will be filled in later once we have P_4_0 and P_4_1.
  // We leave it uninitialized until then.

  // legendre.P_up and legendre.P_up_derivative depend on P_nm_initial and P_nm_derivative_initial so they need to be filled in after.
  // n = 3 is the first row that uses the recurrence relation to calculate the associated legendre functions.
  memcpy(legendre_state->P_up, (double[2][2]){
    // P_2^0                           P_1^0
    { legendre_state->P_nm_initial[3], legendre_state->P_nm_initial[1] },
    // P_2^1                           P_1^1
    { legendre_state->P_nm_initial[4], legendre_state->P_nm_initial[2] }
  }, sizeof(legendre_state->P_up));

  memcpy(legendre_state->P_up_derivative, (double[2][2]){
    { legendre_state->P_nm_derivative_initial[3], legendre_state->P_nm_derivative_initial[1] },
    { legendre_state->P_nm_derivative_initial[4], legendre_state->P_nm_derivative_initial[2] }
  }, sizeof(legendre_state->P_up_derivative));
}
/*
 * m must be either 0 or 1 since we only propagate down in the first two columns.
 */
static void associated_legendre_propagate_down(
  unsigned const n, unsigned const m, 
  associated_legendre_state_t* const legendre, 
  double const cos_theta, double const sin_theta
) {
  assert(m == 0 || m == 1);
  double const A = (double)(2*n - 1)/sqrtf((double)(n*n - m*m));
  double const B = sqrtf((double)((n - 1)*(n - 1) - m)/(double)(n*n - m*m));
  
  double const P_nm = A*cos_theta*legendre->P_up[m][0] - B*legendre->P_up[m][1];
  double const P_nm_derivative = A*(cos_theta*legendre->P_up_derivative[m][0] - sin_theta*legendre->P_up[m][0]) - B*legendre->P_up_derivative[m][1];

  legendre->P_up[m][1] = legendre->P_up[m][0];
  legendre->P_up[m][0] = P_nm;

  legendre->P_up_derivative[m][1] = legendre->P_up_derivative[m][0];
  legendre->P_up_derivative[m][0] = P_nm_derivative;
}
static void associated_legendre_propagate_right(
  unsigned const n, unsigned const m, 
  associated_legendre_state_t* const legendre, 
  double const cos_theta, double const sin_theta
) {
  double const C = 2.*(double)(m - 1)/sqrtf((double)((n - m + 1)*(n + m)));
  double const D = sqrtf((m == 2 ? 2 : 1)*(double)((n - m + 2)*(n + m - 1))/(double)((n + m)*(n - m + 1)));

  double const P_nm = C*cos_theta/sin_theta*legendre->P_left[0] - D*legendre->P_left[1];
  double const P_nm_derivative = C/sin_theta*(cos_theta*legendre->P_left_derivative[0] - legendre->P_left[0]/sin_theta) - D*legendre->P_left_derivative[1];

  legendre->P_left[1] = legendre->P_left[0];
  legendre->P_left[0] = P_nm;

  legendre->P_left_derivative[1] = legendre->P_left_derivative[0];
  legendre->P_left_derivative[0] = P_nm_derivative;
}

typedef struct {
  double P_nm;
  double P_nm_derivative;
} associated_legendre_pair_t;

static associated_legendre_pair_t next_associated_legendre(
  unsigned const n, unsigned const m, unsigned const P_index, 
  associated_legendre_state_t* const legendre_state, 
  double const cos_theta, double const sin_theta
) {
  // If this is true then we use recurrence relations to calculate P_n^m. 
  // Otherwise we use legendre_state->P_nm_initial. (similarly for the derivatives)
  if (n >= 3) {
    if (m == 0 || m == 1) {
      associated_legendre_propagate_down(n, m, legendre_state, cos_theta, sin_theta);
      legendre_state->P_left[1 - m] = legendre_state->P_up[m][0];
      legendre_state->P_left_derivative[1 - m] = legendre_state->P_up_derivative[m][0];
      return (associated_legendre_pair_t){
        .P_nm = legendre_state->P_up[m][0],
        .P_nm_derivative = legendre_state->P_up_derivative[m][0]
      };
    }
    else {
      associated_legendre_propagate_right(n, m, legendre_state, cos_theta, sin_theta);
      return (associated_legendre_pair_t){
        .P_nm = legendre_state->P_left[0],
        .P_nm_derivative = legendre_state->P_left_derivative[0]
      };
    }
  }
  else {
    // If we want to use less memory we could have a switch statement over n with switch statements over m 
    // inside each case and calculating the values here instead of having the arrays P_nm_initial and P_nm_derivative_initial.
    // Then we would fill in P_up and P_up_derivative here also. 

    assert(P_index < sizeof(legendre_state->P_nm_initial)/sizeof(legendre_state->P_nm_initial[0]));
    // The array of 'initial' (as in those that the recurrence relations start from) 
    // associated legendre functions are indexed the same as the g_n_m coefficients.
    return (associated_legendre_pair_t){
      .P_nm = legendre_state->P_nm_initial[P_index],
      .P_nm_derivative = legendre_state->P_nm_derivative_initial[P_index]
    };
  }
}

magnetic_field_vector_t calculate_model_geomagnetic_field(double const latitude, double longitude, double const altitude, double const decimal_year) {
  longitude *= deg_to_rad;
  // theta is the co-latitude.
  double const cos_theta = cosf((90. - latitude)*deg_to_rad);
  double const sin_theta = sqrtf(1 - cos_theta*cos_theta);
  double const r = igrf_earth_radius + altitude;

  // Initialize everything except for the P_nm_left and P_nm_left_derivative, those need to be filled in once we have propagated down for the first time.
  associated_legendre_state_t legendre_state;
  init_associated_legendre_state(&legendre_state, cos_theta, sin_theta);

  magnetic_field_vector_t field = {0};
  unsigned g_index = 0;
  unsigned h_index = 0;
  for (unsigned n = 1; n <= MAX_HARMONIC_DEGREE; ++n) {
    for (unsigned m = 0; m <= n; ++m) {
      double const g_nm = g_coefficients[g_index] + (decimal_year - 2025.)*g_coefficients_rate[g_index];
      ++g_index;

      double h_nm = 0.;
      if (m != 0) {
        h_nm = h_coefficients[h_index] + (decimal_year - 2025.)*h_coefficients_rate[h_index];
        ++h_index;
      }

      associated_legendre_pair_t const new_legendre = next_associated_legendre(n, m, g_index, &legendre_state, cos_theta, sin_theta);
      printf("n = %u m = %u g_index = %u Calculated: %.5g\n", n, m, g_index, new_legendre.P_nm);
      // Convenient testing using the C++ std::assoc_legendre function. A few other modifications in the code are required to build as C++.
#ifdef CPP_ASSOC_LEGENDRE_TESTING
      double const correct = std::assoc_legendre(n, m, cos_theta)*std::sqrt((m == 0 ? 1 : 2)*static_cast<double>(std::tgamma(n - m + 1))/static_cast<double>(std::tgamma(n + m + 1)));
      printf("n = %u m = %u g_index = %u Calculated: %.5g\t Correct: %.5g\n", n, m, g_index, new_legendre.P_nm, correct);
#endif

      field.up += (n + 1)*powf(igrf_earth_radius/r, n + 2)*(g_nm*cosf(m*longitude) + h_nm*sinf(m*longitude))*new_legendre.P_nm;
      field.north += powf(igrf_earth_radius/r, n + 2)*(g_nm*cosf(m*longitude) + h_nm*sinf(m*longitude))*new_legendre.P_nm_derivative;
      field.east -= powf(igrf_earth_radius/r, n + 2)/sin_theta*(h_nm*m*cosf(m*longitude) - g_nm*m*sinf(m*longitude))*new_legendre.P_nm;
    }
  }
  
  return field;
}

