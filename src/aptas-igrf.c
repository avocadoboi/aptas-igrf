#include "aptas-igrf.h"

#include "math_constants.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define MAX_HARMONIC_DEGREE 13

#define G_COEFFICIENT_COUNT MAX_HARMONIC_DEGREE*(MAX_HARMONIC_DEGREE + 3)/2
#define H_COEFFICIENT_COUNT MAX_HARMONIC_DEGREE*(MAX_HARMONIC_DEGREE + 1)/2

static float g_coefficients[G_COEFFICIENT_COUNT];
static float g_coefficients_rate[G_COEFFICIENT_COUNT];

static float h_coefficients[H_COEFFICIENT_COUNT];
static float h_coefficients_rate[H_COEFFICIENT_COUNT];

IGRFError load_IGRF_coefficients(char const* const file_name) {
  FILE* const file = fopen(file_name, "rt");
  if (file == NULL) {
    return IGRFError_FileNotFound;
  }
  int g_index = 0;
  int h_index = 0;
  for (int n = 1; n <= MAX_HARMONIC_DEGREE; ++n) {
    for (int m = 0; m <= n; ++m) {
      if (fscanf(file, "%f%f", g_coefficients + g_index, g_coefficients_rate + g_index) != 2) {
        fclose(file);
        return IGRFError_UnexpectedFileFormat;
      }
      // printf("g_%i^%i: %f %f\n", n, m, g_coefficients[g_index], g_coefficients_rate[g_index]);
      ++g_index;
      if (m != 0) {
        if (fscanf(file, "%f%f", h_coefficients + h_index, h_coefficients_rate + h_index) != 2) {
          fclose(file);
          return IGRFError_UnexpectedFileFormat;
        }
        // printf("h_%i^%i: %f %f\n", n, m, h_coefficients[h_index], h_coefficients_rate[h_index]);
        ++h_index;
      }
    }
  }
  return IGRFError_Success;
}

typedef struct {
  // Starting values for the recursive formulas for the associated legendre functions.
  float P_nm_initial[1+2+3];
  float P_nm_derivative_initial[1+2+3];
  /*
  {
    { P_{n-1}^0, P_{n-2}^0 },
    { P_{n-1}^1, P_{n-2}^1 },
  }
  */
  float P_up[2][2];
  // Derivative of the legendre functions above, with respect to mu = cos(theta).
  float P_up_derivative[2][2];
  /*
  { P_n^{m-1}, P_n^{m-2} }
  */
  float P_left[2];
  float P_left_derivative[2];
} associated_legendre_state_t;

static void init_associated_legendre_state(associated_legendre_state_t* const legendre_state, float const cos_theta, float const sin_theta) {
  // Initialize the associated legendre state struct, containing the initial values that the recurrence relations
  // start from as well as arrays containing P_n^{m-1} and P_n^{m-2} (P_left) that are updated for each recurrence step
  // and similarly for the downward propagation in the two leftmost columns (m = 0 and m = 1).
  memcpy(legendre_state->P_nm_initial, (float[1+2+3]){
    // 0
    1,
    // 1                              2
    cos_theta,                        sin_theta,
    // 3                              4                           5
    0.5f*(3*cos_theta*cos_theta - 1), sqrt_3*cos_theta*sin_theta, 0.5f*sqrt_3*sin_theta*sin_theta
  }, sizeof(legendre_state->P_nm_initial));

  memcpy(legendre_state->P_nm_derivative_initial, (float[1+2+3]){
    // 0
    0,
    // 1                    2
    -sin_theta,             cos_theta,
    // 3                    4                                                   5
    -3*cos_theta*sin_theta, sqrt_3*(cos_theta*cos_theta - sin_theta*sin_theta), sqrt_3*sin_theta*cos_theta
  }, sizeof(legendre_state->P_nm_derivative_initial));

  // P_left and P_left_derivative will be filled in later once we have P_4_0 and P_4_1.
  // We leave it uninitialized until then.

  // legendre.P_up and legendre.P_up_derivative depend on P_nm_initial and P_nm_derivative_initial so they need to be filled in after.
  // n = 3 is the first row that uses the recurrence relation to calculate the associated legendre functions.
  memcpy(legendre_state->P_up, (float[2][2]){
    // P_2^0                    P_1^0
    { legendre_state->P_nm_initial[3], legendre_state->P_nm_initial[1] },
    // P_2^1                    P_1^1
    { legendre_state->P_nm_initial[4], legendre_state->P_nm_initial[2] }
  }, sizeof(legendre_state->P_up));

  memcpy(legendre_state->P_up_derivative, (float[2][2]){
    { legendre_state->P_nm_derivative_initial[3], legendre_state->P_nm_derivative_initial[1] },
    { legendre_state->P_nm_derivative_initial[4], legendre_state->P_nm_derivative_initial[2] }
  }, sizeof(legendre_state->P_up_derivative));
}
/*
 * m must be either 0 or 1 since we only propagate down in the first two columns.
 */
static void associated_legendre_propagate_down(
  int const n, int const m, 
  associated_legendre_state_t* const legendre, 
  float const cos_theta, float const sin_theta
) {
  assert(m == 0 || m == 1);
  float const A = (float)(2*n - 1)/sqrtf((float)(n^2 - m^2));
  float const B = sqrtf((float)((n - 1)^2 - m)/(float)(n^2 - m^2));
  
  float const P_nm = A*cos_theta*legendre->P_up[m][0] - B*legendre->P_up[m][1];
  float const P_nm_derivative = A*(cos_theta*legendre->P_up_derivative[m][0] - sin_theta*legendre->P_up[m][0]) - B*legendre->P_up_derivative[m][1];

  legendre->P_up[m][1] = legendre->P_up[m][0];
  legendre->P_up[m][0] = P_nm;

  legendre->P_up_derivative[m][1] = legendre->P_up_derivative[m][0];
  legendre->P_up_derivative[m][0] = P_nm_derivative;
}
static void associated_legendre_propagate_right(
  int const n, int const m, 
  associated_legendre_state_t* const legendre, 
  float const cos_theta, float const sin_theta
) {
  float const C = 2.f*(float)(m - 1)/sqrtf((float)((n - m + 1)*(n + m)));
  float const D = sqrtf((float)((n - m + 2)*(n + m - 1))/(float)((n + m)*(n - m + 1)));

  float const P_nm = C*cos_theta/sin_theta*legendre->P_left[0] - D*legendre->P_left[1];
  float const P_nm_derivative = C/sin_theta*(cos_theta*legendre->P_left_derivative[0] - legendre->P_left[0]) - D*legendre->P_left_derivative[1];

  legendre->P_left[1] = legendre->P_left[0];
  legendre->P_left[0] = P_nm;

  legendre->P_left_derivative[1] = legendre->P_left_derivative[0];
  legendre->P_left_derivative[0] = P_nm_derivative;
}

typedef struct {
  float P_nm;
  float P_nm_derivative;
} associated_legendre_pair_t;

static associated_legendre_pair_t next_associated_legendre(
  int const n, int const m, int const g_index, 
  associated_legendre_state_t* const legendre_state, 
  float const cos_theta, float const sin_theta
) {
  // If this is true then we use recurrence relations to calculate P_n^m. 
  // Otherwise we use legendre_state->P_nm_initial. (similarly for the derivatives)
  if (n >= 3) {
    if (m == 0 || m == 1) {
      associated_legendre_propagate_down(n, m, legendre_state, cos_theta, sin_theta);
      legendre_state->P_left[1 - m] = legendre_state->P_up[m][0];
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
    assert(g_index < sizeof(legendre_state->P_nm_initial)/sizeof(legendre_state->P_nm_initial[0]));
    // The array of 'initial' (as in those that the recurrence relations start from) 
    // associated legendre functions are indexed the same as the g_n_m coefficients.
    return (associated_legendre_pair_t){
      .P_nm = legendre_state->P_nm_initial[g_index],
      .P_nm_derivative = legendre_state->P_nm_derivative_initial[g_index]
    };
  }
}

magnetic_field_vector_t calculate_model_geomagnetic_field(float const latitude, float const longitude, float const altitude, float const decimal_year) {
  // theta is the co-latitude.
  float const cos_theta = cosf((90.f - latitude)*deg_to_rad);
  float const sin_theta = sqrtf(1 - cos_theta*cos_theta);
  float const r = igrf_earth_radius + altitude;

  // Initialize everything except for the P_nm_left and P_nm_left_derivative, those need to be filled in once we have propagated down for the first time.
  associated_legendre_state_t legendre_state;
  init_associated_legendre_state(&legendre_state, cos_theta, sin_theta);

  magnetic_field_vector_t field = {0};
  int g_index = 0;
  int h_index = 0;
  for (int n = 1; n <= MAX_HARMONIC_DEGREE; ++n) {
    float const north_factor = 1/r*igrf_earth_radius*powf(igrf_earth_radius/r, (float)(n + 1));
    for (int m = 0; m <= n; ++m) {
      float const g_nm = g_coefficients[g_index] + (decimal_year - 2025.f)*g_coefficients_rate[g_index];
      ++g_index;

      float h_nm = 0.f;
      if (m != 0) {
        h_nm = h_coefficients[h_index] + (decimal_year - 2025.f)*h_coefficients_rate[h_index];
        ++h_index;
      }

      associated_legendre_pair_t const new_legendre = next_associated_legendre(n, m, g_index, &legendre_state, cos_theta, sin_theta);
      
    }
  }
  
  return field;
}

