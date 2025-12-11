#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wuninitialized"

#include "aptas-igrf.h"

#include "math-constants.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef CPP_ASSOC_LEGENDRE_TESTING
#include <array>
#include <cmath>
#endif

// The maximum harmonic degree in the IGRF model.
#define N 13

#define G_COEFFICIENT_COUNT N*(N + 3)/2
#define H_COEFFICIENT_COUNT N*(N + 1)/2

// Stored column-by column starting from m = 0 with n increasing for each m.
static float g_coefficients[G_COEFFICIENT_COUNT];
static float g_coefficients_rate[G_COEFFICIENT_COUNT];

static float h_coefficients[H_COEFFICIENT_COUNT];
static float h_coefficients_rate[H_COEFFICIENT_COUNT];

IGRFError load_IGRF_coefficients(char const* const file_name) {
  FILE* const file = fopen(file_name, "rt");
  if (file == NULL) {
    return IGRFError_FileNotFound;
  }
  for (int n = 1; n <= N; ++n) {
    for (int m = 0; m <= n; ++m) {
      // See equation (??) in the accompanying PDF.
      int const g_index = n - 1 + m*N - m*(m - 1)/2; // Integer division here is ok, m*(m-1) is always even.

      if (fscanf(file, "%f%f", g_coefficients + g_index, g_coefficients_rate + g_index) != 2) {
        fclose(file);
        return IGRFError_UnexpectedFileFormat;
      }
      // printf("g_%i^%i: %f %f\n", n, m, g_coefficients[g_index], g_coefficients_rate[g_index]);
      if (m != 0) {
        int const h_index = n - 1 - N + m*N - m*(m - 1)/2;

        if (fscanf(file, "%f%f", h_coefficients + h_index, h_coefficients_rate + h_index) != 2) {
          fclose(file);
          return IGRFError_UnexpectedFileFormat;
        }
        // printf("h_%i^%i: %f %f\n", n, m, h_coefficients[h_index], h_coefficients_rate[h_index]);
      }
    }
  }
  return IGRFError_Success;
}

typedef struct {
  float P_up[2];
  float P_up_derivative[2];

  float P_nn;
  float P_nn_derivative;
} associated_legendre_state_t;

typedef struct {
  float P_nm;
  float P_nm_derivative;
} associated_legendre_pair_t;

static associated_legendre_pair_t next_associated_legendre(
  unsigned const n, unsigned const m, associated_legendre_state_t* const legendre, 
  float const cos_theta, float const sin_theta
) {
  if (n == 1) {
    if (m == 0) {
      legendre->P_up[0] = 1;
      legendre->P_up_derivative[0] = 0;
    }
    else {
      legendre->P_nn = 1;
      legendre->P_nn_derivative = 0;
    }
  }

  // Diagonal recurrence
  if (n == m) {
    float const A = sqrtf((n == 1 ? 2 : 1)*(1 - 1./(float)(2*n)));
    // P_nn is for saving the value until we get to P_{n+1}^{n+1} after being finished looping through the current column,
    // P_up is for saving the value until we get to the very next index, P_{n+1}^n.
    legendre->P_nn_derivative = A*(sin_theta*legendre->P_nn_derivative + cos_theta*legendre->P_nn);
    legendre->P_nn = A*sin_theta*legendre->P_nn;
    legendre->P_up_derivative[0] = legendre->P_nn_derivative;
    legendre->P_up[0] = legendre->P_nn;
    return (associated_legendre_pair_t){
      .P_nm = legendre->P_nn,
      .P_nm_derivative = legendre->P_nn_derivative
    };
  }

  // Otherwise, downwards recurrence
  float const A = (float)(2*n - 1)/sqrtf((float)(n*n - m*m));

  // If there is only one (nonzero) P_n^m (i.e. P_{n-1}^m) above the current one 
  if (m == n - 1) {
    legendre->P_up[1] = legendre->P_up[0];
    legendre->P_up[0] = A*cos_theta*legendre->P_up[1];
    legendre->P_up_derivative[1] = legendre->P_up_derivative[0];
    legendre->P_up_derivative[0] = A*(cos_theta*legendre->P_up_derivative[1] - sin_theta*legendre->P_up[1]);
    return (associated_legendre_pair_t){
      .P_nm = legendre->P_up[0],
      .P_nm_derivative = legendre->P_up_derivative[0]
    };
  }

  // Otherwise we have to use both P_{n-1}^m and P_{n-2}^m
  float const B = sqrtf((float)((n - 1)*(n - 1) - m*m)/(float)(n*n - m*m));
  float const P_nm = A*cos_theta*legendre->P_up[0] - B*legendre->P_up[1];
  float const P_nm_derivative = A*(cos_theta*legendre->P_up_derivative[0] - sin_theta*legendre->P_up[0]) - B*legendre->P_up_derivative[1];

  legendre->P_up[1] = legendre->P_up[0];
  legendre->P_up[0] = P_nm;

  legendre->P_up_derivative[1] = legendre->P_up_derivative[0];
  legendre->P_up_derivative[0] = P_nm_derivative;

  return (associated_legendre_pair_t){
    .P_nm = P_nm,
    .P_nm_derivative = P_nm_derivative
  };
}

magnetic_field_vector_t calculate_model_geomagnetic_field(float latitude, float longitude, float const altitude, float const decimal_year) {
  // If we are exactly at a pole, decrease the latitude slightly such that
  // the returned vector is consistent with the given longitude.
  float const max_latitude = 89.95;
  latitude = fminf(max_latitude, fmaxf(-max_latitude, latitude));

  longitude *= deg_to_rad;
  // theta is the co-latitude.
  float const cos_theta = sinf(latitude*deg_to_rad);
  float const sin_theta = sqrtf(1. - cos_theta*cos_theta);
  float const r = igrf_earth_radius + altitude;
  float const year_factor = decimal_year - 2025.;

  associated_legendre_state_t legendre_state = {0};

  magnetic_field_vector_t field = {0};
  unsigned g_index = 0;
  unsigned h_index = 0;
  for (unsigned m = 0; m <= N; ++m) {
    float const cos_mphi = cosf(m*longitude);
    float const sin_mphi = sinf(m*longitude);
    for (unsigned n = m; n <= N; ++n) {
      if (n == 0) continue;

      float const g_nm = g_coefficients[g_index] + year_factor*g_coefficients_rate[g_index];
      ++g_index;

      float h_nm = 0.;
      if (m != 0) {
        h_nm = h_coefficients[h_index] + year_factor*h_coefficients_rate[h_index];
        ++h_index;
      }

      associated_legendre_pair_t const new_legendre = next_associated_legendre(n, m, &legendre_state, cos_theta, sin_theta);
      // Convenient testing using the C++ std::assoc_legendre function. 
#ifdef CPP_ASSOC_LEGENDRE_TESTING
      float const correct = std::assoc_legendre(n, m, cos_theta)*std::sqrtf((m == 0 ? 1 : 2)*static_cast<float>(std::tgamma(n - m + 1))/static_cast<float>(std::tgamma(n + m + 1)));
      printf("n = %u, m = %u; \t%.5g %.5g %.5g\n", n, m, new_legendre.P_nm, correct, new_legendre.P_nm - correct);
#endif
      float const radial_factor = powf(igrf_earth_radius/r, n + 2);

      field.up += (n + 1)*radial_factor*(g_nm*cos_mphi + h_nm*sin_mphi)*new_legendre.P_nm;
      field.north += radial_factor*(g_nm*cos_mphi + h_nm*sin_mphi)*new_legendre.P_nm_derivative;
      field.east += radial_factor*(g_nm*m*sin_mphi - h_nm*m*cos_mphi)*new_legendre.P_nm/sin_theta;
    }
  }
  
  return field;
}

