#include "aptas-igrf.h"

#include "math_constants.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

// #define CPP_ASSOC_LEGENDRE_TESTING

#ifdef CPP_ASSOC_LEGENDRE_TESTING
#include <array>
#include <cmath>
#endif

// The maximum harmonic degree in the IGRF model.
#define N 13

#define G_COEFFICIENT_COUNT N*(N + 3)/2
#define H_COEFFICIENT_COUNT N*(N + 1)/2

// Stored column-by column starting from m = 0 with n increasing for each m.
static double g_coefficients[G_COEFFICIENT_COUNT];
static double g_coefficients_rate[G_COEFFICIENT_COUNT];

static double h_coefficients[H_COEFFICIENT_COUNT];
static double h_coefficients_rate[H_COEFFICIENT_COUNT];

IGRFError load_IGRF_coefficients(char const* const file_name) {
  FILE* const file = fopen(file_name, "rt");
  if (file == NULL) {
    return IGRFError_FileNotFound;
  }
  for (int n = 1; n <= N; ++n) {
    for (int m = 0; m <= n; ++m) {
      // See equation (??) in the accompanying PDF.
      int const g_index = n - 1 + m*N - m*(m - 1)/2; // Integer division here is ok, m*(m-1) is always even.

      if (fscanf(file, "%lf%lf", g_coefficients + g_index, g_coefficients_rate + g_index) != 2) {
        fclose(file);
        return IGRFError_UnexpectedFileFormat;
      }
      // printf("g_%i^%i: %lf %lf\n", n, m, g_coefficients[g_index], g_coefficients_rate[g_index]);
      if (m != 0) {
        int const h_index = n - 1 - N + m*N - m*(m - 1)/2;

        if (fscanf(file, "%lf%lf", h_coefficients + h_index, h_coefficients_rate + h_index) != 2) {
          fclose(file);
          return IGRFError_UnexpectedFileFormat;
        }
        // printf("h_%i^%i: %lf %lf\n", n, m, h_coefficients[h_index], h_coefficients_rate[h_index]);
      }
    }
  }
  return IGRFError_Success;
}

typedef struct {
  double P_up[2];
  double P_up_derivative[2];

  double P_nn;
  double P_nn_derivative;
} associated_legendre_state_t;

typedef struct {
  double P_nm;
  double P_nm_derivative;
} associated_legendre_pair_t;

static associated_legendre_pair_t next_associated_legendre(
  unsigned const n, unsigned const m, associated_legendre_state_t* const legendre, 
  double const cos_theta, double const sin_theta
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
    double const A = sqrt((n == 1 ? 2 : 1)*(1 - 1./(double)(2*n)));
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
  double const A = (double)(2*n - 1)/sqrt((double)(n*n - m*m));

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
  double const B = sqrt((double)((n - 1)*(n - 1) - m*m)/(double)(n*n - m*m));
  double const P_nm = A*cos_theta*legendre->P_up[0] - B*legendre->P_up[1];
  double const P_nm_derivative = A*(cos_theta*legendre->P_up_derivative[0] - sin_theta*legendre->P_up[0]) - B*legendre->P_up_derivative[1];

  legendre->P_up[1] = legendre->P_up[0];
  legendre->P_up[0] = P_nm;

  legendre->P_up_derivative[1] = legendre->P_up_derivative[0];
  legendre->P_up_derivative[0] = P_nm_derivative;

  return (associated_legendre_pair_t){
    .P_nm = P_nm,
    .P_nm_derivative = P_nm_derivative
  };
}

magnetic_field_vector_t calculate_model_geomagnetic_field(double const latitude, double longitude, double const altitude, double const decimal_year) {
  longitude *= deg_to_rad;
  // theta is the co-latitude.
  double const cos_theta = cos((90. - latitude)*deg_to_rad);
  double const sin_theta = sqrt(1. - cos_theta*cos_theta);
  double const r = igrf_earth_radius + altitude;
  double const year_factor = decimal_year - 2025.;

  // Does not need to be initialized now. 
  associated_legendre_state_t legendre_state;

  magnetic_field_vector_t field = {0};
  unsigned g_index = 0;
  unsigned h_index = 0;
  for (unsigned m = 0; m <= N; ++m) {
    double const cos_mphi = cos(m*longitude);
    double const sin_mphi = sin(m*longitude);
    for (unsigned n = m; n <= N; ++n) {
      if (n == 0) continue;

      double const g_nm = g_coefficients[g_index] + year_factor*g_coefficients_rate[g_index];
      ++g_index;

      double h_nm = 0.;
      if (m != 0) {
        h_nm = h_coefficients[h_index] + year_factor*h_coefficients_rate[h_index];
        ++h_index;
      }

      associated_legendre_pair_t const new_legendre = next_associated_legendre(n, m, &legendre_state, cos_theta, sin_theta);
      // Convenient testing using the C++ std::assoc_legendre function. 
#ifdef CPP_ASSOC_LEGENDRE_TESTING
      double const correct = std::assoc_legendre(n, m, cos_theta)*std::sqrt((m == 0 ? 1 : 2)*static_cast<double>(std::tgamma(n - m + 1))/static_cast<double>(std::tgamma(n + m + 1)));
      printf("n = %u, m = %u; \t%.5g %.5g %.5g\n", n, m, new_legendre.P_nm, correct, new_legendre.P_nm - correct);
#endif
      double const radial_factor = pow(igrf_earth_radius/r, n + 2);

      field.up += (n + 1)*radial_factor*(g_nm*cos_mphi + h_nm*sin_mphi)*new_legendre.P_nm;
      field.north += radial_factor*(g_nm*cos_mphi + h_nm*sin_mphi)*new_legendre.P_nm_derivative;
      // East/west doesn't make sense if we're right at a pole and we get division by 0.
      if (sin_theta > 1e-8) {
        field.east += radial_factor*(g_nm*m*sin_mphi - h_nm*m*cos_mphi)*new_legendre.P_nm/sin_theta;
      }
    }
  }
  
  return field;
}

