#include "aptas-igrf.h"

#include "math_constants.h"

#include <math.h>
#include <stdio.h>

#define MAX_HARMONIC_DEGREE 13

#define G_COEFFICIENT_COUNT MAX_HARMONIC_DEGREE*(MAX_HARMONIC_DEGREE + 3)/2
#define H_COEFFICIENT_COUNT MAX_HARMONIC_DEGREE*(MAX_HARMONIC_DEGREE + 1)/2

static float g_coefficients[G_COEFFICIENT_COUNT];
static float g_coefficients_rate[G_COEFFICIENT_COUNT];
static float extrapolated_g_coefficients[G_COEFFICIENT_COUNT];

static float h_coefficients[H_COEFFICIENT_COUNT];
static float h_coefficients_rate[H_COEFFICIENT_COUNT];
static float extrapolated_h_coefficients[H_COEFFICIENT_COUNT];

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


// void extrapolate_IGRF_coefficients(float const year) {
//   int g_index = 0;
//   int h_index = 0;
//   for (int n = 1; n <= MAX_HARMONIC_DEGREE; ++n) {
//     for (int m = 0; m <= n; ++m) {
//       extrapolated_g_coefficients[g_index] = g_coefficients[g_index] + (year - 2025.f)*g_coefficients_rate[g_index];
//       ++g_index;
//       if (m != 0) {
//         extrapolated_h_coefficients[h_index] = h_coefficients[h_index] + (year - 2025.f)*h_coefficients_rate[h_index];
//         ++h_index;
//       }
//     }
//   }
// }

static float associated_legendre_propagate_right(int const n, int const m, float const P_n_m1, float const P_n_m2, float const cos_theta, float const sin_theta) {
  float const A = 2.f*(float)(m - 1)/sqrtf((float)((n - m + 1)*(n + m)));
  float const B = sqrtf((float)((n - m + 2)*(n + m - 1))/(float)((n + m)*(n - m + 1)));
  return A*cos_theta/sin_theta*P_n_m1 - B*P_n_m2;
}
static float associated_legendre_propagate_down(int const n, int const m, float const P_n1_m, float const P_n2_m, float const cos_theta) {
  float const A = (float)(2*n - 1)/sqrtf((float)(n^2 - m^2));
  float const B = sqrtf((float)((n - 1)^2 - m)/(float)(n^2 - m^2));
  return A*cos_theta*P_n1_m - B*P_n2_m;
}

magnetic_field_vector_t calculate_model_geomagnetic_field(float const latitude, float const longitude, float const altitude, float const decimal_year) {
  magnetic_field_vector_t field = {0};
  float const cos_theta = cosf((90.f - latitude)*deg_to_rad);
  float const sin_theta = sqrtf(1 - cos_theta*cos_theta);
  float const r = igrf_earth_radius + altitude;

  float const P_0_0 = 1;
  float const P_1_0 = cos_theta;
  float const P_1_1 = sin_theta;
  float const P_2_1 = sqrt_3*cos_theta*sin_theta;
  // float P_n1_m = 

  int g_index = 0;
  int h_index = 0;
  for (int n = 1; n <= MAX_HARMONIC_DEGREE; ++n) {
    for (int m = 0; m <= n; ++m) {
      float const g_nm = g_coefficients[g_index] + (decimal_year - 2025.f)*g_coefficients_rate[g_index];
      ++g_index;

      field.north = 1/r*igrf_earth_radius;
      if (m != 0) {
        float const h_nm = h_coefficients[h_index] + (decimal_year - 2025.f)*h_coefficients_rate[h_index];
        ++h_index;
      }
    }
  }
  
  return field;
}


// enum IGRFError load_IGRF_coefficients(char const* const file_name, double const extrapolated_year) {
//   FILE* const stream = fopen(file_name, "rt");
//   if (stream == NULL) {
//     return IGRFError_FileNotFound;
//   }
//   else {
//     int ii = 0;
    // fseek(stream, strec, SEEK_SET);
    // for (int nn = 1; nn <= MAX_HARMONIC_DEGREE; ++nn) {
      // for (int mm = 0; mm <= nn; ++mm) {
        // char in_buffer[buffer_size];
        // int m, n;

        // char irat[9];
        // int line_num;
        // double trash;
        // double g, hh;
        // if (iflag == 1) {
        //   fgets(in_buffer, MAXREAD, stream);
        //   sscanf(in_buffer, "%d%d%lg%lg%lg%lg%s%d",
        //           &n, &m, &g, &hh, &trash, &trash, irat, &line_num);
        // }
        // else {
        //   fgets(in_buffer, MAXREAD, stream);
        //   sscanf(in_buffer, "%d%d%lg%lg%lg%lg%s%d",
        //           &n, &m, &trash, &trash, &g, &hh, irat, &line_num);
        // }
        // if ((nn != n) || (mm != m)) {
        //   ios = -2;
        //   fclose(stream);
        //   return(ios);
        // }
        // ii = ii + 1;
        // switch(gh) {
        //   case 1:  gh1[ii] = g;
        //     break;
        //   case 2:  gh2[ii] = g;
        //     break;
        //   default: printf("\nError in subroutine getshc");
        //     break;
        // }
        // if (m != 0) {
        //     ii = ii+ 1;
        //     switch(gh) {
        //       case 1:  gh1[ii] = hh;
        //         break;
        //       case 2:  gh2[ii] = hh;
        //         break;
        //       default: printf("\nError in subroutine getshc");
        //         break;
        //     }
        // }
  //     }
  //   }
  // }
  // fclose(stream);
  // return(ios);
// }

// vec3_t calculate_model_geomagnetic_field(float latitude, float longitude, float altitude) {
//   double const earths_radius = 6371.2;
//   double bb, cc, dd;
//   double sd;
//   double cd;
//   double rr;
//   double fm,fn;
//   double sl[14];
//   double cl[14];
//   double p[119];
//   double q[119];
//   int ii,j,k,l,m,n;
//   int npq;
//   double argument;
//   double power;
//   double const a2 = 40680631.59;            /* WGS84 */
//   double const b2 = 40408299.98;            /* WGS84 */
//   argument = latitude * deg_to_rad;
//   double slat = sin(argument);
//
//   double aa;
//   if (90. - latitude < 0.001) {
//     aa = 89.999; //  300 ft. from North pole 
//   }
//   else if (90. + latitude < 0.001) {
//     aa = -89.999; //  300 ft. from South pole 
//   }
//   else {
//     aa = latitude;
//   }
//
//   argument = aa * deg_to_rad;
//   double clat = cos(argument);
//   argument = longitude * deg_to_rad;
//   sl[1] = sin(argument);
//   cl[1] = cos(argument);
//   switch(gh) {
//     case 3:  x = 0;
//       y = 0;
//       z = 0;
//       break;
//     case 4:  xtemp = 0;
//       ytemp = 0;
//       ztemp = 0;
//       break;
//     default: printf("\nError in subroutine shval3");
//       break;
//   }
//   sd = 0.0;
//   cd = 1.0;
//   l = 1;
//   n = 0;
//   m = 1;
//   npq = (nmax * (nmax + 3)) / 2;
//   double ratio = earths_radius / altitude;
//   argument = 3.0;
//   aa = sqrt(argument);
//   p[1] = 2. * slat;
//   p[2] = 2. * clat;
//   p[3] = 4.5 * slat * slat - 1.5;
//   p[4] = 3. * aa * clat * slat;
//   q[1] = -clat;
//   q[2] = slat;
//   q[3] = -3. * clat * slat;
//   q[4] = aa * (slat * slat - clat * clat);
//   for (k = 1; k <= npq; ++k) {
//       if (n < m) {
//           m = 0;
//           n = n + 1;
//           argument = ratio;
//           power =  n + 2;
//           rr = pow(argument,power);
//           fn = n;
//       }
//       fm = m;
//       if (k >= 5) {
//           if (m == n) {
//               argument = (1. - 0.5/fm);
//               aa = sqrt(argument);
//               j = k - n - 1;
//               p[k] = (1. + 1.0/fm) * aa * clat * p[j];
//               q[k] = aa * (clat * q[j] + slat/fm * p[j]);
//               sl[m] = sl[m-1] * cl[1] + cl[m-1] * sl[1];
//               cl[m] = cl[m-1] * cl[1] - sl[m-1] * sl[1];
//           }
//           else {
//               argument = fn*fn - fm*fm;
//               aa = sqrt(argument);
//               argument = ((fn - 1.0)*(fn-1.0)) - (fm * fm);
//               bb = sqrt(argument)/aa;
//               cc = (2. * fn - 1.0)/aa;
//               ii = k - n;
//               j = k - 2 * n + 1;
//               p[k] = (fn + 1.0) * (cc * slat/fn * p[ii] - bb/(fn - 1.0) * p[j]);
//               q[k] = cc * (slat * q[ii] - clat/fn * p[ii]) - bb * q[j];
//           }
//       }
//       switch(gh) {
//         case 3:  aa = rr * gha[l];
//           break;
//         case 4:  aa = rr * ghb[l];
//           break;
//         default: printf("\nError in subroutine shval3");
//           break;
//       }
//       if (m == 0) {
//           switch(gh) {
//             case 3:  x = x + aa * q[k];
//               z = z - aa * p[k];
//               break;
//             case 4:  xtemp = xtemp + aa * q[k];
//               ztemp = ztemp - aa * p[k];
//               break;
//             default: printf("\nError in subroutine shval3");
//               break;
//           }
//           l = l + 1;
//       }
//       else {
//           switch(gh) {
//             case 3:  bb = rr * gha[l+1];
//               cc = aa * cl[m] + bb * sl[m];
//               x = x + cc * q[k];
//               z = z - cc * p[k];
//               if (clat > 0) {
//                   y = y + (aa * sl[m] - bb * cl[m]) *
//                     fm * p[k]/((fn + 1.0) * clat);
//               }
//               else {
//                   y = y + (aa * sl[m] - bb * cl[m]) * q[k] * slat;
//               }
//               l = l + 2;
//               break;
//             case 4:  bb = rr * ghb[l+1];
//               cc = aa * cl[m] + bb * sl[m];
//               xtemp = xtemp + cc * q[k];
//               ztemp = ztemp - cc * p[k];
//               if (clat > 0) {
//                   ytemp = ytemp + (aa * sl[m] - bb * cl[m]) *
//                     fm * p[k]/((fn + 1.0) * clat);
//               }
//               else {
//                   ytemp = ytemp + (aa * sl[m] - bb * cl[m]) *
//                     q[k] * slat;
//               }
//               l = l + 2;
//               break;
//             default: printf("\nError in subroutine shval3");
//               break;
//           }
//       }
//       m = m + 1;
//   }
//   if (iext != 0) {
//       aa = ext2 * cl[1] + ext3 * sl[1];
//       switch(gh) {
//         case 3:   x = x - ext1 * clat + aa * slat;
//           y = y + ext2 * sl[1] - ext3 * cl[1];
//           z = z + ext1 * slat + aa * clat;
//           break;
//         case 4:   xtemp = xtemp - ext1 * clat + aa * slat;
//           ytemp = ytemp + ext2 * sl[1] - ext3 * cl[1];
//           ztemp = ztemp + ext1 * slat + aa * clat;
//           break;
//         default:  printf("\nError in subroutine shval3");
//           break;
//       }
//   }
//   switch(gh) {
//     case 3:
//       aa = x;
//       x = x * cd + z * sd;
//       z = z * cd - aa * sd;
//       break;
//     case 4:
//       aa = xtemp;
//       xtemp = xtemp * cd + ztemp * sd;
//       ztemp = ztemp * cd - aa * sd;
//       break;
//     default:
//       printf("\nError in subroutine shval3");
//       break;
//   }
// }
