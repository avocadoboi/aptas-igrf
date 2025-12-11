#include <aptas-igrf.h>
#include <geomag70.h>

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define __USE_MISC
#include <sys/time.h>
#include <sys/resource.h>

double get_time()
{
  struct timeval t;
  struct timezone tzp;
  gettimeofday(&t, &tzp);
  return t.tv_sec + t.tv_usec*1e-6;
}

static float const equal_epsilon = 0.5;

void load_coefficients(void) {
  printf("Loading IGRF14_minified...\n");
  switch (load_IGRF_coefficients("data/IGRF14_minified")) {
  case IGRFError_UnexpectedFileFormat:
    printf("Error: Unexpected file format.\n");
    exit(EXIT_FAILURE);
  case IGRFError_FileNotFound:
    printf("Error: File not found.\n");
    exit(EXIT_FAILURE);
  case IGRFError_Success:
    printf("Succeeded reading IGRF coefficients.\n");
    break;
  default:
    printf("Unknown error code from load_IGRF_coefficients.\n");
    exit(EXIT_FAILURE);
  }
  getshc("old-implementation/IGRF14.COF", 13);
  extrapsh(2025, 2025, 13, 13);
}

void old_performance(void) {
  double const time_before = get_time();
  for (float latitude = -90.; latitude <= 90.; latitude += 90./40) {
    for (float longitude = 0; longitude < 360.; longitude += 360./40) {
      for (float altitude = 0; altitude < 800; altitude += 800./40) {
        float old_field[3] = {0};
        shval3(latitude, longitude, igrf_earth_radius + altitude, 13, old_field);
      }
    }
  }
  double const time_after = get_time();
  printf("Time for old implementation: %.5g s\n", time_after - time_before);
}
void new_performance(void) {
  double const time_before = get_time();
  for (float latitude = -90.; latitude <= 90.; latitude += 90./40) {
    for (float longitude = 0; longitude < 360.; longitude += 360./40) {
      for (float altitude = 0; altitude < 800; altitude += 800./40) {
        magnetic_field_vector_t const new_field = calculate_model_geomagnetic_field(latitude, longitude, altitude, 2025);
      }
    }
  }
  double const time_after = get_time();
  printf("Time for new implementation: %.5g s\n", time_after - time_before);
}

void compare_new_old(void) {
  for (float latitude = -90.; latitude <= 90.; latitude += 90./5) {
    for (float longitude = 0; longitude < 360.; longitude += 360./5) {
      for (float altitude = 0; altitude < 800; altitude += 800./5) {
        float old_field[3] = {0};
        shval3(latitude, longitude, igrf_earth_radius + altitude, 13, old_field);
        magnetic_field_vector_t const new_field = calculate_model_geomagnetic_field(latitude, longitude, altitude, 2025);
        printf("new - old: (%.5g, %.5g, %.5g) nT\n", new_field.east - old_field[1], new_field.north + old_field[0], new_field.up - old_field[2]);
      }
    }
  }
}

void check_against_test_set(void) {
  FILE* const file = fopen("test/data/test_set", "rt");
  if (file == NULL) {
    printf("Failed reading the test set file.\n");
    exit(EXIT_FAILURE);
  }
  while (true) {
    float latitude, longitude, altitude, year;
    magnetic_field_vector_t expected_field;
    if (fscanf(file, "%f%f%f%f%f%f%f", 
               &latitude, &longitude, &altitude, &year, 
               &expected_field.east, &expected_field.north, &expected_field.up) != 7) {
      if (feof(file)) {
        printf("Test succeeded!\n");
        return;
      }
      fclose(file);
      printf("Error occurred while reading the test set.\n");
      exit(EXIT_FAILURE);
    }


    // printf("lat %.5g lon %.5g alt %.5g km yr %.5g\n", latitude, longitude, altitude, year);
    magnetic_field_vector_t const calculated_field = calculate_model_geomagnetic_field(latitude, longitude, altitude, year);
    // printf("ppigrf: (%.5g, %.5g, %.5g) nT\n", expected_field.east, expected_field.north, expected_field.up);
    // printf("new: (%.5g, %.5g, %.5g) nT\n", calculated_field.east, calculated_field.north, calculated_field.up);
    // printf("old: (%.5g, %.5g, %.5g) nT\n", old_field[1], -old_field[0], old_field[2]);
    // printf("new - old: (%.5g, %.5g, %.5g) nT\n", calculated_field.east - old_field[1], calculated_field.north + old_field[0], calculated_field.up - old_field[2]);

    float const err_east = fabs(calculated_field.east - expected_field.east);
    float const err_north = fabs(calculated_field.north - expected_field.north);
    float const err_up = fabs(calculated_field.up - expected_field.up);
    if (err_east > equal_epsilon || err_north > equal_epsilon || err_up > equal_epsilon) {
      printf("(%.3g, %.3g, %.3g)\n", err_east, err_north, err_up);
      printf("TEST FAILED: too large error.\n");
      exit(EXIT_FAILURE);
    }
  }
}

int main(void) {
  load_coefficients();
  old_performance();
  new_performance();
  old_performance();
  new_performance();
  old_performance();
  new_performance();
  old_performance();
  new_performance();
  old_performance();
  new_performance();
  old_performance();
  new_performance();
  old_performance();
  new_performance();
  old_performance();

  // compare_new_old();
  check_against_test_set();

  printf("Testing poles.\n");
  magnetic_field_vector_t const field_0 = calculate_model_geomagnetic_field(-90, 180., 0., 2025.);
  printf("East: %.5g\tNorth: %.5g\tUp: %.5g\n", field_0.east, field_0.north, field_0.up);
  magnetic_field_vector_t const field_1 = calculate_model_geomagnetic_field(90, 180., 0., 2025.);
  printf("East: %.5g\tNorth: %.5g\tUp: %.5g\n", field_1.east, field_1.north, field_1.up);
}
