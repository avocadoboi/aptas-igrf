#include <aptas-igrf.h>

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double const equal_epsilon = 1e-3;

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
}

void check_against_test_set(void) {
  FILE* const file = fopen("test/data/test_set", "rt");
  if (file == NULL) {
    printf("Failed reading the test set file.\n");
    exit(EXIT_FAILURE);
  }
  while (true) {
    double latitude, longitude, altitude, year;
    magnetic_field_vector_t expected_field;
    if (fscanf(file, "%lf%lf%lf%lf%lf%lf%lf", 
               &latitude, &longitude, &altitude, &year, 
               &expected_field.east, &expected_field.north, &expected_field.up) != 7) {
      if (feof(file)) {
        return;
      }
      fclose(file);
      printf("Error occurred while reading the test set.\n");
      exit(EXIT_FAILURE);
    }
    
    printf("%lf %lf %lf %lf\n", latitude, longitude, altitude, year);
    printf("%lf %lf %lf\n", expected_field.east, expected_field.north, expected_field.up);
    magnetic_field_vector_t const calculated_field = calculate_model_geomagnetic_field(latitude, longitude, altitude, year);
    printf("%lf %lf %lf\n", calculated_field.east, calculated_field.north, calculated_field.up);

    printf(
      "%.3g\t%.3g\t%.3g\n\n", 
      fabs(calculated_field.east - expected_field.east),
      fabs(calculated_field.north - expected_field.north),
      fabs(calculated_field.up - expected_field.up)
    );

  }
}

int main(void) {
  load_coefficients();

  // check_against_test_set();
  magnetic_field_vector_t const field = calculate_model_geomagnetic_field(85.f, 0., 0., 2025.);
  printf("East: %.5g\tNorth: %.5g\tUp: %.5g\n", field.east, field.north, field.up);
}
