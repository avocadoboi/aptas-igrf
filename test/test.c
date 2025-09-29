#include <aptas-igrf.h>

#include <stdio.h>
#include <stdlib.h>

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

int main(void) {
  load_coefficients();

  calculate_model_geomagnetic_field(27.f, 0, 0, 2025.5f);
}
