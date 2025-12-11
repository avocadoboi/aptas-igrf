/*

APTAS adaptation of the geomag70 code found at https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field .
It contains functions to read and use spherical harmonic coefficients from an IGRF (International Geomagnetic Reference Field) model file.

*/

#ifndef APTAS_GEOMAG_H
#define APTAS_GEOMAG_H

static int const IGRF_base_year = 2025;

typedef enum {
	IGRFError_Success,
	IGRFError_FileNotFound,
  IGRFError_UnexpectedFileFormat
} IGRFError;

typedef struct {
  float east;
  float north;
  float up;
} magnetic_field_vector_t;

static float const igrf_earth_radius = 6371.2; // [km]

/**
 * Reads Schmidt semi-normalized spherical harmonic coefficients from a specified IGRF model file into an array
 * \param file_name Null-terminated file name string.
 * \note This function must be called before using calculate_model_geomagnetic_field.
 */
IGRFError load_IGRF_coefficients(char const* file_name);

/**
 * Calculates the geomagnetic field at an ECEF point and a date given as a decimal year, based on the IGRF model.
 * \param latitude North latitude in degrees.
 * \param longitude East longitude in degrees.
 * \param altitude Given in kilometres.
 * \param extrapolated_year The date, given in decimal years, to extrapolate the model to.
 * \return Geomagnetic field vector.
 * \note load_IGRF_coefficients must be called before trying to use this function.
 */
magnetic_field_vector_t calculate_model_geomagnetic_field(float latitude, float longitude, float altitude, float decimal_year);

#endif
