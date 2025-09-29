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
  double east;
  double north;
  double up;
} magnetic_field_vector_t;

static double const igrf_earth_radius = 6371.2; // [km]

/**
 * Reads Schmidt quasi-normal internal spherical harmonic coefficients from a specified IGRF model file into an array
 * and extrapolates the coefficients from the date of the base model (which is set in the constant IGRF_base_year) to a given date.
 * \param file_name Null-terminated file name string.
 * \param extrapolated_year The date, given in decimal years, to extrapolate the model to.
 * \param coordinate_system Whether to use geodetic or geocentric coordinate system.
 * \note This function must be called before using calculate_model_geomagnetic_field.
 */
IGRFError load_IGRF_coefficients(char const* file_name);


// void extrapolate_IGRF_coefficients(double year);

/**
 * Calculates the geomagnetic field at an ECEF point based on the IGRF model.
 * \param latitude North latitude in degrees.
 * \param longitude East longitude in degrees.
 * \param altitude Given in metres.
 * \return Geomagnetic field vector where x is eastward, y is northward, and z is vertically upwards from Earth.
 * \note load_IGRF_coefficients must be called before trying to use this function.
 */
magnetic_field_vector_t calculate_model_geomagnetic_field(double latitude, double longitude, double altitude, double decimal_year);

#endif
