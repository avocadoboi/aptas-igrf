#ifndef APTAS_ADCS_GEOMAG70
#define APTAS_ADCS_GEOMAG70

#define PI 3.141592654
#define RAD2DEG (180.0 / PI)

#define RECL 81

#define MAXINBUFF RECL + 14

/** Max size of in buffer **/

#define MAXREAD MAXINBUFF - 2
/** Max to read 2 less than total size (just to be safe) **/

// #define MAXMOD 30
/** Max number of models in a file **/

#define PATH MAXREAD
/** Max path and filename length **/

#define MAXDEG 13
#define MAXCOEFF (MAXDEG * (MAXDEG + 2) + 1) /* index starts with 1!, (from old Fortran?) */

int getshc(char file[PATH], int nmax_of_gh);
int extrapsh(float date, float dte1, int nmax1, int nmax2);
int shval3(float flat, float flon, float elev, int nmax, float *B);

#endif
