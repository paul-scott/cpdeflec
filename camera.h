#ifndef INC_CAMERA_H
#define INC_CAMERA_H

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>

#include "commath.h"

extern double campos[3];

void initcamera();
void freecamera();
int objpixsize(double objsize, double objdist);
int locatecam(double *dots, double *dotsep, double distguess);
void findpix(const double *vec, int *pix);
void finddirpix(const int x, const int y, double *dir);

#endif /* INC_CAMERA_H */
