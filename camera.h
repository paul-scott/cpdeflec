#ifndef INC_CAMERA_H
#define INC_CAMERA_H

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>

#include "commath.h"

extern float campos[3];

void initcamera();
void freecamera();
int objpixsize(float objsize, float objdist);
int locatecam(float *dots, float *dotsep, float distguess);
void findpix(const float *vec, int *pix);
void finddirpix(const int x, const int y, float *dir);

#endif /* INC_CAMERA_H */
