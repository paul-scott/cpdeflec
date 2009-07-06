#ifndef INC_CAMERA_H
#define INC_CAMERA_H

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "commath.h"

void initcamera();
void freecamera();
int objpixsize(float objsize, float objdist);
int locatecam(float *dots, float *dotsep, float distguess);
void findpix(float *vec, int *pix);
void finddirpix(int *pix, float *dir);

#endif /* INC_CAMERA_H */
