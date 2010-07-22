#ifndef INC_PATTERN_H
#define INC_PATTERN_H

#include <stdio.h>
#include <stdlib.h>
#include <tiffio.h>
#include <math.h>

#include "commath.h"

extern double segsize;

void initpattern(const char *relfn);
void freepattern();
double getdist(const uint32 *pix, const double pdist, const int orien);
void transpattvec(double *vec);

#endif /* INC_PATTERN_H */
