#ifndef INC_PATTERN_H
#define INC_PATTERN_H

#include <stdio.h>
#include <stdlib.h>
#include <tiffio.h>
#include <math.h>

#include "commath.h"

extern float segsize;

void initpattern(const char *relfn);
void freepattern();
float getdist(const uint32 *pix, const float pdist, const int orien);
void transpattvec(float *vec);

#endif /* INC_PATTERN_H */
