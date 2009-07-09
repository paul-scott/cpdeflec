#ifndef INC_PATTERN_H
#define INC_PATTERN_H

#include <tiffio.h>
#include <math.h>

#include "commath.h"

extern float segsize;

void initpattern();
float getdist(const uint32 *pix, const float pdist);

#endif /* INC_PATTERN_H */
