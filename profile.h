#ifndef INC_PROFILE_H
#define INC_PROFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <tiffio.h>

#include "commath.h"
#include "camera.h"
#include "pattern.h"
#include "fitting.h"

/* STRUCTURES
 * **********
 */

/* FUNCTIONS
 * *********
 */
int *intalloc(size_t rs, size_t cs);
void set2darray(int *marray, int r, int c, int cs, int val);
void solveprofile(char *imfnh, char *imfnv, int *idots, char *outfn,
		char *calbfn);

#endif /* INC_PROFILE_H */
