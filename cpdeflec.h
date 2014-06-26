/*
 *  Copyright (C) Paul Scott 2011
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INC_CPDEFLEC_H
#define INC_CPDEFLEC_H

#include "commath.h"
#include "camera.h"
#include "pattern.h"
#include "fitting.h"

typedef struct {
	double dotsize; // Physical size of reference dots in mm
	double dotsep[3]; // Distance between TL-BL, BL-TR, TR-TL
	double corns[4][3]; // Mirror boundary corners TL, BL, BR, TR
	double distguess; // Guess of distance from camera to mirror
	double startdepth; // Height of starting point
} Mirror;

void solveprofile(Camera camera, Pattern pattern, Mirror mirror,
		const char *imfnh, const char *imfnv, const int *idots,
		const char *outfn);

#endif /* INC_CPDEFLEC_H */
