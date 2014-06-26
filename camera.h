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

#ifndef INC_CAMERA_H
#define INC_CAMERA_H

typedef struct {
	double pos[3]; // Camera position
	double *trans; // Transform
	double *itrans; // Inverse transform

	// Get the principal values from photogrammetry calibration of camera.
	int dims[2]; // Pixel dimensions
	double pxsize; // Size of a pixel in mm
	double prdist; // Pricipal distance of camera lens in mm
	double soptc[2]; // Sensor optical centre in mm
	double rdisto[2]; // Radial distortion parameters k3, k5
} Camera;

void initcamera(Camera *c);
void freecamera(Camera *c);
int objpixsize(const Camera *c, double objsize, double objdist);
void locatecam(Camera *c, const double *dots, const double *dotsep,
		double distguess);
void pospix(const Camera *c, const double *vec, int *pix);
void pixdir(const Camera *c, int x, int y, double *dir);

#endif /* INC_CAMERA_H */
