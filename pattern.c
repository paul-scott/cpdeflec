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

#include "pattern.h"

/* Pattern parameters
 * ******************
 * Get locating values from photogrammetry. Position is horizontal edge of
 * horiz rotated pattern, and vertical edge of vert rotated pattern.
 */
// NOTE NEED TO FIND BETTER VALUE.
static double patpos[3] = {635.77408,-364.31168,-2334.4864}; // Pos of
static double patxoff = 4.6407583;
static double patyoff = 3.2714379;
// pattern corner
// NOTE MIGHT WANT TO WORK ON NON ORTHOGONAL VECTORS.
static double pattrans[3][3] = {{0.8542358,0.0018653,-0.5198823},
	{0.0015481,0.99998,0.0061316},
	{0.5198834,-0.0060427,0.8542159}}; // Coordinate system translation
double segsize = 47.7; // Width of repeating pattern segment
static int relbins = 274;
static double *rel;

double gethue(const uint32 *rgb);
double hueshift(const double hue, const int orien);

void initpattern(const char *relfn) {
	rel = malloc(relbins*sizeof(*rel));
	const size_t buffsize = 50;
	char line[buffsize];
	int lcount = 0;

	// Need to drag pattern relation in from file.
	FILE *relfile = fopen(relfn, "r");
	if (relfile==NULL) {
		printf("Failed to open hue relation file, exiting...\n");
		exit(1);
	}
	while (fgets(line, buffsize, relfile) != NULL) {
		if (lcount >= relbins) break;
		sscanf(line, "%f", (rel+lcount));
		lcount++;
	}
	fclose(relfile);
}

void freepattern() {
	free(rel);
	rel = NULL;
}

void transpattvec(double *vec) {
	double temp[3];
	// Need to set z component to zero. 
	*(vec) = *(vec) + patxoff;
	*(vec+1) = *(vec+1) + patyoff;
	*(vec+2) = 0.0;
	//printf("pat: %f, %f, %f ", *vec, *(vec+1), *(vec+2));
	matxvec((double *) pattrans, vec, temp);
	*(vec) = temp[0] + patpos[0];
	*(vec+1) = temp[1] + patpos[1];
	*(vec+2) = temp[2] + patpos[2];
	//printf("glo: %f, %f, %f\n", *vec, *(vec+1), *(vec+2));
}

double getdist(const uint32 *rgb, const double pdist, const int orien) {
	// Might like to think of not changing pdist if we get a grey pixel.
	double shift = hueshift(gethue(rgb), orien);
	double pshift = fmod(pdist, segsize);
	//printf("s, ps, %f, %f\n", shift, pshift);

	if (shift < (pshift-0.5*segsize)) {
		return (pdist - pshift + segsize + shift);
	} else if (shift > (pshift+0.5*segsize)) {
		return (pdist - pshift - segsize + shift);
	} else {
		return (pdist - pshift + shift);
	}
}

double gethue(const uint32 *rgb) {
	double R = (double) TIFFGetR(*rgb)/255.0;
	double G = (double) TIFFGetG(*rgb)/255.0;
	double B = (double) TIFFGetB(*rgb)/255.0;

	double maxc = fmax(R,fmax(G,B));
	double minc = fmin(R,fmin(G,B));
	// double L = (maxc + minc)/2.0;
	double hue = 0.0;
	
	if (maxc == minc) {
		hue = 0.0;
		//printf("Grey pixel.\n");
	} else if (R == maxc) {
		hue = fmod(1.0+(G-B)/(6.0*(maxc-minc)),1.0);
	} else if (G == maxc) {
		hue = (2.0 + (B-R)/(maxc-minc))/6.0;
	} else {
		hue = (4.0 + (R-G)/(maxc-minc))/6.0;
	}

	//printf("hue: %f, %f, %f, %f ", hue, R*255.0, G*255.0, B*255.0);
	return hue;
}

double hueshift(const double hue, const int orien) {
	// MIGHT NEED TO CHECK....
	int indice = (int) round(hue*relbins);

	if (indice == relbins) indice = 0;
	// Need to reverse relation depending on orientation.
	if (orien) {
		return segsize*(*(rel+indice));
	} else {
		return segsize*(1.0-(*(rel+indice)));
	}
}
