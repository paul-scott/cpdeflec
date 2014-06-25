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
#include "commath.h"

#include <stdlib.h>
#include <stdio.h>
#include <tiffio.h>
#include <math.h>

double gethue(const uint32_t *rgb);
double hueshift(const Pattern *p, const double hue, const int orien);

void initpattern(Pattern *p, const char *relfn)
{
	p->rel = malloc(p->relbins*sizeof(*p->rel));
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
		if (lcount >= p->relbins) break;
		sscanf(line, "%lf", (p->rel+lcount));
		lcount++;
	}
	fclose(relfile);
}

void freepattern(Pattern *p)
{
	free(p->rel);
	p->rel = NULL;
}

void transpattvec(const Pattern *p, double *vec)
{
	double temp[3];
	// Need to set z component to zero. 
	*(vec) = *(vec) + p->xoff;
	*(vec+1) = *(vec+1) + p->yoff;
	*(vec+2) = 0.0;
	//printf("pat: %f, %f, %f ", *vec, *(vec+1), *(vec+2));
	matxvec((double *) p->trans, vec, temp);
	*(vec) = temp[0] + p->pos[0];
	*(vec+1) = temp[1] + p->pos[1];
	*(vec+2) = temp[2] + p->pos[2];
	//printf("glo: %f, %f, %f\n", *vec, *(vec+1), *(vec+2));
}

double getdist(const Pattern *p, const uint32_t *rgb, const double pdist,
		const int orien)
{
	// Might like to think of not changing pdist if we get a grey pixel.
	double shift = hueshift(p, gethue(rgb), orien);
	double pshift = fmod(pdist, p->segsize);
	//printf("s, ps, %f, %f\n", shift, pshift);

	if (shift < (pshift-0.5*p->segsize)) {
		return (pdist - pshift + p->segsize + shift);
	} else if (shift > (pshift+0.5*p->segsize)) {
		return (pdist - pshift - p->segsize + shift);
	} else {
		return (pdist - pshift + shift);
	}
}

double gethue(const uint32_t *rgb)
{
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

double hueshift(const Pattern *p, const double hue, const int orien)
{
	// MIGHT NEED TO CHECK....
	int indice = (int) round(hue*p->relbins);

	if (indice == p->relbins) indice = 0;
	// Need to reverse relation depending on orientation.
	if (orien) {
		return p->segsize*(*(p->rel+indice));
	} else {
		return p->segsize*(1.0-(*(p->rel+indice)));
	}
}
