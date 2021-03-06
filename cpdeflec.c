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

#include "cpdeflec.h"

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <tiffio.h>

void centroid(const uint32_t *im, uint32_t w, uint32_t h, const int *idot,
		double *dot, int dotpx);
void calcpattvecs(const Pattern *pat, const uint32_t *im, const int iw,
		double *vecs, const int bw, const int *ipix, const int *yb,
		const int *offs, const double idist, const int orien);
void transformpatt(const Pattern *pat, double *vecs, const int bw, const int* yb);
double extrapolate(const Camera *cam, const double *ppos, const double *psnrm,
		const double *ch);
void solvesurface(const Camera *cam, double *vecs, double *poss, const int bw,
		const int *ipix, const int *yb, const int *offs, const double idepth);
void surfnorm(const double *vin, const double *vre, double *snrm);
void slopeerror(const double *vecs, const double *poss, const int bw,
		const int *yb, const int shape, const double *params, double *serr,
		double *sestats);
void pnts_line(const int *p1, const int *p2, Polynom *pol, double *coeffs);
void end_ybounds(const Polynom *l1, const Polynom *l2, const Polynom *l3,
		const Polynom *l4, const int *cornt, const int *cornb, int *ybnd);

void solveprofile(Camera camera, Pattern pattern, Mirror mirror,
		const char *imfnh, const char *imfnv, const int *idots,
		const char *outfn)
{
	TIFF *image; // TIFF file pointer
	uint32_t wid, hei; // Image width and height
	uint32_t *imdata; // Pixel data loaded from TIFF image, y major
	int dotpx; // Approximate width of ref dot in pixels
	double dots[3*2]; // Sub-pixel ref dot locations
	int pxcorns[4*2];

	initcamera(&camera);
	initpattern(&pattern);

	image = TIFFOpen(imfnh, "r");
	if (image == NULL) {
		printf("Failed to open horizontal image.\n");
		exit(1);
	}

	// Load in image data.
	TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &wid);
	TIFFGetField(image, TIFFTAG_IMAGELENGTH, &hei);
	printf("Image dimensions: %d, %d\n", wid, hei);

	imdata = (uint32_t *) _TIFFmalloc(wid*hei*sizeof(uint32_t)); 
	if (imdata == NULL) {
		printf("Failed to allocate memory for horizontal image.\n");
		exit(1);
	}
	if (TIFFReadRGBAImageOriented(image, wid, hei, imdata, ORIENTATION_TOPLEFT,
			0) == 0) {
		printf("Failed to read horizontal image.\n");
		exit(1);
	}
	printf("Horizontal image loaded\n");

	TIFFClose(image);
	image = NULL;

	// Locating camera position
	// ************************

	// Find centre of dots.
	dotpx = objpixsize(&camera, mirror.dotsize, mirror.distguess);
	for (int i=0; i<3; i++) {
		centroid(imdata, wid, hei, (idots+i*2), (dots+i*2), dotpx);
	}

	printf("Centroided dots:\n");
	printf("%f, %f\n%f, %f\n%f, %f\n", dots[0], dots[1], dots[2], dots[3],
			dots[4], dots[5]);

	// Locating position of camera.
	locatecam(&camera, dots, mirror.dotsep, mirror.distguess);
	printf("Camera successfully located.\n");

	printf("Bounding corners:\n");
	// Finding pixels that correspond to mirror bounding corners. 
	for (int i=0; i<4; i++) {
		pospix(&camera, mirror.corns[i], (pxcorns+i*2));
		printf("%d, %d\n", pxcorns[i*2], pxcorns[i*2+1]);
	}

	// Identifying region of interest and finding bounds
	// *************************************************

	// Determining max and min bounds.
	int xpixmax = imax(pxcorns[0], imax(pxcorns[2],
					imax(pxcorns[4], pxcorns[6])));
	int ypixmax = imax(pxcorns[1], imax(pxcorns[3],
					imax(pxcorns[5], pxcorns[7])));
	int xpixmin = imin(pxcorns[0], imin(pxcorns[2],
					imin(pxcorns[4], pxcorns[6])));
	int ypixmin = imin(pxcorns[1], imin(pxcorns[3],
					imin(pxcorns[5], pxcorns[7])));

	int xlen = xpixmax - xpixmin + 1;
	int ylen = ypixmax - ypixmin + 1;
	int offset[2] = {xpixmin, ypixmin};

	// Subtracting offset from pxcorns.
	// Need to use offset for indices of camera and imdata.
	for (int i=0; i<4; i++) {
		pxcorns[i*2+0] = pxcorns[i*2+0] - offset[0];
		pxcorns[i*2+1] = pxcorns[i*2+1] - offset[1];
	}

	printf("Size of area: %d, %d\n", xlen, ylen);
	printf("Offset: %d, %d\n", offset[0], offset[1]);
	printf("Starting pixel value: %d, %d, %d\n",
			TIFFGetR(imdata[offset[1]*wid+offset[0]]),
			TIFFGetG(imdata[offset[1]*wid+offset[0]]),
			TIFFGetB(imdata[offset[1]*wid+offset[0]]));

	double tcoeffs[2];
	double bcoeffs[2];
	double lcoeffs[2];
	double rcoeffs[2];
	Polynom tline = {1, tcoeffs};
	Polynom bline = {1, bcoeffs};
	Polynom lline = {1, lcoeffs};
	Polynom rline = {1, rcoeffs};
	int ybound[2*xlen]; // Alternating top and bottom y bounds (C99 VLA)

	// Finding polynomial to represent the 4 bounding lines
	pnts_line(pxcorns+3*2, pxcorns+0*2, &tline, tcoeffs);
	pnts_line(pxcorns+2*2, pxcorns+1*2, &bline, bcoeffs);
	pnts_line(pxcorns+1*2, pxcorns+0*2, &lline, lcoeffs);
	pnts_line(pxcorns+3*2, pxcorns+2*2, &rline, rcoeffs);

	// Calculating y bounds for ends of region
	end_ybounds(&tline, &lline, &lline, &bline, pxcorns+0*2, pxcorns+1*2, ybound);
	end_ybounds(&rline, &bline, &tline, &rline, pxcorns+3*2, pxcorns+2*2, ybound);

	// Calculating y bounds for centre part of region of interest
	int startmid = imax(pxcorns[0], pxcorns[1*2]);
	int endmid = imin(pxcorns[2*2], pxcorns[3*2]);
	for (int x=startmid+1; x<endmid; x++) {
		ybound[x*2+0] = (int) round(polyget(&tline, (double) x));
		ybound[x*2+1] = (int) round(polyget(&bline, (double) x));
	}

	printf("ybound at 0: %d, %d\n", *ybound, *(ybound+1));

	// Determine which part of pattern each pixel receives
	// ***************************************************

	// Set up array of vectors to hold pattern positions. Later on it holds
	// normal vectors.
	double *vecs = malloc(xlen*ylen*3*sizeof(*vecs));
	if (vecs == NULL) {
		printf("Cannot allocate memory for vecs, exiting...\n");
		exit(1);
	}
	// Initialising array to -1 which is used as a check in calcpattvecs.
	for (int i=0; i<(xlen*ylen*3); i++)
		vecs[i] = -1.0;
	 
	// Start looking in first segment.
	calcpattvecs(&pattern, imdata, wid, vecs, xlen, pxcorns, ybound, offset,
			0.5*pattern.segsize, 1);
	printf("Information extracted from horizontal image.\n");

	// Loading vertical image.
	image = TIFFOpen(imfnv, "r");
	if (image == NULL) {
		printf("Failed to open vertical image.\n");
		exit(1);
	}

	// Load in image data.
	if (TIFFReadRGBAImageOriented(image, wid, hei, imdata, ORIENTATION_TOPLEFT,
			0) == 0) {
		printf("Failed to read vertical image.\n");
		exit(1);
	}
	printf("Vertical image loaded\n");

	// Start looking in first segment.
	calcpattvecs(&pattern, imdata, wid, vecs, xlen, pxcorns, ybound, offset,
			0.5*pattern.segsize, 0);
	printf("Information extracted from vertical image.\n");

	_TIFFfree(imdata);
	imdata = NULL;

	TIFFClose(image);
	image = NULL;

	// Piece together pattern vector components and transform into global
	// coordinates.
	transformpatt(&pattern, vecs, xlen, ybound);
	printf("Pattern vectors translated to global coords.\n");

	// Solve for mirror
	// ****************

	// Set up array of vectors to hold mirror surface positions.
	double *poss = malloc(xlen*ylen*3*sizeof(*poss));
	if (poss == NULL) {
		printf("Cannot allocate memory for poss, exiting...\n");
		exit(1);
	}

	// Initialising array to -1 which is used as a check in solvesurface.
	for (int i=0; i<(xlen*ylen*3); i++)
		poss[i] = -1.0;

	// Testing different solution starts.
	/*
	int initpixi[2];
	//initpixi[0] = (*(pxcorns+6) - *(pxcorns))/2 + *(pxcorns);
	initpixi[0] = imax(*(pxcorns+2), *(pxcorns));
	initpixi[1] = (*(pxcorns+3) - *(pxcorns+1))/2 + *(pxcorns+1);
	printf("Initial pixel: %d, %d\n", initpixi[0], initpixi[1]);
	*/

	solvesurface(&camera, vecs, poss, xlen, pxcorns, ybound, offset,
			mirror.startdepth);
	printf("Surface solved.\n");

	// Centring data points about axis.
	double xmean = 0.0;
	double ymean = 0.0;
	double zmax = -100.0;
	int meancount = 0;
	for (int px=0; px<xlen; px++) {
		for (int py=ybound[px*2+0]; py<=ybound[px*2+1]; py++) {
			xmean = xmean + poss[(py*xlen+px)*3+0];
			ymean = ymean + poss[(py*xlen+px)*3+1];
			if (zmax < poss[(py*xlen+px)*3+2])
				zmax = poss[(py*xlen+px)*3+2];
			meancount++;
		}
	}
	xmean = xmean/meancount;
	ymean = ymean/meancount;
	printf("xmean, ymean, zmax: %f, %f, %f\n", xmean, ymean, zmax);
	for (int px=0; px<xlen; px++) {
		for (int py=ybound[px*2+0]; py<=ybound[px*2+1]; py++) {
			poss[(py*xlen+px)*3+0] = poss[(py*xlen+px)*3+0] - xmean;
			poss[(py*xlen+px)*3+1] = poss[(py*xlen+px)*3+1] - ymean;
			poss[(py*xlen+px)*3+2] = poss[(py*xlen+px)*3+2] - zmax;
		}
	}

	// Find best fitting curves
	// ************************

	// Finding best fitting sphere.
	Errparams fitpars = {poss, xlen, ybound, NULL};
	gsl_vector *sphfitvars = gsl_vector_alloc(4);
	// Setting starting point.
	gsl_vector_set(sphfitvars, 0, 0.0); // x0
	gsl_vector_set(sphfitvars, 1, 0.0); // y0
	gsl_vector_set(sphfitvars, 2, 0.0); // z0 shift (radius accounted)
	gsl_vector_set(sphfitvars, 3, 30000.0); // Radius

	// Step size for first trial.
	gsl_vector *sphstep = gsl_vector_alloc(4);
	gsl_vector_set(sphstep, 0, 500.0);
	gsl_vector_set(sphstep, 1, 500.0);
	gsl_vector_set(sphstep, 2, 10.0);
	gsl_vector_set(sphstep, 3, 5000.0);

	// Fitting sphere to data.
	minerror(&sphesqerr, &fitpars, 4, sphfitvars, sphstep);

	gsl_vector_free(sphstep);
	sphstep = NULL;

	// Saving parameters into array.
	double spvars[4];
	for (int i=0; i<4; i++)
		spvars[i] = gsl_vector_get(sphfitvars, i);

	gsl_vector_free(sphfitvars);
	sphfitvars = NULL;

	printf("Best fit sphere:\n");
	printf("%f, %f, %f, %f\n", spvars[0], spvars[1], spvars[2], spvars[3]);

	// Find best fitting paraboloid.
	gsl_vector *parfitvars = gsl_vector_alloc(8);
	gsl_vector_set(parfitvars, 0, 0.0); // x0
	gsl_vector_set(parfitvars, 1, 0.0); // y0
	gsl_vector_set(parfitvars, 2, 0.0); // z0 shift
	gsl_vector_set(parfitvars, 3, 0.0); // xrot
	gsl_vector_set(parfitvars, 4, 0.0); // yrot
	gsl_vector_set(parfitvars, 5, 0.0); // zrot
	gsl_vector_set(parfitvars, 6, spvars[3]/2.0); // f1
	gsl_vector_set(parfitvars, 7, spvars[3]/2.0); // f2

	// Step size for first trial.
	gsl_vector *parstep = gsl_vector_alloc(8);
	gsl_vector_set(parstep, 0, 500.0);
	gsl_vector_set(parstep, 1, 500.0);
	gsl_vector_set(parstep, 2, 20.0);
	gsl_vector_set(parstep, 3, 5.0*PI/180.0);
	gsl_vector_set(parstep, 4, 5.0*PI/180.0);
	gsl_vector_set(parstep, 5, 90.0*PI/180.0);
	gsl_vector_set(parstep, 6, 5000.0);
	gsl_vector_set(parstep, 7, 5000.0);

	// Fitting paraboloid to data.
	minerror(&parabsqerr, &fitpars, 8, parfitvars, parstep);

	gsl_vector_free(parstep);
	parstep = NULL;

	// Saving parameters into array.
	double pavars[8];
	for (int i=0; i<8; i++)
		pavars[i] = gsl_vector_get(parfitvars, i);

	gsl_vector_free(parfitvars);
	parfitvars = NULL;

	printf("Best fit paraboloid:\n");
	printf("%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", pavars[0], pavars[1],
			pavars[2], pavars[3], pavars[4], pavars[5], pavars[6], pavars[7]);
	
	// Calculating slope errors relative to surfaces
	// *********************************************

	// Setting up arrays for slope error data.
	double *spserr = malloc(xlen*ylen*2*sizeof(*spserr)); // Sphere s/e
	if (spserr == NULL) {
		printf("Cannot allocate memory for spherical s/e, exiting...\n");
		exit(1);
	}
	double *paserr = malloc(xlen*ylen*2*sizeof(*paserr)); // Paraboloid s/e
	if (paserr == NULL) {
		printf("Cannot allocate memory for paraboloidal s/e, exiting...\n");
		exit(1);
	}

	double spsest[4]; // Sphere slope error statistics
	double pasest[4]; // Paraboloid slope error statistics
	// Slope errors for sphere and paraboloid.
	slopeerror(vecs, poss, xlen, ybound, 0, spvars, spserr, spsest);
	printf("Sphere slope errors: %f, %f, %f, %f\n", spsest[0], spsest[1],
			spsest[2], spsest[3]); 
	slopeerror(vecs, poss, xlen, ybound, 1, pavars, paserr, pasest);
	printf("Paraboloid slope errors: %f, %f, %f, %f\n", pasest[0], pasest[1],
			pasest[2], pasest[3]); 

	// Saving panel results to results file.
	FILE *resultfile = fopen("../images/results.out", "a");
	if (resultfile== NULL) {
		printf("Error opening save file\n");	
		exit(1);
	}
	// fn, sphere RoC, parab f1, parab f2, parab zrot, sphere sigx,
	// sphere sigy, parab sigx, parab sigy.
	fprintf(resultfile, "%s,%10.4e,%10.4e,%10.4e,%10.4e,%10.4e,%10.4e,%10.4e,%10.4e\n",
			outfn, spvars[3], pavars[6], pavars[7], pavars[5], spsest[2],
			spsest[3], pasest[2], pasest[3]);
	fclose(resultfile);
	resultfile = NULL;

	// Saving spherical slope errors.
	char sseext[] = ".sse";
	char *ssefn = malloc((strlen(sseext)+strlen(outfn)+1)*sizeof(*ssefn));
	strcpy(ssefn, outfn);
	strcat(ssefn, sseext);
	FILE *ssavefile = fopen(ssefn, "w");
	if (ssavefile == NULL) {
		printf("Error opening save file\n");	
		exit(1);
	}
	for (int x=0; x<xlen; x=x+xlen/40) {
		for (int y=ybound[x*2+0]; y<=ybound[x*2+1]; y=y+xlen/40) {
			fprintf(ssavefile, "%10.3e,%10.3e\n", spserr[(y*xlen+x)*2+0],
				spserr[(y*xlen+x)*2+1]);
		}
	}
	fclose(ssavefile);
 
	free(ssefn);
	ssefn = NULL;

	// Saving paraboloid slope errors.
	char pseext[] = ".pse";
	char *psefn = malloc((strlen(pseext)+strlen(outfn)+1)*sizeof(*psefn));
	strcpy(psefn, outfn);
	strcat(psefn, pseext);
	ssavefile = fopen(psefn, "w");
	if (ssavefile == NULL) {
		printf("Error opening save file\n");	
		exit(1);
	}
	for (int x=0; x<xlen; x=x+xlen/40) {
		for (int y=ybound[x*2+0]; y<=ybound[x*2+1]; y=y+xlen/40) {
			fprintf(ssavefile, "%10.3e,%10.3e\n", paserr[(y*xlen+x)*2+0],
				paserr[(y*xlen+x)*2+1]);
		}
	}
	fclose(ssavefile);
	ssavefile = NULL;
 
	free(psefn);
	psefn = NULL;

	free(spserr);
	spserr = NULL;
	free(paserr);
	paserr = NULL;

	// Saving to file selection of points.
	char posext[] = ".pos";
	char *posfn = malloc((strlen(posext)+strlen(outfn)+1)*sizeof(*posfn));
	strcpy(posfn, outfn);
	strcat(posfn, posext);
	FILE *savefile = fopen(posfn, "w");
	if (savefile == NULL) {
		printf("Error opening save file\n");	
		exit(1);
	}
	for (int x=0; x<xlen; x=x+xlen/100) {
		for (int y=ybound[x*2+0]; y<=ybound[x*2+1]; y=y+xlen/100) {
			fprintf(savefile, "%10.3e,%10.3e,%10.3e\n", poss[(y*xlen+x)*3+0],
				poss[(y*xlen+x)*3+1], poss[(y*xlen+x)*3+2]);
		}
	}
	fclose(savefile);
	free(posfn);
	posfn = NULL;

	// Saving to file a lessor selection of points.
	char poslext[] = ".posl";
	char *poslfn = malloc((strlen(poslext)+strlen(outfn)+1)*sizeof(*poslfn));
	strcpy(poslfn, outfn);
	strcat(poslfn, poslext);
	savefile = fopen(poslfn, "w");
	if (savefile == NULL) {
		printf("Error opening save file\n");	
		exit(1);
	}
	for (int x=0; x<xlen; x=x+xlen/40) {
		for (int y=ybound[x*2+0]; y<=ybound[x*2+1]; y=y+xlen/40) {
			fprintf(savefile, "%10.3e,%10.3e,%10.3e\n", poss[(y*xlen+x)*3+0],
				poss[(y*xlen+x)*3+1], poss[(y*xlen+x)*3+2]);
		}
	}
	fclose(savefile);
	savefile = NULL;
	free(poslfn);
	poslfn = NULL;
	
	free(vecs);
	vecs = NULL;
	free(poss);
	poss = NULL;

	freecamera(&camera);
	freepattern(&pattern);
}

// Determines which part of pattern is reflected for each pixel.
void calcpattvecs(const Pattern *pat, const uint32_t *im, const int iw,
		double *vecs, const int bw, const int *ipix, const int *yb,
		const int *offs, const double idist, const int orien)
{
	// Orientation 1 is horizontal, 0 is vertical.

	// Working out first point.
	vecs[(ipix[1]*bw+ipix[0])*3+orien] = getdist(pat,
			im+(ipix[1]+offs[1])*iw+ipix[0]+offs[0], idist, orien);
	// Run down first column, shouldn't need to go up since at corner.
	for (int y=ipix[1]+1; y<=yb[ipix[0]*2+1]; y++) {
		// Using previously calculated distance.
		double pdist = vecs[((y-1)*bw+ipix[0])*3+orien];
		vecs[(y*bw+ipix[0])*3+orien] = getdist(pat,
				im+(y+offs[1])*iw+ipix[0]+ offs[0], pdist, orien);
	}

	// Work our way left.
	for (int x=ipix[0]-1; x>=0; x--) {
		// Find a good starting point.
		int iy = imin(yb[x*2+1], imax(ipix[1], yb[x*2]));
		// Determining distance for starting point.
		double pdist = vecs[(iy*bw+x+1)*3+orien];
		vecs[(iy*bw+x)*3+orien] = getdist(pat,
				im+(iy+offs[1])*iw+x+offs[0], pdist, orien);
		// Shouldn't need to move up.

		// Working down.
		for (int y=iy+1; y<=yb[x*2+1]; y++) {
			int count = 1;
			// Averaging from up to 3 surrounding distances.
			pdist = vecs[((y-1)*bw+x)*3+orien];
			double dcorn = vecs[((y-1)*bw+x+1)*3+orien];
			double dnext = vecs[(y*bw+x+1)*3+orien];
			if (dcorn != -1.0) {
				pdist = pdist + dcorn;
				count++;
			}
			if (dnext != -1.0) {
				pdist = pdist + dnext;
				count++;
			}
			pdist = pdist/count; // Normalizing

			vecs[(y*bw+x)*3+orien] = getdist(pat,
					im+(y+offs[1])*iw+x+offs[0], pdist, orien);
		}
	}

	// Work our way right.
	for (int x=ipix[0]+1; x<bw; x++) {
		// Find a good starting point.
		int iy = imin(yb[x*2+1], imax(ipix[1], yb[x*2]));
		// Determining distance for starting point.
		double pdist = vecs[(iy*bw+x-1)*3+orien];
		vecs[(iy*bw+x)*3+orien] = getdist(pat,
				im+(iy+offs[1])*iw+x+offs[0], pdist, orien);
		// Working up.
		for (int y=iy-1; y>=yb[x*2]; y--) {
			int count = 1;
			// Averaging from up to 3 surrounding distances.
			pdist = vecs[((y+1)*bw+x)*3+orien];
			double dcorn = vecs[((y+1)*bw+x-1)*3+orien];
			double dnext = vecs[(y*bw+x-1)*3+orien];
			if (dcorn != -1.0) {
				pdist = pdist + dcorn;
				count++;
			}
			if (dnext != -1.0) {
				pdist = pdist + dnext;
				count++;
			}
			pdist = pdist/count; // Normalizing

			vecs[(y*bw+x)*3+orien] = getdist(pat,
					im+(y+offs[1])*iw+x+offs[0], pdist, orien);
		}
		// Working down.
		for (int y=iy+1; y<=yb[x*2+1]; y++) {
			int count = 1;
			// Averaging from up to 3 surrounding distances.
			pdist = vecs[((y-1)*bw+x)*3+orien];
			double dcorn = vecs[((y-1)*bw+x-1)*3+orien];
			double dnext = vecs[(y*bw+x-1)*3+orien];
			if (dcorn != -1.0) {
				pdist = pdist + dcorn;
				count++;
			}
			if (dnext != -1.0) {
				pdist = pdist + dnext;
				count++;
			}
			pdist = pdist/count; // Normalizing

			vecs[(y*bw+x)*3+orien] = getdist(pat,
					(im+(y+offs[1])*iw+x+offs[0]), pdist, orien);
		}
	}
}

void transformpatt(const Pattern *pat, double *vecs, const int bw, const int* yb)
{
	for (int x=0; x<bw; x++) {
		for (int y=*(yb+x*2); y<=*(yb+x*2+1); y++) {
			transpattvec(pat, (vecs+(y*bw+x)*3));
		}
	}
}

void solvesurface(const Camera *cam, double *vecs, double *poss, const int bw,
		const int *ipix, const int *yb, const int *offs, const double idepth)
{
	// After position and norm estimated, could redo extrapolation taking,
	// into account the estimated norm.
	double ch[3]; // Direction from camera to point
	double ph[3]; // Direction from pattern to point
	// Calculating direction from camera to point.
	pixdir(cam, *(ipix)+*(offs), *(ipix+1)+*(offs+1), ch);
	// Estimating position of point.
	*(poss+((*(ipix+1))*bw+*ipix)*3) = ch[0]*(idepth-cam->pos[2])/ch[2] +
		cam->pos[0];
	*(poss+((*(ipix+1))*bw+*ipix)*3+1) = ch[1]*(idepth-cam->pos[2])/ch[2] +
		cam->pos[1];
	*(poss+((*(ipix+1))*bw+*ipix)*3+2) = idepth;

	ph[0] = *(vecs+((*(ipix+1))*bw+*ipix)*3) -
		*(poss+((*(ipix+1))*bw+*ipix)*3);
	ph[1] = *(vecs+((*(ipix+1))*bw+*ipix)*3+1) -
		*(poss+((*(ipix+1))*bw+*ipix)*3+1);
	ph[2] = *(vecs+((*(ipix+1))*bw+*ipix)*3+2) -
		*(poss+((*(ipix+1))*bw+*ipix)*3+2);
	scale(1.0/norm(ph), ph);

	surfnorm(ch, ph, (vecs+((*(ipix+1))*bw+*ipix)*3));
	
	// Working up starting column.
	for (int y=*(ipix+1)-1; y>=*(yb+(*ipix)*2); y--) {
		pixdir(cam, *(ipix)+*(offs), y+*(offs+1), ch);
		double dist = extrapolate(cam, (poss+((y+1)*bw+*ipix)*3),
				(vecs+((y+1)*bw+*ipix)*3), ch);
		*(poss+(y*bw+*ipix)*3) = dist*ch[0]+ cam->pos[0];
		*(poss+(y*bw+*ipix)*3+1) = dist*ch[1]+ cam->pos[1];
		*(poss+(y*bw+*ipix)*3+2) = dist*ch[2]+ cam->pos[2];
		ph[0] = *(vecs+(y*bw+*ipix)*3) -
			*(poss+(y*bw+*ipix)*3);
		ph[1] = *(vecs+(y*bw+*ipix)*3+1) -
			*(poss+(y*bw+*ipix)*3+1);
		ph[2] = *(vecs+(y*bw+*ipix)*3+2) -
			*(poss+(y*bw+*ipix)*3+2);
		scale(1.0/norm(ph), ph);
		surfnorm(ch, ph, (vecs+(y*bw+*ipix)*3));
	}

	// Working way down starting column.
	for (int y=*(ipix+1)+1; y<=*(yb+(*ipix)*2+1); y++) {
		pixdir(cam, *(ipix)+*(offs), y+*(offs+1), ch);
		double dist = extrapolate(cam, (poss+((y-1)*bw+*ipix)*3),
				(vecs+((y-1)*bw+*ipix)*3), ch);
		*(poss+(y*bw+*ipix)*3) = dist*ch[0]+ cam->pos[0];
		*(poss+(y*bw+*ipix)*3+1) = dist*ch[1]+ cam->pos[1];
		*(poss+(y*bw+*ipix)*3+2) = dist*ch[2]+ cam->pos[2];
		ph[0] = *(vecs+(y*bw+*ipix)*3) -
			*(poss+(y*bw+*ipix)*3);
		ph[1] = *(vecs+(y*bw+*ipix)*3+1) -
			*(poss+(y*bw+*ipix)*3+1);
		ph[2] = *(vecs+(y*bw+*ipix)*3+2) -
			*(poss+(y*bw+*ipix)*3+2);
		scale(1.0/norm(ph), ph);
		surfnorm(ch, ph, (vecs+(y*bw+*ipix)*3));
	}

	// Work our way left.
	for (int x=*ipix-1; x>=0; x--) {
		// Find a good starting point.
		int iy = imin(*(yb+x*2+1),imax(*(ipix+1),*(yb+x*2)));
		pixdir(cam, x+*(offs), iy+*(offs+1), ch);
		double dist = extrapolate(cam, (poss+(iy*bw+x+1)*3),
				(vecs+(iy*bw+x+1)*3), ch);
		*(poss+(iy*bw+x)*3) = dist*ch[0]+ cam->pos[0];
		*(poss+(iy*bw+x)*3+1) = dist*ch[1]+ cam->pos[1];
		*(poss+(iy*bw+x)*3+2) = dist*ch[2]+ cam->pos[2];
		ph[0] = *(vecs+(iy*bw+x)*3) -
			*(poss+(iy*bw+x)*3);
		ph[1] = *(vecs+(iy*bw+x)*3+1) -
			*(poss+(iy*bw+x)*3+1);
		ph[2] = *(vecs+(iy*bw+x)*3+2) -
			*(poss+(iy*bw+x)*3+2);
		scale(1.0/norm(ph), ph);
		surfnorm(ch, ph, (vecs+(iy*bw+x)*3));
		
		// Working up.
		for (int y=iy-1; y>=*(yb+x*2); y--) {
			int count = 1;
			pixdir(cam, x+*(offs), y+*(offs+1), ch);
			dist = extrapolate(cam, (poss+((y+1)*bw+x)*3),
					(vecs+((y+1)*bw+x)*3), ch);
			if (*(poss+((y+1)*bw+x+1)*3) != -1.0) {
				dist = dist + extrapolate(cam, (poss+((y+1)*bw+x+1)*3),
							(vecs+((y+1)*bw+x+1)*3), ch);
				count++;
			}
			if (*(poss+(y*bw+x+1)*3) != -1.0) {
				dist = dist + extrapolate(cam, (poss+(y*bw+x+1)*3),
							(vecs+(y*bw+x+1)*3), ch);
				count++;
			}
			dist = dist/count;
			*(poss+(y*bw+x)*3) = dist*ch[0]+ cam->pos[0];
			*(poss+(y*bw+x)*3+1) = dist*ch[1]+ cam->pos[1];
			*(poss+(y*bw+x)*3+2) = dist*ch[2]+ cam->pos[2];
			ph[0] = *(vecs+(y*bw+x)*3) -
				*(poss+(y*bw+x)*3);
			ph[1] = *(vecs+(y*bw+x)*3+1) -
				*(poss+(y*bw+x)*3+1);
			ph[2] = *(vecs+(y*bw+x)*3+2) -
				*(poss+(y*bw+x)*3+2);
			scale(1.0/norm(ph), ph);
			surfnorm(ch, ph, (vecs+(y*bw+x)*3));
		}

		// Working down.
		for (int y=iy+1; y<=*(yb+x*2+1); y++) {
			int count = 1;
			pixdir(cam, x+*(offs), y+*(offs+1), ch);
			dist = extrapolate(cam, (poss+((y-1)*bw+x)*3),
					(vecs+((y-1)*bw+x)*3), ch);
			if (*(poss+((y-1)*bw+x+1)*3) != -1.0) {
				dist = dist + extrapolate(cam, (poss+((y-1)*bw+x+1)*3),
							(vecs+((y-1)*bw+x+1)*3), ch);
				count++;
			}
			if (*(poss+(y*bw+x+1)*3) != -1.0) {
				dist = dist + extrapolate(cam, (poss+(y*bw+x+1)*3),
							(vecs+(y*bw+x+1)*3), ch);
				count++;
			}
			dist = dist/count;
			*(poss+(y*bw+x)*3) = dist*ch[0]+ cam->pos[0];
			*(poss+(y*bw+x)*3+1) = dist*ch[1]+ cam->pos[1];
			*(poss+(y*bw+x)*3+2) = dist*ch[2]+ cam->pos[2];
			ph[0] = *(vecs+(y*bw+x)*3) -
				*(poss+(y*bw+x)*3);
			ph[1] = *(vecs+(y*bw+x)*3+1) -
				*(poss+(y*bw+x)*3+1);
			ph[2] = *(vecs+(y*bw+x)*3+2) -
				*(poss+(y*bw+x)*3+2);
			scale(1.0/norm(ph), ph);
			surfnorm(ch, ph, (vecs+(y*bw+x)*3));
		}
	}

	// Work our way right.
	for (int x=*ipix+1; x<bw; x++) {
		// Find a good starting point.
		int iy = imin(*(yb+x*2+1),imax(*(ipix+1),*(yb+x*2)));
		pixdir(cam, x+*(offs), iy+*(offs+1), ch);
		double dist = extrapolate(cam, (poss+(iy*bw+x-1)*3),
				(vecs+(iy*bw+x-1)*3), ch);
		*(poss+(iy*bw+x)*3) = dist*ch[0]+ cam->pos[0];
		*(poss+(iy*bw+x)*3+1) = dist*ch[1]+ cam->pos[1];
		*(poss+(iy*bw+x)*3+2) = dist*ch[2]+ cam->pos[2];
		ph[0] = *(vecs+(iy*bw+x)*3) -
			*(poss+(iy*bw+x)*3);
		ph[1] = *(vecs+(iy*bw+x)*3+1) -
			*(poss+(iy*bw+x)*3+1);
		ph[2] = *(vecs+(iy*bw+x)*3+2) -
			*(poss+(iy*bw+x)*3+2);
		scale(1.0/norm(ph), ph);
		surfnorm(ch, ph, (vecs+(iy*bw+x)*3));

		// Working up.
		for (int y=iy-1; y>=*(yb+x*2); y--) {
			int count = 1;
			pixdir(cam, x+*(offs), y+*(offs+1), ch);
			dist = extrapolate(cam, (poss+((y+1)*bw+x)*3),
					(vecs+((y+1)*bw+x)*3), ch);
			if (*(poss+((y+1)*bw+x-1)*3) != -1.0) {
				dist = dist + extrapolate(cam, (poss+((y+1)*bw+x-1)*3),
							(vecs+((y+1)*bw+x-1)*3), ch);
				count++;
			}
			if (*(poss+(y*bw+x-1)*3) != -1.0) {
				dist = dist + extrapolate(cam, (poss+(y*bw+x-1)*3),
							(vecs+(y*bw+x-1)*3), ch);
				count++;
			}
			dist = dist/count;
			*(poss+(y*bw+x)*3) = dist*ch[0]+ cam->pos[0];
			*(poss+(y*bw+x)*3+1) = dist*ch[1]+ cam->pos[1];
			*(poss+(y*bw+x)*3+2) = dist*ch[2]+ cam->pos[2];
			ph[0] = *(vecs+(y*bw+x)*3) -
				*(poss+(y*bw+x)*3);
			ph[1] = *(vecs+(y*bw+x)*3+1) -
				*(poss+(y*bw+x)*3+1);
			ph[2] = *(vecs+(y*bw+x)*3+2) -
				*(poss+(y*bw+x)*3+2);
			scale(1.0/norm(ph), ph);
			surfnorm(ch, ph, (vecs+(y*bw+x)*3));
		}
		// Working down.
		for (int y=iy+1; y<=*(yb+x*2+1); y++) {
			int count = 1;
			pixdir(cam, x+*(offs), y+*(offs+1), ch);
			dist = extrapolate(cam, (poss+((y-1)*bw+x)*3),
					(vecs+((y-1)*bw+x)*3), ch);
			if (*(poss+((y-1)*bw+x-1)*3) != -1.0) {
				dist = dist + extrapolate(cam, (poss+((y-1)*bw+x-1)*3),
							(vecs+((y-1)*bw+x-1)*3), ch);
				count++;
			}
			if (*(poss+(y*bw+x-1)*3) != -1.0) {
				dist = dist + extrapolate(cam, (poss+(y*bw+x-1)*3),
							(vecs+(y*bw+x-1)*3), ch);
				count++;
			}
			dist = dist/count;
			*(poss+(y*bw+x)*3) = dist*ch[0]+ cam->pos[0];
			*(poss+(y*bw+x)*3+1) = dist*ch[1]+ cam->pos[1];
			*(poss+(y*bw+x)*3+2) = dist*ch[2]+ cam->pos[2];
			ph[0] = *(vecs+(y*bw+x)*3) -
				*(poss+(y*bw+x)*3);
			ph[1] = *(vecs+(y*bw+x)*3+1) -
				*(poss+(y*bw+x)*3+1);
			ph[2] = *(vecs+(y*bw+x)*3+2) -
				*(poss+(y*bw+x)*3+2);
			scale(1.0/norm(ph), ph);
			surfnorm(ch, ph, (vecs+(y*bw+x)*3));
		}
	}
}

void surfnorm(const double *vin, const double *vre, double *snrm)
{
	// vin - incident ray, vre - reflected ray.
	snrm[0] = vre[0] - vin[0];
	snrm[1] = vre[1] - vin[1];
	snrm[2] = vre[2] - vin[2];
	scale(1.0/norm(snrm), snrm);
}

double extrapolate(const Camera *cam, const double *ppos, const double *psnrm,
		const double *ch)
{
	double pc[3] = {ppos[0]-cam->pos[0], ppos[1]-cam->pos[1],
		ppos[2]-cam->pos[2]};
	return dot(pc, psnrm)/dot(ch, psnrm);
}

void slopeerror(const double *vecs, const double *poss, const int bw,
		const int *yb, const int shape, const double *params, double *serr,
		double *sestats)
{
	// Calculating slope errors relative to sphere.
	// Shape 0 is sphere, shape 1 is elliptic paraboloid.
	// Parameters are in terms of acting on raw data.
	// For sphere: params = xshift yshift zshift rad.
	// For paraboloid: params = xshift yshift zshift xrot yrot zrot f1 f2.
	size_t pointc = 0;
	double inorm[3]; // Ideal normal
	double mux = 0.0; // Slope error mean
	double muy = 0.0; // Slope error mean
	double locx[3];
	double locy[3];
	double xaxis[3] = {1.0, 0.0, 0.0};
	double tpos1[3];
	double tpos2[3];
	double tnrm[3];
	double rotm[9];
	double inrotm[9];
	if (shape == 1) {
		rotmatxyz(params[3], params[4], params[5], rotm);
		rotmatzyx(-params[5], -params[4], -params[3], inrotm);
	}
	for (int px=0; px<bw; px++) {
		for (int py=yb[px*2]; py<=yb[px*2+1]; py++) {
			// Work out ideal normal for x and y coords of given pixel.
			if (shape == 0) {
				sphereslope(poss[(py*bw+px)*3] + params[0],
						poss[(py*bw+px)*3+1] + params[1], params[3], inorm);
			} else if (shape == 1) {
				// Data shifted before being rotated.
				tpos1[0] = poss[(py*bw+px)*3+0] + params[0];
				tpos1[1] = poss[(py*bw+px)*3+1] + params[1];
				tpos1[2] = poss[(py*bw+px)*3+2] + params[2];
				// Data rotated. z in tpos2 not used past here.
				matxvec(rotm, tpos1, tpos2);
				parabslope(tpos2[0], tpos2[1], params[6], params[7], tnrm);
				// Rotate inorm so that it matches up with data normals.
				matxvec(inrotm, tnrm, inorm);
			}

			// Take cross product of inorm with global x axis to produce locy.
			// Note that y axis ends up pointing in roughly opposite direction
			// to global y.
			cross(inorm, xaxis, locy);
			cross(locy, inorm, locx);
			
			// Working out local x and y components of surface norm.
			double epsx = dot(locx, vecs+(py*bw+px)*3);
			double epsy = dot(locy, vecs+(py*bw+px)*3);

			// For now just save in vecs.
			serr[(py*bw+px)*2+0] = epsx;
			serr[(py*bw+px)*2+1] = epsy;

			mux = mux + epsx;
			muy = muy + epsy;
			pointc++;
		}
	}

	mux = mux/pointc;
	muy = muy/pointc;
	double sigx = 0.0;
	double sigy = 0.0;
	for (int px=0; px<bw; px++) {
		for (int py=yb[px*2]; py<=yb[px*2+1]; py++) {
			sigx = sigx + pow(serr[(py*bw+px)*2+0]-mux, 2.0);
			sigy = sigy + pow(serr[(py*bw+px)*2+1]-muy, 2.0);
		}
	}
	sigx = sqrt(sigx/(pointc-1));
	sigy = sqrt(sigy/(pointc-1));
	sestats[0] = mux;
	sestats[1] = muy;
	sestats[2] = sigx;
	sestats[3] = sigy;
}

void centroid(const uint32_t *im, uint32_t w, uint32_t h, const int *idot,
		double *dot, int dotpx)
{

	double *region = NULL; // Region of pixels about dot guess
	int regw, regh; // Region width and height

	// Finding the region boundaries which are 2 times dot size either side.
	// Define x along width and y along height.
	int lowx = idot[0] - 2*dotpx;
	int highx = idot[0] + 2*dotpx;
	int lowy = idot[1] - 2*dotpx;
	int highy = idot[1] + 2*dotpx;
	double backgr = 0.0;
	double weight = 0.0;
	double xmean = 0.0;
	double ymean = 0.0;
	
	// Check if region bounds are outside image bounds and correct.
	if (lowx < 0) lowx = 0;
	if (highx >= w) highx = w - 1;
	if (lowy < 0) lowy = 0;
	if (highy >= h) highy = h - 1;
	// Region includes end points.
	regw = highx-lowx+1;
	regh = highy-lowy+1;

 	region = malloc(regw*regh*sizeof(*region));
	if (region == NULL) {
		printf("Failled to allocate memory for centroiding, exiting...\n");
		exit(1);
	}

	// Each pixel in the region is set to its greyscale value.
	for (int rx=0; rx<regw; rx++) {
		for (int ry=0; ry<regh; ry++) {
			double greypix;
			uint32_t pix = im[(ry+lowy)*w+rx+lowx];
			greypix = TIFFGetR(pix) + TIFFGetG(pix) + TIFFGetB(pix);
			greypix = 765.0 - greypix;
			region[ry*regw+rx] = greypix;
		}
	}

	// Find the background intensity by looking at boundaries.
	for (int rx=0; rx<regw; rx++) {
		backgr = backgr + region[rx] + region[(regh-1)*regw+rx];
	}
	for (int ry=0; ry<regh; ry++) {
		backgr = backgr + region[ry*regw] + region[ry*regw+regh-1];
	}
	backgr = backgr/(2.0*regw + 2.0*regh);

	// Subtract background intensity and find weighting.
	for (int rx=0; rx<regw; rx++) {
		for (int ry=0; ry<regh; ry++) {
			double greypix = region[ry*regw+rx] - backgr;
			if (greypix < 0)
				greypix = 0.0;
			region[ry*regw+rx] = greypix;
			weight = weight + greypix;
		}
	}

	// Find x and y weighted means.
	for (int rx=0; rx<regw; rx++) {
		for (int ry=0; ry<regh; ry++) {
			xmean = xmean + rx*region[ry*regw+rx]/weight;
			ymean = ymean + ry*region[ry*regw+rx]/weight;
		}
	}

	free(region);
	region = NULL;

	dot[0] = xmean + lowx;
	dot[1] = ymean + lowy;
}

// Find line that intercepts two points
void pnts_line(const int *p1, const int *p2, Polynom *pol, double *coeffs)
{
	// First checking to make sure that we don't get infinite slope.
	// If so then polynomial set to constant zero (it shouldn't get called).
	if (p1[0] == p2[0]) {
		pol->degree = 0;
		coeffs[0] = 0.0;
	} else {
		coeffs[1] = ((double) (p1[1] - p2[1]))/(p1[0] - p2[0]);
		coeffs[0] = (double) (p2[1] - coeffs[1]*p2[0]);
	}
}

/*
 * Find upper and lower y bounds for ends of region of interest.
 */
void end_ybounds(const Polynom *l1, const Polynom *l2, const Polynom *l3,
		const Polynom *l4, const int *cornt, const int *cornb, int *ybnd)
{
	// If corners have same x value then we have a vertical line.
	if (cornt[0] == cornb[0]) {
		ybnd[cornt[0]*2+0] = cornt[1];
		ybnd[cornt[0]*2+1] = cornb[2];
	} else {
		// At most one of the two loops will execute.
		for (int x=cornt[0]; x<=cornb[0]; x++) {
			ybnd[x*2+0] = (int) round(polyget(l1, (double) x));
			ybnd[x*2+1] = (int) round(polyget(l2, (double) x));
		}
		for (int x=cornb[0]; x<=cornt[0]; x++) {
			ybnd[x*2+0] = (int) round(polyget(l3, (double) x));
			ybnd[x*2+1] = (int) round(polyget(l4, (double) x));
		}
	}
}

