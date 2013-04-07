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

/* Mirror parameters
 * *****************
 */
static double dotsize = 7.0; // Physiscal size of reference of dots in mm
static double dotsep[3] = {1192.8986,1745.6891,1271.5304}; // Distance
// between TL-BL, BL-TR, TR-TL
static double corns[4][3] = {{61.010941,17.819132,-24.028367},
	{81.010941,1167.8191,-24.028367},
	{1211.0109,1167.8191,-24.028367},
	{1211.0109,17.819132,-24.028367}}; // Mirror boundary corners:
// TL, BL, BR, TR
static double camdistguess = 3000.0; // Guess of distance from mirror to camera
static double startdepth = -24.028367; // Height of starting point.
static char *relfn = "data/huerel.csv";

/* LOCAL FUNCTIONS
 * ***************
 * Declared here instead of header.
 */
void solveprofile(char *imfnh, char *imfnv, int *idots, char *outfn);
void centroid(uint32 *im, uint32 w, uint32 h, int *idot, double *dot,
		int dotpx);
void calcpattvecs(const uint32 *im, const int iw, double *vecs, const int bw,
		const int *ipix, const int *yb, const int *offs, const double idist,
		const int orien);
void transformpatt(double *vecs, const int bw, const int* yb);
double extrapolate(const double *ppos, const double *psnrm, const double *ch);
void solvesurface(double *vecs, double *poss, const int bw, const int *ipix,
		const int *yb, const int *offs, const double idepth);
void surfnorm(const double *vin, const double *vre, double *snrm);
void slopeerror(const double *vecs, const double *poss, const int bw,
		const int *yb, const int shape, const double *params, double *serr,
		double *sestats);

void solveprofile(char *imfnh, char *imfnv, int *idots, char *outfn) {
	TIFF *image; // TIFF file pointer
	uint32 wid, hei; // Image width and height
	uint32 *imdata; // Pixel data loaded from TIFF image, y major
	int dotpx; // Approximate width of ref dot in pixels
	double *dots = malloc(3*2*sizeof(*dots)); // Sub-pixel ref dot locations
	int *pxcorns = malloc(4*2*sizeof(*pxcorns));

	//printf("%d, %d\n%d, %d\n%d, %d\n", *(idots), *(idots+1), *(idots+2),
	//	*(idots+3), *(idots+4), *(idots+5));

	// Might want to move but initcamera needs to be run.
	initcamera();
	initpattern(relfn);

	image = TIFFOpen(imfnh, "r");
	if (image == NULL) {
		printf("Failed to open horizontal image.\n");
		exit(1);
	}

	// Load in image data.
	TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &wid);
	TIFFGetField(image, TIFFTAG_IMAGELENGTH, &hei);
	printf("Image dimensions: %d, %d\n", wid, hei);

	imdata = (uint32 *) _TIFFmalloc(wid*hei*sizeof(uint32)); 
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

	// Find centre of dots.
	dotpx = objpixsize(dotsize, camdistguess);
	for (int i=0; i<3; i++) {
		centroid(imdata, wid, hei, (idots+i*2), (dots+i*2), dotpx);
	}

	printf("Centroided dots:\n");
	printf("%f, %f\n%f, %f\n%f, %f\n", *(dots), *(dots+1), *(dots+2),
		*(dots+3), *(dots+4), *(dots+5));

	// Locating position of camera.
	if (locatecam(dots, dotsep, camdistguess) != 0) {
		printf("Failed to locate camera, exiting...\n");
		exit(1);
	}
	printf("Camera successfully located.\n");

	printf("Bouding corners:\n");
	// Finding pixels that correspond to mirror bounding corners. 
	for (int i=0; i<4; i++) {
		findpix(corns[i], (pxcorns+i*2));
		printf("%d, %d\n", *(pxcorns+i*2), *(pxcorns+i*2+1));
	}

	//double tempdirec[3];
	//finddirpix(2129, 2918, tempdirec);
	//printf("%f, %f, %f\n", tempdirec[0], tempdirec[1], tempdirec[2]);

	free(dots);
	dots = NULL;

	// Determining max and min bounds.
	int xpixmax = imax(*(pxcorns),imax(*(pxcorns+2),
					imax(*(pxcorns+4),*(pxcorns+6))));
	int ypixmax = imax(*(pxcorns+1),imax(*(pxcorns+3),
					imax(*(pxcorns+5),*(pxcorns+7))));
	int xpixmin = imin(*(pxcorns),imin(*(pxcorns+2),
					imin(*(pxcorns+4),*(pxcorns+6))));
	int ypixmin = imin(*(pxcorns+1),imin(*(pxcorns+3),
					imin(*(pxcorns+5),*(pxcorns+7))));

	int xlen = xpixmax - xpixmin + 1;
	int ylen = ypixmax - ypixmin + 1;
	int offset[2] = {xpixmin, ypixmin};

	// Subtracting offset from pxcorns.
	// Need to use offset for indices of camera and imdata.
	for (int i=0; i<4; i++) {
		*(pxcorns+i*2) = *(pxcorns+i*2) - *(offset);
		*(pxcorns+i*2+1) = *(pxcorns+i*2+1) - *(offset+1);
	}

	printf("Size of area: %d, %d\n", xlen, ylen);
	printf("Offset: %d, %d\n", offset[0], offset[1]);
	printf("Starting pixel value: %d, %d, %d\n",
			TIFFGetR(*(imdata+(*(offset+1))*wid+*offset)),
			TIFFGetG(*(imdata+(*(offset+1))*wid+*offset)),
			TIFFGetB(*(imdata+(*(offset+1))*wid+*offset)));

	double *tcoeffs = malloc(2*sizeof(*tcoeffs));
	double *bcoeffs = malloc(2*sizeof(*bcoeffs));
	double *lcoeffs = malloc(2*sizeof(*lcoeffs));
	double *rcoeffs = malloc(2*sizeof(*rcoeffs));
	Polynom tline = {1, tcoeffs};
	Polynom bline = {1, bcoeffs};
	Polynom lline = {1, lcoeffs};
	Polynom rline = {1, rcoeffs};
	int *ybound = malloc(2*xlen*sizeof(*ybound));

	// Finding slope and offset for top bounding line.
	// First checking to make sure that we don't get infinite slope.
	// If so then polynomial set to constant zero (it shouldn't get called).
	if (*(pxcorns) == *(pxcorns+3*2)) {
		tline.degree = 0;
		*tcoeffs = 0.0;
	} else {
		*(tcoeffs+1) = ((double) (*(pxcorns+3*2+1) - *(pxcorns+1)))/
							(*(pxcorns+3*2) - *(pxcorns));
		*(tcoeffs) = (double) (*(pxcorns+1) -
							(*(tcoeffs+1))*(*(pxcorns)));
		//printf("coeffs: %f, %f\n", *(tline.coeffs), *(tline.coeffs+1));
	}

	// Finding slope and offset for bottom bounding line.
	// First checking to make sure that we don't get infinite slope.
	// If so then polynomial set to constant zero (it shouldn't get called).
	if (*(pxcorns+2*2) == *(pxcorns+1*2)) {
		bline.degree = 0;
		*bcoeffs = 0.0;
	} else {
		*(bcoeffs+1) = ((double) (*(pxcorns+2*2+1) - *(pxcorns+1*2+1)))/
							(*(pxcorns+2*2) - *(pxcorns+1*2));
		*(bcoeffs) = (double) (*(pxcorns+1*2+1) -
							(*(bcoeffs+1))*(*(pxcorns+1*2)));
		//printf("coeffs: %f, %f\n", *(bline.coeffs), *(bline.coeffs+1));
	}

	// Finding slope and offset for left bounding line.
	// First checking to make sure that we don't get infinite slope.
	// If so then polynomial set to constant zero (it shouldn't get called).
	if (*(pxcorns) == *(pxcorns+1*2)) {
		lline.degree = 0;
		*lcoeffs = 0.0;
	} else {
		*(lcoeffs+1) = ((double) (*(pxcorns+1*2+1) - *(pxcorns+1)))/
							(*(pxcorns+1*2) - *(pxcorns));
		*(lcoeffs) = (double) (*(pxcorns+1) -
							(*(lcoeffs+1))*(*(pxcorns)));
		//printf("coeffs: %f, %f\n", *(lline.coeffs), *(lline.coeffs+1));
	}

	// Finding slope and offset for right bounding line.
	// First checking to make sure that we don't get infinite slope.
	// If so then polynomial set to constant zero (it shouldn't get called).
	if (*(pxcorns+2*2) == *(pxcorns+3*2)) {
		rline.degree = 0;
		*rcoeffs = 0.0;
	} else {
		*(rcoeffs+1) = ((double) (*(pxcorns+3*2+1) - *(pxcorns+2*2+1)))/
							(*(pxcorns+3*2) - *(pxcorns+2*2));
		*(rcoeffs) = (double) (*(pxcorns+2*2+1) -
							(*(rcoeffs+1))*(*(pxcorns+2*2)));
		//printf("coeffs: %f, %f\n", *(rline.coeffs), *(rline.coeffs+1));
	}

	// Calculating bounds.
	if (*(pxcorns) == *(pxcorns+1*2)) {
		*(ybound) = *(pxcorns+1);
		*(ybound+1) = *(pxcorns+1*2+1);
	} else {
		// At most one of the two loops will execute.
		for (int x=*(pxcorns); x<=(*(pxcorns+1*2)); x++) {
			*(ybound+x*2) = (int) round(polyget(&tline,(double) x));
			*(ybound+x*2+1) = (int) round(polyget(&lline,(double) x));
		}
		for (int x=*(pxcorns+1*2); x<=(*(pxcorns)); x++) {
			*(ybound+x*2) = (int) round(polyget(&lline,(double) x));
			*(ybound+x*2+1) = (int) round(polyget(&bline,(double) x));
		}
	}

	if (*(pxcorns+2*2) == *(pxcorns+3*2)) {
		*(ybound+(xlen-1)*2) = *(pxcorns+3*2+1);
		*(ybound+(xlen-1)*2+1) = *(pxcorns+2*2+1);
	} else {
		// At most one of the two loops will execute.
		for (int x=*(pxcorns+3*2); x<=(*(pxcorns+2*2)); x++) {
			*(ybound+x*2) = (int) round(polyget(&rline,(double) x));
			*(ybound+x*2+1) = (int) round(polyget(&bline,(double) x));
		}
		for (int x=*(pxcorns+2*2); x<=(*(pxcorns+3*2)); x++) {
			*(ybound+x*2) = (int) round(polyget(&tline,(double) x));
			*(ybound+x*2+1) = (int) round(polyget(&rline,(double) x));
		}
	}

	int startmid = imax(*(pxcorns),*(pxcorns+1*2));
	int endmid = imin(*(pxcorns+2*2),*(pxcorns+3*2));
	for (int x=startmid+1; x<endmid; x++) {
		*(ybound+x*2) = (int) round(polyget(&tline,(double) x));
		*(ybound+x*2+1) = (int) round(polyget(&bline,(double) x));
	}

	printf("ybound at 0: %d, %d\n", *ybound, *(ybound+1));

	free(tcoeffs);
	tcoeffs = NULL;
	free(bcoeffs);
	bcoeffs = NULL;
	free(lcoeffs);
	lcoeffs = NULL;
	free(rcoeffs);
	rcoeffs = NULL;

	// Set up array of vectors to hold pattern positions. Later on it holds
	// normal vectors.
	double *vecs = malloc(xlen*ylen*3*sizeof(*vecs));
	if (vecs == NULL) {
		printf("Cannot allocate memory for vecs, exiting...\n");
		exit(1);
	}
	// Initialising array to -1 which is used as a check in calcpattvecs.
	for (int i=0; i<(xlen*ylen*3); i++) {
		*(vecs+i) = -1.0;
	}
	 
	// Start looking in first segment.
	calcpattvecs(imdata, wid, vecs, xlen, pxcorns, ybound, offset,
			0.5*segsize, 1);
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
	calcpattvecs(imdata, wid, vecs, xlen, pxcorns, ybound, offset,
			0.5*segsize, 0);
	printf("Information extracted from vertical image.\n");

	_TIFFfree(imdata);
	imdata = NULL;

	// Need to now find best fitting curve.
	// And fit ideal curve.

	// Piece together pattern vector components and transform into global
	// coordinates.
	transformpatt(vecs, xlen, ybound);
	printf("Pattern vectors translated to global coords.\n");

	// Set up array of vectors to hold mirror surface positions.
	double *poss = malloc(xlen*ylen*3*sizeof(*poss));
	if (poss == NULL) {
		printf("Cannot allocate memory for poss, exiting...\n");
		exit(1);
	}
	// Initialising array to -1 which is used as a check in solvesurface.
	for (int i=0; i<(xlen*ylen*3); i++) {
		*(poss+i) = -1.0;
	}

	// Testing different solution starts.
	/*
	int initpixi[2];
	//initpixi[0] = (*(pxcorns+6) - *(pxcorns))/2 + *(pxcorns);
	initpixi[0] = imax(*(pxcorns+2), *(pxcorns));
	initpixi[1] = (*(pxcorns+3) - *(pxcorns+1))/2 + *(pxcorns+1);
	printf("Initial pixel: %d, %d\n", initpixi[0], initpixi[1]);
	*/

	solvesurface(vecs, poss, xlen, pxcorns, ybound, offset, startdepth);
	printf("Surface solved.\n");

	// Centring data points about axis.
	double xmean = 0.0;
	double ymean = 0.0;
	double zmax = -100.0;
	int meancount = 0;
	for (int px=0; px<xlen; px++) {
		for (int py=*(ybound+px*2); py<=*(ybound+px*2+1); py++) {
			xmean = xmean + *(poss+(py*xlen+px)*3);
			ymean = ymean + *(poss+(py*xlen+px)*3+1);
			if (zmax < *(poss+(py*xlen+px)*3+2)) {
				zmax = *(poss+(py*xlen+px)*3+2);
			}
			meancount++;
		}
	}
	xmean = xmean/meancount;
	ymean = ymean/meancount;
	printf("xmean, ymean, zmax: %f, %f, %f\n", xmean, ymean, zmax);
	for (int px=0; px<xlen; px++) {
		for (int py=*(ybound+px*2); py<=*(ybound+px*2+1); py++) {
			*(poss+(py*xlen+px)*3) = *(poss+(py*xlen+px)*3) - xmean;
			*(poss+(py*xlen+px)*3+1) = *(poss+(py*xlen+px)*3+1) - ymean;
			*(poss+(py*xlen+px)*3+2) = *(poss+(py*xlen+px)*3+2) - zmax;
		}
	}

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
	double *spvars = malloc(4*sizeof(*spvars));
	*(spvars) = gsl_vector_get(sphfitvars, 0);
	*(spvars+1) = gsl_vector_get(sphfitvars, 1);
	*(spvars+2) = gsl_vector_get(sphfitvars, 2);
	*(spvars+3) = gsl_vector_get(sphfitvars, 3);

	gsl_vector_free(sphfitvars);
	sphfitvars = NULL;

	printf("Best fit sphere:\n");
	printf("%f, %f, %f, %f\n", *(spvars), *(spvars+1), *(spvars+2),
			*(spvars+3));

	// Need to now find best fitting paraboloid.
	gsl_vector *parfitvars = gsl_vector_alloc(8);
	gsl_vector_set(parfitvars, 0, 0.0); // x0
	gsl_vector_set(parfitvars, 1, 0.0); // y0
	gsl_vector_set(parfitvars, 2, 0.0); // z0 shift
	gsl_vector_set(parfitvars, 3, 0.0); // xrot
	gsl_vector_set(parfitvars, 4, 0.0); // yrot
	gsl_vector_set(parfitvars, 5, 0.0); // zrot
	gsl_vector_set(parfitvars, 6, *(spvars+3)/2.0); // f1
	gsl_vector_set(parfitvars, 7, *(spvars+3)/2.0); // f2

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
	printf("done\n");

	// Saving parameters into array.
	double *pavars = malloc(8*sizeof(*pavars));
	*(pavars) = gsl_vector_get(parfitvars, 0);
	*(pavars+1) = gsl_vector_get(parfitvars, 1);
	*(pavars+2) = gsl_vector_get(parfitvars, 2);
	*(pavars+3) = gsl_vector_get(parfitvars, 3);
	*(pavars+4) = gsl_vector_get(parfitvars, 4);
	*(pavars+5) = gsl_vector_get(parfitvars, 5);
	*(pavars+6) = gsl_vector_get(parfitvars, 6);
	*(pavars+7) = gsl_vector_get(parfitvars, 7);

	gsl_vector_free(parfitvars);
	parfitvars = NULL;

	printf("Best fit paraboloid:\n");
	printf("%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", *(pavars), *(pavars+1),
			*(pavars+2), *(pavars+3), *(pavars+4), *(pavars+5), *(pavars+6),
			*(pavars+7));
	
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
			outfn, *(spvars+3), *(pavars+6), *(pavars+7), *(pavars+5),
			spsest[2], spsest[3], pasest[2], pasest[3]);
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
		for (int y=*(ybound+x*2); y<=*(ybound+x*2+1); y=y+xlen/40) {
			fprintf(ssavefile, "%10.3e,%10.3e\n", *(spserr+(y*xlen+x)*2),
				*(spserr+(y*xlen+x)*2+1));
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
		for (int y=*(ybound+x*2); y<=*(ybound+x*2+1); y=y+xlen/40) {
			fprintf(ssavefile, "%10.3e,%10.3e\n", *(paserr+(y*xlen+x)*2),
				*(paserr+(y*xlen+x)*2+1));
		}
	}
	fclose(ssavefile);
	ssavefile = NULL;
 
	free(psefn);
	psefn = NULL;

	free(spvars);
	spvars = NULL;
	free(pavars);
	pavars = NULL;

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
		for (int y=*(ybound+x*2); y<=*(ybound+x*2+1); y=y+xlen/100) {
			fprintf(savefile, "%10.3e,%10.3e,%10.3e\n", *(poss+(y*xlen+x)*3),
				*(poss+(y*xlen+x)*3+1), *(poss+(y*xlen+x)*3+2));
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
		for (int y=*(ybound+x*2); y<=*(ybound+x*2+1); y=y+xlen/40) {
			fprintf(savefile, "%10.3e,%10.3e,%10.3e\n", *(poss+(y*xlen+x)*3),
				*(poss+(y*xlen+x)*3+1), *(poss+(y*xlen+x)*3+2));
		}
	}
	fclose(savefile);
	savefile = NULL;
	free(poslfn);
	poslfn = NULL;
	
	free(pxcorns);
	pxcorns = NULL;
	free(ybound);
	ybound = NULL;

	free(vecs);
	vecs = NULL;
	free(poss);
	poss = NULL;

	freecamera();
	freepattern();
}

// Determines which part of pattern is reflected for each pixel.
void calcpattvecs(const uint32 *im, const int iw, double *vecs, const int bw,
		const int *ipix, const int *yb, const int *offs, const double idist,
		const int orien) {
	// Orientation 1 is horizontal, 0 is vertical.

	// Working out first point.
	*(vecs+((*(ipix+1))*bw+*ipix)*3+orien) = getdist((im+(*(ipix+1)+
									*(offs+1))*iw+*(ipix)+*offs), idist, orien);
	// Run down first column, shouldn't need to go up since at corner.
	for (int y=*(ipix+1)+1; y<=*(yb+(*ipix)*2+1); y++) {
		// Using previously calculated distance.
		double pdist = *(vecs+((y-1)*bw+*ipix)*3+orien);
		*(vecs+(y*bw+*ipix)*3+orien) = getdist((im+(y+*(offs+1))*iw+*(ipix)+
									*offs), pdist, orien);
	}

	// Work our way left.
	for (int x=*ipix-1; x>=0; x--) {
		// Find a good starting point.
		int iy = imin(*(yb+x*2+1),imax(*(ipix+1),*(yb+x*2)));
		// Determining distance for starting point.
		double pdist = *(vecs+(iy*bw+x+1)*3+orien);
		*(vecs+(iy*bw+x)*3+orien) = getdist((im+(iy+*(offs+1))*iw+x+*offs),
										pdist, orien);
		// Shouldn't need to move up.

		// Working down.
		for (int y=iy+1; y<=*(yb+x*2+1); y++) {
			int count = 1;
			// Averaging from up to 3 surrounding distances.
			pdist = *(vecs+((y-1)*bw+x)*3+orien);
			double dcorn = *(vecs+((y-1)*bw+x+1)*3+orien);
			double dnext = *(vecs+(y*bw+x+1)*3+orien);
			if (dcorn != -1.0) {
				pdist = pdist + dcorn;
				count++;
			}
			if (dnext != -1.0) {
				pdist = pdist + dnext;
				count++;
			}
			pdist = pdist/count; // Normalizing

			*(vecs+(y*bw+x)*3+orien) = getdist((im+(y+*(offs+1))*iw+x+*offs),
										pdist, orien);
		}
	}

	// Work our way right.
	for (int x=*ipix+1; x<bw; x++) {
		// Find a good starting point.
		int iy = imin(*(yb+x*2+1),imax(*(ipix+1),*(yb+x*2)));
		// Determining distance for starting point.
		double pdist = *(vecs+(iy*bw+x-1)*3+orien);
		*(vecs+(iy*bw+x)*3+orien) = getdist((im+(iy+*(offs+1))*iw+x+*offs),
										pdist, orien);
		// Working up.
		for (int y=iy-1; y>=*(yb+x*2); y--) {
			int count = 1;
			// Averaging from up to 3 surrounding distances.
			pdist = *(vecs+((y+1)*bw+x)*3+orien);
			double dcorn = *(vecs+((y+1)*bw+x-1)*3+orien);
			double dnext = *(vecs+(y*bw+x-1)*3+orien);
			if (dcorn != -1.0) {
				pdist = pdist + dcorn;
				count++;
			}
			if (dnext != -1.0) {
				pdist = pdist + dnext;
				count++;
			}
			pdist = pdist/count; // Normalizing

			*(vecs+(y*bw+x)*3+orien) = getdist((im+(y+*(offs+1))*iw+x+*offs),
										pdist, orien);
		}
		// Working down.
		for (int y=iy+1; y<=*(yb+x*2+1); y++) {
			int count = 1;
			// Averaging from up to 3 surrounding distances.
			pdist = *(vecs+((y-1)*bw+x)*3+orien);
			double dcorn = *(vecs+((y-1)*bw+x-1)*3+orien);
			double dnext = *(vecs+(y*bw+x-1)*3+orien);
			if (dcorn != -1.0) {
				pdist = pdist + dcorn;
				count++;
			}
			if (dnext != -1.0) {
				pdist = pdist + dnext;
				count++;
			}
			pdist = pdist/count; // Normalizing

			*(vecs+(y*bw+x)*3+orien) = getdist((im+(y+*(offs+1))*iw+x+*offs),
										pdist, orien);
		}
	}
}

void transformpatt(double *vecs, const int bw, const int* yb) {
	for (int x=0; x<bw; x++) {
		for (int y=*(yb+x*2); y<=*(yb+x*2+1); y++) {
			transpattvec((vecs+(y*bw+x)*3));
		}
	}
}

void solvesurface(double *vecs, double *poss, const int bw, const int *ipix,
		const int *yb, const int *offs, const double idepth) {
	// After position and norm estimated, could redo extrapolation taking,
	// into account the estimated norm.
	double ch[3]; // Direction from camera to point
	double ph[3]; // Direction from pattern to point
	// Calculating direction from camera to point.
	finddirpix(*(ipix)+*(offs), *(ipix+1)+*(offs+1), ch);
	// Estimating position of point.
	*(poss+((*(ipix+1))*bw+*ipix)*3) = ch[0]*(idepth-campos[2])/ch[2] +
		campos[0];
	*(poss+((*(ipix+1))*bw+*ipix)*3+1) = ch[1]*(idepth-campos[2])/ch[2] +
		campos[1];
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
		finddirpix(*(ipix)+*(offs), y+*(offs+1), ch);
		double dist = extrapolate((poss+((y+1)*bw+*ipix)*3),
				(vecs+((y+1)*bw+*ipix)*3), ch);
		*(poss+(y*bw+*ipix)*3) = dist*ch[0]+ campos[0];
		*(poss+(y*bw+*ipix)*3+1) = dist*ch[1]+ campos[1];
		*(poss+(y*bw+*ipix)*3+2) = dist*ch[2]+ campos[2];
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
		finddirpix(*(ipix)+*(offs), y+*(offs+1), ch);
		double dist = extrapolate((poss+((y-1)*bw+*ipix)*3),
				(vecs+((y-1)*bw+*ipix)*3), ch);
		*(poss+(y*bw+*ipix)*3) = dist*ch[0]+ campos[0];
		*(poss+(y*bw+*ipix)*3+1) = dist*ch[1]+ campos[1];
		*(poss+(y*bw+*ipix)*3+2) = dist*ch[2]+ campos[2];
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
		finddirpix(x+*(offs), iy+*(offs+1), ch);
		double dist = extrapolate((poss+(iy*bw+x+1)*3),
				(vecs+(iy*bw+x+1)*3), ch);
		*(poss+(iy*bw+x)*3) = dist*ch[0]+ campos[0];
		*(poss+(iy*bw+x)*3+1) = dist*ch[1]+ campos[1];
		*(poss+(iy*bw+x)*3+2) = dist*ch[2]+ campos[2];
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
			finddirpix(x+*(offs), y+*(offs+1), ch);
			dist = extrapolate((poss+((y+1)*bw+x)*3),
					(vecs+((y+1)*bw+x)*3), ch);
			if (*(poss+((y+1)*bw+x+1)*3) != -1.0) {
				dist = dist + extrapolate((poss+((y+1)*bw+x+1)*3),
							(vecs+((y+1)*bw+x+1)*3), ch);
				count++;
			}
			if (*(poss+(y*bw+x+1)*3) != -1.0) {
				dist = dist + extrapolate((poss+(y*bw+x+1)*3),
							(vecs+(y*bw+x+1)*3), ch);
				count++;
			}
			dist = dist/count;
			*(poss+(y*bw+x)*3) = dist*ch[0]+ campos[0];
			*(poss+(y*bw+x)*3+1) = dist*ch[1]+ campos[1];
			*(poss+(y*bw+x)*3+2) = dist*ch[2]+ campos[2];
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
			finddirpix(x+*(offs), y+*(offs+1), ch);
			dist = extrapolate((poss+((y-1)*bw+x)*3),
					(vecs+((y-1)*bw+x)*3), ch);
			if (*(poss+((y-1)*bw+x+1)*3) != -1.0) {
				dist = dist + extrapolate((poss+((y-1)*bw+x+1)*3),
							(vecs+((y-1)*bw+x+1)*3), ch);
				count++;
			}
			if (*(poss+(y*bw+x+1)*3) != -1.0) {
				dist = dist + extrapolate((poss+(y*bw+x+1)*3),
							(vecs+(y*bw+x+1)*3), ch);
				count++;
			}
			dist = dist/count;
			*(poss+(y*bw+x)*3) = dist*ch[0]+ campos[0];
			*(poss+(y*bw+x)*3+1) = dist*ch[1]+ campos[1];
			*(poss+(y*bw+x)*3+2) = dist*ch[2]+ campos[2];
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
		finddirpix(x+*(offs), iy+*(offs+1), ch);
		double dist = extrapolate((poss+(iy*bw+x-1)*3),
				(vecs+(iy*bw+x-1)*3), ch);
		*(poss+(iy*bw+x)*3) = dist*ch[0]+ campos[0];
		*(poss+(iy*bw+x)*3+1) = dist*ch[1]+ campos[1];
		*(poss+(iy*bw+x)*3+2) = dist*ch[2]+ campos[2];
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
			finddirpix(x+*(offs), y+*(offs+1), ch);
			dist = extrapolate((poss+((y+1)*bw+x)*3),
					(vecs+((y+1)*bw+x)*3), ch);
			if (*(poss+((y+1)*bw+x-1)*3) != -1.0) {
				dist = dist + extrapolate((poss+((y+1)*bw+x-1)*3),
							(vecs+((y+1)*bw+x-1)*3), ch);
				count++;
			}
			if (*(poss+(y*bw+x-1)*3) != -1.0) {
				dist = dist + extrapolate((poss+(y*bw+x-1)*3),
							(vecs+(y*bw+x-1)*3), ch);
				count++;
			}
			dist = dist/count;
			*(poss+(y*bw+x)*3) = dist*ch[0]+ campos[0];
			*(poss+(y*bw+x)*3+1) = dist*ch[1]+ campos[1];
			*(poss+(y*bw+x)*3+2) = dist*ch[2]+ campos[2];
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
			finddirpix(x+*(offs), y+*(offs+1), ch);
			dist = extrapolate((poss+((y-1)*bw+x)*3),
					(vecs+((y-1)*bw+x)*3), ch);
			if (*(poss+((y-1)*bw+x-1)*3) != -1.0) {
				dist = dist + extrapolate((poss+((y-1)*bw+x-1)*3),
							(vecs+((y-1)*bw+x-1)*3), ch);
				count++;
			}
			if (*(poss+(y*bw+x-1)*3) != -1.0) {
				dist = dist + extrapolate((poss+(y*bw+x-1)*3),
							(vecs+(y*bw+x-1)*3), ch);
				count++;
			}
			dist = dist/count;
			*(poss+(y*bw+x)*3) = dist*ch[0]+ campos[0];
			*(poss+(y*bw+x)*3+1) = dist*ch[1]+ campos[1];
			*(poss+(y*bw+x)*3+2) = dist*ch[2]+ campos[2];
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

void surfnorm(const double *vin, const double *vre, double *snrm) {
	// vin - incident ray, vre - reflected ray.
	*(snrm) = *(vre) - *(vin);
	*(snrm+1) = *(vre+1) - *(vin+1);
	*(snrm+2) = *(vre+2) - *(vin+2);
	scale(1.0/norm(snrm), snrm);
}

double extrapolate(const double *ppos, const double *psnrm, const double *ch) {
	double pc[3] = {*ppos-campos[0],*(ppos+1)-campos[1],*(ppos+2)-campos[2]};
	return dot(pc,psnrm)/dot(ch,psnrm);
}

void slopeerror(const double *vecs, const double *poss, const int bw,
		const int *yb, const int shape, const double *params, double *serr,
		double *sestats) {
	// Calculating slope errors relative to sphere.
	// Shape 0 is sphere, shape 1 is elliptic paraboloid.
	// Parameters are in terms of acting on raw data.
	// For sphere: params = xshift yshift zshift rad.
	// For paraboloid: params = xshift yshift zshift xrot yrot zrot f1 f2.
	size_t pointc = 0;
	double *inorm = malloc(3*sizeof(*inorm)); // Ideal normal
	double mux = 0.0; // Slope error mean
	double muy = 0.0; // Slope error mean
	double *locx = malloc(3*sizeof(*locx));
	double *locy = malloc(3*sizeof(*locy));
	double xaxis[3] = {1.0, 0.0, 0.0};
	double *tpos1 = malloc(3*sizeof(*tpos1));
	double *tpos2 = malloc(3*sizeof(*tpos2));
	double *tnrm = malloc(3*sizeof(*tnrm));
	double *rotm = malloc(9*sizeof(*rotm));
	double *inrotm = malloc(9*sizeof(*rotm));
	if (shape == 1) {
		rotmatxyz(*(params+3), *(params+4), *(params+5), rotm);
		rotmatzyx(-*(params+5), -*(params+4), -*(params+3), inrotm);
	}
	for (int px=0; px<bw; px++) {
		for (int py=*(yb+px*2); py<=*(yb+px*2+1); py++) {
			// Work out ideal normal for x and y coords of given pixel.
			if (shape == 0) {
				sphereslope(*(poss+(py*bw+px)*3) + *(params),
						*(poss+(py*bw+px)*3+1) + *(params+1), *(params+3),
						inorm);
			} else if (shape == 1) {
				// Data shifted before being rotated.
				*(tpos1) = *(poss+(py*bw+px)*3) + *(params);
				*(tpos1+1) = *(poss+(py*bw+px)*3+1) + *(params+1);
				*(tpos1+2) = *(poss+(py*bw+px)*3+2) + *(params+2);
				// Data rotated. z in tpos2 not used past here.
				matxvec(rotm, tpos1, tpos2);
				parabslope(*(tpos2), *(tpos2+1), *(params+6), *(params+7),
						tnrm);
				// Rotate inorm so that it matches up with data normals.
				matxvec(inrotm, tnrm, inorm);
			}

			// Take cross product of inorm with global x axis to produce locy.
			// Note that y axis ends up pointing in roughly opposite direction
			// to global y.
			cross(inorm, xaxis, locy);
			cross(locy, inorm, locx);
			
			// Working out local x and y components of surface norm.
			double epsx = dot(locx, (vecs+(py*bw+px)*3));
			double epsy = dot(locy, (vecs+(py*bw+px)*3));

			// For now just save in vecs.
			*(serr+(py*bw+px)*2) = epsx;
			*(serr+(py*bw+px)*2+1) = epsy;

			mux = mux + epsx;
			muy = muy + epsy;
			pointc++;
		}
	}
	free(inorm);
	inorm = NULL;
	free(locx);
	locx = NULL;
	free(locy);
	locy = NULL;
	free(tpos1);
	tpos1 = NULL;
	free(tpos2);
	tpos2 = NULL;
	free(tnrm);
	tnrm = NULL;
	free(rotm);
	rotm = NULL;
	free(inrotm);
	inrotm = NULL;

	mux = mux/pointc;
	muy = muy/pointc;
	double sigx = 0.0;
	double sigy = 0.0;
	for (int px=0; px<bw; px++) {
		for (int py=*(yb+px*2); py<=*(yb+px*2+1); py++) {
			sigx = sigx + pow((*(serr+(py*bw+px)*2)-mux), 2.0);
			sigy = sigy + pow((*(serr+(py*bw+px)*2+1)-muy), 2.0);
		}
	}
	sigx = sqrt(sigx/(pointc-1));
	sigy = sqrt(sigy/(pointc-1));
	*(sestats) = mux;
	*(sestats+1) = muy;
	*(sestats+2) = sigx;
	*(sestats+3) = sigy;
}

void centroid(uint32 *im, uint32 w, uint32 h, int *idot, double *dot,
		int dotpx) {

	double *region; // Region of pixels about dot guess
	int regw, regh; // Region width and height

	// Finding the region boundaries which are 2 times dot size either side.
	// Define x along width and y along height.
	int lowx = *(idot) - 2*dotpx;
	int highx = *(idot) + 2*dotpx;
	int lowy = *(idot+1) - 2*dotpx;
	int highy = *(idot+1) + 2*dotpx;
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
			uint32 pix = *(im+(ry+lowy)*w+rx+lowx);
			greypix = TIFFGetR(pix) + TIFFGetG(pix) + TIFFGetB(pix);
			greypix = 765.0 - greypix;
			*(region+ry*regw+rx) = greypix;
		}
	}

	// Find the background intensity by looking at boundaries.
	for (int rx=0; rx<regw; rx++) {
		backgr = backgr + *(region+rx) + *(region+(regh-1)*regw+rx);
	}
	for (int ry=0; ry<regh; ry++) {
		backgr = backgr + *(region+ry*regw) + *(region+ry*regw+regh-1);
	}
	backgr = backgr/(2.0*regw + 2.0*regh);

	// Subtract background intensity and find weighting.
	for (int rx=0; rx<regw; rx++) {
		for (int ry=0; ry<regh; ry++) {
			double greypix = *(region+ry*regw+rx) - backgr;
			if (greypix < 0) greypix = 0.0;
			*(region+ry*regw+rx) = greypix;
			weight = weight + greypix;
		}
	}

	// Find x and y weighted means.
	for (int rx=0; rx<regw; rx++) {
		for (int ry=0; ry<regh; ry++) {
			xmean = xmean + rx*(*(region+ry*regw+rx))/weight;
			ymean = ymean + ry*(*(region+ry*regw+rx))/weight;
		}
	}

	free(region);
	region = NULL;

	*(dot) = xmean + lowx;
	*(dot+1) = ymean + lowy;
}

int main(int argc, char *argv[]) {
	printf("Reading in arguments...\n");

	if (argc != 10) {
		printf("Incorrect number of arguments.\n\nprofile.o himage vimage p1x p1y p2x p2y p3x p3y outfn\n\n");
		exit(1);
	}

	int *idots = malloc(3*2*sizeof(*idots));

	*(idots) = atoi(argv[3]);
	*(idots+1) = atoi(argv[4]);
	*(idots+2) = atoi(argv[5]);
	*(idots+3) = atoi(argv[6]);
	*(idots+4) = atoi(argv[7]);
	*(idots+5) = atoi(argv[8]);

	solveprofile(argv[1], argv[2], idots, argv[9]);
	free(idots);
	idots = NULL;

	return 0;
}
