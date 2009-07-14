#include "profile.h"

/* Mirror parameters
 * *****************
 */
static float dotsize = 7.0f; // Physiscal size of reference of dots in mm
static float dotsep[3] = {1192.8986f,1745.6891f,1271.5304f}; // Distance
// between TL-BL, BL-TR, TR-TL
static float corns[4][3] = {{61.010941f,17.819132f,-24.028367f},
	{61.010941f,1167.8191f,-24.028367f},
	{1211.0109f,1167.8191f,-24.028367f},
	{1211.0109f,17.819132f,-24.028367f}}; // Mirror boundary corners:
// TL, BL, BR, TR
static float camdistguess = 3000.0f; // Guess of distance from mirror to camera
static float startdepth = -24.028367f; // Height of starting point.
static char *relfn = "huerel.csv";

/* LOCAL FUNCTIONS
 * ***************
 * Declared here instead of header.
 */
void extractcalib(FILE *file);
void centroid(uint32 *im, uint32 w, uint32 h, int *idot, float *dot, int dotpx);
void calcpattvecs(const uint32 *im, const int iw, float *vecs, const int bw,
		const int *ipix, const int *yb, const int *offs, const float idist,
		const int orien);
void transformpatt(float *vecs, const int bw, const int* yb);
float extrapolate(const float *ppos, const float *psnrm, const float *ch);
void solvesurface(float *vecs, float *poss, const int bw, const int *ipix,
		const int *yb, const int *offs, const float idepth);
void surfnorm(const float *vin, const float *vre, float *snrm);

/* Functions for creating pointers within python
 * In python for (3,2) array use:
 * ptr = profile.intalloc(3, 2)
 * To set y val of first dot to 200:
 * profile.set2darray(ptr, 0, 1, 2, 200);
 */
int *intalloc(size_t rs, size_t cs) {
	return malloc(rs*cs*sizeof(int));
}

void set2darray(int *marray, int r, int c, int cs, int val) {
	*(marray + r*cs + c) = val;
}

// MIGHT WANT TO CONSIDER CHECKING EACH MALLOC FOR ERROR.
void solveprofile(char *imfnh, char *imfnv, int *idots, char *outfn,
		char *calbfn) {
	TIFF *image; // TIFF file pointer
	uint32 wid, hei; // Image width and height
	uint32 *imdata; // Pixel data loaded from TIFF image, y major
	int R1, R2, R3;
	int dotpx; // Approximate width of ref dot in pixels
	float *dots = malloc(3*2*sizeof(*dots)); // Sub-pixel ref dot locations
	int *pxcorns = malloc(4*2*sizeof(*pxcorns));

	printf("%s\n", imfnh);
	printf("%d, %d\n%d, %d\n%d, %d\n", *(idots), *(idots+1), *(idots+2),
		*(idots+3), *(idots+4), *(idots+5));
	// Load in and extract parameters from calibration file if supplied.
	// Need to check if "" in python gives NULL pointer.
	if (outfn != NULL) {
		FILE *file = fopen(calbfn, "r");
		if (file != NULL) extractcalib(file);
		else printf("Error loading calibration file. Using defaults...\n");	
	} else printf("No calibration file given. Using defaults...\n");	

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
	printf("%d, %d\n", wid, hei);

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

	// Testing.
	R1 = TIFFGetR(*(imdata));
	R2 = TIFFGetR(*(imdata+1));
	R3 = TIFFGetR(*(imdata+wid));
	printf("R, R+1, R+wid: %d, %d, %d\n", R1, R2, R3);

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

	// Finding pixels that correspond to mirror bounding corners. 
	for (int i=0; i<4; i++) {
		findpix(corns[i], (pxcorns+i*2));
		printf("%d, %d\n", *(pxcorns+i*2), *(pxcorns+i*2+1));
	}

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
	printf("StartPix %d, %d, %d\n",
			TIFFGetR(*(imdata+(*(offset+1))*wid+*offset)),
			TIFFGetG(*(imdata+(*(offset+1))*wid+*offset)),
			TIFFGetB(*(imdata+(*(offset+1))*wid+*offset)));

	float *tcoeffs = malloc(2*sizeof(*tcoeffs));
	float *bcoeffs = malloc(2*sizeof(*bcoeffs));
	float *lcoeffs = malloc(2*sizeof(*lcoeffs));
	float *rcoeffs = malloc(2*sizeof(*rcoeffs));
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
		*tcoeffs = 0.0f;
	} else {
		*(tcoeffs+1) = ((float) (*(pxcorns+3*2+1) - *(pxcorns+1)))/
							(*(pxcorns+3*2) - *(pxcorns));
		*(tcoeffs) = (float) (*(pxcorns+1) -
							(*(tcoeffs+1))*(*(pxcorns)));
		printf("coeffs: %f, %f\n", *(tline.coeffs), *(tline.coeffs+1));
	}

	// Finding slope and offset for bottom bounding line.
	// First checking to make sure that we don't get infinite slope.
	// If so then polynomial set to constant zero (it shouldn't get called).
	if (*(pxcorns+2*2) == *(pxcorns+1*2)) {
		bline.degree = 0;
		*bcoeffs = 0.0f;
	} else {
		*(bcoeffs+1) = ((float) (*(pxcorns+2*2+1) - *(pxcorns+1*2+1)))/
							(*(pxcorns+2*2) - *(pxcorns+1*2));
		*(bcoeffs) = (float) (*(pxcorns+1*2+1) -
							(*(bcoeffs+1))*(*(pxcorns+1*2)));
		printf("coeffs: %f, %f\n", *(bline.coeffs), *(bline.coeffs+1));
	}

	// Finding slope and offset for left bounding line.
	// First checking to make sure that we don't get infinite slope.
	// If so then polynomial set to constant zero (it shouldn't get called).
	if (*(pxcorns) == *(pxcorns+1*2)) {
		lline.degree = 0;
		*lcoeffs = 0.0f;
	} else {
		*(lcoeffs+1) = ((float) (*(pxcorns+1*2+1) - *(pxcorns+1)))/
							(*(pxcorns+1*2) - *(pxcorns));
		*(lcoeffs) = (float) (*(pxcorns+1) -
							(*(lcoeffs+1))*(*(pxcorns)));
		printf("coeffs: %f, %f\n", *(lline.coeffs), *(lline.coeffs+1));
	}

	// Finding slope and offset for right bounding line.
	// First checking to make sure that we don't get infinite slope.
	// If so then polynomial set to constant zero (it shouldn't get called).
	if (*(pxcorns+2*2) == *(pxcorns+3*2)) {
		rline.degree = 0;
		*rcoeffs = 0.0f;
	} else {
		*(rcoeffs+1) = ((float) (*(pxcorns+3*2+1) - *(pxcorns+2*2+1)))/
							(*(pxcorns+3*2) - *(pxcorns+2*2));
		*(rcoeffs) = (float) (*(pxcorns+2*2+1) -
							(*(rcoeffs+1))*(*(pxcorns+2*2)));
		printf("coeffs: %f, %f\n", *(rline.coeffs), *(rline.coeffs+1));
	}

	// Calculating bounds.
	if (*(pxcorns) == *(pxcorns+1*2)) {
		*(ybound) = *(pxcorns+1);
		*(ybound+1) = *(pxcorns+1*2+1);
	} else {
		// At most one of the two loops will execute.
		for (int x=*(pxcorns); x<=(*(pxcorns+1*2)); x++) {
			*(ybound+x*2) = (int) roundf(polyget(&tline,(float) x));
			*(ybound+x*2+1) = (int) roundf(polyget(&lline,(float) x));
		}
		for (int x=*(pxcorns+1*2); x<=(*(pxcorns)); x++) {
			*(ybound+x*2) = (int) roundf(polyget(&lline,(float) x));
			*(ybound+x*2+1) = (int) roundf(polyget(&bline,(float) x));
		}
	}

	if (*(pxcorns+2*2) == *(pxcorns+3*2)) {
		*(ybound+(xlen-1)*2) = *(pxcorns+3*2+1);
		*(ybound+(xlen-1)*2+1) = *(pxcorns+2*2+1);
	} else {
		// At most one of the two loops will execute.
		for (int x=*(pxcorns+3*2); x<=(*(pxcorns+2*2)); x++) {
			*(ybound+x*2) = (int) roundf(polyget(&rline,(float) x));
			*(ybound+x*2+1) = (int) roundf(polyget(&bline,(float) x));
		}
		for (int x=*(pxcorns+2*2); x<=(*(pxcorns+3*2)); x++) {
			*(ybound+x*2) = (int) roundf(polyget(&tline,(float) x));
			*(ybound+x*2+1) = (int) roundf(polyget(&rline,(float) x));
		}
	}

	int startmid = imax(*(pxcorns),*(pxcorns+1*2));
	int endmid = imin(*(pxcorns+2*2),*(pxcorns+3*2));
	for (int x=startmid+1; x<endmid; x++) {
		*(ybound+x*2) = (int) roundf(polyget(&tline,(float) x));
		*(ybound+x*2+1) = (int) roundf(polyget(&bline,(float) x));
	}

	printf("ybound at 0:, %d, %d\n", *ybound, *(ybound+1));

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
	float *vecs = malloc(xlen*ylen*3*sizeof(*vecs));
	if (vecs == NULL) {
		printf("Cannot allocate memory for vecs, exiting...\n");
		exit(1);
	}
	// Initialising array to -1 which is used as a check in calcpattvecs.
	for (int i=0; i<(xlen*ylen*3); i++) {
		*(vecs+i) = -1.0f;
	}
	 
	// Start looking in first segment.
	calcpattvecs(imdata, wid, vecs, xlen, pxcorns, ybound, offset,
			0.5f*segsize, 1);
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
			0.5f*segsize, 0);
	printf("Information extracted from vertical image.\n");

	_TIFFfree(imdata);
	imdata = NULL;

	// Printing out a selection of pattvecs.
	int noselec = 5;
	for (int y=0; y<ylen; y=y+ylen/noselec) {
		for (int x=0; x<xlen; x=x+xlen/noselec) {
			printf("(%f,%f,%f), ", *(vecs+(y*xlen+x)*3),
					*(vecs+(y*xlen+x)*3+1), *(vecs+(y*xlen+x)*3+2));
		}
		printf("\n");
	}

	// Need to now find best fitting curve.
	// And fit ideal curve.

	// Piece together pattern vector components and transform into global
	// coordinates.
	transformpatt(vecs, xlen, ybound);
	printf("Pattern vectors translated to global coords.\n");

	// Set up array of vectors to hold pattern positions. Later on it holds
	// normal vectors.
	float *poss = malloc(xlen*ylen*3*sizeof(*poss));
	if (poss == NULL) {
		printf("Cannot allocate memory for poss, exiting...\n");
		exit(1);
	}
	// Initialising array to -1 which is used as a check in solvesurface.
	for (int i=0; i<(xlen*ylen*3); i++) {
		*(poss+i) = -1.0f;
	}

	solvesurface(vecs, poss, xlen, pxcorns, ybound, offset, startdepth);
	printf("Surface solved.\n");

	// Printing out a selection of points.
	for (int y=0; y<ylen; y=y+ylen/20) {
		printf("%f, %f, %f\n", *(poss+(y*xlen+500)*3),
				*(poss+(y*xlen+500)*3+1), *(poss+(y*xlen+500)*3+2));
	}

	// Need to now find best fitting curve.
	// And fit ideal curve.

	// Finding best fitting sphere.
	Errparams sphfitpars = {&sphere, poss, xlen, ybound, NULL};
	gsl_vector *sphfitvars = gsl_vector_alloc(4);
	// Setting starting point.
	gsl_vector_set(sphfitvars, 0, 500.0); // x0
	gsl_vector_set(sphfitvars, 1, 500.0); // y0
	gsl_vector_set(sphfitvars, 2, -10.0); // z0 shift (radius accounted)
	gsl_vector_set(sphfitvars, 3, 30000.0); // Radius

	// Step size for first trial.
	gsl_vector *sphstep = gsl_vector_alloc(4);
	gsl_vector_set(sphstep, 0, 500.0);
	gsl_vector_set(sphstep, 1, 500.0);
	gsl_vector_set(sphstep, 2, 100.0);
	gsl_vector_set(sphstep, 3, 5000.0);

	// Fitting ideal curve to data.
	minerror(&sphfitpars, 4, sphfitvars, sphstep);

	gsl_vector_free(sphstep);
	sphstep = NULL;

	gsl_vector_free(sphfitvars);
	sphfitvars = NULL;

	// Then calculate slope errors against best fitting curve and ideal curve.

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

void extractcalib(FILE *file) {

}

// Determines which part of pattern is reflected for each pixel.
void calcpattvecs(const uint32 *im, const int iw, float *vecs, const int bw,
		const int *ipix, const int *yb, const int *offs, const float idist,
		const int orien) {
	// Orientation 1 is horizontal, 0 is vertical.

	// Working out first point.
	*(vecs+((*(ipix+1))*bw+*ipix)*3+orien) = getdist((im+(*(ipix+1)+
									*(offs+1))*iw+*(ipix)+*offs), idist, orien);
	// Run down first column, shouldn't need to go up since at corner.
	for (int y=*(ipix+1)+1; y<=*(yb+(*ipix)*2+1); y++) {
		// Using previously calculated distance.
		float pdist = *(vecs+((y-1)*bw+*ipix)*3+orien);
		*(vecs+(y*bw+*ipix)*3+orien) = getdist((im+(y+*(offs+1))*iw+*(ipix)+
									*offs), pdist, orien);
	}

	// Work our way left.
	for (int x=*ipix-1; x>=0; x--) {
		// Find a good starting point.
		int iy = imin(*(yb+x*2+1),imax(*(ipix+1),*(yb+x*2)));
		// Determining distance for starting point.
		float pdist = *(vecs+(iy*bw+x+1)*3+orien);
		*(vecs+(iy*bw+x)*3+orien) = getdist((im+(iy+*(offs+1))*iw+x+*offs),
										pdist, orien);
		// Shouldn't need to move up.

		// Working down.
		for (int y=iy+1; y<=*(yb+x*2+1); y++) {
			int count = 1;
			// Averaging from up to 3 surrounding distances.
			pdist = *(vecs+((y-1)*bw+x)*3+orien);
			float dcorn = *(vecs+((y-1)*bw+x+1)*3+orien);
			float dnext = *(vecs+(y*bw+x+1)*3+orien);
			if (dcorn != -1.0f) {
				pdist = pdist + dcorn;
				count++;
			}
			if (dnext != -1.0f) {
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
		float pdist = *(vecs+(iy*bw+x-1)*3+orien);
		*(vecs+(iy*bw+x)*3+orien) = getdist((im+(iy+*(offs+1))*iw+x+*offs),
										pdist, orien);
		// Working up.
		for (int y=iy-1; y>=*(yb+x*2); y--) {
			int count = 1;
			// Averaging from up to 3 surrounding distances.
			pdist = *(vecs+((y+1)*bw+x)*3+orien);
			float dcorn = *(vecs+((y+1)*bw+x-1)*3+orien);
			float dnext = *(vecs+(y*bw+x-1)*3+orien);
			if (dcorn != -1.0f) {
				pdist = pdist + dcorn;
				count++;
			}
			if (dnext != -1.0f) {
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
			float dcorn = *(vecs+((y-1)*bw+x-1)*3+orien);
			float dnext = *(vecs+(y*bw+x-1)*3+orien);
			if (dcorn != -1.0f) {
				pdist = pdist + dcorn;
				count++;
			}
			if (dnext != -1.0f) {
				pdist = pdist + dnext;
				count++;
			}
			pdist = pdist/count; // Normalizing

			*(vecs+(y*bw+x)*3+orien) = getdist((im+(y+*(offs+1))*iw+x+*offs),
										pdist, orien);
		}
	}
}

void transformpatt(float *vecs, const int bw, const int* yb) {
	for (int x=0; x<bw; x++) {
		for (int y=*(yb+x*2); y<=*(yb+x*2+1); y++) {
			transpattvec((vecs+(y*bw+x)*3));
		}
	}
}

void solvesurface(float *vecs, float *poss, const int bw, const int *ipix,
		const int *yb, const int *offs, const float idepth) {
	// After position and norm estimated, could redo extrapolation taking,
	// into account the estimated norm.
	float ch[3]; // Direction from camera to point
	float ph[3]; // Direction from pattern to point
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
	fscale(1.0f/fnorm(ph), ph);

	surfnorm(ch, ph, (vecs+((*(ipix+1))*bw+*ipix)*3));

	// Working way down starting column.
	for (int y=*(ipix+1)+1; y<=*(yb+(*ipix)*2+1); y++) {
		finddirpix(*(ipix)+*(offs), y+*(offs+1), ch);
		float dist = extrapolate((poss+((y-1)*bw+*ipix)*3),
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
		fscale(1.0f/fnorm(ph), ph);
		surfnorm(ch, ph, (vecs+(y*bw+*ipix)*3));
	}

	// Work our way left.
	for (int x=*ipix-1; x>=0; x--) {
		// Find a good starting point.
		int iy = imin(*(yb+x*2+1),imax(*(ipix+1),*(yb+x*2)));
		finddirpix(x+*(offs), iy+*(offs+1), ch);
		float dist = extrapolate((poss+(iy*bw+x+1)*3),
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
		fscale(1.0f/fnorm(ph), ph);
		surfnorm(ch, ph, (vecs+(iy*bw+x)*3));

		// Working down.
		for (int y=iy+1; y<=*(yb+x*2+1); y++) {
			int count = 1;
			finddirpix(x+*(offs), y+*(offs+1), ch);
			dist = extrapolate((poss+((y-1)*bw+x)*3),
					(vecs+((y-1)*bw+x)*3), ch);
			if (*(poss+((y-1)*bw+x+1)*3) != -1.0f) {
				dist = dist + extrapolate((poss+((y-1)*bw+x+1)*3),
							(vecs+((y-1)*bw+x+1)*3), ch);
				count++;
			}
			if (*(poss+(y*bw+x+1)*3) != -1.0f) {
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
			fscale(1.0f/fnorm(ph), ph);
			surfnorm(ch, ph, (vecs+(y*bw+x)*3));
		}
	}

	// Work our way right.
	for (int x=*ipix+1; x<bw; x++) {
		// Find a good starting point.
		int iy = imin(*(yb+x*2+1),imax(*(ipix+1),*(yb+x*2)));
		finddirpix(x+*(offs), iy+*(offs+1), ch);
		float dist = extrapolate((poss+(iy*bw+x-1)*3),
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
		fscale(1.0f/fnorm(ph), ph);
		surfnorm(ch, ph, (vecs+(iy*bw+x)*3));

		// Working up.
		for (int y=iy-1; y>=*(yb+x*2); y--) {
			int count = 1;
			finddirpix(x+*(offs), y+*(offs+1), ch);
			dist = extrapolate((poss+((y+1)*bw+x)*3),
					(vecs+((y+1)*bw+x)*3), ch);
			if (*(poss+((y+1)*bw+x-1)*3) != -1.0f) {
				dist = dist + extrapolate((poss+((y+1)*bw+x-1)*3),
							(vecs+((y+1)*bw+x-1)*3), ch);
				count++;
			}
			if (*(poss+(y*bw+x-1)*3) != -1.0f) {
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
			fscale(1.0f/fnorm(ph), ph);
			surfnorm(ch, ph, (vecs+(y*bw+x)*3));
		}
		// Working down.
		for (int y=iy+1; y<=*(yb+x*2+1); y++) {
			int count = 1;
			finddirpix(x+*(offs), y+*(offs+1), ch);
			dist = extrapolate((poss+((y-1)*bw+x)*3),
					(vecs+((y-1)*bw+x)*3), ch);
			if (*(poss+((y-1)*bw+x-1)*3) != -1.0f) {
				dist = dist + extrapolate((poss+((y-1)*bw+x-1)*3),
							(vecs+((y-1)*bw+x-1)*3), ch);
				count++;
			}
			if (*(poss+(y*bw+x-1)*3) != -1.0f) {
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
			fscale(1.0f/fnorm(ph), ph);
			surfnorm(ch, ph, (vecs+(y*bw+x)*3));
		}
	}
}

void surfnorm(const float *vin, const float *vre, float *snrm) {
	// vin - incident ray, vre - reflected ray.
	*(snrm) = *(vre) - *(vin);
	*(snrm+1) = *(vre+1) - *(vin+1);
	*(snrm+2) = *(vre+2) - *(vin+2);
	fscale(1.0f/fnorm(snrm), snrm);
}

float extrapolate(const float *ppos, const float *psnrm, const float *ch) {
	float pc[3] = {*ppos-campos[0],*(ppos+1)-campos[1],*(ppos+2)-campos[2]};
	return fdot(pc,psnrm)/fdot(ch,psnrm);
}

void centroid(uint32 *im, uint32 w, uint32 h, int *idot, float *dot,
		int dotpx) {

	float *region; // Region of pixels about dot guess
	int regw, regh; // Region width and height

	// Finding the region boundaries which are 2 times dot size either side.
	// Define x along width and y along height.
	int lowx = *(idot) - 2*dotpx;
	int highx = *(idot) + 2*dotpx;
	int lowy = *(idot+1) - 2*dotpx;
	int highy = *(idot+1) + 2*dotpx;
	float backgr = 0.0f;
	float weight = 0.0f;
	float xmean = 0.0f;
	float ymean = 0.0f;
	
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
			float greypix;
			uint32 pix = *(im+(ry+lowy)*w+rx+lowx);
			greypix = TIFFGetR(pix) + TIFFGetG(pix) + TIFFGetB(pix);
			greypix = 765.0f - greypix;
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
	backgr = backgr/(2.0f*regw + 2.0f*regh);

	// Subtract background intensity and find weighting.
	for (int rx=0; rx<regw; rx++) {
		for (int ry=0; ry<regh; ry++) {
			float greypix = *(region+ry*regw+rx) - backgr;
			if (greypix < 0) greypix = 0.0f;
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

int main() {
	// Might want to make this a callable program instead of interfacing with
	// python.
	float tick = 50.0f;
	float boom = 100.0f;
	int *idots = malloc(3*2*sizeof(idots));
	printf("Starting...\n");
	printf("%f, %f\n", tick-boom, boom-tick);

	*(idots) = 1131;
	*(idots+1) = 671;
	*(idots+2) = 1122;
	*(idots+3) = 2342;
	*(idots+4) = 2623;
	*(idots+5) = 817;
	solveprofile("./130000-h.tiff", "./130000-v.tiff", idots, "out", "cal");
	free(idots);

	return 0;
}
