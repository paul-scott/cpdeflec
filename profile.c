#include "profile.h"

/* Mirror parameters
 * *****************
 */
static float dotsize = 12.0f; // Physiscal size of reference of dots in mm
static float dotsep[3] = {1200.0f,1697.056f,1200.0f}; // Distance between
// TL-BL, BL-TR, TR-TL
static float mirredge = 1175.0f;
static float corns[4][3] = {{35.0f,35.0f,50.0f},{35.0f,1035.0f,50.0f},
	{1035.0f,1035.0f,50.0f},{1035.0f,35.0f,50.0f}}; // Mirror boundary corners:
// TL, BL, BR, TR
static float camdistguess = 2500.0f; // Guess of distance from mirror to camera

/* LOCAL FUNCTIONS
 * ***************
 * Declared here instead of header.
 */
void extractcalib(FILE *file);
void centroid(uint32 *im, uint32 w, uint32 h, int *idot, float *dot, int dotpx);
void calcpattvecs(const uint32 *im, const int iw, float *vecs, const int bw,
		const int *ipix, const int *yb, const int *offs, const float idist,
		const int orien);

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
	initpattern();

	image = TIFFOpen(imfnh, "r");
	if (image == NULL) {
		printf("Failed to load horizontal image.\n");
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
	// Need to use offset for indices of for camera and imdata.
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
			0.5f*segsize, 0);
	printf("Information extracted from horizontal image.\n");

	free(pxcorns);
	pxcorns = NULL;

	_TIFFfree(imdata);
	imdata = NULL;

	free(ybound);
	ybound = NULL;

	free(vecs);
	vecs = NULL;

	freecamera();
}

void extractcalib(FILE *file) {

}

// Determines which part of pattern is reflected for each pixel.
void calcpattvecs(const uint32 *im, const int iw, float *vecs, const int bw,
		const int *ipix, const int *yb, const int *offs, const float idist,
		const int orien) {
	// Orientation 0 is horizontal, 1 is vertical.

	float dist = getdist((im+(*(ipix+1)+*(offs+1))*iw+*(ipix)+*offs), idist);
	// Run down first column, shouldn't need to go up since at corner.
	for (int y=*(ipix+1)+1; y<=*(yb+(*ipix)*2+1); y++) {
		*(vecs+(y*bw+*ipix)*3+orien) = getdist((im+(y+*(offs+1))*iw+*(ipix)+
									*offs), *(vecs+((y-1)*bw+*ipix)*3+orien));
		printf("%f, %d\n", *(vecs+(y*bw+*ipix)*3+orien), y); 
	}

	// Work our way left.

	/*

	for x in xrange(ipix[0]-1,-1,-1):
		# Find a value between the bounds.
		ystart = min(max(ipix[1],yb[x,0]),yb[x,1])
		pdist = vecs[x+1,ystart,orien]
		# Search at and above the starting point.
		for y in xrange(ystart,yb[x,0]-1,-1):
			pdist = patt.getdist(image[x,y], pdist)
			vecs[x,y,orien] = pdist
	*/
}

void centroid(uint32 *im, uint32 w, uint32 h, int *idot, float *dot,
		int dotpx) {

	float *region; // Region of pixels about dot guess
	int regw, regh; // Region width and height

	// Finding the region boundaries which are 3 times dot size either side.
	// Define x along width and y along height.
	int lowx = *(idot) - 3*dotpx;
	int highx = *(idot) + 3*dotpx;
	int lowy = *(idot+1) - 3*dotpx;
	int highy = *(idot+1) + 3*dotpx;
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
	float tick = 50.0f;
	float boom = 100.0f;
	int *idots = malloc(3*2*sizeof(idots));
	printf("Starting...\n");
	printf("%f, %f\n", tick-boom, boom-tick);

	*(idots) = 934;
	*(idots+1) = 332;
	*(idots+2) = 942;
	*(idots+3) = 2542;
	*(idots+4) = 3143;
	*(idots+5) = 335;
	solveprofile("./120199-h.tiff", "file2", idots, "out", "cal");
	free(idots);

	return 0;
}
