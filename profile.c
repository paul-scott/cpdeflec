#include "profile.h"

/* Mirror parameters
 * *****************
 */
static float dotsize = 12.0f; // Physiscal size of reference of dots in mm
static float dotsep[3] = {1200.0f,1697.056f,1200.0f}; // Distance between
// TL-BL, BL-TR, TR-TL
static float mirredge = 1175.0f;
static float corns[4][3] = {{35.0f,35.0f,50.0f},{35.0f,1035.0f,50.0f},
	{25.0f,25.0f,0.0f},{25.0f,25.0f,0.0f}}; // Mirror boundary corners:
// TL, BL, BR, TR
static float camdistguess = 2500.0f; // Guess of distance from mirror to camera

/* LOCAL FUNCTIONS
 * ***************
 * Declared here instead of header.
 */
void extractcalib(FILE *file);
void centroid(uint32 *im, uint32 w, uint32 h, int *idot, float *dot, int dotpx);

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

void solveprofile(char *imfnh, char *imfnv, int *idots, char *outfn,
		char *calbfn) {
	TIFF *image; // TIFF file pointer
	uint32 w, h; // Image width and height
	uint32 *imdata; // Pixel data loaded from TIFF image, y major
	int R1, R2, G1, G2, B1, B2, A;
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

	image = TIFFOpen(imfnh, "r");
	if (image == NULL) {
		printf("Failed to load horizontal image.\n");
		exit(1);
	}

	// Load in image data.
	TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &w);
	TIFFGetField(image, TIFFTAG_IMAGELENGTH, &h);
	printf("%d, %d\n", w, h);

	imdata = (uint32 *) _TIFFmalloc(w*h*sizeof(uint32)); 
	if (imdata == NULL) {
		printf("Failed to allocate memory for horizontal image.\n");
		exit(1);
	}
	if (TIFFReadRGBAImageOriented(image, w, h, imdata, ORIENTATION_TOPLEFT,
			0) == 0) {
		printf("Failed to read horizontal image.\n");
		exit(1);
	}
	printf("Horizontal image loaded\n");

	// TIFFGetR is a macro.
	R1 = TIFFGetR(*imdata);
	R2 = TIFFGetR(*(imdata+w));
	G1 = TIFFGetG(*imdata);
	G2 = TIFFGetG(*(imdata+w));
	B1 = TIFFGetB(*imdata);
	B2 = TIFFGetB(*(imdata+w));
	A = TIFFGetA(*imdata);
	printf("%d,%d,%d,%d,%d,%d,%d\n",A,R1,G1,B1,R2,G2,B2);

	// Find centre of dots.
	dotpx = objpixsize(dotsize, camdistguess);
	for (int i=0; i<3; i++) {
		centroid(imdata, w, h, (idots+i*2), (dots+i*2), dotpx);
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
		getpix((float *) (corns+i*2), (pxcorns+i*2));
	}
	
	free(dots);
	dots = NULL;

	free(pxcorns);
	pxcorns = NULL;

	_TIFFfree(imdata);
	imdata = NULL;

	TIFFClose(image);
	image = NULL;

}

void extractcalib(FILE *file) {

	initcamera();
	initpattern();
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
