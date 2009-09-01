#include "camera.h"

/* Local Structures
 */
typedef struct {
	float *dsep; // Need to change to double??
	float *ph;
} Funcpars;

/* Camera parameters
 * *****************
 * Get the principal values from photogrammetry calibration of camera.
 * Need to sort out coord system of principal point from photogrammetry.
 * For now assume that prpoint is position where optical axis strikes array
 * relative to centre of array. And that increase in direction related to px
 * corresponds in increase in ax dir.
 */
// Can probably pull w and h in from profile.c:
static int camdims[2] = {4288,2848}; // Width and height of CCD array in pixels
static float pxsize = 0.00554f; // Size of a pixel in mm
static float prdist = 20.53f; // Pricipal distance of camera lens in mm
static float prpoint[2] = {0.1065f,0.2374f}; // Pricipal point in mm
// Take negative of vms y component.

float campos[3] = {0.0f,0.0f,0.0f};
static float *camtrans;
static float *caminvtrans;

void print_state(size_t iter, gsl_multiroot_fsolver *s);
void finddirpixcam(const float *pix, float *dir);
int func(const gsl_vector *p, void *params, gsl_vector *f);

void initcamera() {
	camtrans = malloc(3*3*sizeof(*camtrans));
	caminvtrans = malloc(3*3*sizeof(*caminvtrans));
}

void freecamera() {
	free(camtrans);
	camtrans = NULL;
	free(caminvtrans);
	caminvtrans = NULL;
}

int objpixsize(float objsize, float objdist) {
	return (int) ceilf((objsize*prdist/objdist)/pxsize);
}

int locatecam(float *dots, float *dotsep, float distguess) {
	/* Assumes first point is origin and second is global y direction from
	 * origin. Global x lies in plane of points 1-2-3.
	 */
	float *gxh = malloc(3*sizeof(*gxh));
	float *gyh = malloc(3*sizeof(*gyh));
	float *gzh = malloc(3*sizeof(*gzh));
	float *fdists = malloc(3*sizeof(*fdists));

	gsl_multiroot_fsolver *solver;
	gsl_multiroot_function multif;
	Funcpars params;
	gsl_vector *dists = gsl_vector_alloc(3);
	float *ph = malloc(3*3*sizeof(*ph));
	int status;
	size_t iter = 0;

	gsl_vector_set_all(dists, (double) distguess);

	// Find directions from pixel values.
	for (int i=0; i<3; i++) finddirpixcam((dots+i*2), (ph+i*3));

	params.dsep = dotsep;
	params.ph = ph;
	// Set up multiroot function.
	multif.f = &func; // System of 3 functions address
	multif.n = 3; // Number of functions
	multif.params = &params; // Parameters to feed func

	// Create solver.
	// Other solver options are _dnewton, _broyden, _hybrid.
	solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, 3);
	if (solver == NULL) {
		printf("Cannot allocate memory for camera locator, exiting...");
		exit(1);
	}
	// Give solver functions and starting vector.
	gsl_multiroot_fsolver_set(solver, &multif, dists);

	do {
		status = gsl_multiroot_fsolver_iterate(solver);

		if (status == GSL_ENOPROG) {
			printf("Solver stuck, exiting...\n");
			exit(1);
		} else if (status == GSL_EBADFUNC) {
			printf("Bad function solver evaluation, exiting...\n");
			exit(1);
		}

		print_state(iter, solver);
 
		status = gsl_multiroot_test_residual(solver->f, 1e-7);

		iter++;
	} while (status == GSL_CONTINUE && iter < 1000);

	*(fdists) = (float) gsl_vector_get(solver->x, 0);
	*(fdists+1) = (float) gsl_vector_get(solver->x, 1);
	*(fdists+2) = (float) gsl_vector_get(solver->x, 2);

	gsl_multiroot_fsolver_free(solver);
	solver = NULL;
	
	// Find unit vectors of global coordinate system described in camera
	// array system.
	// Set dot direction vectors to positions.
	fscale(*(fdists), ph);
	fscale(*(fdists+1), ph+1*3);
	fscale(*(fdists+2), ph+2*3);
	
	// Global y is from dot 0 to dot 1.
	*(gyh) = *(ph+1*3) - *(ph);
	*(gyh+1) = *(ph+1*3+1) - *(ph+1);
	*(gyh+2) = *(ph+1*3+2) - *(ph+2);
	fscale(1.0f/fnorm(gyh), gyh);

	// Appoximate x axis.
	*(gxh) = *(ph+2*3) - *(ph);
	*(gxh+1) = *(ph+2*3+1) - *(ph+1);
	*(gxh+2) = *(ph+2*3+2) - *(ph+2);

	// Find z axis.
	fcross(gxh, gyh, gzh);
	fscale(1.0f/fnorm(gzh), gzh);
	
	// Find x axis.
	fcross(gyh, gzh, gxh);

	// Set camera transform values.
	*(camtrans) = *(gxh);
	*(camtrans+1) = *(gxh+1);
	*(camtrans+2) = *(gxh+2);
	*(camtrans+1*3) = *(gyh);
	*(camtrans+1*3+1) = *(gyh+1);
	*(camtrans+1*3+2) = *(gyh+2);
	*(camtrans+2*3) = *(gzh);
	*(camtrans+2*3+1) = *(gzh+1);
	*(camtrans+2*3+2) = *(gzh+2);

	// Set first ph vector to cam position.
	fscale(-1.0f, ph); // Only uses first vector in ph
	fmatxvec(camtrans, ph, campos); // Only uses first vector in ph

	printf("campos: %f, %f, %f\n", *campos, *(campos+1), *(campos+2));

	// Calculate inverse camera transform.
	gsl_matrix *ctran = gsl_matrix_alloc(3,3);
	gsl_matrix *cinvtran = gsl_matrix_alloc(3,3);
	gsl_permutation *perm = gsl_permutation_alloc(3);
	int signum = 0;
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			gsl_matrix_set(ctran, i, j, (double) *(camtrans+i*3+j));
			//printf("%f, ", *(camtrans+i*3+j));
		}
	}
	gsl_linalg_LU_decomp(ctran, perm, &signum);
	gsl_linalg_LU_invert(ctran, perm, cinvtran);
	
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			*(caminvtrans+i*3+j) = (float) gsl_matrix_get(cinvtran, i, j);
			//printf("%f, ", *(caminvtrans+i*3+j));
		}
	}

	/* USED FOR CHECKING INVERSE
	fmatxvec(caminvtrans, campos, (ph+3)); // Only uses first vector in ph
	printf("orig: %f, %f, %f\n", *ph, *(ph+1), *(ph+2));
	printf("transed: %f, %f, %f\n", *(ph+3), *(ph+3+1), *(ph+3+2));
	*/

	gsl_matrix_free(ctran);
	ctran = NULL;
	gsl_matrix_free(cinvtran);
	cinvtran = NULL;
	gsl_permutation_free(perm);
	perm = NULL;

	gsl_vector_free(dists);
	dists = NULL;

	free(ph);
	ph = NULL;

	free(gxh);
	gxh = NULL;
	free(gyh);
	gyh = NULL;
	free(gzh);
	gzh = NULL;
	free(fdists);
	fdists = NULL;

	return 0; // Return zero if successful.
}

// Finds pixel that intercepts position vector.
void findpix(const float *vec, int *pix) {
	float *pos = malloc(3*sizeof(*pos)); // Position from camera
	float *dir = malloc(3*sizeof(*dir));
	// Changing origin.
	*(pos) = *(vec) - *(campos);
	*(pos+1) = *(vec+1) - *(campos+1);
	*(pos+2) = *(vec+2) - *(campos+2);
	// Transforming to camera coord system.
	fmatxvec(caminvtrans, pos, dir);
	free(pos);
	pos = NULL;

	// Position on array from optical centre.
	// Might want to round here...
	fscale(-prdist/(*(dir+2)), dir);
	printf("%f, %f, %f\n", *(dir), *(dir+1), *(dir+2));
	*(pix) = (int) ((*(dir) + prpoint[0])/pxsize + camdims[0]/2.0f - 0.5f);
	*(pix+1) = (int) ((*(dir+1) + prpoint[1])/pxsize + camdims[1]/2.0f - 0.5f);

	free(dir);
	dir = NULL;
}

// Calculates global direction of pixel from camera.
void finddirpix(const int x, const int y, float *dir) {
	// Similar to finddirpixcam except it works on integers and transforms
	// the direction vector to global coords.
	float temp[3];
	temp[0] = -(pxsize*(x - *(camdims)/2.0f + 0.5f) - *(prpoint));
	temp[1] = -(pxsize*(y - *(camdims+1)/2.0f + 0.5f) - *(prpoint+1));
	temp[2] = prdist;
	// Normalizing vector.
	fscale(1.0f/fnorm(temp),temp);

	fmatxvec(camtrans, temp, dir);
}

// Calculates camera coord direction of pixel from camera. Note the use of
// float to get sub-pixel accuracy for locating.
void finddirpixcam(const float *pix, float *dir) {
	/* Define camera ax opposite to px and ay opposite to py. az is along
	 * optical axis, therefore (x,y,z)=(0,0,0) optical centre. As position
	 * of pixel in array increases, image pixel value increases.
	 * Signs used here to get the direction right. Reverse signs and we
	 * get position of pixel in the array from optical centre.
	 */
	*(dir) = -(pxsize*(*(pix) - *(camdims)/2.0f + 0.5f) - *(prpoint));
	*(dir+1) = -(pxsize*(*(pix+1) - *(camdims+1)/2.0f + 0.5f) - *(prpoint+1));
	*(dir+2) = prdist;
	// Normalizing vector.
	fscale(1.0f/fnorm(dir),dir);
}

int func(const gsl_vector *p, void *params, gsl_vector *f) {
	Funcpars *par = (Funcpars *) params;
	const double p0 = gsl_vector_get(p,0);
	const double p1 = gsl_vector_get(p,1);
	const double p2 = gsl_vector_get(p,2);
	// NEED TO CHECK WHAT HAPPENS WHEN WE MULTIPLY A FLOAT WITH DOUBLE.
	double f1 = p0*p0-2*p0*p1*fdot((par->ph),(par->ph+1*3))+p1*p1-
			(*(par->dsep))*(*(par->dsep));
	double f2 = p1*p1-2*p1*p2*fdot((par->ph+1*3),(par->ph+2*3))+p2*p2-
			(*(par->dsep+1))*(*(par->dsep+1));
	double f3 = p2*p2-2*p2*p0*fdot((par->ph+2*3),(par->ph))+p0*p0-
			(*(par->dsep+2))*(*(par->dsep+2));

	gsl_vector_set(f, 0, f1);
	gsl_vector_set(f, 1, f2);
	gsl_vector_set(f, 2, f3);

	return GSL_SUCCESS;
}

void print_state(size_t iter, gsl_multiroot_fsolver *s) {
	// Prints current state of solver.
	printf("iter = %3u x = %.3f %.3f %.3f f(x) = %.3e %.3e %.3e\n",
			iter,
			gsl_vector_get(s->x, 0),
			gsl_vector_get(s->x, 1),
			gsl_vector_get(s->x, 2),
			gsl_vector_get(s->f, 0),
			gsl_vector_get(s->f, 1),
			gsl_vector_get(s->f, 2));
}

