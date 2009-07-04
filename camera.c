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
static float pixsize = 0.00554f; // Size of a pixel in mm
static float prdist = 20.52f; // Pricipal distance of camera lens in mm
static float prpoint[2] = {0.1065f,-0.2374f}; // Pricipal point in mm


void initcamera() {

}

int objpixsize(float objsize, float objdist) {
	return (int) ceilf((objsize*prdist/objdist)/pixsize);
}

int locatecam(float *dots, float *dotsep, float distguess) {
	/* Assumes first point is origin and second is global y direction from
	 * origin. Global x lies in plane of points 1-2-3.
	 */
	gsl_multiroot_folver *solver;
	gsl_multiroot_function multif;
	Funcpars params;
	gsl_vector *dists = gsl_vector_alloc(3);
	float *ph = malloc(3*3*sizeof(*ph));

	gsl_vector_set_all(dists, (double) distguess);

	// Find directions from pixel values.
	for (int i=0; i<3; i++) finddirpixcam((dots+i*2), (ph+i*3));

	params.dsep = dotsep;
	params.ph = ph;
	// Set up multiroot function.
	multif.f = &func; // System of 3 functions adress
	multif.n = 3; // Number of functions
	multif.params = &params; // Parameters to feed func

	// Create solver.
	// Other solver options are _dnewton, _broyden, _hybrids.
	solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 3);
	if (solver == NULL) {
		printf("Cannot allocate memory for camera locator, exiting...");
		exit(1);
	}
	// Give solver functions and starting vector.
	gsl_multiroot_fsolver_set(solver, &multif, dists);
	

	gsl_vector_free(dists);
	dists = NULL;

	free(ph);
	ph = NULL;

	gsl_multiroot_fsolver_free(solver);
	solver = NULL;
	/*
		# Find unit vectors of global coordinate system described in camera
		# array system.
		gyh = (dists[1]*dirs[1]-dists[0]*dirs[0])
		gyh = gyh/norm(gyh)
		gzh = cross((dists[2]*dirs[2]-dists[0]*dirs[0]), gyh)
		gzh = gzh/norm(gzh)
		gxh = cross(gyh, gzh)
		self.camtrans = array([gxh,gyh,gzh])
		self.campos = dot(self.camtrans,-dists[0]*dirs[0])
*/
	return 0; // Return zero if successful.
}

// Finds pixel that intercepts position vector.
void findpix(float *vec, int *pix) {

}

// Calculates global direction of pixel from camera.
void finddirpix(int *pix, float *dir) {

	// Note that finddirpixcam requires a float...
	// Might consider making another function that only operates on integers.
}

// Calculates camera coord direction of pixel from camera. Note the use of
// float to get sub-pixel accuracy for locating.
void finddirpixcam(float *pix, float *dir) {
	/* Define camera ax opposite to px and ay opposite to py. az is along
	 * optical axis, therefore (x,y,z)=(0,0,0) optical centre. As position
	 * of pixel in array increases, image pixel value increases.
	 * Signs used here to get the direction right. Reverse signs and we
	 * get position of pixel in the array from optical centre.
	 */
	float mag;
	*(dir) = -(pxsize*(*(pix) - *(dim)/2.0f + 0.5f) + *(prpnt));
	*(dir+1) = -(pxsize*(*(pix+1) - *(dim+1)/2.0f + 0.5f) + *(prpnt+1));
	*(dir+2) = prdist;
	// Normalizing vector.
	fscale(1.0f/fnorm(dir),dir);
}

int func(gsl_vector *p, void *params, gsl_vector *f) {
	Funcpars *par = (Funcpars *) params;
	const double p0 = gsl_vector_get(p,0);
	const double p1 = gsl_vector_get(p,1);
	const double p2 = gsl_vector_get(p,2);
	// NEED TO CHECK WHAT HAPPENS WHEN WE MULTIPLY A FLOAT WITH DOUBLE.
	// NEED TO TAKE ABS STILL. fabsf for float.
	double f1 = fabs(p0*p0-2*p0*p1*fdot((par->ph),(par->ph+1*3))+p1*p1-
			(*(par->dsep))*(*(par->dsep)));
	double f2 = 
	double f3 = 

	gsl_vector_set(f, 0, f1);
	gsl_vector_set(f, 1, f2);
	gsl_vector_set(f, 2, f3);

	return GSL_SUCCESS;
}
/*
	def system3(p, ph, l):
		out = [(norm(p[0]*ph[0]-p[1]*ph[1])-l[0])**2]
		out.append((norm(p[1]*ph[1]-p[2]*ph[2])-l[1])**2)
		out.append((norm(p[0]*ph[0]-p[2]*ph[2])-l[2])**2)
		print('err:  ')
		print(out)
		print(p)
		return out
		*/
