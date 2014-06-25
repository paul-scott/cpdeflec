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

#include "camera.h"
#include "commath.h"

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>

typedef struct {
	const double *dsep;
	const double *ph;
} Funcpars;

void print_state(size_t iter, gsl_multiroot_fsolver *s);
void finddirpixcam(const Camera *c, const double *pix, double *dir);
int func(const gsl_vector *p, void *params, gsl_vector *f);

void initcamera(Camera *c)
{
	c->trans = malloc(3*3*sizeof(*c->trans));
	c->itrans = malloc(3*3*sizeof(*c->itrans));
	c->pos[0] = 0.;
	c->pos[1] = 0.;
	c->pos[2] = 0.;
}

void freecamera(Camera *c)
{
	free(c->trans);
	c->trans = NULL;
	free(c->itrans);
	c->itrans = NULL;
}

int objpixsize(const Camera *c, double objsize, double objdist)
{
	return (int) ceil((objsize*c->prdist/objdist)/c->pxsize);
}

void locatecam(Camera *c, const double *dots, const double *dotsep,
		double distguess)
{
	/* Assumes first point is origin and second is global y direction from
	 * origin. Global x lies in plane of points 1-2-3.
	 */
	double gxh[3];
	double gyh[3];
	double gzh[3];
	double fdists[3];
	double ph[9];

	gsl_multiroot_fsolver *solver;
	gsl_multiroot_function multif;
	Funcpars params;
	gsl_vector *dists = gsl_vector_alloc(3);

	gsl_vector_set_all(dists, distguess);

	// Find directions from pixel values.
	for (size_t i=0; i<3; i++)
		finddirpixcam(c, (dots+i*2), (ph+i*3));

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

	int status;
	size_t iter = 0;
	do {
		status = gsl_multiroot_fsolver_iterate(solver);

		if (status == GSL_ENOPROG) {
			printf("Solver stuck, exiting...\n");
			exit(1);
		} else if (status == GSL_EBADFUNC) {
			printf("Bad function solver evaluation, exiting...\n");
			exit(1);
		}

		//print_state(iter, solver);
		status = gsl_multiroot_test_residual(solver->f, 1e-7);
		iter++;
	} while (status == GSL_CONTINUE && iter < 1000);

	fdists[0] = gsl_vector_get(solver->x, 0);
	fdists[1] = gsl_vector_get(solver->x, 1);// +15;
	fdists[2] = gsl_vector_get(solver->x, 2);//+200;

	gsl_multiroot_fsolver_free(solver);
	solver = NULL;
	
	// Find unit vectors of global coordinate system described in camera
	// array system.
	// Set dot direction vectors to positions.
	scale(fdists[0], ph+0*3);
	scale(fdists[1], ph+1*3);
	scale(fdists[2], ph+2*3);
	
	// Global y is from dot 0 to dot 1.
	gyh[0] = ph[3] - ph[0];
	gyh[1] = ph[4] - ph[1];
	gyh[2] = ph[5] - ph[2];
	scale(1.0/norm(gyh), gyh);

	// Appoximate x axis.
	gxh[0] = ph[6] - ph[0];
	gxh[1] = ph[7] - ph[1];
	gxh[2] = ph[8] - ph[2];

	// Find z axis.
	cross(gxh, gyh, gzh);
	scale(1.0/norm(gzh), gzh);
	
	// Find x axis.
	cross(gyh, gzh, gxh);

	// Set camera transform values.
	c->trans[0] = gxh[0];
	c->trans[1] = gxh[1];
	c->trans[2] = gxh[2];
	c->trans[3] = gyh[0];
	c->trans[4] = gyh[1];
	c->trans[5] = gyh[2];
	c->trans[6] = gzh[0];
	c->trans[7] = gzh[1];
	c->trans[8] = gzh[2];

	// Set first ph vector to cam position.
	scale(-1.0, ph); // Only uses first vector in ph
	matxvec(c->trans, ph, c->pos); // Only uses first vector in ph

	printf("Camera position: %f, %f, %f\n", c->pos[0], c->pos[1], c->pos[2]);

	// Calculate inverse camera transform.
	gsl_matrix *ctran = gsl_matrix_alloc(3,3);
	gsl_matrix *cinvtran = gsl_matrix_alloc(3,3);
	gsl_permutation *perm = gsl_permutation_alloc(3);
	int signum = 0;
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			gsl_matrix_set(ctran, i, j, c->trans[i*3+j]);

	gsl_linalg_LU_decomp(ctran, perm, &signum);
	gsl_linalg_LU_invert(ctran, perm, cinvtran);
	
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			c->itrans[i*3+j] = gsl_matrix_get(cinvtran, i, j);

	/* USED FOR CHECKING INVERSE
	matxvec(c->itrans, c->pos, (ph+3)); // Only uses first vector in ph
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
}

/*
 * It seems like the radial distortions are being incorrectly applied.
 * It is likely, given the context where the distortion values came from,
 * that they are from going from the pixel frame of reference to the apparent
 * pixel frame of reference.  This would suggest that the signs should be
 * positive for the finddir functions.  For the findpix function an inverse of
 * this distortion would need to be take (just taking negatives is likely
 * wrong).  However findpix might not need to be accurate as it only seems
 * to be used for finding the pixels that bound the mirror corners.
 * Actually I take that back, by playing round with raderr we can see it has
 * a big impact on the solution.
 * Should go back to original document to investigate.
 */

// Finds pixel that intercepts position vector.
void findpix(const Camera *c, const double *vec, int *pix)
{
	double pos[3]; // Position from camera
	double dir[3];
	// Changing origin.
	pos[0] = vec[0] - c->pos[0];
	pos[1] = vec[1] - c->pos[1];
	pos[2] = vec[2] - c->pos[2];

	scale(1.0/norm(pos), pos);
	//printf("%f, %f, %f\n", *(pos), *(pos+1), *(pos+2));
	// Transforming to camera coord system.
	matxvec(c->itrans, pos, dir);

	// Position on array from optical centre.
	// Might want to round here...
	// Scalling direction vector.
	scale(-c->prdist/dir[2], dir);

	// Applying radial distortion.
	double radius = sqrt(dir[0]*dir[0] + dir[1]*dir[1]);
	double raderr = 1.0 + c->rdisto[0]*pow(radius, 2.0) +
			c->rdisto[1]*pow(radius, 4.0);

	dir[0] *= raderr;
	dir[1] *= raderr;

	pix[0] = (int) ((dir[0] + c->soptc[0])/c->pxsize + c->dims[0]/2.0 - 0.5);
	pix[1] = (int) ((dir[1] + c->soptc[1])/c->pxsize + c->dims[1]/2.0 - 0.5);
}

// Calculates global direction of pixel from camera.
void finddirpix(const Camera *c, int x, int y, double *dir)
{
	// Similar to finddirpixcam except it works on integers and transforms
	// the direction vector to global coords.
	double temp[3];
	temp[0] = -(c->pxsize*(x - c->dims[0]/2.0 + 0.5) - c->soptc[0]);
	temp[1] = -(c->pxsize*(y - c->dims[1]/2.0 + 0.5) - c->soptc[1]);
	temp[2] = c->prdist;

	// Correcting for radial lens distortions.
	double radius = sqrt(temp[0]*temp[0] + temp[1]*temp[1]);
	double raderr = 1.0 - c->rdisto[0]*pow(radius, 2.0) -
			c->rdisto[1]*pow(radius, 4.0);

	temp[0] *= raderr;
	temp[1] *= raderr;

	// Normalizing vector.
	scale(1.0/norm(temp), temp);

	matxvec(c->trans, temp, dir);
}

// Calculates camera coord direction of pixel from camera. Note the use of
// double to get sub-pixel accuracy for locating.
void finddirpixcam(const Camera *c, const double *pix, double *dir)
{
	/* Define camera ax opposite to px and ay opposite to py. az is along
	 * optical axis, therefore (x,y,z)=(0,0,0) optical centre. As position
	 * of pixel in array increases, image pixel value increases.
	 * Signs used here to get the direction right. Reverse signs and we
	 * get position of pixel in the array from optical centre.
	 */
	dir[0] = -(c->pxsize*(pix[0] - c->dims[0]/2.0 + 0.5) - c->soptc[0]);
	dir[1] = -(c->pxsize*(pix[1] - c->dims[1]/2.0 + 0.5) - c->soptc[1]);
	dir[2] = c->prdist;

	// Correcting for radial lens distortions.
	double radius = sqrt(dir[0]*dir[0] + dir[1]*dir[1]);
	double raderr = 1.0 - c->rdisto[0]*pow(radius, 2.0) -
			c->rdisto[1]*pow(radius, 4.0);

	dir[0] *= raderr;
	dir[1] *= raderr;

	// Normalizing vector.
	scale(1.0/norm(dir), dir);
}

int func(const gsl_vector *p, void *params, gsl_vector *f)
{
	Funcpars *par = (Funcpars *) params;
	const double p0 = gsl_vector_get(p, 0);
	const double p1 = gsl_vector_get(p, 1);
	const double p2 = gsl_vector_get(p, 2);

	double f1 = sqrt(p0*p0 - 2*p0*p1*dot(par->ph, par->ph+1*3) + p1*p1) -
			par->dsep[0];
	double f2 = sqrt(p1*p1 - 2*p1*p2*dot(par->ph+1*3, par->ph+2*3) + p2*p2)-
			par->dsep[1];
	double f3 = sqrt(p2*p2 - 2*p2*p0*dot(par->ph+2*3, par->ph) + p0*p0)-
			par->dsep[2];

	gsl_vector_set(f, 0, f1);
	gsl_vector_set(f, 1, f2);
	gsl_vector_set(f, 2, f3);

	return GSL_SUCCESS;
}

void print_state(size_t iter, gsl_multiroot_fsolver *s)
{
	// Prints current state of solver.
	printf("iter = %3lu x = %.3f %.3f %.3f f(x) = %.3e %.3e %.3e\n",
			iter,
			gsl_vector_get(s->x, 0),
			gsl_vector_get(s->x, 1),
			gsl_vector_get(s->x, 2),
			gsl_vector_get(s->f, 0),
			gsl_vector_get(s->f, 1),
			gsl_vector_get(s->f, 2));
}
