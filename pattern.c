#include "pattern.h"

/* Pattern parameters
 * ******************
 * Get locating values from photogrammetry. Position is horizontal edge of
 * horiz rotated pattern, and vertical edge of vert rotated pattern.
 */
// NOTE NEED TO FIND BETTER VALUE.
static float patpos[3] = {635.77408f,-364.31168f,-2334.4864f}; // Pos of
static float patxoff = 4.6407583;
static float patyoff = 3.2714379;
// pattern corner
// NOTE MIGHT WANT TO WORK ON NON ORTHOGONAL VECTORS.
static float pattrans[3][3] = {{0.8542358f,0.0018653f,-0.5198823f},
	{0.0015481f,0.99998f,0.0061316f},
	{0.5198834f,-0.0060427f,0.8542159f}}; // Coordinate system translation
float segsize = 47.7f; // Width of repeating pattern segment
static int relbins = 274;
static float *rel;

float gethue(const uint32 *rgb);
float hueshift(const float hue, const int orien);

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

void transpattvec(float *vec) {
	float temp[3];
	// Need to set z component to zero. 
	*(vec) = *(vec) + patxoff;
	*(vec+1) = *(vec+1) + patyoff;
	*(vec+2) = 0.0f;
	//printf("pat: %f, %f, %f ", *vec, *(vec+1), *(vec+2));
	fmatxvec((float *) pattrans, vec, temp);
	*(vec) = temp[0] + patpos[0];
	*(vec+1) = temp[1] + patpos[1];
	*(vec+2) = temp[2] + patpos[2];
	//printf("glo: %f, %f, %f\n", *vec, *(vec+1), *(vec+2));
}

float getdist(const uint32 *rgb, const float pdist, const int orien) {
	// Might like to think of not changing pdist if we get a grey pixel.
	float shift = hueshift(gethue(rgb), orien);
	float pshift = fmodf(pdist, segsize);
	//printf("s, ps, %f, %f\n", shift, pshift);

	if (shift < (pshift-0.5f*segsize)) {
		return (pdist - pshift + segsize + shift);
	} else if (shift > (pshift+0.5f*segsize)) {
		return (pdist - pshift - segsize + shift);
	} else {
		return (pdist - pshift + shift);
	}
}

float gethue(const uint32 *rgb) {
	float R = (float) TIFFGetR(*rgb)/255.0f;
	float G = (float) TIFFGetG(*rgb)/255.0f;
	float B = (float) TIFFGetB(*rgb)/255.0f;

	float maxc = fmaxf(R,fmaxf(G,B));
	float minc = fminf(R,fminf(G,B));
	// float L = (maxc + minc)/2.0f;
	float hue = 0.0f;
	
	if (maxc == minc) {
		hue = 0.0f;
		//printf("Grey pixel.\n");
	} else if (R == maxc) {
		hue = fmodf(1.0f+(G-B)/(6.0f*(maxc-minc)),1.0f);
	} else if (G == maxc) {
		hue = (2.0f + (B-R)/(maxc-minc))/6.0f;
	} else {
		hue = (4.0f + (R-G)/(maxc-minc))/6.0f;
	}

	//printf("hue: %f, %f, %f, %f ", hue, R*255.0f, G*255.0f, B*255.0f);
	return hue;
}

float hueshift(const float hue, const int orien) {
	// MIGHT NEED TO CHECK....
	int indice = (int) roundf(hue*relbins);

	if (indice == relbins) indice = 0;
	// Need to reverse relation depending on orientation.
	if (orien) {
		return segsize*(*(rel+indice));
	} else {
		return segsize*(1.0f-(*(rel+indice)));
	}
}
