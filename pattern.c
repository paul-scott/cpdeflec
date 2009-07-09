#include "pattern.h"

/* Pattern parameters
 * ******************
 * Get locating values from photogrammetry. Position is horizontal edge of
 * horiz rotated pattern, and vertical edge of vert rotated pattern.
 */
static float patpos[3] = {2000.0f,-500.0f,-3000.0f}; // Pos of pattern corner
static float pattrans[3][3] = {{1.0f,0.0f,0.0f},{0.0f,1.0f,0.0f},
	{0.0f,0.0f,1.0f}}; // Coordinate system translation
float segsize = 50.0f; // Width of repeating pattern segment
static float coeffs[2] = {0.0f,1.0f}; // Coeffs of hue/distance relation
static int relbins = 10;
static float *rel;
// TEMPORARY UNTIL CALIB FILE USED.
static float temprelbins[10] = {0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,
									0.9f};

float gethue(const uint32 *rgb);
float hueshift(float hue);

void initpattern() {
	// Need to drag pattern relation in from file.
	// For now use dummy data.
	rel = temprelbins;
}

float getdist(const uint32 *rgb, const float pdist) {
	// Might like to think of not changing pdist if we get a grey pixel.
	float shift = hueshift(gethue(rgb));
	float pshift = fmodf(pdist, segsize);
	printf("s, ps, %f, %f\n", shift, pshift);

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
		printf("Grey pixel.\n");
	} else if (R == maxc) {
		hue = fmodf(1.0f+(G-B)/(6.0f*(maxc-minc)),1.0f);
	} else if (G == maxc) {
		hue = (2.0f + (B-R)/(maxc-minc))/6.0f;
	} else {
		hue = (4.0f + (R-G)/(maxc-minc))/6.0f;
	}

	printf("hue: %f, %f, %f, %f ", hue, R*255.0f, G*255.0f, B*255.0f);
	return hue;
}

float hueshift(float hue) {
	// MIGHT NEED TO CHECK....
	int indice = (int) roundf(hue*relbins);
	if (indice == relbins) indice = 0;
	return segsize*(*(rel+indice));
}
