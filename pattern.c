#include "pattern.h"

/* Pattern parameters
 * ******************
 * Get locating values from photogrammetry. Position is horizontal edge of
 * horiz rotated pattern, and vertical edge of vert rotated pattern.
 */
static float patpos[3] = {2000.0f,-500.0f,-3000.0f}; // Pos of pattern corner
static float pattrans[3][3] = {{1.0f,0.0f,0.0f},{0.0f,1.0f,0.0f},
	{0.0f,0.0f,1.0f}}; // Coordinate system translation
static float segsize = 50.0f; // Width of repeating pattern segment
static float coeffs[2] = {0.0f,1.0f}; // Coeffs of hue/distance relation
static Polynom pathuerela = {2,{0.0f,1.0f}};


void initpattern() {

}
