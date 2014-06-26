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

#ifndef INC_PATTERN_H
#define INC_PATTERN_H

#include <stdint.h>

// Get locating values from photogrammetry. Position is horizontal edge of
// horiz rotated pattern, and vertical edge of vert rotated pattern.
typedef struct {
	double *rel; // Pattern colour to distance relation
	double relbins; // Number of bins in relation
	double pos[3]; // Pattern position
	double xoff; // Pattern start offset
	double yoff;
	double trans[3][3]; // Transform
	double segsize; // Width of repeating pattern segment
	char *relfn; // Filename for pattern hue relation
} Pattern;

void initpattern(Pattern *p);
void freepattern(Pattern *p);
double getdist(const Pattern *p, const uint32_t *pix, double pdist, int orien);
void transpattvec(const Pattern *p, double *vec);

#endif /* INC_PATTERN_H */
