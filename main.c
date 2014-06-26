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

#include <stdlib.h>
#include <stdio.h>

#include "cpdeflec.h"

int main(int argc, char *argv[])
{
	printf("Reading in arguments...\n");

	if (argc != 10) {
		printf("Incorrect number of arguments.\n\n"
				"profile.o himage vimage p1x p1y p2x p2y p3x p3y outfn\n\n");
		exit(1);
	}

	int idots[3*2];

	idots[0] = atoi(argv[3]);
	idots[1] = atoi(argv[4]);
	idots[2] = atoi(argv[5]);
	idots[3] = atoi(argv[6]);
	idots[4] = atoi(argv[7]);
	idots[5] = atoi(argv[8]);

	Camera camera;
	Pattern pattern;
	Mirror mirror;

	// Setup camera parameters
	camera.dims[0] = 4288;
	camera.dims[1] = 2848;
	camera.pxsize = 0.00554;
	camera.prdist = 20.53;
	camera.soptc[0] = 0.1065;
	camera.soptc[1] = 0.2374;
	camera.rdisto[0] = -2.6652e-4;
	camera.rdisto[1] = 5.3876e-7;

	// Setup pattern parameters
	pattern.relbins = 274;
	pattern.pos[0] = 635.77408;
	pattern.pos[1] = -364.31168;
	pattern.pos[2] = -2334.4864;
	pattern.xoff = 4.6407583;
	pattern.yoff = 3.2714379;
	pattern.trans[0][0] = 0.8542358;
	pattern.trans[0][1] = 0.0018653;
	pattern.trans[0][2] = -0.5198823;
	pattern.trans[1][0] = 0.0015481;
	pattern.trans[1][1] = 0.99998;
	pattern.trans[1][2] = 0.0061316;
	pattern.trans[2][0] = 0.5198834;
	pattern.trans[2][1] = -0.0060427;
	pattern.trans[2][2] = 0.8542159;
	pattern.segsize = 47.7;
	pattern.relfn = "data/huerel.csv";

	// Setup mirror parameters
	mirror.dotsize = 7.0;
	mirror.dotsep[0] = 1192.8986;
	mirror.dotsep[1] = 1745.6891;
	mirror.dotsep[2] = 1271.5304;
	mirror.corns[0][0] = 61.010941;
	mirror.corns[0][1] = 17.819132;
	mirror.corns[0][2] = -24.028367;
	mirror.corns[1][0] = 81.010941;
	mirror.corns[1][1] = 1167.8191;
	mirror.corns[1][2] = -24.028367;
	mirror.corns[2][0] = 1211.0109;
	mirror.corns[2][1] = 1167.8191;
	mirror.corns[2][2] = -24.028367;
	mirror.corns[3][0] = 1211.0109;
	mirror.corns[3][1] = 17.819132;
	mirror.corns[3][2] = -24.028367;
	mirror.distguess = 3000.0;
	mirror.startdepth = -24.028367;

	solveprofile(camera, pattern, mirror, argv[1], argv[2], idots, argv[9]);

	return 0;
}
