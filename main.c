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

	int *idots = malloc(3*2*sizeof(*idots));

	*(idots) = atoi(argv[3]);
	*(idots+1) = atoi(argv[4]);
	*(idots+2) = atoi(argv[5]);
	*(idots+3) = atoi(argv[6]);
	*(idots+4) = atoi(argv[7]);
	*(idots+5) = atoi(argv[8]);

	solveprofile(argv[1], argv[2], idots, argv[9]);
	free(idots);
	idots = NULL;

	return 0;
}
