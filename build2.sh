#! /bin/bash

gcc -lm -ltiff -lgsl -lgslcblas -o profile.o profile.c commath.c camera.c pattern.c fitting.c -std=c99 -ltiff
