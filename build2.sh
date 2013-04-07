#! /bin/bash

gcc -o cpdeflec cpdeflec.c commath.c camera.c pattern.c fitting.c -O2 -std=c99 -lm -ltiff -lgsl -lgslcblas 
