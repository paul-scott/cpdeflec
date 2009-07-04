#! /bin/bash

swig -python profile.i
gcc -c -fPIC profile.c commath.c camera.c pattern.c profile_wrap.c -I/usr/include/python2.6/ -std=c99
ld -shared -ltiff -lm -lgsl -lgslcblas profile.o commath.o camera.o pattern.o profile_wrap.o -o _profile.so -ltiff

rm profile.o commath.o camera.o pattern.o profile_wrap.o profile_wrap.c
