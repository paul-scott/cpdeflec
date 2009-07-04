%module profile
%{
/* Could use #include "profile.h" instead if we have header file, otherwise
 * need to use extern since space is already carved out in profile.c.
 * Functions don't need extern since they are counted as extern by default.
 */
#include "profile.h"
%}

void solveprofile(char *imfnh, char *imfnv, int *idots, char *outfn,
    char *calbfn);
int *intalloc(size_t rs, size_t cs);
void free(void *__ptr); 
void set2darray(int *marray, int r, int c, int cs, int val);
