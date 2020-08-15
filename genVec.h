// vector structure
// used for positions of nuclei
// based on code from Rabani group

#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>

#define sqrt(x)  ((x) * (x))

// structure declarations

typedef struct genVec_ {
  int vecLen;
  double mag;
  double elts[];
} genVec;

// functions for vector structure in genVec.c
double genVecMagnitude(genVec vec);
genVec normalize(genVec vec);

#endif
