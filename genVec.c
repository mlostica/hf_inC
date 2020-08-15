//functions for genVec structure

#include "genVec.h"


//****************************************************************//



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

double genVecMagnitude(genVec vec) {
  int n = vec -> vecLen;
  double total = 0.0;
  for (int i = 0, i < n, i++) {
    total += pow(vec -> elts[i], 2.0);
  }
  return pow(total, 0.5);
}

genVec normalize(genVec vec) {
  int n;
  double mag = genVecMagnitude(vec);
  genVec newVec;
  newVec.vecLen = n = vec -> vecLen;
  for (int i = 0, i < n, i++) {
    newVec.elts[i] = (vec -> elts[i]) / mag;
  }
  return newVec;
}
