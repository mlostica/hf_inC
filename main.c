//sto3g test
#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "basis.h"

int main(){
  orb *o;
  o = (orb *)malloc(sizeof(orb) + sizeof(prim)*3);
  o -> orbLen = 3;

  printf("%i \n", o -> orbLen);

  o -> orbPrims[0].ex = 0.3425250914E+01;
  o -> orbPrims[0].co = 0.1543289673;
  o -> orbPrims[0].pN = normPrim(o -> orbPrims[0]);
  printf("Normalization const for first prim %f\n", o -> orbPrims[0].pN);

  o -> orbPrims[1].ex = 0.6239137298;
  o -> orbPrims[1].co = 0.5353281423;
  o -> orbPrims[1].pN = normPrim(o -> orbPrims[1]);
  printf("Normalization const for first prim %f\n", o -> orbPrims[1].pN);

  o -> orbPrims[2].ex= 0.1688554040;
  o -> orbPrims[2].co = 0.4446345422;
  o -> orbPrims[2].pN = normPrim(o -> orbPrims[2]);
  printf("Normalization const for first prim %f\n", o -> orbPrims[2].pN);

  o -> oN = normOrbPrims( o -> orbPrims, o -> orbLen);
  //printf("%f\n", o.orbPrims[1].co);
  printf("Normalization for full gto %f \n", o -> oN);
  return 0;
}
