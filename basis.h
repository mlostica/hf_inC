//*******************************************//
// basis strcture
// for s gaussian-type orbitals
//*******************************************//

#include "vector.h"
#include <math.h>
#include <stddef.h>

//*******************************************//
// structure declarations
//*******************************************//


//////////////////
// primitive gto
//////////////////
#ifndef PRIM_H
#define PRIM_H

typedef struct prim_{
  
  // ctr, where gto is ctrd
  vector ctr;

  // ex, exponent
  // co, coefficient
  // pN, primitive normalization const
  double ex, co, pN;
  
} prim;

//////////////////
// prim funcs
//////////////////

double normPrim(prim p);

#endif


//////////////////
// gto
//////////////////
#ifndef ORB_H
#define ORB_H 

typedef struct orb_ {

  // orbLen, number of prims in orb
  int orbLen;

  // oN, orb normalization const
  double oN;

  // prims, array of primitive gtos
  prim orbPrims[];
  
} orb;

//////////////////
// orb funcs
//////////////////
double normOrbPrims(prim o[], int orbLen);

#endif

//////////////////
// nucleus strut
//////////////////
#ifndef NUC_H
#define NUC_H 

typedef struct nuc_ {

  // orbLen, number of prims in orb
  vector coords;

  // z, charge of nucleus 
  double z;

} nuc;

#endif
