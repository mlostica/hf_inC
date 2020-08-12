//*******************************************//
// functions for prim and orb strutures
//*******************************************//

#include "basis.h"
#include "vector.h"
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#define pi 3.14159265358979323846

//////////////////
// prim funcs
//////////////////

double normPrim(prim p){

  return pow((pi/(2.0*p.ex)),-0.75);
}

//////////////////
// orb funcs
//////////////////

double normOrbPrims(prim o[], int orbLen){

  // initialize indices and total to be returned
  int p, q;
  double total = 0.0;
  //int n = sizeof(o)/sizeof(o[0]);
  //int n = sizeof(o);
  //int n = sizeof(o[0]);
  //int n = orbLen;

  vector vecDiff;
  double distSq, pre1, pre2, exponent;

  for( p = 0 ; p < orbLen ; p++){
    for( q = 0 ; q < orbLen ; q++){
      
      // get square distances btw prim orbs
      vecDiff = retSubtractedVectors(o[p].ctr, o[q].ctr);
      distSq = retDotProduct(vecDiff, vecDiff);

      // compute factors
      pre1 = o[p].co* o[q].co* o[p].pN * o[q].pN;
      pre2 = pow((pi/(o[p].ex + o[q].ex)),1.5);
      exponent = exp(-1.0*((o[p].ex*o[q].ex)/(o[p].ex+o[q].ex))*distSq);

      // increment total
      total += (pre1 * pre2 * exponent);
    }
  }

  return(pow(total,-0.5));
}
