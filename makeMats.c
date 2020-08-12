//*******************************************//
// functions for computing matrix elements
//*******************************************//

#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "vector.h"
#include "basis.h"
#include "matrixElts.h"
#include "makeMats.h"

#define pi 3.14159265358979323846

//*******************************************//
// matrix construction
//*******************************************//


//////////////////
// overlap matrix S
//////////////////

double makeS(orb basis[], int blen);
  // basis[], set of total gtos
  // blen, number of gtos

  // initialize indices
  // and array to be returned
  int i,j;
  double s[blen][blen];

  // loop over all primitive orbitals
  for( i = 0 ; i < blen ; i++){
    for( j = 0 ; j < blen ; j++){

      s[i][j] = sElt(basis[i].oN, basis[j].oN
        , basis[i].orbPrims, basis[j].orbPrims
        , basis[i].orbLen, basis[j].orbLen);
    }
  }
  return s;
}

