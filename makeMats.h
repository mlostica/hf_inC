//*******************************************//
// functions for assembling matrices
//*******************************************//

#include "vector.h"
#include <math.h>
#include <stddef.h>
#include "basis.h"
#include "matrixElts.h"

//*******************************************//
// matrix construction
//*******************************************//


//////////////////
// overlap matrix S
//////////////////

double makeS(orb basis[], int blen);

//////////////////
// kinetic energy matrix T
//////////////////

double makeT(orb basis[], int blen);

//////////////////
// nuclear potential energy matrix Vn
//////////////////

double makeVn(orb basis[], int blen, nuc nuclei[], int nlen);

//////////////////
// yeet eri
//////////////////

double makeERI(orb basis[], int blen);

//////////////////
// make core Hamiltonian
//////////////////

double makeHcore(orb basis[], int blen, nuc nuclei[], int nlen);
