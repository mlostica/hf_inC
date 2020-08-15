//sto3g test for h2 molecule
#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "basis.h"
#include "matrixElts.h"
#include "makeMats.h"

#define plen 3
#define blen 2
#define nat 2
#define hDist 1.4
// plen, number of prim orbs
// blen, number of total gtos
// nat, number of atoms
// hDist, distance btw hydrogens

int main(){

  // initialize nuclei
  // eventually these should be read in 
  nuc *hA, *hB;
  hA = (nuc *)malloc(sizeof(nuc) + sizeof(double)*5);
  hB = (nuc *)malloc(sizeof(nuc) + sizeof(double)*5);

  hA -> coords = newVec(0.0, 0.0, 0.0);
  hB -> coords = newVec(0.0, 0.0, 1.4);

  hA -> z = 1.0;
  hB -> z = 1.0;

  nuc *nuclei[nat];
  nuclei[0] = hA;
  nuclei[1] = hB;

  // gto for hydrogens A and B
  orb *oA, *oB;
  oA = (orb *)malloc(sizeof(orb) + sizeof(prim)*plen + sizeof(vector)*plen
    + sizeof(double)*plen*7);
  oB = (orb *)malloc(sizeof(orb) + sizeof(prim)*plen + sizeof(vector)*plen
    + sizeof(double)*plen*7);
  oA -> orbLen = plen;
  oB -> orbLen = plen;

  // initialize primitive orbs for each full gto
  // hard coded lol
  // eventually these should be read in 
  oA -> orbPrims[0].ex = 0.3425250914E+01;
  oA -> orbPrims[0].co = 0.1543289673;

  oB -> orbPrims[0].ex = 0.3425250914E+01;
  oB -> orbPrims[0].co = 0.1543289673;

  oA -> orbPrims[1].ex = 0.6239137298;
  oA -> orbPrims[1].co = 0.5353281423;

  oB -> orbPrims[1].ex = 0.6239137298;
  oB -> orbPrims[1].co = 0.5353281423;

  oA -> orbPrims[2].ex = 0.1688554040;
  oA -> orbPrims[2].co = 0.4446345422;

  oB -> orbPrims[2].ex = 0.1688554040;
  oB -> orbPrims[2].co = 0.4446345422;

  // assign positions for primitive orbitals
  // and normalize them
  // also initialize inidices
  int i,j,k,l;
  for (i = 0; i < plen; i++) {
    
    oA -> orbPrims[i].ctr = newVec(0.0, 0.0, 0.0);
    oB -> orbPrims[i].ctr = newVec(0.0, 0.0, 1.4);

    oA -> orbPrims[i].pN = normPrim(oA -> orbPrims[i]);
    oB -> orbPrims[i].pN = normPrim(oB -> orbPrims[i]);
  }

  // normalize total gtos
  oA -> oN = normOrbPrims( oA -> orbPrims, oA -> orbLen);
  oB -> oN = normOrbPrims( oB -> orbPrims, oB -> orbLen);

  // assign positions for orbitals
  // also initialize inidices
  //int i,j,k,l;
  //for (i = 0; i < plen; i++) {
  //  
  //  oA -> orbPrims[i].ctr = newVec(0.0, 0.0, 0.0);
  //  oB -> orbPrims[i].ctr = newVec(0.0, 0.0, 1.4);
  //}

  // make array of orbs
  orb *orbs[blen];
  orbs[0] = oA;
  orbs[1] = oB;

  // initialize tensors
  double **sMat, **tMat, **vnMat, **hCore, ****eri;

  sMat = (double **)malloc(sizeof(double *)*blen);
  tMat = (double **)malloc(sizeof(double *)*blen);
  vnMat = (double **)malloc(sizeof(double *)*blen);
  hCore = (double **)malloc(sizeof(double *)*blen);

  // compute matrices for overlap and hCore
  for (i = 0; i < blen; i++) {

    sMat[i] = (double *)malloc(sizeof(double)*blen);
    tMat[i] = (double *)malloc(sizeof(double)*blen);
    vnMat[i] = (double *)malloc(sizeof(double)*blen);
    hCore[i] = (double *)malloc(sizeof(double)*blen);

    for (j = 0; j <= i ; j++) {

      // overlap 
      sMat[i][j] = sMat[j][i] = sElt(orbs[i] -> oN, orbs[j] -> oN
        , orbs[i] -> orbPrims, orbs[j] -> orbPrims
        , orbs[i] -> orbLen, orbs[j] -> orbLen);

      // kinetic energy
      tMat[i][j] = tMat[j][i] = tElt(orbs[i] -> oN, orbs[j] -> oN
        , orbs[i] -> orbPrims, orbs[j] -> orbPrims
        , orbs[i] -> orbLen, orbs[j] -> orbLen);

      for (k = 0; k < nat; k++) {

        // electron-nuclear attraction
        vnMat[i][j] += vnElt(orbs[i] -> oN, orbs[j] -> oN
          , orbs[i] -> orbPrims, orbs[j] -> orbPrims
          , orbs[i] -> orbLen, orbs[j] -> orbLen
          , *nuclei[k]);

      vnMat[j][i] = vnMat[i][j];

      // core hamiltonian elts 
      hCore[i][j] = hCore[j][i] = tMat[i][j] + vnMat[i][j]; 
      }
    }
  }

  // compute eri tensor
  eri = (double ****)malloc(sizeof(double ***)*blen);
  for (i = 0; i < blen; i++) {

    eri[i] = (double ***)malloc(sizeof(double **)*blen);

    for (j = 0; j < blen; j++) {

      eri[i][j] = (double **)malloc(sizeof(double *)*blen);

      for (k = 0; k < blen; k++) {

        eri[i][j][k] = (double *)malloc(sizeof(double )*blen);

        for (l = 0; l < blen; l++) {
          
          eri[i][j][k][l] = eriElt(orbs[i] -> oN, orbs[j] -> oN, orbs[k] -> oN, orbs[l] -> oN
            , orbs[i] -> orbPrims, orbs[j] -> orbPrims, orbs[k] -> orbPrims, orbs[l] -> orbPrims
            , orbs[i] -> orbLen, orbs[j] -> orbLen, orbs[k] -> orbLen, orbs[l] -> orbLen);
        }
      }
    }
  }

  // verify eri
  printf("\n");
  printf("\n");
  printf("eri[0][0][0][0]  = %f\n", eri[0][0][0][0]);
  printf("eri[0][0][1][1]  = %f\n", eri[0][0][1][1]);
  printf("eri[1][0][0][0]  = %f\n", eri[1][0][0][0]);
  printf("eri[1][0][1][0]  = %f\n", eri[1][0][1][0]);
  printf("\n");
  printf("\n");


  // print matrices
  printf("sMat\n");
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {
      printf("%f ", sMat[i][j]);
    }
    printf("\n");
  }
  printf("\n\n");

  printf("hCore\n");
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {
      printf("%f ", hCore[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  printf("\n");

  // free memory
  for (i = 0; i < blen; i++) {
    free(orbs[i]);
    free(sMat[i]);
    free(tMat[i]);
    free(vnMat[i]);
    free(hCore[i]);

    for (j = 0; j < blen; j++) {

      for (k = 0; k < blen; k++) {
        free(eri[i][j][k]);
      }
      free(eri[i][j]);
    }
    free(eri[i]);
  }
  for (i = 0; i < nat; i++) {
    free(nuclei[i]);
  }
  free(sMat);
  free(tMat);
  free(vnMat);
  free(hCore);
  free(eri);

  return 0;
}
