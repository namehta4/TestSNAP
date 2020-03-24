/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Aidan Thompson, Christian Trott, SNL
------------------------------------------------------------------------- */

#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <omp.h>
#include "sna.h"
#include "memory.h"
#include <iostream>
/* ----------------------------------------------------------------------

   this implementation is based on the method outlined
   in Bartok[1], using formulae from VMK[2].

   for the Clebsch-Gordan coefficients, we
   convert the VMK half-integral labels
   a, b, c, alpha, beta, gamma
   to array offsets j1, j2, j, m1, m2, m
   using the following relations:

   j1 = 2*a
   j2 = 2*b
   j =  2*c

   m1 = alpha+a      2*alpha = 2*m1 - j1
   m2 = beta+b    or 2*beta = 2*m2 - j2
   m =  gamma+c      2*gamma = 2*m - j

   in this way:

   -a <= alpha <= a
   -b <= beta <= b
   -c <= gamma <= c

   becomes:

   0 <= m1 <= j1
   0 <= m2 <= j2
   0 <= m <= j

   and the requirement that
   a+b+c be integral implies that
   j1+j2+j must be even.
   The requirement that:

   gamma = alpha+beta

   becomes:

   2*m - j = 2*m1 - j1 + 2*m2 - j2

   Similarly, for the Wigner U-functions U(J,m,m') we
   convert the half-integral labels J,m,m' to
   array offsets j,ma,mb:

   j = 2*J
   ma = J+m
   mb = J+m'

   so that:

   0 <= j <= 2*Jmax
   0 <= ma, mb <= j.

   For the bispectrum components B(J1,J2,J) we convert to:

   j1 = 2*J1
   j2 = 2*J2
   j = 2*J

   and the requirement:

   |J1-J2| <= J <= J1+J2, for j1+j2+j integral

   becomes:

   |j1-j2| <= j <= j1+j2, for j1+j2+j even integer

   or

   j = |j1-j2|, |j1-j2|+2,...,j1+j2-2,j1+j2

   [1] Albert Bartok-Partay, "Gaussian Approximation..."
   Doctoral Thesis, Cambrindge University, (2009)

   [2] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii,
   "Quantum Theory of Angular Momentum," World Scientific (1988)

------------------------------------------------------------------------- */

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))
static const SNADOUBLE MY_PI  = 3.14159265358979323846; // pi

SNA::SNA(int nnatoms, int nnbor, Memory* memory_in, SNADOUBLE rfac0_in,
         int twojmax_in,
         SNADOUBLE rmin0_in, int switch_flag_in, int /*bzero_flag_in*/)
{
  wself = 1.0;
  
  num_atoms = nnatoms;
  num_nbors = nnbor;
  memory = memory_in;
  rfac0 = rfac0_in;
  rmin0 = rmin0_in;
  switch_flag = switch_flag_in;

  twojmax = twojmax_in;

  ncoeff = compute_ncoeff();
  nmax = 0;

  create_twojmax_arrays();
}

/* ---------------------------------------------------------------------- */

SNA::~SNA()
{
}

/* ---------------------------------------------------------------------- */

void SNA::omp_offload_init()
{
#pragma omp target enter data map(to:this[0:1])
#pragma omp target enter data map(to: \
     this->rij.dptr[0:this->rij.size], \
     this->rcutij.dptr[0:this->rcutij.size], \
     this->wj.dptr[0:this->wj.size], \
     this->beta.dptr[0:this->beta.size], \
     this->fact.dptr[0:this->fact.size], \
     this->ulist_parity.dptr[0:this->ulist_parity.size], \
     this->idxz_trd.dptr[0:this->idxz_trd.size], \
     this->idxb_trd.dptr[0:this->idxb_trd.size], \
     this->idxc_trd.dptr[0:this->idxc_trd.size], \
     this->idxcg_block.dptr[0:this->idxcg_block.size], \
     this->idxdu_block.dptr[0:this->idxdu_block.size]) nowait

#pragma omp target enter data map(alloc: \
     this->idxb_block.dptr[0:this->idxb_block.size], \
     this->idxzbeta.dptr[0:this->idxzbeta.size], \
     this->idxz.dptr[0:this->idxz.size], \
     this->cglist.dptr[0:this->cglist.size], \
     this->rootpqarray.dptr[0:this->rootpqarray.size], \
     this->rootpqparityarray.dptr[0:this->rootpqparityarray.size], \
     this->ulisttot.dptr[0:this->ulisttot.size]) nowait

#pragma omp target enter data map(to: \
     this->ulist.dptr[0:this->ulist.size], \
     this->dulist.dptr[0:this->dulist.size], \
     this->ylist.dptr[0:this->ylist.size], \
     this->dedr.dptr[0:this->dedr.size]) nowait
}

/* ---------------------------------------------------------------------- */

void SNA::omp_offload_update()
{
#pragma omp target update to \
  (this->rij.dptr[0:this->rij.size], \
   this->rcutij.dptr[0:this->rcutij.size], \
   this->wj.dptr[0:this->wj.size])
}

/* ---------------------------------------------------------------------- */


void SNA::build_indexlist()
{
  int jdim = twojmax + 1;
  int j1j2tot = jdim*(jdim+1)/2;
  int idxb_count, idxz_count;
  
  
  // index list for duarray, yarray
  // only include left half
  // NOTE: idxdu indicates lefthalf only
  //       idxu indicates both halves

  // index list for beta and B

#if defined(openmp_version)
#pragma omp target teams distribute parallel for private(idxb_count) nowait
#endif
  for(int pxc = 0; pxc < j1j2tot; pxc++){
      int j1 = int((-1 + sqrt(1+8*pxc))/2);
      int j2 = pxc - j1*(j1+1)/2;
      idxb_count = idxb_trd(j1,j2);
      for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {
        if (j < j1) continue;
        idxb_block(j1,j2,j) = idxb_count; 
        idxb_count++;
      }
  }
#pragma omp target exit data map(from: this->idxb_trd.dptr[0:this->idxb_trd.size]) nowait
  
  // index list for zlist

#if defined(openmp_version)
#pragma omp target teams distribute parallel for private(idxz_count) nowait
#endif
  for(int pxc = 0; pxc < j1j2tot; pxc++){
      int j1 = int((-1 + sqrt(1+8*pxc))/2);
      int j2 = pxc - j1*(j1+1)/2;
      idxz_count = idxz_trd(j1,j2);
      for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {
        // find right beta[jjb] entries
        // multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3
        // this should not be computed here
        SNADOUBLE betaj; 
        if (j >= j1) {
          const int jjb = idxb_block(j1,j2,j);
          if (j1 == j) {
            if (j2 == j) betaj = 3*beta(jjb+1);
            else betaj = 2*beta(jjb+1);
          } else betaj = beta(jjb+1); 
        } else if (j >= j2) {
          const int jjb = idxb_block(j,j2,j1);
          if (j2 == j) betaj = 2*beta(jjb+1)*(j1+1)/(j+1.0);
          else betaj = beta(jjb+1)*(j1+1)/(j+1.0);
        } else {
          const int jjb = idxb_block(j2,j,j1);
          betaj = beta(jjb+1)*(j1+1)/(j+1.0); 
        }
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            idxzbeta(idxz_count) = betaj;
            idxz_count++;
	  }
      }
  }



#if defined(openmp_version)
#pragma omp target teams distribute parallel for private(idxz_count) nowait
#endif
  for(int pxc = 0; pxc < j1j2tot; pxc++){
      int j1 = int((-1 + sqrt(1+8*pxc))/2);
      int j2 = pxc - j1*(j1+1)/2;
      idxz_count = idxz_trd(j1,j2);
      for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            idxz(idxz_count,0) = j1;
            idxz(idxz_count,1) = j2;
            idxz(idxz_count,2) = j;

            int ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
            idxz(idxz_count,3) = ma1min;
            idxz(idxz_count,4) = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
            idxz(idxz_count,5) = MIN(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;

            int mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
            idxz(idxz_count,6) = mb1min;
            idxz(idxz_count,7) = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
            idxz(idxz_count,8) = MIN(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;

            idxz_count++;
          }
      }
  }
#pragma omp target exit data map(from: this->idxz_trd.dptr[0:this->idxz_trd.size]) nowait

}

/* ---------------------------------------------------------------------- */

void SNA::init()
{
#if defined(openmp_version)
  omp_offload_init();
#endif
  init_clebsch_gordan();
  init_rootpqarray();
  build_indexlist();
}


/* ---------------------------------------------------------------------- */


void SNA::grow_rij(int newnmax)
{
  if(newnmax <= nmax) return;

  nmax = newnmax;

  rij.resize(num_atoms,num_nbors,3);
  inside.resize(num_atoms,num_nbors);
  wj.resize(num_atoms,num_nbors);
  rcutij.resize(num_atoms,num_nbors);
  dedr.resize(num_atoms,num_nbors,3);
}

/* ----------------------------------------------------------------------
   compute Ui by summing over neighbors j
------------------------------------------------------------------------- */

void SNA::compute_ui()
{
  zero_uarraytot();
  addself_uarraytot(wself);

#if defined(openmp_version)
#pragma omp target teams distribute parallel for collapse(2)
#endif
  for (int natom = 0; natom < num_atoms; natom++) {
    for(int nbor = 0; nbor < num_nbors; nbor++) {
      SNADOUBLE  x = rij(natom,nbor,0);
      SNADOUBLE  y = rij(natom,nbor,1);
      SNADOUBLE  z = rij(natom,nbor,2);
      SNADOUBLE  rsq = x * x + y * y + z * z;
      SNADOUBLE  r = sqrt(rsq);
  
      SNADOUBLE  theta0 = (r - rmin0) * rfac0 * MY_PI / (rcutij(natom,nbor) - rmin0);
      SNADOUBLE  z0 = r / tan(theta0);
  
      SNADOUBLE  r0inv = 1.0 / sqrt(r * r + z0 * z0);
      SNADOUBLE  a_r = r0inv * z0;
      SNADOUBLE  a_i = -r0inv * z;
      SNADOUBLE  b_r = r0inv * y;
      SNADOUBLE  b_i = -r0inv * x;
      
      compute_uarray(natom, nbor, a_r, a_i, b_r, b_i);
      
      SNADOUBLE  sfac = compute_sfac(r, rcutij(natom,nbor));
      sfac *= wj(natom,nbor);

      add_uarraytot(natom,nbor,sfac);
    }
  }
}

/* ----------------------------------------------------------------------
   compute Yi from Ui without storing Zi, looping over zlist
------------------------------------------------------------------------- */

void SNA::compute_yi()
{
#if defined(openmp_version)
#pragma omp target teams distribute parallel for collapse(2)
#endif
  for (int natom = 0; natom < num_atoms; natom++) {
    for(int jjdu = 0; jjdu < idxdu_max; jjdu++) {
      ylist(natom,jjdu) = {0.0, 0.0};
    }
  }

#if defined(openmp_version)
#pragma omp target teams distribute parallel for collapse(2)
#endif
  for (int natom = 0; natom < num_atoms; natom++) {
    for(int jjz = 0; jjz < idxz_max; jjz++) {
      const int j1 = idxz(jjz,0);
      const int j2 = idxz(jjz,1);
      const int j = idxz(jjz,2);
      const int ma1min = idxz(jjz,3);
      const int ma2max = idxz(jjz,4);
      const int na = idxz(jjz,5);
      const int mb1min = idxz(jjz,6);
      const int mb2max = idxz(jjz,7);
      const int nb = idxz(jjz,8);
  
      const SNADOUBLE betaj = idxzbeta(jjz);
  
      int mb = (2 * (mb1min+mb2max) - j1 - j2 + j) / 2;
      int ma = (2 * (ma1min+ma2max) - j1 - j2 + j) / 2;
      const int jjdu = idxdu_block(j) + (j+1)*mb + ma;
  
      int jju1 = (j1*(j1+1)*(2*j1+1)/6) + (j1+1)*mb1min;
      int jju2 = (j2*(j2+1)*(2*j2+1)/6) + (j2+1)*mb2max;
      int icgb = mb1min*(j2+1) + mb2max + idxcg_block(j1,j2,j);
      int icga = ma1min*(j2+1) + ma2max + idxcg_block(j1,j2,j);
  
      SNADOUBLE ztmp_r = 0.0;
      SNADOUBLE ztmp_i = 0.0;
  
      // loop over columns of u1 and corresponding
      // columns of u2 satisfying Clebsch-Gordan constraint 
      //      2*mb-j = 2*mb1-j1 + 2*mb2-j2
  
      for(int ib = 0; ib < nb; ib++) {
  
        SNADOUBLE suma1_r = 0.0;
        SNADOUBLE suma1_i = 0.0;
            
        // loop over elements of row u1[mb1] and corresponding elements 
        // of row u2[mb2] satisfying Clebsch-Gordan constraint 
        //      2*ma-j = 2*ma1-j1 + 2*ma2-j
        for(int ia = 0; ia < na; ia++)
        {
          int ma1 = ma1min+ia;
          int ma2 = ma2max-ia;
          int inca = icga + j2*ia;
          suma1_r += cglist(inca) *
            (ulisttot(natom,jju1+ma1).re * ulisttot(natom,jju2+ma2).re
              - ulisttot(natom,jju1+ma1).im * ulisttot(natom,jju2+ma2).im);
          suma1_i += cglist(inca) *
            (ulisttot(natom,jju1+ma1).re * ulisttot(natom,jju2+ma2).im
             + ulisttot(natom,jju1+ma1).im * ulisttot(natom,jju2+ma2).re);
        }
  
  
        ztmp_r += cglist(icgb) * suma1_r;
        ztmp_i += cglist(icgb) * suma1_i;
        jju1 += j1+1;
        jju2 -= j2+1;
        icgb += j2;
      } // end loop over ib
  
      // apply z(j1,j2,j,ma,mb) to unique element of y(j)
  
#if defined(openmp_version)
#pragma omp atomic
#endif
      ylist(natom,jjdu).re += betaj*ztmp_r;
#if defined(openmp_version)
#pragma omp atomic
#endif
      ylist(natom,jjdu).im += betaj*ztmp_i;
  
    } // end jjz loop
  }
}



/* ----------------------------------------------------------------------
   calculate derivative of Ui w.r.t. atom j
------------------------------------------------------------------------- */

void SNA::compute_duidrj()
{
  compute_duarray();
}


/* ----------------------------------------------------------------------
   compute dEidRj
------------------------------------------------------------------------- */

void SNA::compute_deidrj()
{
#if defined(openmp_version)
#pragma omp target teams distribute parallel for collapse(2)
#endif
  for(int natom = 0; natom < num_atoms; natom++) {
    for (int nbor = 0; nbor < num_nbors; nbor++) {
      for(int k = 0; k < 3; k++)
        dedr(natom,nbor,k) = 0.0;
      for(int j = 0; j <= twojmax; j++) {
        int jjdu = idxdu_block(j);
    
        for(int mb = 0; 2*mb < j; mb++)
          for(int ma = 0; ma <= j; ma++) {
            for(int k = 0; k < 3; k++)
              dedr(natom,nbor,k) +=
                dulist(natom,nbor,jjdu,k).re * ylist(natom,jjdu).re +
                dulist(natom,nbor,jjdu,k).im * ylist(natom,jjdu).im;
            jjdu++;
          } //end loop over ma mb
    
        // For j even, handle middle column
    
        if (j%2 == 0) {
          int mb = j/2;
          for(int ma = 0; ma < mb; ma++) {
            for(int k = 0; k < 3; k++)
              dedr(natom,nbor,k) +=
                dulist(natom,nbor,jjdu,k).re * ylist(natom,jjdu).re +
                dulist(natom,nbor,jjdu,k).im * ylist(natom,jjdu).im;
            jjdu++;
          }
    
          for(int k = 0; k < 3; k++)
            dedr(natom,nbor,k) += 
              (dulist(natom,nbor,jjdu,k).re * ylist(natom,jjdu).re +
               dulist(natom,nbor,jjdu,k).im * ylist(natom,jjdu).im)*0.5;
          jjdu++;
        } // end if jeven
      } // end loop over j
    }
  }

#if defined(openmp_version)
#pragma omp target teams distribute parallel for collapse(3)
#endif
  for(int natom = 0; natom < num_atoms; natom++)
    for (int nbor = 0; nbor < num_nbors; nbor++)
      for(int k = 0; k < 3; k++)
        dedr(natom,nbor,k) *= 2.0;

#pragma omp target update from(this->dedr.dptr[0:this->dedr.size])
}



/* ---------------------------------------------------------------------- */

void SNA::zero_uarraytot()
{
#if defined(openmp_version)
#pragma omp target teams distribute parallel for collapse(2)
#endif
  for (int natom = 0; natom < num_atoms; natom++)
    for (int jju = 0; jju < idxu_max; jju++)
        ulisttot(natom,jju) = {0.0,0.0};
}

/* ---------------------------------------------------------------------- */

void SNA::addself_uarraytot(SNADOUBLE wself_in)
{
#if defined(openmp_version)
#pragma omp target teams distribute parallel for
#endif
  for (int natom = 0; natom < num_atoms; natom++) {
    for (int j = 0; j <= twojmax; j++) {
      int jju = j*(j+1)*(2*j+1)/6;
      for (int ma = 0; ma <= j; ma++) {
        ulisttot(natom,jju) = {wself_in, 0.0};
        jju += j+2;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   add Wigner U-functions for one neighbor to the total
------------------------------------------------------------------------- */

void SNA::add_uarraytot(int natom, int nbor, SNADOUBLE sfac)
{

  for (int j = 0; j <= twojmax; j++) {
    int jju = j*(j+1)*(2*j+1)/6;
    for (int mb = 0; mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
#if defined(openmp_version)
#pragma omp atomic
#endif
        ulisttot(natom,jju).re += sfac * ulist(natom,nbor,jju).re;
#if defined(openmp_version)
#pragma omp atomic
#endif
        ulisttot(natom,jju).im += sfac * ulist(natom,nbor,jju).im;
        jju++;
      }
  }
}

/* ----------------------------------------------------------------------
   compute Wigner U-functions for one neighbor
------------------------------------------------------------------------- */

void SNA::compute_uarray(int natom, int nbor, SNADOUBLE a_r, SNADOUBLE a_i, SNADOUBLE b_r,
                         SNADOUBLE b_i)
{
  SNADOUBLE rootpq;

  // compute Cayley-Klein parameters for unit quaternion

  // VMK Section 4.8.2
  ulist(natom, nbor, 0) = {1.0,0.0};
 
  for (int j = 1; j <= twojmax; j++) {
    int maxxb1 = int((j+1)/2)*j;
    for (int xb = 0; xb < maxxb1; xb++){
      int ma = xb%j;
      int mb = int(xb/j);
      int jju = (j*(j+1)*(2*j+1)/6)+(xb+mb);
      int jjup = ((j-1)*j*(2*j-1)/6)+xb;
      
      if(ma == 0)ulist(natom,nbor,jju) = {0.0,0.0};

      SNADOUBLE up_r = ulist(natom,nbor,jjup).re;
      SNADOUBLE up_i = ulist(natom,nbor,jjup).im;

      rootpq = rootpqarray(j - ma,j - mb);
      ulist(natom,nbor,jju).re += rootpq * (a_r * up_r + a_i * up_i);
      ulist(natom,nbor,jju).im += rootpq * (a_r * up_i - a_i * up_r);

      rootpq = rootpqarray(ma + 1,j - mb);
      ulist(natom,nbor,jju+1).re = -rootpq * (b_r * up_r + b_i * up_i);
      ulist(natom,nbor,jju+1).im = -rootpq * (b_r * up_i - b_i * up_r);
    }


    if (j%2 == 0) {
      int mb = j/2;
      int jju = (j*(j+1)*(2*j+1)/6)+maxxb1+int(j/2);
      int jjup = ((j-1)*j*(2*j-1)/6)+maxxb1-1;

      ulist(natom,nbor,jju) = {0.0,0.0};
      for (int ma = 0; ma < j; ma++) {
	SNADOUBLE ud_r = ulist(natom,nbor,jjup-ma).re;
	SNADOUBLE ud_i = ulist(natom,nbor,jjup-ma).im;
        
	rootpq = ulist_parity(jjup-ma)*rootpqarray(j - ma,j - mb);
        ulist(natom,nbor,jju+ma).re += rootpq * (a_r * ud_r + a_i * -ud_i);
        ulist(natom,nbor,jju+ma).im += rootpq * (a_r * -ud_i - a_i * ud_r);

        rootpq = ulist_parity(jjup-ma)*rootpqarray(ma + 1,j - mb);
        ulist(natom,nbor,jju+1+ma).re = -rootpq * (b_r * ud_r + b_i * -ud_i);
        ulist(natom,nbor,jju+1+ma).im = -rootpq * (b_r * -ud_i - b_i * ud_r);
      }
    }

    int jju = j*(j+1)*(2*j+1)/6;
    int jjup = jju+(j+1)*(j+1)-1;
    for (int mb = 0; 2*mb < j; mb++) {
      for (int ma = 0; ma <= j; ma++) {
        ulist(natom,nbor,jjup) = {ulist_parity(jju)*ulist(natom,nbor,jju).re,ulist_parity(jju)*-ulist(natom,nbor,jju).im};
        jju++;
        jjup--;
      }
    }
  }
}


/* ----------------------------------------------------------------------
   Compute derivatives of Wigner U-functions for one neighbor
   see comments in compute_uarray()
------------------------------------------------------------------------- */

void SNA::compute_duarray()
{
#if defined(openmp_version)
#pragma omp target teams distribute parallel for collapse(2)
#endif
  for (int natom = 0; natom < num_atoms; natom++) {
    for (int nbor = 0; nbor < num_nbors; nbor++) {
      SNADOUBLE rcut = rcutij(natom,nbor);
      SNADOUBLE wxj = wj(natom,nbor);
    
      SNADOUBLE x = rij(natom,nbor,0);
      SNADOUBLE y = rij(natom,nbor,1);
      SNADOUBLE z = rij(natom,nbor,2);
      SNADOUBLE rsq = x * x + y * y + z * z;
      SNADOUBLE r = sqrt(rsq);
      SNADOUBLE rscale0 = rfac0 * MY_PI / (rcutij(natom,nbor) - rmin0);
      SNADOUBLE theta0 = (r - rmin0) * rscale0;
      SNADOUBLE cs = cos(theta0);
      SNADOUBLE sn = sin(theta0);
      SNADOUBLE z0 = r * cs / sn;
      SNADOUBLE dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;
    
      SNADOUBLE r0inv;
      SNADOUBLE a_r, a_i, b_r, b_i;
      SNADOUBLE da_r[3], da_i[3], db_r[3], db_i[3];
      SNADOUBLE dz0[3], dr0inv[3], dr0invdr;
      SNADOUBLE rootpq;
    
      SNADOUBLE rinv = 1.0 / r;
      SNADOUBLE ux = x * rinv;
      SNADOUBLE uy = y * rinv;
      SNADOUBLE uz = z * rinv;
      
      int jju,jjdu,jjup,jjdup;

      r0inv = 1.0 / sqrt(r * r + z0 * z0);
      a_r = z0 * r0inv;
      a_i = -z * r0inv;
      b_r = y * r0inv;
      b_i = -x * r0inv;
    
      dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);
    
      dr0inv[0] = dr0invdr * ux;
      dr0inv[1] = dr0invdr * uy;
      dr0inv[2] = dr0invdr * uz;
    
      dz0[0] = dz0dr * ux;
      dz0[1] = dz0dr * uy;
      dz0[2] = dz0dr * uz;
    
      for (int k = 0; k < 3; k++) {
        da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
        da_i[k] = -z * dr0inv[k];
      }
      da_i[2] += -r0inv;
      for (int k = 0; k < 3; k++) {
        db_r[k] = y * dr0inv[k];
        db_i[k] = -x * dr0inv[k];
      }
      db_i[0] += -r0inv;
      db_r[1] += r0inv;
    
      ulist(natom,nbor,0) = {1.0,0.0};
      for (int k = 0; k < 3; k++)
    	  dulist(natom,nbor,0,k) = {0.0,0.0};
     

      for (int j = 1; j <= twojmax; j++) {
    	int maxxb1 = int((j+1)/2)*j;
        for (int xb = 0; xb < maxxb1; xb++){
          int ma = xb%j;
          int mb = int(xb/j);
          jju = (j*(j+1)*(2*j+1)/6)+(xb+mb);
          jjdu = idxdu_block(j)+(xb+mb);
          jjup = ((j-1)*j*(2*j-1)/6)+xb;
          jjdup = idxdu_block(j-1)+xb;
          
          if(ma == 0){
            ulist(natom,nbor,jju) = {0.0,0.0};
            for(int k = 0; k < 3; k++)
              dulist(natom,nbor,jjdu,k) = {0.0,0.0};
          }
           
          rootpq = rootpqarray(j - ma,j - mb);
          ulist(natom,nbor,jju).re += rootpq *
                                 (a_r *  ulist(natom,nbor,jjup).re +
                                  a_i *  ulist(natom,nbor,jjup).im);
          ulist(natom,nbor,jju).im += rootpq *
                                 (a_r *  ulist(natom,nbor,jjup).im -
                                  a_i *  ulist(natom,nbor,jjup).re);
            
          for (int k = 0; k < 3; k++) {
            dulist(natom,nbor,jjdu,k).re +=
              rootpq * (da_r[k] * ulist(natom,nbor,jjup).re +
                        da_i[k] * ulist(natom,nbor,jjup).im +
                        a_r * dulist(natom,nbor,jjdup,k).re +
                        a_i * dulist(natom,nbor,jjdup,k).im);
            dulist(natom,nbor,jjdu,k).im +=
              rootpq * (da_r[k] * ulist(natom,nbor,jjup).im -
                        da_i[k] * ulist(natom,nbor,jjup).re +
                        a_r * dulist(natom,nbor,jjdup,k).im -
        		a_i * dulist(natom,nbor,jjdup,k).re);
          }
        
        
          rootpq = rootpqarray(ma + 1,j - mb);
          ulist(natom,nbor,jju+1).re =
            -rootpq * (b_r *  ulist(natom,nbor,jjup).re +
                       b_i *  ulist(natom,nbor,jjup).im);
          ulist(natom,nbor,jju+1).im =
            -rootpq * (b_r *  ulist(natom,nbor,jjup).im -
                       b_i *  ulist(natom,nbor,jjup).re);
        
          for (int k = 0; k < 3; k++) {
            dulist(natom,nbor,jjdu+1,k).re =
              -rootpq * (db_r[k] * ulist(natom,nbor,jjup).re +
                         db_i[k] * ulist(natom,nbor,jjup).im +
                         b_r * dulist(natom,nbor,jjdup,k).re +
                         b_i * dulist(natom,nbor,jjdup,k).im);
            dulist(natom,nbor,jjdu+1,k).im =
              -rootpq * (db_r[k] * ulist(natom,nbor,jjup).im -
                         db_i[k] * ulist(natom,nbor,jjup).re +
                         b_r * dulist(natom,nbor,jjdup,k).im -
                         b_i * dulist(natom,nbor,jjdup,k).re);
          }
        }


        // handle middle column using inversion symmetry of previous layer
        if (j%2 == 0) {
          int mb = j/2;
	  jjup = ((j-1)*j*(2*j-1)/6)+maxxb1-1;
	  jjdup = idxdu_block(j-1)+maxxb1-1;
          jju = (j*(j+1)*(2*j+1)/6)+maxxb1+int(j/2);
          jjdu = idxdu_block(j)+maxxb1+int(j/2);

	  ulist(natom,nbor,jju) = {0.0,0.0};
	  for (int k = 0; k < 3; k++)
    	      dulist(natom,nbor,jjdu,k) = {0.0,0.0};

	  for (int ma = 0; ma < j; ma++) {
            rootpq = rootpqparityarray(j - ma,j - mb);
            ulist(natom,nbor,jju+ma).re += rootpq *
                                   (a_r *  ulist(natom,nbor,jjup-ma).re +
                                    a_i *  -ulist(natom,nbor,jjup-ma).im);
            ulist(natom,nbor,jju+ma).im += rootpq *
                                   (a_r *  -ulist(natom,nbor,jjup-ma).im -
                                    a_i *  ulist(natom,nbor,jjup-ma).re);
    
            for (int k = 0; k < 3; k++) {
              dulist(natom,nbor,jjdu+ma,k).re +=
                rootpq * (da_r[k] * ulist(natom,nbor,jjup-ma).re +
                          da_i[k] * -ulist(natom,nbor,jjup-ma).im +
                          a_r * dulist(natom,nbor,jjdup-ma,k).re +
                          a_i * -dulist(natom,nbor,jjdup-ma,k).im);
              dulist(natom,nbor,jjdu+ma,k).im +=
                rootpq * (da_r[k] * -ulist(natom,nbor,jjup-ma).im -
                          da_i[k] * ulist(natom,nbor,jjup-ma).re +
                          a_r * -dulist(natom,nbor,jjdup-ma,k).im -
                          a_i * dulist(natom,nbor,jjdup-ma,k).re);
            }
    
            rootpq = -rootpqparityarray(ma + 1,j - mb);
            ulist(natom,nbor,jju+1+ma).re =
              -rootpq * (b_r *  ulist(natom,nbor,jjup-ma).re +
                         b_i *  -ulist(natom,nbor,jjup-ma).im);
            ulist(natom,nbor,jju+1+ma).im =
              -rootpq * (b_r *  -ulist(natom,nbor,jjup-ma).im -
                         b_i *  ulist(natom,nbor,jjup-ma).re);
    
            for (int k = 0; k < 3; k++) {
              dulist(natom,nbor,jjdu+1+ma,k).re =
                -rootpq * (db_r[k] * ulist(natom,nbor,jjup-ma).re +
                           db_i[k] * -ulist(natom,nbor,jjup-ma).im +
                           b_r * dulist(natom,nbor,jjdup-ma,k).re +
                           b_i * -dulist(natom,nbor,jjdup-ma,k).im);
              dulist(natom,nbor,jjdu+1+ma,k).im =
                -rootpq * (db_r[k] * -ulist(natom,nbor,jjup-ma).im -
                           db_i[k] * ulist(natom,nbor,jjup-ma).re +
                           b_r * -dulist(natom,nbor,jjdup-ma,k).im -
                           b_i * dulist(natom,nbor,jjdup-ma,k).re);
            }
          }
        }
      }

      
      SNADOUBLE sfac = compute_sfac(r, rcut);
      SNADOUBLE dsfac = compute_dsfac(r, rcut);
    
      sfac *= wxj;
      dsfac *= wxj;
      for (int j = 0; j <= twojmax; j++) {
        int jju = j*(j+1)*(2*j+1)/6;
        int jjdu = idxdu_block(j);
	int maxxb3 = (int(j/2)+1)*(j+1);
        for (int xb = 0; xb < maxxb3; xb++){
            dulist(natom,nbor,jjdu+xb,0).re = dsfac * ulist(natom,nbor,jju+xb).re * ux +
                                      sfac * dulist(natom,nbor,jjdu+xb,0).re;
            dulist(natom,nbor,jjdu+xb,0).im = dsfac * ulist(natom,nbor,jju+xb).im * ux +
                                      sfac * dulist(natom,nbor,jjdu+xb,0).im;
            dulist(natom,nbor,jjdu+xb,1).re = dsfac * ulist(natom,nbor,jju+xb).re * uy +
                                      sfac * dulist(natom,nbor,jjdu+xb,1).re;
            dulist(natom,nbor,jjdu+xb,1).im = dsfac * ulist(natom,nbor,jju+xb).im * uy +
                                      sfac * dulist(natom,nbor,jjdu+xb,1).im;
            dulist(natom,nbor,jjdu+xb,2).re = dsfac * ulist(natom,nbor,jju+xb).re * uz +
                                      sfac * dulist(natom,nbor,jjdu+xb,2).re;
            dulist(natom,nbor,jjdu+xb,2).im = dsfac * ulist(natom,nbor,jju+xb).im * uz +
                                      sfac * dulist(natom,nbor,jjdu+xb,2).im;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of arrays
------------------------------------------------------------------------- */

SNADOUBLE SNA::memory_usage()
{
  int jdimpq = twojmax + 2;
  int jdim = twojmax + 1;
  SNADOUBLE bytes;
  bytes = ncoeff * sizeof(SNADOUBLE);                       // coeff

  bytes += jdimpq*jdimpq * sizeof(SNADOUBLE);               // pqarray
  bytes += idxcg_max * sizeof(SNADOUBLE);                   // cglist
  bytes += jdim * jdim * jdim * sizeof(int);                // idxcg_block

  bytes += idxu_max * sizeof(SNADOUBLE) * 2;                // ulist
  bytes += idxu_max * sizeof(SNADOUBLE) * 2;                // ulisttot
  bytes += idxdu_max * 3 * sizeof(SNADOUBLE) * 2;           // dulist
  bytes += idxu_max * sizeof(int);                          // ulist_parity
  bytes += jdim * sizeof(int);                              // idxu_block
  bytes += jdim * sizeof(int);                              // idxdu_block

  bytes += idxz_max * 9 * sizeof(int);                      // idxz
  bytes += idxz_max * sizeof(SNADOUBLE);                    // idxzbeta
  bytes += jdim * jdim * jdim * sizeof(int);                // idxz_block

  bytes += idxdu_max * sizeof(SNADOUBLE) * 2;               // ylist
  bytes += idxb_max * 3 * sizeof(int);                      // idxb

  bytes += jdim * jdim * jdim * sizeof(int);                // idxb_block

  return bytes;
}

/* ---------------------------------------------------------------------- */

void SNA::create_twojmax_arrays()
{
  int jdimpq = twojmax + 2;
  int jdim = twojmax + 1;
  int idxb_count, idxz_count;
  int m, aa2, bb2; 

  idxcg_block.resize(jdim,jdim,jdim);
  idxdu_block.resize(jdim);
  
  idxb_trd.resize(jdim,jdim);
  idxc_trd.resize(jdim,jdim);
  idxz_trd.resize(jdim,jdim);
  idxb_block.resize(jdim,jdim,jdim);
  
  fact.resize(40);
  for(int i=0;i<40;i++)
    fact(i)=factorial(i);
  
  idxb_count = 0;  
  idxz_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++) {
      idxb_trd(j1,j2) = idxb_count;
      idxz_trd(j1,j2) = idxz_count;
      for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {
        if (j >= j1) idxb_count++;
        for (int mb = 0; 2*mb <= j; mb++) 
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;
      }
    }
  idxb_max = idxb_count;
  idxz_max = idxz_count;
  
  
  
  int idxcg_count = 0;
  int idycg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++){
      idxc_trd(j1,j2) = idycg_count;
      for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {
        idxcg_block(j1,j2,j) = idxcg_count; 
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++){
            idxcg_count++;
            aa2 = 2 * m1 - j1;
            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;
            if(m < 0 || m > j) {
              idycg_count++;
              continue;
            }
            idycg_count++;
	  }
      }
    }
  idxcg_max = idxcg_count;
  idxu_max = jdim*(jdim+1)*(2*jdim+1)/6;
  
     
  
  int idxdu_count = 0;
  for(int j = 0; j <= twojmax; j++) {
    idxdu_block(j) = idxdu_count;
    for(int mb = 0; 2*mb <= j; mb++) 
      for(int ma = 0; ma <= j; ma++)
	 idxdu_count++; 
  }
  idxdu_max = idxdu_count;

  
  // parity list for uarray inversion symmetry
  // parity +1: u[ma-j][mb-j] = +Conj([u[ma][mb])
  // parity -1: u[ma-j][mb-j] = -Conj([u[ma][mb])
  ulist_parity.resize(idxu_max);
  int idxu_count = 0;
  for(int j = 0; j <= twojmax; j++) {
    int mbpar = 1;
    for(int mb = 0; mb <= j; mb++) {
      int mapar = mbpar;
      for(int ma = 0; ma <= j; ma++) {
        ulist_parity(idxu_count) = mapar;
        mapar = -mapar;
        idxu_count++;
      }
      mbpar = -mbpar;
    }
  }


  cglist.resize(idxcg_max);
  rootpqarray.resize(jdimpq,jdimpq);
  rootpqparityarray.resize(jdimpq,jdimpq);
  ulist.resize(num_atoms,num_nbors,idxu_max);
  ulisttot.resize(num_atoms,idxu_max);
  dulist.resize(num_atoms,num_nbors,idxdu_max,3);
  ylist.resize(num_atoms,idxdu_max);
  

  idxz.resize(idxz_max,9);
  idxzbeta.resize(idxz_max);
}

/* ---------------------------------------------------------------------- */

void SNA::destroy_twojmax_arrays()
{
}

/* ----------------------------------------------------------------------
   factorial n, wrapper for precomputed table
------------------------------------------------------------------------- */

SNADOUBLE SNA::factorial(int n)
{
  if (n < 0 || n > nmaxfactorial) {
    char str[128];
    //printf("Invalid argument to factorial %d", n);
    exit(1);
  }

  return nfac_table[n];
}

/* ----------------------------------------------------------------------
   factorial n table, size SNA::nmaxfactorial+1
------------------------------------------------------------------------- */

const SNADOUBLE SNA::nfac_table[] = {
  1,
  1,
  2,
  6,
  24,
  120,
  720,
  5040,
  40320,
  362880,
  3628800,
  39916800,
  479001600,
  6227020800,
  87178291200,
  1307674368000,
  20922789888000,
  355687428096000,
  6.402373705728e+15,
  1.21645100408832e+17,
  2.43290200817664e+18,
  5.10909421717094e+19,
  1.12400072777761e+21,
  2.5852016738885e+22,
  6.20448401733239e+23,
  1.5511210043331e+25,
  4.03291461126606e+26,
  1.08888694504184e+28,
  3.04888344611714e+29,
  8.8417619937397e+30,
  2.65252859812191e+32,
  8.22283865417792e+33,
  2.63130836933694e+35,
  8.68331761881189e+36,
  2.95232799039604e+38,
  1.03331479663861e+40,
  3.71993326789901e+41,
  1.37637530912263e+43,
  5.23022617466601e+44,
  2.03978820811974e+46,  // nmaxfactorial = 39
};

/* ----------------------------------------------------------------------
   the function delta given by VMK Eq. 8.2(1)
------------------------------------------------------------------------- */

SNADOUBLE SNA::deltacg(int j1, int j2, int j)
{
  SNADOUBLE sfaccg = fact((j1 + j2 + j) / 2 + 1);
  return sqrt(fact((j1 + j2 - j) / 2) *
              fact((j1 - j2 + j) / 2) *
              fact((-j1 + j2 + j) / 2) / sfaccg);
}

/* ----------------------------------------------------------------------
   assign Clebsch-Gordan coefficients using
   the quasi-binomial formula VMK 8.2.1(3)
------------------------------------------------------------------------- */

void SNA::init_clebsch_gordan()
{
  SNADOUBLE sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac,idxcg_count; 
  int jdim = twojmax+1;
  int j1j2tot = jdim*(jdim+1)/2;
  
#if defined(openmp_version)
#pragma omp target teams distribute parallel for private(idxcg_count,sum)
#endif
  for(int pxc = 0; pxc < j1j2tot; pxc++){
      int j1 = int((-1 + sqrt(1+8*pxc))/2);
      int j2 = pxc - j1*(j1+1)/2;
      idxcg_count = idxc_trd(j1,j2);
      for(int j = abs(j1 - j2); j <= MIN(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;
          for (int m2 = 0; m2 <= j2; m2++) {

            // -c <= cc <= c
            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;
            if(m < 0 || m > j) {
              cglist(idxcg_count) = 0.0;
              idxcg_count++;
              continue;
            }

            sum = 0.0;
            for (int z = MAX(0, MAX(-(j - j2 + aa2)
                                    / 2, -(j - j1 - bb2) / 2));
                 z <= MIN((j1 + j2 - j) / 2,
                          MIN((j1 - aa2) / 2, (j2 + bb2) / 2));
                 z++) {
              ifac = z % 2 ? -1 : 1;
              sum += ifac /
                (fact(z) *
                 fact((j1 + j2 - j) / 2 - z) *
                 fact((j1 - aa2) / 2 - z) *
                 fact((j2 + bb2) / 2 - z) *
                 fact((j - j2 + aa2) / 2 + z) *
                 fact((j - j1 - bb2) / 2 + z));
            }
           
            cc2 = 2 * m - j;
            dcg = deltacg(j1, j2, j);
            sfaccg = sqrt(fact((j1 + aa2) / 2) *
                          fact((j1 - aa2) / 2) *
                          fact((j2 + bb2) / 2) *
                          fact((j2 - bb2) / 2) *
                          fact((j  + cc2) / 2) *
                          fact((j  - cc2) / 2) *
                          (j + 1));
            
            cglist(idxcg_count) = sum * dcg * sfaccg;
            idxcg_count++;
          }
        }
      }
  }
#pragma omp target exit data map(from:this->fact.dptr[0:this->fact.size], \
	       	this->idxc_trd.dptr[0:this->idxc_trd.size]) nowait
}

/* ----------------------------------------------------------------------
   pre-compute table of sqrt[p/m2], p, q = 1,twojmax
   the p = 0, q = 0 entries are allocated and skipped for convenience.
   a second table is computed with +/-1 parity factor
------------------------------------------------------------------------- */

void SNA::init_rootpqarray()
{
#if defined(openmp_version)
#pragma omp target teams distribute parallel for collapse(2)
#endif
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++){
      rootpqarray(p,q) = sqrt(static_cast<SNADOUBLE>(p)/q);
      rootpqparityarray(p,q) = pow(-1,(p+q))*rootpqarray(p,q);
    }
}

/* ---------------------------------------------------------------------- */

int SNA::compute_ncoeff()
{
  int ncount;

  ncount = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = abs(j1 - j2);
           j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) ncount++;

  return ncount;
}

/* ---------------------------------------------------------------------- */

SNADOUBLE SNA::compute_sfac(SNADOUBLE r, SNADOUBLE rcut)
{ 
  if (switch_flag == 0) return 1.0;
  if (switch_flag == 1) {
    if(r <= rmin0) return 1.0;
    else if(r > rcut) return 0.0;
    else {
      SNADOUBLE rcutfac = MY_PI / (rcut - rmin0);
      return 0.5 * (cos((r - rmin0) * rcutfac) + 1.0);
    }
  }
  return 0.0;
}

/* ---------------------------------------------------------------------- */

SNADOUBLE SNA::compute_dsfac(SNADOUBLE r, SNADOUBLE rcut)
{
  if (switch_flag == 0) return 0.0;
  if (switch_flag == 1) {
    if(r <= rmin0) return 0.0;
    else if(r > rcut) return 0.0;
    else {
      SNADOUBLE rcutfac = MY_PI / (rcut - rmin0);
      return -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
    }
  }
  return 0.0;
}

