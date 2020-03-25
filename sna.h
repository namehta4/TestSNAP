/* -*- c++ -*- -------------------------------------------------------------
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

#ifndef LMP_SNA_H
#define LMP_SNA_H

#include "arrayMDgpu.h"

typedef double SNADOUBLE;
typedef float SNAFLOAT;

using SNAcomplex = struct{SNADOUBLE re, im;};

struct SNA_BINDICES {
  int j1, j2, j;
};

class SNA {

public:
  SNA(int, int, class Memory*, SNADOUBLE, int, SNADOUBLE, int, int);
  ~SNA();
  void build_indexlist();
  void init();
  SNADOUBLE memory_usage();

  int ncoeff;
  int num_atoms, num_nbors;
 

  void omp_offload_init();
  void omp_offload_update();

  // functions for bispectrum coefficients
  void compute_ui();
  void compute_yi();

  // functions for derivatives
  void compute_duidrj();
  void compute_duarray();
  void compute_deidrj();
  Array3D<SNADOUBLE> dedr;

  // per sna class instance for OMP use
  Array3D<SNADOUBLE> rij;
  Array2D<int> inside;
  Array2D<SNADOUBLE> wj;
  Array2D<SNADOUBLE> rcutij;
  int nmax;

  void grow_rij(int);
  Array1D<SNA_BINDICES> idxb;

  Array1D<double> beta; 
//  double* beta;
private:
  Memory* memory;
  SNADOUBLE rmin0, rfac0;

  // use indexlist instead of loops, constructor generates these
  Array2D<int> idxz;
  Array1D<SNADOUBLE> idxzbeta;
  int idxcg_max, idxu_max, idxdu_max, idxz_max, idxb_max;

  // data for bispectrum coefficients
  int twojmax;
  Array3D<int> idxcg_block;
  Array2D<SNADOUBLE> rootpqarray;
  Array2D<SNADOUBLE> rootpqparityarray;
  Array1D<SNADOUBLE> cglist;

  Array2D<SNAcomplex> ulisttot;
  Array3D<SNAcomplex> ulist;
  Array1D<int> idxu_block;
  Array1D<int> ulist_parity;
  Array1D<int> idxdu_block;
  
  Array3D<int> idxb_block;

  Array2D<int> idxb_trd;
  Array2D<int> idxc_trd;
  Array2D<int> idxz_trd;
  
  Array1D<SNADOUBLE> fact;

  // derivatives of data
  Array4D<SNAcomplex> dulist;
  Array2D<SNAcomplex> ylist;

  static const int nmaxfactorial = 167;
  static const SNADOUBLE nfac_table[];
  SNADOUBLE factorial(int);

  void create_twojmax_arrays();
  void destroy_twojmax_arrays();
  void init_clebsch_gordan();
  void init_rootpqarray();
  void zero_uarraytot();
  void addself_uarraytot(SNADOUBLE);
  void add_uarraytot(int, int, SNADOUBLE);
  void compute_uarray(int, int, SNADOUBLE, SNADOUBLE, SNADOUBLE,
                      SNADOUBLE);
  SNADOUBLE deltacg(int, int, int);
  int compute_ncoeff();
  SNADOUBLE compute_sfac(SNADOUBLE, SNADOUBLE);
  SNADOUBLE compute_dsfac(SNADOUBLE, SNADOUBLE);

  // Sets the style for the switching function
  // 0 = none
  // 1 = cosine

  int switch_flag;

  // self-weight

  SNADOUBLE wself;

};

#endif

