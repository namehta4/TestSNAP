// ----------------------------------------------------------------------
// Copyright (2019) Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000
// with Sandia Corporation, the U.S. Government
// retains certain rights in this software. This
// software is distributed under the Zero Clause
// BSD License
//
// TestSNAP - A prototype for the SNAP force kernel
// Version 0.0.2
// Main changes: Y array trick, memory compaction
//
// Original author: Aidan P. Thompson, athomps@sandia.gov
// http://www.cs.sandia.gov/~athomps, Sandia National Laboratories
//
// Additional authors:
// Sarah Anderson
// Rahul Gayatri
// Steve Plimpton
// Christian Trott
//
// Collaborators:
// Stan Moore
// Evan Weinberg
// Nick Lubbers
// Mitch Wood
//
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   Contributing authors: Aidan Thompson, Christian Trott, SNL
------------------------------------------------------------------------- */

#ifndef LMP_SNA_H
#define LMP_SNA_H

// Array class
#include "arrayMDcpu.h"
//using SNADOUBLE = double ;
//using SNAFLOAT = float ;
//using SNAcomplex = struct {SNADOUBLE re, im ;};

class SNA {

  int num_atoms, num_nbor;

public:
  SNA(int, int, SNADOUBLE, int, SNADOUBLE, int, int, const double*);
  ~SNA();
  void build_indexlist(const double*);
  void init();
  SNADOUBLE memory_usage();

  int ncoeff;

  // functions for bispectrum coefficients

  void compute_ui(int);

  void compute_yi(int,SNADOUBLE*);

  // functions for derivatives

  void compute_duidrj(int, int);
  void compute_deidrj(SNADOUBLE*);

  // per sna class instance for OMP use

  Array3D<SNADOUBLE> rij;
  Array2D<int> inside;
  Array2D<SNADOUBLE> wj;
  Array2D<SNADOUBLE> rcutij;
  int nmax;

  void grow_rij(int);

private:
  SNADOUBLE rmin0, rfac0;

  // use indexlist instead of loops, constructor generates these

  Array2D<int> idxz_j1j2j;
  Array2D<int> idxz_ma;
  Array2D<int> idxz_mb;
  Array2D<int> idxb;
  int idxcg_max, idxu_max, idxz_max, idxb_max;
  int idxz_j1j2j_max, idxz_ma_max, idxz_mb_max;

  // data for bispectrum coefficients

  int twojmax;
  Array2D<SNADOUBLE> rootpqarray;
  Array1D<SNADOUBLE> cglist;
  Array3D<int> idxcg_block;

  Array1D<SNAcomplex> ulisttot;
  Array1D<SNAcomplex> ulist;
  Array1D<int> idxu_block;

  Array3D<int> idxb_block;

  // derivatives of data

  Array2D<SNAcomplex> dulist;
  Array1D<SNAcomplex> ylist;

  static const int nmaxfactorial = 167;
  static const SNADOUBLE nfac_table[];
  SNADOUBLE factorial(int);

  void create_twojmax_arrays();
  void destroy_twojmax_arrays();
  void init_clebsch_gordan();
  void init_rootpqarray();
  void zero_uarraytot(int);
  void addself_uarraytot(int,SNADOUBLE);
  void add_uarraytot(int, SNADOUBLE, SNADOUBLE, SNADOUBLE);
  void compute_uarray(int, SNADOUBLE, SNADOUBLE, SNADOUBLE,
                      SNADOUBLE, SNADOUBLE);
  SNADOUBLE deltacg(int, int, int);
  int compute_ncoeff();
  void compute_duarray(int,SNADOUBLE, SNADOUBLE, SNADOUBLE,
                       SNADOUBLE, SNADOUBLE, SNADOUBLE, SNADOUBLE, SNADOUBLE);
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

