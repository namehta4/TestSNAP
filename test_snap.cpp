/* ----------------------------------------------------------------------
   Copyright (2018) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.

   Author: Aidan P. Thompson
------------------------------------------------------------------------- */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono> 
#include <cmath>
#include <omp.h>
#include <sys/time.h>
#include "sna.h"
#include "memory.h"
#include "arrayMDcpu.h"
#include "test_snap.h"

#if REFDATA_TWOJ==8
#include "refdata_2J8_W.h"
#elif REFDATA_TWOJ==14
#include "refdata_2J14_W.h"
#elif REFDATA_TWOJ==2
#include "refdata_2J2_W.h"
#elif REFDATA_TWOJ==4
#include "refdata_2J4_W.h"
#endif

/* ----------------------------------------------------------------------
    Vars to record timings of individual routines
------------------------------------------------------------------------- */
static double elapsed_ui = 0.0,  
            elapsed_zi = 0.0,
            elapsed_yi = 0.0,
            elapsed_duidrj = 0.0,
            elapsed_deidrj = 0.0;
            
/* ----------------------------------------------------------------------
  Elapsed Time
------------------------------------------------------------------------- */
inline double elapsedTime(timeval start_time, timeval end_time)
{
    return ((end_time.tv_sec - start_time.tv_sec) +1e-6*(end_time.tv_usec - start_time.tv_usec));
}   

 /* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
  Compute error
------------------------------------------------------------------------- */

inline void init_forces()
{
  // initialize all forces to zero
  for (int j = 0; j < ntotal; j++) {
    f(j,0) = 0.0;
    f(j,1) = 0.0;
    f(j,2) = 0.0;
  }
}

/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
  Read reference data
------------------------------------------------------------------------- */

    // generate neighbors, dummy values
inline void read_data()
{
  int jt = 0;
  int jjt = 0;
  for (int natom = 0; natom < nlocal; natom++) {
    for (int nbor = 0; nbor < ninside; nbor++) {
      snaptr->rij(natom,nbor,0) = refdata.rij[jt++];
      snaptr->rij(natom,nbor,1) = refdata.rij[jt++];
      snaptr->rij(natom,nbor,2) = refdata.rij[jt++];
      snaptr->inside(natom,nbor) = refdata.jlist[jjt++];
      snaptr->wj(natom,nbor) = 1.0;
      snaptr->rcutij(natom,nbor) = rcutfac;
    }
  }
}

/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
  Compute error
------------------------------------------------------------------------- */

inline void compute_forces(SNA* snaptr)
{
  // Fij = dEi/dRj = -dEi/dRi => add to Fi, subtract from Fj
  for (int natom = 0; natom < nlocal; natom++) {
    for (int nbor = 0; nbor < ninside; nbor++) {
      int j = snaptr->inside(natom,nbor);
      f(natom,0) += snaptr->dedr(natom,nbor,0); 
      f(natom,1) += snaptr->dedr(natom,nbor,1); 
      f(natom,2) += snaptr->dedr(natom,nbor,2); 
      f(j,0) -= snaptr->dedr(natom,nbor,0); 
      f(j,1) -= snaptr->dedr(natom,nbor,1); 
      f(j,2) -= snaptr->dedr(natom,nbor,2); 
    } // loop over neighbor forces
  } // loop over atoms
}

/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
  Compute error
------------------------------------------------------------------------- */

inline void compute_error(SNA *snaptr)
{
  int jt = 0;
  for (int j = 0; j < ntotal; j++) {
    double ferrx = f(j,0)-refdata.fj[jt++];
    double ferry = f(j,1)-refdata.fj[jt++];
    double ferrz = f(j,2)-refdata.fj[jt++];
    sumsqferr += ferrx*ferrx + ferry*ferry + ferrz*ferrz;
  }
}

/* ---------------------------------------------------------------------- */


int main(int argc, char* argv[])
{
    int tid=0, numThreads=0, numTeams=0;
#if defined(openmp_version)
    printf("*******OpenMP TARGET version *******\n");
#pragma omp target map(tofrom: numTeams, numThreads)
#pragma omp teams shared(numTeams) private(tid)
   {
        tid = omp_get_team_num();
        if(tid == 0)
        {
            numTeams = omp_get_num_teams();
#pragma omp parallel
            {
                int ttid = omp_get_thread_num();
                if(ttid == 0)
                    numThreads = omp_get_num_threads();
            }
        }
    }
    printf("Number of OpenMP teams = %d\n", numTeams);
    printf("Number of OpenMP Device Threads = %d\n", numThreads);
#else
    printf("*******SEQ version *******\n");
#endif

  // process command line options
  options(argc, argv);

  // initialize data structures
  init();

    // loop over steps
  auto start = myclock::now();
  for (int istep = 0; istep < nsteps; istep++) {
    
    // evaluate force kernel
    compute();
  }
  
  auto stop = myclock::now();
  myduration elapsed = stop - start;

  printf("-----------------------\n");
  printf("Summary of TestSNAP run\n");
  printf("-----------------------\n");
  printf("natoms = %d \n",nlocal);
  printf("nghostatoms = %d \n",nghost);
  printf("nsteps = %d \n",nsteps);
  printf("nneighs = %d \n",ninside);
  printf("twojmax = %d \n",twojmax);
  printf("duration = %g [sec]\n",elapsed.count());
  printf("step time = %g [sec/step]\n",elapsed.count()/nsteps);
  printf("grind time = %g [msec/atom-step]\n",1000.0*elapsed.count()/(nlocal*nsteps));
  printf("RMS |Fj| deviation %g [eV/A]\n",sqrt(sumsqferr/(ntotal*nsteps)));

  printf("\n Individual routine timings\n");
  printf("compute_ui = %f\n", elapsed_ui);
  printf("compute_yi = %f\n", elapsed_yi);
  printf("compute_duidrj = %f\n", elapsed_duidrj);
  printf("compute_deidrj = %f\n", elapsed_deidrj);  
}

/* ----------------------------------------------------------------------
   Allocate memory and initialize data structures
------------------------------------------------------------------------- */

void options(int argc, char* argv[]) {

  for (int i = 1; i < argc; i++) {

    if ( (strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0) ) {
      printf("TestSNAP 1.0 (stand-alone SNAP force kernel)\n\n");
      printf("The following optional command-line switches override default values\n");
      printf("-ns, --nsteps <val>: set the number of force calls to val (default 1)\n");
      exit(0);
    } else if ( (strcmp(argv[i], "-ns") == 0) || (strcmp(argv[i], "--nsteps") == 0) ) {
      nsteps = atoi(argv[++i]);
    } else {
      printf("ERROR: Unknown command line argument: %s\n",argv[i]);
      exit(1);
    }
  }
}

/* ----------------------------------------------------------------------
   Allocate memory and initialize data structures
------------------------------------------------------------------------- */

void init() {

  // initialize SNAP model using reference data

  ninside = refdata.ninside;
  ncoeff = refdata.ncoeff;
  nlocal = refdata.nlocal;
  nghost = refdata.nghost;
  ntotal = nlocal+nghost;
  twojmax = refdata.twojmax;
  rcutfac = refdata.rcutfac;

  // allocate SNA object
  memory = new Memory();

  // omit beta0 from beta vector
  snaptr = new SNA(nlocal,ninside,memory,rfac0,twojmax,
                   rmin0,switchflag,bzeroflag);
  int tmp = snaptr->ncoeff;
  if (tmp != ncoeff) {
    printf("ERROR: ncoeff from SNA does not match reference data\n");
    exit(1);
  }
  
  snaptr->beta.resize(ncoeff+1);
  for (int icoeff = 0; icoeff < ncoeff+1; icoeff++){
    snaptr->beta(icoeff) = refdata.coeff[icoeff];
  }

  f.resize(ntotal,3);
  snaptr->num_atoms = nlocal;
  snaptr->num_nbors = ninside;
  snaptr->grow_rij(ninside);

  // initialize SNA object
  snaptr->init();

  // initialize error tally
  sumsqferr = 0.0;
}

/* ----------------------------------------------------------------------
   Calculate forces on all local atoms 
------------------------------------------------------------------------- */

void compute() {
  timeval startTimer, endTimer;
  
  init_forces();
  read_data();
  
#if defined(openmp_version)
  snaptr->omp_offload_update();
#endif

  // compute Ui for all atoms (loops over all atoms)
  gettimeofday(&startTimer, NULL);
  snaptr->compute_ui();
  gettimeofday(&endTimer, NULL);
  elapsed_ui += elapsedTime(startTimer, endTimer);

  // compute Yi for all atoms (loops over all atoms)
  gettimeofday(&startTimer, NULL);
  snaptr->compute_yi();
  gettimeofday(&endTimer, NULL);
  elapsed_yi += elapsedTime(startTimer, endTimer);


  // compute dUi/drj by looping over neighbors for all atoms within cutoff
  gettimeofday(&startTimer, NULL);
  snaptr->compute_duidrj();
  gettimeofday(&endTimer, NULL);
  elapsed_duidrj += elapsedTime(startTimer, endTimer);
      
  // compute dBi/drj by looping over neighbors for all atoms within cutoff
  gettimeofday(&startTimer, NULL);
  snaptr->compute_deidrj();
  gettimeofday(&endTimer, NULL);
  elapsed_deidrj += elapsedTime(startTimer, endTimer);
       
  compute_forces(snaptr); 
  compute_error(snaptr);
}

