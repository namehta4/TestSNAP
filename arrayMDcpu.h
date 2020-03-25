#ifndef LMP_ARRAYMD_H
#define LMP_ARRAYMD_H

/* ----------------------------------------------------------------
   Class library implementing simple array structures similar to Fortran,
   except since this is c++ the fast running index are the rightmost
   ( For this "CPU" version. )

   These are the necessary set for the snap demo code, and may not be
   completely defined for all combinations of sub-dimensional references
   and copy constructors.

   2D through 6D arrays are supported,
       ArrahxD<type> newarray;
   The default constructor initializes bounds, but not memory
   For a 3D example, to allocate memory use newarray.resize( n1,n2,n3)

   The '()' operator is overlaaded here to access array elements.

   Copy constructors are provided for openmp private/firstprivate.
   In this case, the data pointer must be allocated with resize or
   assignment as above.

   When using the GPU, the data pointers must be set up another way,
   presumably with "#pragma enter data map(alloc..." and a deep copy
   of the pointer. The "rebound" routine is a copy of "resize" without
   allocation.

   Sarah Anderson, Cray Inc.
   ----------------------------------------------------------------*/

// version that prints
#ifdef DEBUG
#include "sna.h"
#include <strings.h>
#include <openacc.h>

/* ----------------------------------------------------------------
   Debug utility - just prints array (better be small),
   no longer does checksumming
   ----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

class Checksum {
public:
  Checksum() {}
  ~Checksum() {}

  void add( SNAreal *data, int num_data ) {
    SNAreal v;
    for(int i=0; i<num_data; ++i)  {
      v= data[i];      fprintf(stderr,"%10.6g  ",v);
    }
  }
  void add( SNAcomplex *data, int num_data ) {
    SNAreal re,im;
    for(int i=0; i<num_data; ++i)  {
      re= data[i].real; im= data[i].imag;
      fprintf(stderr,"%10.6g  ",re);
      fprintf(stderr,"%10.6g  ",im);
    }
  }
};

#if 0
// test, demonstration of ussage
int main(int argc, char* argv[]) {
  srand48( 6516458945L + 1 );
  Checksum *check = new Checksum();
#define N 100000
  double *array = new double[N];
  for(int i=0; i<N; ++i) array[i]= drand48();
  check->add( array,N );
  delete check;
}
#endif

#endif

/*----------------------------------------------------------------
Here the rightmost array index varies most rapidly, as in C
  ----------------------------------------------------------------*/

template<typename T> class Array1D
{
 public:
  unsigned n1;
  unsigned size;
  T * dptr;

  inline T& operator() (unsigned i1) { return dptr[i1]; }

  Array1D() { n1=0; size=0; dptr=NULL; }
  Array1D(const Array1D &p) { n1=p.n1;size=0; dptr=p.dptr; }
  Array1D(int in1) { n1=in1; size=n1; dptr=(T*)malloc(size*sizeof(T)); }
  ~Array1D() { if (size && dptr) free(dptr); }

  void resize(unsigned in1) {
    if (size && dptr) free(dptr);
    n1=in1; size=n1;
    dptr= (T*)malloc(size*sizeof(T));
  }
  unsigned getSize() { return size * sizeof(T); }  // NB: in bytes
};

template<typename T> class Array2D
{
 public:
  unsigned n1, n2, b1;
  unsigned size;
  T * dptr;

  inline T& operator() (unsigned i1, unsigned i2) { return dptr[i2+(n2*i1)]; }

  Array2D() { n1=n2= 0; size=0; dptr=NULL; }
  Array2D(const Array2D &p) { n1=p.n1; n2=p.n2; size=0; dptr=p.dptr; }
  Array2D(int in1, int in2) { n1=in1; n2=in2; size=n1*n2; dptr=(T*)malloc(size*sizeof(T)); }
  ~Array2D() { if (size && dptr) free(dptr); }

  void resize(unsigned in1,unsigned in2) {
    if (size && dptr) free(dptr);
    n1=in1; n2=in2; size=n1*n2;
    dptr= (T*)malloc(size*sizeof(T));
  }
  void rebound(unsigned in1,unsigned in2) {
    n1=in1; n2=in2; size=n1*n2; dptr=NULL;
  }
  void setBase(unsigned i1,unsigned i2) { b1= i1*n2 + i2; }
  inline T& operator() (unsigned i2) { return dptr[b1+i2]; }

  void setSize(unsigned in1, unsigned in2) { n1 = in1; n2 = in2; size = n1*n2; }  // NB: in bytes
  unsigned getSize() { return size * sizeof(T); }  // NB: in bytes
};

/* ----------------------------------------------------------------
For classes larger than 2D, a short cut access arrangment is made for
repeated use of invariant left-hand-side indices. For example, if

a( ii,jj,k ) is repeatedly used in loops only over "k", ahead of
the k-loop set the base with
   a.setBase(ii,jj)

and then in the loop refer to
   a(k)

!!!! Note this _REQUIRES_ declaring the class "a" as
"firstprivate" so that parallel copies can have different base offsets
----------------------------------------------------------------*/

template<typename T> class Array3D
{
 private:
  unsigned f2,f1, b1,b2;
 public:
  unsigned n1,n2,n3;
  unsigned size;
  T * dptr;
  inline T& operator() (unsigned i1, unsigned i2, unsigned i3)
  { return dptr[i3+i2*f2+i1*f1]; }

  Array3D() { n1=n2=n3=0; size=0; dptr=NULL;  }
  Array3D(const Array3D &p) { n1=p.n1; n2=p.n2; n3=p.n3; size=0; dptr=p.dptr; f2=n3; f1=f2*n2; }
  Array3D(unsigned in1, unsigned in2, unsigned in3) {
    n1=in1; n2=in2; n3=in3; size= n1*n2*n3; dptr= (T*)malloc(size*sizeof(T)); f2=n3; f1=f2*n2;
  }
  ~Array3D() { if ( size && dptr ) free(dptr); }

  void resize(unsigned in1, unsigned in2, unsigned in3) {
    n1=in1; n2=in2; n3=in3; size=n1*n2*n3; f2=n3; f1=f2*n2;
    if ( size && dptr ) free(dptr);
    dptr= (T*)malloc(size*sizeof(T));
  }
  void rebound(unsigned in1, unsigned in2, unsigned in3) {
    n1=in1; n2=in2; n3=in3; size=n1*n2*n3;  dptr=NULL; f2=n3; f1=f2*n2;
  }

  inline void setBase(unsigned i1,unsigned i2) { b2= i1*f1+i2*f2; }
  inline void setBase(unsigned i1) { b1= i1*f1; }
  inline T& operator() (unsigned i2,unsigned i3) { return dptr[b1+i2*f2+i3]; }
  inline T& operator() (unsigned i3) { return dptr[b2+i3]; }

  void setSize(unsigned in1, unsigned in2) { n1 = in1; n2 = in2; size = n1*n2; }  // NB: in bytes
  unsigned getSize() { return size * sizeof(T); }
};

template<typename T> class Array4D
{
 private:
  unsigned f3,f2,f1, b1,b2,b3;
 public:
  unsigned n1,n2,n3,n4;
  unsigned size;
  T * dptr;
  inline T& operator() (unsigned i1, unsigned i2, unsigned i3, unsigned i4)
  { return dptr[i4+i3*f3+i2*f2+i1*f1]; }

  Array4D() { n1=n2=n3=n4=0; size=0; dptr=NULL; f3=n4; f2=f3*n3; f1=f2*n2; }
  Array4D(const Array4D &p) { n1=p.n1; n2=p.n2; n3=p.n3; n4=p.n4; size=0; dptr=p.dptr; f3=n4; f2=f3*n3; f1=f2*n2; }

  Array4D(unsigned in1, unsigned in2, unsigned in3, unsigned in4) {
    n1=in1; n2=in2; n3=in3; n4=in4; size= n1*n2*n3*n4; f3=n4; f2=f3*n3; f1=f2*n2;
    dptr= (T*)malloc(size*sizeof(T));
  }
  ~Array4D() { if ( size && dptr ) free(dptr); }

  void resize(unsigned in1, unsigned in2, unsigned in3, unsigned in4) {
    n1=in1; n2=in2; n3=in3; n4=in4; size= n1*n2*n3*n4; f3=n4; f2=f3*n3; f1=f2*n2;
    if ( size && dptr) free(dptr);
    dptr= (T*)malloc(size*sizeof(T));
  }
  void rebound(unsigned in1, unsigned in2, unsigned in3, unsigned in4) {
    n1=in1; n2=in2; n3=in3; n4=in4; size= n1*n2*n3*n4; dptr=NULL; f3=n4; f2=f3*n3; f1=f2*n2;
  }

  inline void setBase(unsigned i1,unsigned i2, unsigned i3) { b3= i1*f1+i2*f2+i3*f3; }
  inline void setBase(unsigned i1,unsigned i2) { b2= i1*f1+i2*f2; }
  inline void setBase(unsigned i1) { b3= i1*f1; }
  inline T& operator() (unsigned i2, unsigned i3, unsigned i4) { return dptr[b1+i2*f2+i3*f3+i4]; }
  inline T& operator() (unsigned i3, unsigned i4) { return dptr[b2+i3*f3+i4]; }
  inline T& operator() (unsigned i4) { return dptr[b3+i4]; }

  unsigned getSize() { return size * sizeof(T); }
  void setSize(unsigned in1, unsigned in2) { n1 = in1; n2 = in2; size = n1*n2; }  // NB: in bytes
};

template<typename T> class Array5D
{
 private:
  unsigned f4,f3,f2,f1, b1,b2,b3,b4;
 public:
  unsigned n1,n2,n3,n4,n5;
  unsigned size;
  T * dptr;

  // fully dereference with 5d reference
  inline T& operator() (unsigned i1, unsigned i2, unsigned i3, unsigned i4, unsigned i5)
  { return dptr[i5+i4*f4+i3*f3+i2*f2+i1*f1]; }

  Array5D() { n1=n2=n3=n4=n5=0; size=0; dptr=NULL;  }
  Array5D(const Array5D &p) {
    n1=p.n1; n2=p.n2; n3=p.n3; n4=p.n4; n5=p.n5; size=0; dptr=p.dptr;
    f4=n5; f3=f4*n4; f2=f3*n3; f1=f2*n2;
  }
  Array5D(unsigned in1, unsigned in2, unsigned in3, unsigned in4, unsigned in5) {
    n1=in1; n2=in2; n3=in3; n4=in4; n5=in5; size= n1*n2*n3*n4*n5; dptr= (T*)malloc(size*sizeof(T));
    f4=n5; f3=f4*n4; f2=f3*n3; f1=f2*n2;
  }
  ~Array5D() { if ( size && dptr ) free(dptr); }

  void resize(unsigned in1, unsigned in2, unsigned in3, unsigned in4, unsigned in5) {
    n1=in1; n2=in2; n3=in3; n4=in4; n5=in5; size= n1*n2*n3*n4*n5;
    f4=n5; f3=f4*n4; f2=f3*n3; f1=f2*n2;
    if (size && dptr) free(dptr);
    dptr= (T*)malloc(size*sizeof(T));
  }
  void rebound(unsigned in1, unsigned in2, unsigned in3, unsigned in4, unsigned in5) {
    n1=in1; n2=in2; n3=in3; n4=in4; n5=in5; size= n1*n2*n3*n4*n5; dptr=NULL;
    f4=n5; f3=f4*n4; f2=f3*n3; f1=f2*n2;
  }

  inline void setBase(unsigned i1,unsigned i2, unsigned i3, unsigned i4) { b4= i1*f1+i2*f2+i3*f3+i4*f4; }
  inline void setBase(unsigned i1,unsigned i2, unsigned i3) { b3= i1*f1+i2*f2+i3*f3; }
  inline void setBase(unsigned i1,unsigned i2) { b2= i1*f1+i2*f2; }
  inline void setBase(unsigned i1) { b1= i1*f1; }
  inline T& operator() (unsigned i2, unsigned i3, unsigned i4, unsigned i5) { return dptr[b1+i2*f2+i3*f3+i4*f4+i5]; }
  inline T& operator() (unsigned i3, unsigned i4, unsigned i5) { return dptr[b2+i3*f3+i4*f4+i5]; }
  inline T& operator() (unsigned i4, unsigned i5) { return dptr[b3+i4*f4+i5]; }
  inline T& operator() (unsigned i5) { return dptr[b4+i5]; }

  unsigned getSize() { return size * sizeof(T); }
  void setSize(unsigned in1, unsigned in2) { n1 = in1; n2 = in2; size = n1*n2; }  // NB: in bytes
};

template<typename T> class Array6D
{
 private:
  unsigned f1,f2,f3,f4,f5, b6,b5,b4,b3,b2,b1;
 public:
  unsigned n1,n2,n3,n4,n5,n6;
  unsigned size;
  T * dptr;

  // fully dereference with 6d reference
  inline T& operator() (unsigned i1, unsigned i2, unsigned i3, unsigned i4, unsigned i5, unsigned i6)
  { return dptr[i6+i5*f5+i4*f4+i3*f3+i2*f2+i1*f1]; }

  Array6D() { n1=n2=n3=n4=n5=n6=0; size=0; dptr=NULL; }
  Array6D(const Array6D &p) { n1=p.n1; n2=p.n2; n3=p.n3; n4=p.n4; n5=p.n5; n6=p.n6; size=0; dptr=p.dptr;
    f5=n6; f4=f5*n5; f3=f4*n4; f2=f3*n3; f1=f2*n2;
  }
  Array6D(unsigned in1, unsigned in2, unsigned in3, unsigned in4, unsigned in5, unsigned in6) {
    n1=in1; n2=in2; n3=in3; n4=in4; n5=in5; n6=in6; size= n1*n2*n3*n4*n5*n6;
    f5=n6; f4=f5*n5; f3=f4*n4; f2=f3*n3; f1=f2*n2;
    dptr= (T*)malloc(size*sizeof(T));
  }
  ~Array6D() { if ( size && dptr ) free(dptr); }

  void resize(unsigned in1, unsigned in2, unsigned in3, unsigned in4, unsigned in5, unsigned in6) {
    n1=in1; n2=in2; n3=in3; n4=in4; n5=in5; n6=in6; size= n1*n2*n3*n4*n5*n6;
    f5=n6; f4=f5*n5; f3=f4*n4; f2=f3*n3; f1=f2*n2;
    if (size && dptr) free(dptr);
    dptr= (T*)malloc(size*sizeof(T));
  }
  void rebound(unsigned in1, unsigned in2, unsigned in3, unsigned in4, unsigned in5, unsigned in6) {
    n1=in1; n2=in2; n3=in3; n4=in4; n5=in5; n6=in6; size= n1*n2*n3*n4*n5*n6; dptr=NULL;
    f5=n6; f4=f4*n5; f3=f4*n4; f2=f3*n3; f1=f2*n2;
  }

  inline void setBase(unsigned i1,unsigned i2, unsigned i3, unsigned i4,unsigned i5) { b5= i1*f1+i2*f2+i3*f3+i4*f4+i5*f5; }
  inline void setBase(unsigned i1,unsigned i2, unsigned i3, unsigned i4) { b4= i1*f1+i2*f2+i3*f3+i4*f4; }
  inline void setBase(unsigned i1,unsigned i2, unsigned i3) { b3= i1*f1+i2*f2+i3*f3; }
  inline void setBase(unsigned i1,unsigned i2) { b2= i1*f1+i2*f2; }
  inline void setBase(unsigned i1) { b1= i1*f1; }
  inline T& operator() (unsigned i2, unsigned i3, unsigned i4, unsigned i5,unsigned i6) { return dptr[b1+i2*f2+i3*f3+i4*f4+i5*f5+i6]; }
  inline T& operator() (unsigned i3, unsigned i4, unsigned i5, unsigned i6) { return dptr[b2+i3*f3+i4*f4+i5*f5+i6]; }
  inline T& operator() (unsigned i4, unsigned i5, unsigned i6) { return dptr[b3+i4*f4+i5*f5+i6]; }
  inline T& operator() (unsigned i5, unsigned i6) { return dptr[b4+i5*f5+i6]; }
  inline T& operator() (unsigned i6) { return dptr[b5+i6]; }

  unsigned getSize() { return size * sizeof(T); }
};

#endif
