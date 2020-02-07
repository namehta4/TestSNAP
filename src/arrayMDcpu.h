#ifndef LMP_ARRAYMD_H
#define LMP_ARRAYMD_H

/* ----------------------------------------------------------------
   Class library implementing simple array structures similar to Fortran,
   except since this is c++ the fast running index is the rightmost.

   These are the necessary set features needed by the TestSNAP code, and may not be
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

   Sarah Anderson, Cray Inc
   Rahul Gayatri, NERSC.
   ----------------------------------------------------------------*/

/*----------------------------------------------------------------
The rightmost array index varies most rapidly, as in C
  ----------------------------------------------------------------*/

#include <iostream>

using SNADOUBLE = double ;
using SNAFLOAT = float ;
using SNAcomplex = struct {SNADOUBLE re, im ;};


/*
 * Common member of the class
 * n1 = number of elements in 1st dimension
 * n2 = number of elements in 2nd dimension
 * n3 = number of elements in 3rd dimension
 * n4 = number of elements in 4th dimension
 * num_elements - total number of elements (i.e., n1*n2 for 2 dimensional arrays)
 * size - total size of the num_elements in the memory
 * dptr - base pointer of the 1D array in the memory
 */

// 1D array
template<typename T> class Array1D
{
  private:
    unsigned n1, num_elements;
    size_t size;
    T * dptr;

  public:
    inline T& operator() (unsigned i1) { return dptr[i1]; }

    //Default constructor
    Array1D() = default;

    // Copy constructor
    Array1D(const Array1D &p)
    {
      n1 = p.n1;
      num_elements = p.num_elements;
      size = 0;
      dptr = p.dptr;
    }

    // Constructor for creating in1 number of elements.
    Array1D(unsigned in1)
    {
      n1 = num_elements = in1;
      size = n1*sizeof(T);
      dptr = new T[n1];
    }

    // Destructor
    ~Array1D() { if (size && dptr) delete[]dptr; }

    // Resize a declared array to a new dimension
    void resize(unsigned in1)
    {
      if (size && dptr) delete[]dptr;
      n1 = num_elements = in1;
      size = num_elements*sizeof(T);
      dptr = new T[num_elements];
    }

    // Returns the base pointer. This is needed in the SNA.
    T* getBase() {return dptr;}

    // Returns the size of the allocated memory.
    size_t getSize() { return size; }
};

// 2D array
template<typename T> class Array2D
{
  private:
    unsigned n1, n2, num_elements;
    unsigned b1, b2;
    size_t size;
    T * dptr;

  public:
    inline T& operator() (unsigned i1, unsigned i2) { return dptr[i2+(n2*i1)]; }

    //Default constructor
    Array2D() = default;

    // Copy constructor
    Array2D(const Array2D &p)
    {
      n1 = p.n1; n2 = p.n2;
      num_elements = p.num_elements;
      size = 0;
      dptr = p.dptr;
    }

    // Constructor for creating in1 number of elements.
    Array2D(unsigned in1, unsigned in2)
    {
      n1 = in1; n2 = in2;
      num_elements = n1*n2;
      size = num_elements*sizeof(T);
      dptr = new T[num_elements];
    }

    // Destructor
    ~Array2D() { if (size && dptr) delete[]dptr; }

    // Resize a declared array to a new dimension
    void resize(unsigned in1, unsigned in2)
    {
      if (size && dptr) delete[]dptr;
      n1 = in1; n2 = in2;
      num_elements = n1*n2;
      size = num_elements*sizeof(T);
      dptr = new T[num_elements];
    }

    void setBase(unsigned i1)
    {
      b1 = i1*n2;
    }
    inline T& operator() (unsigned i2) { return dptr[i2+b1]; }

    // Returns the size of the allocated memory.
    size_t getSize() { return size; }
};

// 3D array
template<typename T> class Array3D
{
  private:
    unsigned n1, n2, n3, num_elements;
    unsigned f1, f2, b1, b2;
    size_t size;
    T * dptr;

  public:
    inline T& operator() (unsigned i1, unsigned i2, unsigned i3)
    { return dptr[i3+ i2*f2+ i1*f1]; }

    //Default constructor
    Array3D() = default;

    // Copy constructor
    Array3D(const Array3D &p)
    {
      n1 = p.n1; n2 = p.n2; n3 = p.n3;
      num_elements = p.num_elements;
      f1 = n2*n3; f2 = n3;
      size = 0;
      dptr = p.dptr;
    }

    // Constructor for creating in1 number of elements.
    Array3D(unsigned in1, unsigned in2, unsigned in3)
    {
      n1 = in1; n2 = in2; n3 = in3;
      num_elements = n1*n2*n3;
      f1 = n2*n3; f2 = n3;
      size = num_elements*sizeof(T);
      dptr = new T[num_elements];
    }

    // Destructor
    ~Array3D() { if (size && dptr) delete dptr; }

    // Resize a declared array to a new dimension
    void resize(unsigned in1, unsigned in2, unsigned in3)
    {
      if (size && dptr) delete[] dptr;
      n1 = in1; n2 = in2; n3 = in3;
      f1 = n2*n3; f2 = n3;
      num_elements = n1*n2*n3;
      size = num_elements*sizeof(T);
      dptr = new T[num_elements];
    }

    // Set base is a shortcut access to 2nd and 3rd dimension.
    void setBase(unsigned i1)
    {
      b1 = i1*f1;
    }
    inline T& operator() (unsigned i2, unsigned i3)
    { return dptr[i3+ i2*f2+ b1]; }

    void setBase(unsigned i1, unsigned i2)
    {
      b2 = i1*f1 + i2*f2;
    }
    inline T& operator() (unsigned i3)
    { return dptr[i3+ b2]; }

    // Returns the size of the allocated memory.
    size_t getSize() { return size; }
};

// 4D array
template<typename T> class Array4D
{
  private:
    unsigned n1, n2, n3, n4, num_elements;
    unsigned f1, f2, f3, b1, b2, b3;
    size_t size;
    T * dptr;

  public:
    inline T& operator() (unsigned i1, unsigned i2, unsigned i3, unsigned i4)
    { return dptr[i4+ i3*f3+ i2*f2+ i1*f1]; }

    //Default constructor
    Array4D() = default;

    // Copy constructor
    Array4D(const Array4D &p)
    {
      n1 = p.n1; n2 = p.n2; n3 = p.n3;
      num_elements = p.num_elements;
      f1 = n2*n3*n4; f2 = n3*n4; f3 = n3;
      size = 0;
      dptr = p.dptr;
    }

    // Constructor for creating in1 number of elements.
    Array4D(unsigned in1, unsigned in2, unsigned in3, unsigned in4)
    {
      n1 = in1; n2 = in2; n3 = in3, n4 = in4;
      num_elements = n1*n2*n3*n4;
      f1 = n2*n3*n4; f2 = n3*n4; f3 = n3;
      size = num_elements*sizeof(T);
      dptr = new T[num_elements];
    }

    // Destructor
    ~Array4D() { if (size && dptr) delete[]dptr; }

    // Resize a declared array to a new dimension
    void resize(unsigned in1, unsigned in2, unsigned in3, unsigned in4)
    {
      if (size && dptr) delete[]dptr;
      n1 = in1; n2 = in2; n3 = in3, n4 = in4;
      num_elements = n1*n2*n3*n4;
      f1 = n2*n3*n4; f2 = n3*n4; f3 = n3;
      size = num_elements*sizeof(T);
      dptr = new T[num_elements];
    }

    // Set base is a shortcut access to 2nd and 3rd dimension.
    void setBase(unsigned i1)
    {
      b1 = i1*f1;
    }
    inline T& operator() (unsigned i2, unsigned i3, unsigned i4)
    { return dptr[i4 + i3*f3+ i2*f2+ b1]; }

    void setBase(unsigned i1, unsigned i2)
    {
      b2 = i1*f1 + i2*f2;
    }
    inline T& operator() (unsigned i3, unsigned i4)
    { return dptr[i4 + i3*f3 + b2]; }

    // Returns the size of the allocated memory.
    size_t getSize() { return size; }
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
