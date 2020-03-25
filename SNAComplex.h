#ifndef __SNAComplex
#define __SNAComplex

#include <iostream>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <sys/time.h>

#if(__OpenACC)
#include "arrayMDgpu.h"
#else
#include "arrayMDcpu.h"
#endif

#include"SNAComplex.h"

#if(__OpenMP)
#include <omp.h>
#endif

using SNAreal = double;
using SNAcomplex = struct{SNAreal real, imag;};

template<class type>
class SNACustomComplex {
    private:

    public:
        SNAreal real;
        SNAreal imag;

    explicit SNACustomComplex () {
        real = 0.00;
        imag = 0.00;
    }


     explicit SNACustomComplex(const double& a, const double& b) {
        real = a;
        imag = b;
    }

     SNACustomComplex(const SNACustomComplex& src) {
        real = src.x;
        imag = src.y;
    }

     ~SNACustomComplex()
     {
     }

};

#endif
