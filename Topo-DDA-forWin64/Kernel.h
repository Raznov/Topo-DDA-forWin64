#ifndef TOPO_KERNEL_H_
#define TOPO_KERNEL_H_

#include <cuda_runtime_api.h>
#include <cuda.h>
#include "cufft.h"

//kernel wrappers:
void A2As(double* A, cufftDoubleComplex* A00, cufftDoubleComplex* A01, cufftDoubleComplex* A02, cufftDoubleComplex* A11, cufftDoubleComplex* A12, cufftDoubleComplex* A22, int NxFFT, int NyFFT, int NzFFT);
void B2Bs(double* bDev, cufftDoubleComplex* bxDev, cufftDoubleComplex* byDev, cufftDoubleComplex* bzDev, int NxFFT, int NyFFT, int NzFFT);
void Conv(cufftDoubleComplex* Convx, cufftDoubleComplex* Convy, cufftDoubleComplex* Convz, cufftDoubleComplex* A00, cufftDoubleComplex* A01, cufftDoubleComplex* A02, cufftDoubleComplex* A11, cufftDoubleComplex* A12, cufftDoubleComplex* A22, cufftDoubleComplex* bxDev, cufftDoubleComplex* byDev, cufftDoubleComplex* bzDev, int NxFFT, int NyFFT, int NzFFT);
void Conv2B(cufftDoubleComplex* Convx, cufftDoubleComplex* Convy, cufftDoubleComplex* Convz, double* bDev, int NxFFT, int NyFFT, int NzFFT);
void APtoESum(cufftDoubleComplex* A00, cufftDoubleComplex* A01, cufftDoubleComplex* A02, cufftDoubleComplex* A11, cufftDoubleComplex* A12, cufftDoubleComplex* A22, cufftDoubleComplex* PxDev, cufftDoubleComplex* PyDev, cufftDoubleComplex* PzDev, cufftDoubleComplex* ESumxDev, cufftDoubleComplex* ESumyDev, cufftDoubleComplex* ESumzDev, int NxFFT, int NyFFT, int NzFFT, int NxA, int NyA, int NzA, int index1, int index2, int index3, int deduction);

#endif
