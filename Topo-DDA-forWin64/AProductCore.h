#ifndef TOPO_APRODUCTCORE_H_
#define TOPO_APRODUCTCORE_H_

#include <cuda_runtime_api.h>
#include <cuda.h>
#include "cufft.h"

#include "CoreStructure.h"
#include "Kernel.h"
#include "SiCi.h"

class AProductCore {
private:
    //---------------------------------Necessary values for A matrix-----------------------------------------------------
    CoreStructure* CStr;
    double K;
    double lam;

    //FFT related variables;
    double* AHos;                              //A_dicDoubl
    double* ADev;                              // double, but actually double*2 course real and imag are both stored in this double
    cufftDoubleComplex* A00, * A01, * A02, * A11, * A12, * A22; //The components in A for FFT, only in device
    double* bHos;                              //b in Aproduct
    double* bDev;                              // double, but actually double*2 course real and imag are both stored in this double
    cufftDoubleComplex* bxDev, * byDev, * bzDev; //b components for FFT, only in device
    int NxFFT, NyFFT, NzFFT;                   //2*(Nx, Ny, Nz) - 1
    int NFFT;                                  //NxFFT*NyFFT*NzFFT
    cufftDoubleComplex* Convx, * Convy, * Convz; //convolution of A and b on device

    //FFT plan for all
    cufftHandle Plan;

    //---------------------------------------Necessary for periodic A matrix----------------------------------------
    int MAXm;                            //-MAXm<=m<=MAXm
    int MAXn;                             //-MAXn<=n<=MAXn
    double Lm;                         //desplacement vector for one period, Currently should only be in x and y direction; d included do not need to time d
    double Ln;

    //--------------------------------Not necessary for A matrix but should be the same for diff DDAModel using the same A matrix------------------------------

    //VectorXcd diel;                   //real diel after 0~1 corresponds to real numbers
    VectorXcd material;
    //VectorXcd diel_max;                         //corresponds to the previous maximum obj
    double nback;                       //background material refractive index. 1.0 for air.

    //------------------------------For FCD and LDR choice-----------------------
    string AMatrixMethod;
    SiCi* SiCiValue;

public:
    AProductCore(CoreStructure* CStr_, double lam_, VectorXcd material_, double nback_, string AMatrixMethod_);
    AProductCore(CoreStructure* CStr_, double lam_, VectorXcd material_, double nback_, int MAXm_, int MAXn_, double Lm_, double Ln_, string AMatrixMethod_);
    ~AProductCore();
    Matrix3cd A_dic_generator(double x, double y, double z);
    Matrix3cd A_dic_generator(double x, double y, double z, int m, int n);
    VectorXcd Aproduct(VectorXcd& b);                                          //without al*b because al is in DDAModel and can be diff for the same AMatrix
    //void UpdateStr(VectorXd step);                                      //alpha not updated because in DDAModel, do not forget!
    //void output_to_file();
    //void output_to_file(string save_position, int iteration);

    CoreStructure* get_CStr();
    int get_N();
    int get_Nx();
    int get_Ny();
    int get_Nz();
    double get_Lm();
    double get_Ln();
    tuple<list<int>, list<int>, list<int>, list<int>> get_para_info();
    VectorXi* get_R();
    double get_d();
    double get_lam();
    double get_K();
    //VectorXcd* get_diel();        
    VectorXd* get_diel_old();
    VectorXcd* get_material();
    double get_nback();
    //VectorXcd* get_diel_max();                        
    VectorXd* get_diel_old_max();
    Matrix3cd FCD_inter(double x, double y, double z);
    Matrix3cd LDR_inter(double x, double y, double z);
};

#endif