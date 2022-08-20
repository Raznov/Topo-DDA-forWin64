#ifndef PREPROCESSING_H_INCLUDED
#define PREPROCESSING_H_INCLUDED
//#define _USE_MATH_DEFINES
//#include <fstream>
//#include <iostream>
//#include <complex>
//#include <cmath>
//#include <vector>
//#include <chrono>
//#include <omp.h>
//#include <string>
//#include <cstdlib>
//#include <map>
//#include <list>
//#include <tuple>
//#include <set>
//#include <time.h>
//
//#include "Eigen/Dense"
//#include "Eigen/Core"
//#include "Eigen/IterativeLinearSolvers"
//#include "Eigen/Sparse"
//#include "unsupported/Eigen/CXX11/Tensor"
//
////#include "cufftw.h"
//#include <cuda_runtime_api.h>
//#include <cuda.h>
//#include "cufft.h"
//#include "direct.h"
//#include "INIReader.h"
//#include "filterOption.h"
//#include "filterReader.h"
//#include "symReader.h"
//#include "EvoDDAModel.h"
//#include "ObjDDAModel.h"
//
//using namespace std;
//using namespace Eigen;
//using namespace std::chrono;
//
////void task();
//
//
//
//
//
//
//
//
//
//
//
//
////class task {
////    private:
////        string name;
////        map<string, string> findtask;
////    public:
////        void execute();
////};
//
//
//
//
//
//
//
//
//
////Abstract parent class for Obj function.
//
//
////class FOMscattering2D{
////private:
////    VectorXd FOMParameters;
////    double d;
////    int N;
////    VectorXcd* P;
////    VectorXi* R;
////    DDAModel* model;
////    double K;
////    int Paralength;
////    list<Vector3d> n_K_s_l;
////    double ATUC;
////    double E0;
////    Vector3d n_E0;
////    //Matrix3d FconstM;
////    
////
////public:
////    FOMscattering2D(list<double> parameters, DDAModel* model_);
////    list<double> GetVal();                                                        //Return list of far field abs(Es) at specified directions
////    Vector3cd FTUC(Vector3d n_K_s);
////};
////
////class FOMreflect2D {
////private:
////    VectorXd FOMParameters;
////    double d;
////    int N;
////    VectorXcd* P;
////    VectorXi* R;
////    DDAModel* model;
////    double K;
////    int Paralength;
////    list<Vector3d> n_K_s_l;
////    double ATUC;
////    double E0;
////    Vector3d n_E0;
////    //Matrix3d FconstM;
////
////
////public:
////    FOMreflect2D(list<double> parameters, DDAModel* model_);
////    list<double> GetVal();                                                        //Return list of far field abs(Es) at specified directions
////    Vector3cd FTUC(Vector3d n_K_s);
////};
////
////class FOMscattering0D {
////private:
////
////    double d;
////    int N;
////    VectorXcd* P;
////    VectorXi* R;
////    DDAModel* model;
////    double K;
////    double E0;
////    int Paralength;
////    list<Vector3d> n_K_s_l;
////    //double ATUC;
////    //Matrix3d FconstM;
////
////
////public:
////    FOMscattering0D(list<double> parameters, DDAModel* model_);
////    list<double> GetVal();                                                      //Return list of dCsca/dOmega at specified directions
////    double FTUCnsquare(Vector3d n_K_s);
////};
//
///*
//class FOMAbs {
//private:
//
//    double d;
//    int N;
//    VectorXcd* P;
//    VectorXi* R;
//    DDAModel* model;
//    double K;
//    double E0;
//    int Paralength;
//    list<Vector3d> n_K_s_l;
//    //double ATUC;
//    //Matrix3d FconstM;
//
//
//public:
//    FOMAbs(list<double> parameters, DDAModel* model_);
//    list<double> GetVal();                                                      //Return list of dCsca/dOmega at specified directions
//    double FTUCnsquare(Vector3d n_K_s);
//};
//*/
//
//
//
//
//
//
//
//
//
//
////void Evo_Focus(SpacePara* spacepara_tmp, CoreStructure* CStr, DDAModel* TestModel, string save_position, int start_num, int max_evo,
////    int min_num, int max_num, Vector3d lower_bound, Vector3d upper_bound, bool sym  //Parameters for focus generation
////);
////void Evo_single(string save_position, Vector3i bind, Vector3d l, int MAX_ITERATION_EVO, Vector3d move_focus);
////
////
////
////
////void eval_FOM(string name, DDAModel* TestModel, list<double> theta, list<double> phi);
////void eval_FOM_2Dperiod(string name, DDAModel* TestModel, list<double> theta, list<double> phi);

#endif // PREPROCESSING_H_INCLUDED