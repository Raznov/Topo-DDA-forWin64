#ifndef TOPO_OBJ_H_
#define TOPO_OBJ_H_

#include "DDAModel.h"

class ObjDDAModel {
public:
    bool Have_Devx;
    //bool Have_Penalty;
    virtual void SingleResponse(int idx, bool deduction) = 0;
    virtual double GroupResponse() = 0;
    virtual void Reset() = 0;
    virtual double GetVal() = 0;
};


//Child classes for Obj function


class ObjPointEDDAModel : public ObjDDAModel {
private:
    double x;
    double y;
    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    double d;
    int N;
    VectorXcd* P;
    VectorXi* R;
    DDAModel* model;
    //EvoDDAModel* evomodel;
    Vector3cd E_sum;
    Vector3cd E_ext;
public:
    ObjPointEDDAModel(vector<double> parameters, DDAModel* model_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

//class ObjPointListEDDAModel : public ObjDDAModel {
//private:
//    VectorXd x;
//    VectorXd y;
//    VectorXd z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
//    int PNum;
//    double d;
//    int N;
//    VectorXcd* P;
//    VectorXi* R;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    MatrixXcd E_sum;
//    MatrixXcd E_ext;
//public:
//    ObjPointListEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//};
//
//class ObjPointIDDAModel : public ObjDDAModel {
//private:
//    double x;
//    double y;
//    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
//    double d;
//    int N;
//    VectorXcd* P;
//    VectorXi* R;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    Vector3cd E_sum;
//    Vector3cd E_ext;
//public:
//    ObjPointIDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//};
//
//class ObjIntegratedEDDAModel : public ObjDDAModel {
//private:
//    double d;
//    int N;
//    int Nx;
//    int Ny;
//    int Nz;
//    VectorXcd* P;
//    VectorXcd* al;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    VectorXcd E;
//    double E_int;
//    VectorXi* R;
//public:
//    ObjIntegratedEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//};
//
//class ObjMidAvgEDDAModel : public ObjDDAModel {
//private:
//    double d;
//    int N;
//    int Nx;
//    int Ny;
//    int Nz;
//    VectorXcd* P;
//    VectorXcd* al;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    VectorXcd E;
//    double E_avg;
//    double r;                 //radius of the middle regions (for 2D)
//    VectorXi* R;
//    double centerx;
//    double centery;
//public:
//    ObjMidAvgEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//};
//
//class Objscattering2D : public ObjDDAModel {
//private:
//    double d;
//    int N;
//    VectorXcd* P;
//    VectorXi* R;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    double K;
//    double E0;
//    int Paralength;
//    list<Vector3d> n_K_s_l;
//    list<Vector3cd> PSum_l;
//    list<Matrix3d> FconstM_l;
//    double ATUC;
//
//public:
//    Objscattering2D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//    double FTUCnsquareoversinal();
//};
//
//class Objreflect2D : public ObjDDAModel {
//private:
//    double d;
//    int N;
//    VectorXcd* P;
//    VectorXi* R;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    double K;
//    double E0;
//    int Paralength;
//    list<Vector3d> n_K_s_l;
//    list<Vector3cd> PSum_l;
//    list<Matrix3d> FconstM_l;
//    double ATUC;
//
//public:
//    Objreflect2D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//    double FTUCnsquareoversinal();
//};
//
//class Objscattering0D : public ObjDDAModel {
//private:
//
//    double d;
//    int N;
//    VectorXcd* P;
//    VectorXi* R;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    double K;
//    double E0;
//    int Paralength;
//    list<Vector3d> n_K_s_l;
//    list<Vector3cd> PSum_l;
//    list<Matrix3d> FconstM_l;
//    //double ATUC;
//    //Matrix3d FconstM;
//
//
//public:
//    Objscattering0D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//    double FTUCnsquare();
//};
//
//class ObjAbs : public ObjDDAModel {
//private:
//    double d;
//    int N;
//    int Nx;
//    int Ny;
//    int Nz;
//    double K;
//    double K3;
//    double E0;
//    VectorXcd* P;
//    VectorXcd* al;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    VectorXcd E;
//    double Cabs;
//    VectorXi* R;
//public:
//    ObjAbs(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//};
//
//class ObjAbsPartial : public ObjDDAModel {
//private:
//    double d;
//    int N;
//    int Nx;
//    int Ny;
//    int Nz;
//    double K;
//    double K3;
//    double E0;
//    VectorXcd* P;
//    VectorXcd* al;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    VectorXcd E;
//    double Cabs;
//    VectorXi* R;
//    set<int> integralpos;             //Only geometries with pixels in this set will be considered.
//public:
//    ObjAbsPartial(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//};
//
//class ObjAbsPartialzslice : public ObjDDAModel {
//private:
//    double d;
//    int N;
//    int Nx;
//    int Ny;
//    int Nz;
//    double K;
//    double K3;
//    double E0;
//    VectorXcd* P;
//    VectorXcd* al;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    VectorXcd E;
//    double Cabs;
//    VectorXi* R;
//    set<int> integralpos;             //Only geometries with pixels in this set will be considered.
//    set<int> zslices{ 27, 28, 29 };
//public:
//    ObjAbsPartialzslice(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//};
//
//class ObjAbsbyfar : public ObjDDAModel {
//private:
//    double d;
//    int N;
//    VectorXcd* P;
//    VectorXi* R;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    double K;
//    double E0;
//    int Paralength;
//    list<Vector3d> n_K_s_l;
//    list<Vector3cd> PSum_l;
//    list<Matrix3d> FconstM_l;
//    double ATUC;
//
//public:
//    ObjAbsbyfar(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//    double FTUCnsquareoversinal();
//};
//
//
//class ObjIntegrateEPartial : public ObjDDAModel {
//private:
//    double d;
//    int N;
//    int Nx;
//    int Ny;
//    int Nz;
//    double K;
//    double K3;
//    double E0;
//    VectorXcd* P;
//    VectorXcd* al;
//    DDAModel* model;
//    EvoDDAModel* evomodel;
//    VectorXcd E;
//    //double Cabs;
//    double Total;
//    VectorXi* R;
//    set<int> integralpos;             //Only geometries with pixels in this set will be considered.
//public:
//    ObjIntegrateEPartial(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
//    void SingleResponse(int idx, bool deduction);
//    double GroupResponse();
//    double GetVal();
//    void Reset();
//};
#endif
