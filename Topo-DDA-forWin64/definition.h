#ifndef PREPROCESSING_H_INCLUDED
#define PREPROCESSING_H_INCLUDED
#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <chrono>
#include <omp.h>
#include <string>
#include <cstdlib>
#include <map>
#include <list>
#include <tuple>
#include <set>
#include <time.h>

//#include "Eigen/Dense"
//#include "Eigen/Core"
//#include "Eigen/IterativeLinearSolvers"
//#include "Eigen/Sparse"
//#include "unsupported/Eigen/CXX11/Tensor"

//#include "cufftw.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "cufft.h"
#include "direct.h"



using namespace std;
//using namespace Eigen;
using namespace std::chrono;
class Vectorcd;
class Vectord;
class Vectori;


class Vectorcd {
private:
    vector<complex<double>> m_vec;
public:
    Vectorcd();
    Vectorcd(int input_size, double initial = 0.0);
    Vectorcd(vector<complex<double>> input);
    //Copy constructor
    Vectorcd(const Vectorcd& copy);
    //Move constructor
    Vectorcd(Vectorcd&& copy) = default;
    //Assignment operator
    Vectorcd& operator= (const Vectorcd& input);
    Vectorcd& operator= (Vectorcd&& input);

    //For () to be used to get the element in the Vector
    complex<double>& operator() (const int pos);
    const complex<double>& operator() (const int pos) const;

    //Vector * Scalar
    /*Vectorcd operator*(const complex<double>& input);
    Vectorcd operator*(const double& input);
    Vectorcd operator*(const int& input);*/
    //Vector / Scalar
    Vectorcd operator/(const complex<double>& input);
    Vectorcd operator/(const double& input);
    Vectorcd operator/(const int& input);
    //Vector += Vector
    Vectorcd& operator+=(const Vectorcd& rhs);
    Vectorcd& operator+=(const Vectord& rhs);
    Vectorcd& operator+=(const Vectori& rhs);
    //Vector -= Vector
    Vectorcd& operator-=(const Vectorcd& rhs);
    Vectorcd& operator-=(const Vectord& rhs);
    Vectorcd& operator-=(const Vectori& rhs);

    //Vector dot* Vector
    complex<double> dot(Vectorcd& input);
    complex<double> dot(Vectord& input);
    complex<double> dot(Vectori& input);
    void print();
    int size() const;
    double norm();
    complex<double> sum();
    Vectorcd vecpow(int input);
    Vectorcd cwiseAbs();

    //friend Vectorcd operator+(const Vectorcd& x, const Vectorcd& y);
    //friend Vectorcd operator+(const Vectord& x, const Vectorcd& y);
};

class Vectord {
private:
    vector<double> m_vec;
public:
    Vectord();
    Vectord(int input_size, double initial = 0.0);
    Vectord(vector<double> input);
    //Copy constructor
    Vectord(const Vectord& copy);
    //Move constructor
    Vectord(Vectord&& copy) = default;
    //Assignment operator
    Vectord& operator= (const Vectord& input);
    Vectord& operator= (Vectord&& input);
    //For () to be used to get the element in the Vector
    double& operator()(const int pos);
    const double& operator()(const int pos) const;

    //Vector * Scalar
    /*Vectorcd operator*(const complex<double>& input);
    Vectord operator*(const double& input);
    Vectord operator*(const int& input);*/
    //Vector / Scalar
    Vectorcd operator/(const complex<double>& input);
    Vectord operator/(const double& input);
    Vectord operator/(const int& input);
    //Vector += Vector
    Vectord& operator+=(const Vectord& rhs);
    Vectord& operator+=(const Vectori& rhs);
    //Vector -= Vector
    Vectord& operator-=(const Vectord& rhs);
    Vectord& operator-=(const Vectori& rhs);
    //Vector dot* Vector
    complex<double> dot(Vectorcd& input);
    double dot(Vectord& input);
    double dot(Vectori& input);
    void print();
    const int size() const;
    double norm();
    double sum();
    Vectord vecpow(int input);
    Vectord cwiseAbs();

    //friend Vectorcd operator+(const Vectord& x, const Vectorcd& y);
};

class Vectori {
private:
    vector<int> m_vec;
public:
    Vectori();
    Vectori(int input_size, double initial = 0.0);
    Vectori(vector<int> input);
    //Copy constructor
    Vectori(const Vectori& copy);
    //Move constructor
    Vectori(Vectori&& copy) = default;
    //Assignment operator
    Vectori& operator= (const Vectori& input);
    Vectori& operator= (Vectori&& input);
    //For () to be used to get the element in the Vector
    int& operator()(const int pos);
    const int& operator()(const int pos) const;

    //Vector * Scalar
    /*Vectorcd operator*(const complex<double>& input);
    Vectord operator*(const double& input);
    Vectori operator*(const int& input);*/
    //Vector / Scalar
    Vectorcd operator/(const complex<double>& input);
    Vectord operator/(const double& input);
    Vectori operator/(const int& input);
    //Vector += Vector
    Vectori& operator+=(const Vectori& rhs);
    //Vector -= Vector
    Vectori& operator-=(const Vectori& rhs);
    //Vector dot* Vector
    complex<double> dot(Vectorcd& input);
    double dot(Vectord& input);
    int dot(Vectori& input);
    void print();
    int size() const;
    double norm();
    int sum();
    Vectori vecpow(int input);
    Vectori cwiseAbs();
};

//Vectorcd vecadd(Vectorcd x, Vectorcd y);
//Vectorcd vecadd(Vectord x, Vectorcd y);
//Vectorcd vecadd(Vectorcd x, Vectord y);
//Vectord vecadd(Vectord x, Vectord y);
//Vectorcd vecsub(Vectorcd x, Vectorcd y);
//Vectorcd vecsub(Vectord x, Vectorcd y);
//Vectorcd vecsub(Vectorcd x, Vectord y);
//Vectord vecsub(Vectord x, Vectord y);


Vectorcd operator+(const Vectorcd &x, const Vectorcd &y);
Vectorcd operator+(const Vectord& x, const Vectorcd& y);
Vectorcd operator+(const Vectorcd& x, const Vectord& y);
Vectord operator+(const Vectord& x, const Vectord& y);
Vectorcd operator-(const Vectorcd& x, const Vectorcd& y);
Vectorcd operator-(const Vectord& x, const Vectorcd& y);
Vectorcd operator-(const Vectorcd& x, const Vectord& y);
Vectord operator-(const Vectord& x, const Vectord& y);

Vectorcd operator*(const Vectorcd &x, const complex<double>& input);
Vectorcd operator*(const complex<double>& input, const Vectorcd& x);
Vectorcd operator*(const Vectorcd& x, const double& input);
Vectorcd operator*(const double& input, const Vectorcd& x);
Vectorcd operator*(const Vectorcd& x, const int& input);
Vectorcd operator*(const int& input, const Vectorcd& x);

Vectorcd operator*(const Vectord& x, const complex<double>& input);
Vectorcd operator*(const complex<double>& input, const Vectord& x);
Vectord operator*(const Vectord& x, const double& input);
Vectord operator*(const double& input, const Vectord& x);
Vectord operator*(const Vectord& x, const int& input);
Vectord operator*(const int& input, const Vectord& x);

Vectorcd operator*(const Vectori& x, const complex<double>& input);
Vectorcd operator*(const complex<double>& input, const Vectori& x);
Vectord operator*(const Vectori& x, const double& input);
Vectord operator*(const double& input, const Vectori& x);
Vectori operator*(const Vectori& x, const int& input);
Vectori operator*(const int& input, const Vectori& x);


complex<double> vecmean(const Vectorcd& x);
double vecmean(const Vectord& x);
double vecmean(const Vectori& x);
void vectofile(ofstream& fout, Vectorcd object);
void vectofile(ofstream& fout, Vectord object);
void vectofile(ofstream& fout, Vectori object);


class Matrixcd {
private:
    vector<vector<complex<double>>> m_vec;
public:
    Matrixcd();
    Matrixcd(int input_size1, int input_size2);
    //Copy constructor
    Matrixcd(const Matrixcd& copy);
    Matrixcd(Matrixcd&& copy)=default;
    //Assignment operator
    Matrixcd& operator= (const Matrixcd& input);
    Matrixcd& operator= (Matrixcd&& input);
    //+=
    Matrixcd& operator+=(const Matrixcd& rhs);
    //For () to be used to get the element in the Vector
    complex<double>& operator()(const int pos1, const int pos2);
    const complex<double>& operator()(const int pos1, const int pos2) const;
    //Matrix * Scalar
    Matrixcd operator*(const complex<double>& input);
    Matrixcd operator*(const double& input);
    Matrixcd operator*(const int& input);
    //Matrix / Scalar
    Matrixcd operator/(const complex<double>& input);
    Matrixcd operator/(const double& input);
    Matrixcd operator/(const int& input);
    void print();
    Vectori getShape() const;
};

class Matrixd {
private:
    vector<vector<double>> m_vec;
public:
    Matrixd();
    Matrixd(int input_size1, int input_size2);
    //Copy constructor
    Matrixd(const Matrixd& copy);
    Matrixd(Matrixd&& copy)=default;
    //Assignment operator
    Matrixd& operator= (const Matrixd& input);
    Matrixd& operator= (Matrixd&& input);
    //For () to be used to get the element in the Vector
    double& operator()(const int pos1, const int pos2);
    const double& operator()(const int pos1, const int pos2) const;
    //Matrix * Scalar
    Matrixcd operator*(const complex<double>& input);
    Matrixd operator*(const double& input);
    Matrixd operator*(const int& input);
    //Matrix / Scalar
    Matrixcd operator/(const complex<double>& input);
    Matrixd operator/(const double& input);
    Matrixd operator/(const int& input);
    void print();
    Vectori getShape() const;
};

class Matrixi {
private:
    vector<vector<int>> m_vec;
public:
    Matrixi();
    Matrixi(int input_size1, int input_size2);
    //Copy constructor
    Matrixi(const Matrixi& copy);
    Matrixi(Matrixi&& copy) = default;
    //Assignment operator
    Matrixi& operator= (const Matrixi& input);
    Matrixi& operator= (Matrixi&& input);
    //For () to be used to get the element in the Vector
    int& operator()(const int pos1, const int pos2);
    const int& operator()(const int pos1, const int pos2) const;
    //Matrix * Scalar
    Matrixcd operator*(const complex<double>& input);
    Matrixd operator*(const double& input);
    Matrixi operator*(const int& input);
    //Matrix / Scalar
    Matrixcd operator/(const complex<double>& input);
    Matrixd operator/(const double& input);
    Matrixi operator/(const int& input);
    void print();
    Vectori getShape() const;
};

Matrixcd operator+(const Matrixcd& x, const Matrixcd& y);





//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------















//kernel wrappers:
void A2As(double *A, cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, int NxFFT, int NyFFT, int NzFFT);
void B2Bs(double *bDev, cufftDoubleComplex *bxDev, cufftDoubleComplex *byDev, cufftDoubleComplex *bzDev, int NxFFT, int NyFFT, int NzFFT);
void Conv(cufftDoubleComplex *Convx, cufftDoubleComplex *Convy, cufftDoubleComplex *Convz, cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, cufftDoubleComplex *bxDev, cufftDoubleComplex *byDev, cufftDoubleComplex *bzDev, int NxFFT, int NyFFT, int NzFFT);
void Conv2B(cufftDoubleComplex *Convx, cufftDoubleComplex *Convy, cufftDoubleComplex *Convz, double *bDev, int NxFFT, int NyFFT, int NzFFT);
void APtoESum(cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, cufftDoubleComplex *PxDev, cufftDoubleComplex *PyDev, cufftDoubleComplex *PzDev, cufftDoubleComplex *ESumxDev, cufftDoubleComplex *ESumyDev, cufftDoubleComplex *ESumzDev, int NxFFT, int NyFFT, int NzFFT, int NxA, int NyA, int NzA, int index1, int index2, int index3, int deduction);

complex<double> Get_Alpha(double lam, double K, double d, complex<double> diel, Vectord n_E0, Vectord n_K);
complex<double> Get_Alpha_FCD(double lam, double K, double d, complex<double> diel);
complex<double> Get_material(string mat, double wl, string unit);                  //name of mat to get its diel function at certain wavlength              
Vectorcd Get_2_material(string sub, string mat, double wl, string unit);          //a wrapper for Get_material
Vectorcd Get_X_material(list<string> mat_l, double wl, string unit);
double Average(Vectorcd* E, int N, double exponent);
double Get_Max(Vectorcd* E, int N);
double G(Vectorcd* E, int N, double exponent, double E0);
//ArrayXcd FFT(int nx,int ny,int nz,ArrayXcd in,int _direction);
Vectori build_a_bulk(int Nx, int Ny, int Nz);
bool CheckPerp(Vectord v1, Vectord v2);
Vectord nEPerpinXZ(double theta, double phi);
Matrixi find_scope_3_dim(Vectori* x);
Vectord initial_diel_func(string initial_diel, int N);
double initial_diel_func(string initial_diel);
list<double> makelist(double start, double end, double interval);
list<double> makelist(double start, double end, int number);
pair<Vectori, Vectord> InputInitial(string open_position, string model_label);

double exp_update(const double x, const double x_max, const double y_min, const double y_max);
double piecewise_update(const double x, const double x_max, const double y_min, const double y_max);
double linear_update(const double x, const double x_max, const double y_min, const double y_max);

int makedirect(string name);



struct filterinfo {
    int iteration;
    double rfilter;
};

class FilterOption {
    private:
        double beta;
        double ita;
        string beta_type;
        double beta_min;
        double beta_max;
        double rfilter;
        vector<filterinfo> rfilterlist;
        bool fixit;
        int MAX_ITERATION_FIXED;

    public:
        FilterOption(double beta_min_, double beta_max_, double ita_, string beta_type_, vector<filterinfo> rfilterlist_, bool fixit_=false, int MAX_ITERATION_FIXED_=100);
        double get_beta();
        double get_ita();
        double get_rfilter();
        void update_beta(const int iteration, const int Max_iteration);
        double SmoothDensity(double input);
        bool filterchange(int iteration);
};



class Structure{
    private:
        Vectori geometry;
        bool para_cond;
    public:
        //--------------------------------------Dont support dependent para build-------------------------------------------------------
        Structure(Vectori *total_space, Vectori *geometry_, bool para_cond_=false);            //vanilla initialization

        //Sphere
        Structure(Vectori* total_space, double r, Vectord center, bool para_cond_ = false);  //r: actual radius/d. center: actual center/d. In charge of Sphere 
        Structure(Vectori* total_space, double r, double h, Vectord center, bool para_cond_ = false);
        
        //Circle
        //Structure(Vectori *total_space, string initial_diel, double r, Vectori center, Vectori direction, int para_); //build a circle, direction is its normalized direction in Cart. coord.

        //From file
        //Structure(Vectori* total_space, string FileName, int para_);                                                    //Read a structure from txt file
        
        //--------------------------------------Support dependent para build-------------------------------------------------------------
        //Bulk
        Structure(Vectori *total_space, Vectord l, Vectord center, bool para_cond_ = false);    //Ractangular(both 2D and 3D). Para can only = 0&1 
        Structure(Vectori* total_space, Vectord l, Vectord center, Structure* Str, bool para_cond_ = false);
        //Duplicate
        //Structure(Vectori *total_space, Structure *s, Vectori direction, int times, int para_);                     //Initializa a Structure by duplicating a existing structure along a certain direction for several times. Direction is normalized and can only be alone x, y or z.
        //The original structure is not included. original structure + new structure = times * original structure. If set para=2 and original para=1, then depend on origin str as geometry_dep. If para=2 and original para=2, will copy origin geometry_dep.
        //Structure(Vectori* total_space, Structure* s, int dep_way);     //Special copy for setting up 2-fold and 4-fold symmetry dependence in xy plane. Para auto set to 2.
        //The original str at left down corner of the xy plane(left down corner is (0,0)). dep_way=1, 2, 3 corresponds to the other blocks in clock-wise. 


        //-------------------------------------Other member functions--------------------------------------------------------------------
        Vectori* get_geometry();
        void cut(Vectori* big, Vectori* smalll);
        int get_geometry_size();
        //bool sym_or_not();
        //vector<string> get_sym_condition();
        //vector<double> get_sym_axis();
        bool para_or_not();
};



class Space{
    private:
        Vectori *total_space;
        int Nx, Ny, Nz;
        int N;                        //total size of all the geometry inside the list(dipole size which is actual size/3)
        vector<Structure> *ln;
    public:
        Space(Vectori *total_space_, int Nx_, int Ny_, int Nz_, int N_, vector<Structure>*ln_);
        Vectori *get_total_space();
        int get_ln_size();
        tuple<int, int, int, int> get_Ns();
        vector<Structure> *get_ln();
        void show_something_about_Structures() const;
        friend Space operator+(const Space &s1, Structure &s2);
        

};

struct WeightPara {
    double weight;
    int position;
};

class SpacePara {
private:
    Space* space;
    Vectori geometry;                //3N dimension
    Vectori geometryPara;            //N dimension. N=number of dipoles. Each position stores the para index in Vectori Para : 0->Para[0]...
    vector<vector<int>> Paratogeometry;
    Vectord Para;                    //P dimension. P=number of parameters. Same as Para_origin if Filter=False. Filtered and biased para if Filter=True.
    Vectord Para_origin;             //Un-filtered, unbiased para. No use when Filter=False.
    Vectord Para_filtered;           //Filtered but unbiased para. No use when Filter=False.
    Matrixi scope;                   //[[xmin, xmax],[ymin, ymax],[zmin, zmax]]
    Vectori bind;
    Vectori FreeparatoPara;          //Position of free parameters inside Para. dimension<=P. FreeparatoPara[i] is the index of a free parameter inside Para.
    //vector<list<int>> Paratogeometry;  //P dimension. Each position stores a list of corresponding dipole index for parameter for this specific position.
    bool Filter;                      //True for with Filter. False for without Filter. Defualt as False for most initilizations.
    FilterOption* Filterstats;        //Only used when Filter=True
    vector<vector<WeightPara>> FreeWeight;
    vector<int> ParaDividePos;

public:
    //SpacePara(Space* space_, string initial_diel, Vectori geometry_, Vectord diel_); //Can freeze part of the parameter space
    //SpacePara(Space* space_, Vectori bind_, Vectori* geometryPara_, Vectord* Para_, Vectori* FreeparatoPara_);
    //SpacePara(Space* space_, Vectori bind_, string initial_diel); //l, center similar to bulk build in Structure class. Every 'bind' nearby dipoles correspond 
                                                                    //to 1 parameter in this bulk. bind=(2,2,2): 2*2*2; bind=(1,1,3):1*1*3
    //SpacePara(Space* space_, Vectori bind_, string initial_diel1, string initial_diel2); //One constant layer at bottom. One design region on top. divide_pos is the divide in geometry array of the two parts.
    
    //SpacePara(Space* space_, Vectori bind_, string initial_diel, Vectori* geometryPara_);
    //SpacePara(Space* space_, Vectori bind_, string initial_diel_center, string initial_diel_ring, double r, string type);   //ONly for 2d cylinder or spheres. r is raidus/d.

    //SpacePara(Space* space_, Vectori bind_, string initial_diel_background, list<string>* initial_diel_list, list<double>* r_list, list<Vector2d>* center_list);
    //Build 2d cylinders with diel in the list, rest of the diel is the backgroudn diel.

    SpacePara(Space* space_, Vectori bind_, int number, double limitx1, double limitx2, double limity1, double limity2);  //random rect in a region with extruded 2D geometry
    SpacePara(Space* space_, Vectori bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2);  //random rect with 3D
    SpacePara(Space* space_, Vectori bind_, int number, double limitx1, double limitx2, double limity1, double limity2, Vectori* geometryPara_);
    SpacePara(Space* space_, Vectori bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2, Vectori* geometryPara_);

    SpacePara(Space* space_, Vectori bind_, string initial_diel, list<Vectori*> FParaGeometry_, list<Vectori*> BParaGeometry_, list<double> BPara_);
    SpacePara(Space* space_, Vectori bind_, vector<string> initial_diel_list, list<Vectori*> FParaGeometry_, list<Vectori*> BParaGeometry_, list<double> BPara_, bool Filter_, FilterOption* Filterstats_ = NULL, string symmetry = "None", vector<double> symaxis = vector<double>{0.0,0.0}); //Should be exactly the same with the previous one except of Filter
    SpacePara(Space* space_, Vectori bind_, vector<string> initial_diel_list, vector<double> BPara_, bool Filter_ = false, FilterOption* Filterstats_ = NULL, string symmetry = "None", vector<double> symaxis = vector<double>{ 0.0,0.0 });
    SpacePara(Space* space_, Vectori bind_, Vectori* InputGeo, Vectord* Inputdiel, bool Filter_ = false, FilterOption* Filterstats_ = NULL, string symmetry = "None", vector<double> symaxis = vector<double>{ 0.0,0.0 });
    
    void ChangeBind(Vectori bind_);                                  //Change bind number
    void ChangeFilter();
    Vectori cut(Vectori* big, Vectori* smalll);

    Space* get_space();
    Vectori get_geometry();
    Vectori* get_geometryPara();
    Vectord* get_Para();
    Vectord* get_Para_origin();
    Vectord* get_Para_filtered();
    Vectori* get_bind();
    Vectori* get_Free();
    //vector<list<int>>* get_Paratogeometry();
    bool get_Filter();
    FilterOption* get_Filterstats();
    vector<vector<WeightPara>>* get_FreeWeight();
    vector<int>* get_ParaDividePos();
    vector<vector<int>>* get_Paratogeometry();
};

//Abstract parent class for objective function.
class Objective;
class ObjectivePointE;
class ObjectiveSurfaceEExp;

class SiCi {
public:
    int numberSi;
    int numberCi;
    double disSi;
    double disCi;
    Vectord Si;
    Vectord Ci;
    SiCi();
    double get_Si(double x);
    double get_Ci(double y);

};

class CoreStructure {
private:
    //---------------------------------Geometries, not related to wavelength-------------------------------
    Space* space;
    SpacePara* spacepara;
    int N;                        //Number of dipoles
    int Nx;                       //scope of space. Nx*Ny*Nz!=N
    int Ny;
    int Nz;
    double d;
    Vectori R;                      //Position of dipoles. Both R and RResult are unitless, so need to time d to get real number.

    //Vectori RDep;                   //Position of the dependent para points in space of the points in the same position as R
    //list<list<int>> PositionDep;    //First D has the D of para(para=1). Second D is the positions of other points dependent on the para in the 1stD.
    //Vectori PositionPara;          //Position i(in R) of the parameters (3*i=x, 3*i+1=y, 3*i+2=z)
    //list<int> para_nums;
    //list<int> para_starts;
    //list<int> para_dep_nums;
    //list<int> para_dep_starts;
    Vectord diel_old;                //The 0~1 version of diel, 3*N
    Vectord diel_old_max;
public:
    CoreStructure(SpacePara* spacepara_, double d_);
    void UpdateStr(Vectord step, int current_it, int Max_it);
    void UpdateStr(SpacePara* spacepara_);
    void UpdateStrSingle(int idx, double value);
    void output_to_file();
    void output_to_file(string save_position, int iteration, string mode = "normal");

    int get_N();
    int get_Nx();
    int get_Ny();
    int get_Nz();
    Vectori* get_R();
    double get_d();
    SpacePara* get_spacepara();
    //list<list<int>>* get_PositionDep();
    //Vectori* get_PositionPara();
    //list<int>* get_para_nums();
    //list<int>* get_para_starts();
    //list<int>* get_para_dep_nums();
    //list<int>* get_para_dep_starts();
    Vectord* get_diel_old();
    Vectord* get_diel_old_max();

};

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
    
    //Vectorcd diel;                   //real diel after 0~1 corresponds to real numbers
    Vectorcd material;
    //Vectorcd diel_max;                         //corresponds to the previous maximum obj
    double nback;                       //background material refractive index. 1.0 for air.

    //------------------------------For FCD and LDR choice-----------------------
    string AMatrixMethod;
    SiCi* SiCiValue;

public:
    AProductCore(CoreStructure* CStr_, double lam_, Vectorcd material_, double nback_, string AMatrixMethod_);
    AProductCore(CoreStructure* CStr_, double lam_, Vectorcd material_, double nback_, int MAXm_, int MAXn_, double Lm_, double Ln_, string AMatrixMethod_);
    ~AProductCore();
    Matrixcd A_dic_generator(double x, double y, double z);
    Matrixcd A_dic_generator(double x, double y, double z, int m, int n);
    Vectorcd Aproduct(Vectorcd& b);                                          //without al*b because al is in DDAModel and can be diff for the same AMatrix
    //void UpdateStr(Vectord step);                                      //alpha not updated because in DDAModel, do not forget!
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
    Vectori* get_R();
    double get_d();
    double get_lam();
    double get_K();
    //Vectorcd* get_diel();        
    Vectord* get_diel_old();               
    Vectorcd* get_material();
    double get_nback();
    //Vectorcd* get_diel_max();                        
    Vectord* get_diel_old_max();
    Matrixcd FCD_inter(double x, double y, double z);
    Matrixcd LDR_inter(double x, double y, double z);
};

class DDAModel {
private:
    //------------------------------------Get from AProductCore------------------------------ For the same AMatrix, these are all the same
    AProductCore* Core;

    //-----------------------------------Not from AProductCore------------------------------- For the same AMatrix, these can be diff for diff DDAModel
    bool RResultSwitch;               //0(false) for plot only E field on the structure points (fast), 1(true) for using RResult different from R to plot (slow but adjustable).
    Vectori RResult;                //The position matrix for the EResult (where you want to plot the E field)
    double E0; 
    Vectord n_E0;
    Vectord n_K;
    Vectorcd P;
    Vectorcd E;
    Vectorcd Einternal;              //E field on structure points
    Vectorcd EResult;                //E field on designated points
    Vectorcd al;                       // 1 over alpha instead of alpha.
    Vectorcd diel;                     //Real dielectric from diel_old. Needed to calculate the Lorentz factor.
    bool verbose;
    Vectorcd P_max;
    Vectorcd al_max;

    //------------------different for different angles------------------
    int time;
    int ITERATION;
    double Error;


public:
    DDAModel(AProductCore* AProductCore_, Vectord n_K_, double E0_, Vectord n_E0_);
    DDAModel(AProductCore* AProductCore_, Vectord n_K_, double E0_, Vectord n_E0_, Vectori* RResult_);
    void bicgstab(int MAX_ITERATION, double MAX_ERROR);
    void bicgstab(int MAX_ITERATION, double MAX_ERROR, int EVOITERATION);  //FOR DEBUG ONLY. OUTPUT SOME VALUE AT CERTAIN EVO ITERATION.
    void change_E(Vectorcd E_);
    void reset_E();             //reset E to E0                                
    void UpdateAlpha();                                //update alpha according to updated diel in AProductCore.
    void UpdateAlphaSingle(int idx);
    void solve_E();                                                        //update the result E field on each dipole or on a designated space
    void update_E_in_structure();                                          //update the result E field on each dipole 
    Vectorcd Aproductwithalb(Vectorcd& b);                    //add the al*b term on base of AproductCore
    void output_to_file();
    void output_to_file(string save_position, int iteration, int ModelLabel);              //especially used for EvoOptimization
    void output_to_file(string save_position, int iteration);             //For simplify output
    //void output_to_file(string save_position, double wavelength, int iteration);
    void InitializeP(Vectorcd& Initializer);
    Vectorcd* get_P();
    Vectord get_nE0();
    Vectord get_nK();
    double get_E0();
    Vectorcd* get_Einternal();
    AProductCore* get_Core();
    Vectorcd* get_al();
    Vectorcd* get_P_max();
    Vectorcd* get_al_max();
    int get_ITERATION();

    //-----------------From AProductCore-----------------------

    int get_N();
    int get_Nx();
    int get_Ny();
    int get_Nz();
    Vectori* get_R();
    double get_d();
    SpacePara* get_spacepara();
    double get_lam();
    double get_K();
    Vectord* get_diel_old();
    Vectorcd* get_material();
    Vectord* get_diel_old_max();
};

class ObjectiveDDAModel;


class EvoDDAModel {
private:
    double output_time;
    CoreStructure* CStr;
    list<DDAModel*> ModelList;                    //List of DDA models sharing the same AProductCore : "Core"
    int ModelNum;                                 //number of DDA model
    string save_position;

    list<Vectorcd> PforOrigin;
    list<Vectorcd> PforAdjoint;
    list<Vectorcd> PforOriginMax;
    list<Vectorcd> PforAdjointMax;

    list<list<double>*>* ObjectParameters;
    list<double> MajorObjectParameters;
    list<list<double>> MinorObjectParameters;

    list<string>* ObjectFunctionNames;
    string MajorObjectFunctionName;
    double MajorObjectFunctionResult;
    list<string> MinorObjectFunctionNames;
    list<double> MinorObjectFunctionResults;
    list<ObjectiveDDAModel*> ObjList;
    Vectord Originarray;                               //Record the objective function for partial derivative (the value before change)   
    bool HavePenalty;
    double PenaltyFactor;

    double MaxObj;                                //Record the historical maximum obj func
    double PreviousObj;                            //The previous obj
    int CutoffHold;
    Vectord MaxObjarray;                         //the individual objs for each model when the average obj is maximum(not necessaily the maximum individual objs)
    double epsilon_fix;
    double epsilon_tmp;                         //The epsilon used for calculation (can be different from the fixed input epsilon)
    bool HavePathRecord;
    bool HaveOriginHeritage;
    bool HaveAdjointHeritage;
    int Stephold;

    Vectord gradientsquare;                    //cumulative summation of gradients square. Used in Adagrad.
public:
    EvoDDAModel(list<string>* ObjectFunctionNames_, list<list<double>*>* ObjectParameters_, double epsilon_fix_, bool HavePathRecord_, bool HavePenalty_, bool HaveOriginHeritage_, bool HaveAdjointHeritage_, double PenaltyFactor_, string save_position_, CoreStructure* CStr_, list<DDAModel*> ModelList_);
    
    //functions used to calculate partial derivatives                                 
    tuple<Vectord, Vectorcd> devx_and_Adevxp(double epsilon, DDAModel* CurrentModel, ObjectiveDDAModel* objective, double origin);                       //partial derivative of obj to parameter and A to x times p
    tuple<Vectord, Vectorcd> devx_and_Adevxp_tmp(double epsilon, DDAModel* CurrentModel, ObjectiveDDAModel* objective, double origin);
    Vectorcd devp(double epsilon, DDAModel* CurrentModel, ObjectiveDDAModel* objective, double origin);                       //partial derivative of obj to P. Size of P

    void EvoOptimization(int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method, double start_num=0);
    void EvoOptimization(int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method, Vectord* V_, Vectord* S_);
    double CalTheObjForSingleStr(int MAX_ITERATION, double MAX_ERROR, int Name);                    //If you want to calculate the objective for single DDA structure.

    //The objective choosing function:
    ObjectiveDDAModel* ObjectiveFactory(string ObjectName, list<double> ObjectParameters, DDAModel* ObjDDAModel);

    double get_output_time();
    double L1Norm();
    Vectord gradients_filtered(Vectord gradients, int current_it, int Max_it);



};



class ObjectiveDDAModel {
public:
    bool Have_Devx;
    bool Have_Penalty;
    virtual void SingleResponse(int idx, bool deduction) = 0;
    virtual double GroupResponse() = 0;
    virtual void Reset() = 0;
    virtual double GetVal() = 0;
};


//Child classes for objective function


class ObjectivePointEDDAModel : public ObjectiveDDAModel {
private:
    double x;
    double y;
    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    double d;
    int N;
    Vectorcd* P;
    Vectori* R;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vectorcd E_sum;
    Vectorcd E_ext;
public:
    ObjectivePointEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectivePointListEDDAModel : public ObjectiveDDAModel {
private:
    Vectord x;
    Vectord y;
    Vectord z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    int PNum;
    double d;
    int N;
    Vectorcd* P;
    Vectori* R;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Matrixcd E_sum;
    Matrixcd E_ext;
public:
    ObjectivePointListEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectivePointIDDAModel : public ObjectiveDDAModel {
private:
    double x;
    double y;
    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    double d;
    int N;
    Vectorcd* P;
    Vectori* R;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vectorcd E_sum;
    Vectorcd E_ext;
public:
    ObjectivePointIDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveIntegratedEDDAModel : public ObjectiveDDAModel {
private:
    double d;
    int N;
    int Nx;
    int Ny;
    int Nz;
    Vectorcd* P;
    Vectorcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vectorcd E;
    double E_int;
    Vectori* R;
public:
    ObjectiveIntegratedEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveMidAvgEDDAModel : public ObjectiveDDAModel {
private:
    double d;
    int N;
    int Nx;
    int Ny;
    int Nz;
    Vectorcd* P;
    Vectorcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vectorcd E;
    double E_avg;
    double r;                 //radius of the middle regions (for 2D)
    Vectori* R;
    double centerx;
    double centery;
public:
    ObjectiveMidAvgEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class Objectivescattering2D : public ObjectiveDDAModel {
private:
    double d;
    int N;
    Vectorcd* P;
    Vectori* R;
    DDAModel* model;
    EvoDDAModel* evomodel;
    double K;
    double E0;
    int Paralength;
    list<Vectord> n_K_s_l;
    list<Vectorcd> PSum_l;
    list<Matrixd> FconstM_l;
    double ATUC;

public:
    Objectivescattering2D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
    double FTUCnsquareoversinal();
};

class Objectivereflect2D : public ObjectiveDDAModel {
private:
    double d;
    int N;
    Vectorcd* P;
    Vectori* R;
    DDAModel* model;
    EvoDDAModel* evomodel;
    double K;
    double E0;
    int Paralength;
    list<Vectord> n_K_s_l;
    list<Vectorcd> PSum_l;
    list<Matrixd> FconstM_l;
    double ATUC;

public:
    Objectivereflect2D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
    double FTUCnsquareoversinal();
};

class Objectivescattering0D : public ObjectiveDDAModel {
private:

    double d;
    int N;
    Vectorcd* P;
    Vectori* R;
    DDAModel* model;
    EvoDDAModel* evomodel;
    double K;
    double E0;
    int Paralength;
    list<Vectord> n_K_s_l;
    list<Vectorcd> PSum_l;
    list<Matrixd> FconstM_l;
    //double ATUC;
    //Matrixd FconstM;


public:
    Objectivescattering0D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
    double FTUCnsquare();
};

class ObjectiveAbs : public ObjectiveDDAModel {
private:
    double d;
    int N;
    int Nx;
    int Ny;
    int Nz;
    double K;
    double K3;
    double E0;
    Vectorcd* P;
    Vectorcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vectorcd E;
    double Cabs;
    Vectori* R;
public:
    ObjectiveAbs(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveAbsPartial : public ObjectiveDDAModel {
private:
    double d;
    int N;
    int Nx;
    int Ny;
    int Nz;
    double K;
    double K3;
    double E0;
    Vectorcd* P;
    Vectorcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vectorcd E;
    double Cabs;
    Vectori* R;
    set<int> integralpos;             //Only geometries with pixels in this set will be considered.
public:
    ObjectiveAbsPartial(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveAbsPartialzslice : public ObjectiveDDAModel {
private:
    double d;
    int N;
    int Nx;
    int Ny;
    int Nz;
    double K;
    double K3;
    double E0;
    Vectorcd* P;
    Vectorcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vectorcd E;
    double Cabs;
    Vectori* R;
    set<int> integralpos;             //Only geometries with pixels in this set will be considered.
    set<int> zslices{ 27, 28, 29 };
public:
    ObjectiveAbsPartialzslice(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveAbsbyfar : public ObjectiveDDAModel {
private:
    double d;
    int N;
    Vectorcd* P;
    Vectori* R;
    DDAModel* model;
    EvoDDAModel* evomodel;
    double K;
    double E0;
    int Paralength;
    list<Vectord> n_K_s_l;
    list<Vectorcd> PSum_l;
    list<Matrixd> FconstM_l;
    double ATUC;

public:
    ObjectiveAbsbyfar(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
    double FTUCnsquareoversinal();
};


class ObjectiveIntegrateEPartial : public ObjectiveDDAModel {
private:
    double d;
    int N;
    int Nx;
    int Ny;
    int Nz;
    double K;
    double K3;
    double E0;
    Vectorcd* P;
    Vectorcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vectorcd E;
    //double Cabs;
    double Total;
    Vectori* R;
    set<int> integralpos;             //Only geometries with pixels in this set will be considered.
public:
    ObjectiveIntegrateEPartial(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class FOMscattering2D{
private:
    Vectord FOMParameters;
    double d;
    int N;
    Vectorcd* P;
    Vectori* R;
    DDAModel* model;
    double K;
    int Paralength;
    list<Vectord> n_K_s_l;
    double ATUC;
    double E0;
    Vectord n_E0;
    //Matrixd FconstM;
    

public:
    FOMscattering2D(list<double> parameters, DDAModel* model_);
    list<double> GetVal();                                                        //Return list of far field abs(Es) at specified directions
    Vectorcd FTUC(Vectord n_K_s);
};

class FOMreflect2D {
private:
    Vectord FOMParameters;
    double d;
    int N;
    Vectorcd* P;
    Vectori* R;
    DDAModel* model;
    double K;
    int Paralength;
    list<Vectord> n_K_s_l;
    double ATUC;
    double E0;
    Vectord n_E0;
    //Matrixd FconstM;


public:
    FOMreflect2D(list<double> parameters, DDAModel* model_);
    list<double> GetVal();                                                        //Return list of far field abs(Es) at specified directions
    Vectorcd FTUC(Vectord n_K_s);
};

class FOMscattering0D {
private:

    double d;
    int N;
    Vectorcd* P;
    Vectori* R;
    DDAModel* model;
    double K;
    double E0;
    int Paralength;
    list<Vectord> n_K_s_l;
    //double ATUC;
    //Matrixd FconstM;


public:
    FOMscattering0D(list<double> parameters, DDAModel* model_);
    list<double> GetVal();                                                      //Return list of dCsca/dOmega at specified directions
    double FTUCnsquare(Vectord n_K_s);
};

/*
class FOMAbs {
private:

    double d;
    int N;
    Vectorcd* P;
    Vectori* R;
    DDAModel* model;
    double K;
    double E0;
    int Paralength;
    list<Vectord> n_K_s_l;
    //double ATUC;
    //Matrixd FconstM;


public:
    FOMAbs(list<double> parameters, DDAModel* model_);
    list<double> GetVal();                                                      //Return list of dCsca/dOmega at specified directions
    double FTUCnsquare(Vectord n_K_s);
};
*/










void Evo_Focus(SpacePara* spacepara_tmp, CoreStructure* CStr, DDAModel* TestModel, string save_position, int start_num, int max_evo,
    int min_num, int max_num, Vectord lower_bound, Vectord upper_bound, bool sym  //Parameters for focus generation
);
void Evo_single(string save_position, Vectori bind, Vectord l, int MAX_ITERATION_EVO, Vectord move_focus);




void eval_FOM(string name, DDAModel* TestModel, list<double> theta, list<double> phi);
void eval_FOM_2Dperiod(string name, DDAModel* TestModel, list<double> theta, list<double> phi);

#endif // PREPROCESSING_H_INCLUDED