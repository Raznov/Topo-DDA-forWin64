#ifndef TOPO_TOOLS_H_
#define TOPO_TOOLS_H_
#define _USE_MATH_DEFINES

#include <complex>
#include <list>
#include <vector>

#include "Eigen/Core"


using namespace std;
using namespace Eigen;

complex<double> Get_Alpha(double lam, double K, double d, complex<double> diel, Vector3d n_E0, Vector3d n_K);
complex<double> Get_Alpha_FCD(double lam, double K, double d, complex<double> diel);
complex<double> Get_material(string mat, double wl, string unit);                  //name of mat to get its diel function at certain wavlength              
Vector2cd Get_2_material(string sub, string mat, double wl, string unit);          //a wrapper for Get_material
VectorXcd Get_X_material(list<string> mat_l, double wl, string unit);
double Average(VectorXcd* E, int N, double exponent);
double Get_Max(VectorXcd* E, int N);
double G(VectorXcd* E, int N, double exponent, double E0);
//ArrayXcd FFT(int nx,int ny,int nz,ArrayXcd in,int _direction);
VectorXi build_a_bulk(int Nx, int Ny, int Nz);
bool CheckPerp(Vector3d v1, Vector3d v2);
Vector3d nEPerpinXZ(double theta, double phi);
MatrixXi find_scope_3_dim(VectorXi* x);
VectorXd initial_diel_func(string initial_diel, int N);
double initial_diel_func(string initial_diel);
double initial_diel_func(double initial_diel);
list<double> makelist(double start, double end, double interval);
list<double> makelist(double start, double end, int number);
pair<VectorXi, VectorXd> InputInitial(string open_position, string model_label);
list<string> ReadMat(string input);
vector<double> ReadLam(string input);

double exp_update(const double x, const double x_max, const double y_min, const double y_max);
double piecewise_update(const double x, const double x_max, const double y_min, const double y_max);
double linear_update(const double x, const double x_max, const double y_min, const double y_max);

//int makedirect(string name);
vector<string> splitInputStr(string input, string delimiter);
pair<VectorXi, VectorXd> getInputStr(string pathCommonData, string pathPara);
tuple<int, int, int> getInputNs(string pathCommonData);
#endif
