#ifndef TOPO_CORESTRUCTURE_H_
#define TOPO_CORESTRUCTURE_H_

#include "Space.h"
#include "SpacePara.h"

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
    VectorXi R;                      //Position of dipoles. Both R and RResult are unitless, so need to time d to get real number.

    //VectorXi RDep;                   //Position of the dependent para points in space of the points in the same position as R
    //list<list<int>> PositionDep;    //First D has the D of para(para=1). Second D is the positions of other points dependent on the para in the 1stD.
    //VectorXi PositionPara;          //Position i(in R) of the parameters (3*i=x, 3*i+1=y, 3*i+2=z)
    //list<int> para_nums;
    //list<int> para_starts;
    //list<int> para_dep_nums;
    //list<int> para_dep_starts;
    VectorXd diel_old;                //The 0~1 version of diel, 3*N
    VectorXd diel_old_max;
public:
    CoreStructure(SpacePara* spacepara_, double d_);
    void UpdateStr(VectorXd step, int current_it, int Max_it);
    void UpdateStr(SpacePara* spacepara_);
    void UpdateStrSingle(int idx, double value);
    void output_to_file();
    void output_to_file(string save_position, int iteration, string mode = "normal");

    int get_N();
    int get_Nx();
    int get_Ny();
    int get_Nz();
    VectorXi* get_R();
    double get_d();
    SpacePara* get_spacepara();
    //list<list<int>>* get_PositionDep();
    //VectorXi* get_PositionPara();
    //list<int>* get_para_nums();
    //list<int>* get_para_starts();
    //list<int>* get_para_dep_nums();
    //list<int>* get_para_dep_starts();
    VectorXd* get_diel_old();
    VectorXd* get_diel_old_max();

};

#endif