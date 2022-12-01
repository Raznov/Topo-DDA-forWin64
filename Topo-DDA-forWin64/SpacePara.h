#ifndef TOPO_SPACEPARA_H_
#define TOPO_SPACEPARA_H_

#include <list>

#include "Space.h"
#include "filterOption.h"

struct WeightPara {
    double weight;
    int position;
};

class SpacePara {
private:
    Space* space;
    VectorXi geometry;                //3N dimension
    VectorXi geometryPara;            //N dimension. N=number of dipoles. Each position stores the para index in VectorXi Para : 0->Para[0]...
    vector<vector<int>> Paratogeometry;
    VectorXd Para;                    //P dimension. P=number of parameters. Same as Para_origin if Filter=False. Filtered and biased para if Filter=True.
    VectorXd Para_origin;             //Un-filtered, unbiased para. No use when Filter=False.
    VectorXd Para_filtered;           //Filtered but unbiased para. No use when Filter=False.
    MatrixXi scope;                   //[[xmin, xmax],[ymin, ymax],[zmin, zmax]]
    Vector3i bind;
    VectorXi FreeparatoPara;          //Position of free parameters inside Para. dimension<=P. FreeparatoPara[i] is the index of a free parameter inside Para.
    //vector<list<int>> Paratogeometry;  //P dimension. Each position stores a list of corresponding dipole index for parameter for this specific position.
    bool Filter;                      //True for with Filter. False for without Filter. Defualt as False for most initilizations.
    FilterOption* Filterstats;        //Only used when Filter=True
    vector<vector<WeightPara>> FreeWeight;
    vector<int> ParaDividePos;
    bool Periodic;
    int Lx;
    int Ly;

public:
    //SpacePara(Space* space_, string initial_diel, VectorXi geometry_, VectorXd diel_); //Can freeze part of the parameter space
    //SpacePara(Space* space_, Vector3i bind_, VectorXi* geometryPara_, VectorXd* Para_, VectorXi* FreeparatoPara_);
    //SpacePara(Space* space_, Vector3i bind_, string initial_diel); //l, center similar to bulk build in Structure class. Every 'bind' nearby dipoles correspond 
                                                                    //to 1 parameter in this bulk. bind=(2,2,2): 2*2*2; bind=(1,1,3):1*1*3
    //SpacePara(Space* space_, Vector3i bind_, string initial_diel1, string initial_diel2); //One constant layer at bottom. One design region on top. divide_pos is the divide in geometry array of the two parts.

    //SpacePara(Space* space_, Vector3i bind_, string initial_diel, VectorXi* geometryPara_);
    SpacePara(Space* space_, Vector3i bind_, string initial_diel_center, string initial_diel_ring, double r, string type);   //ONly for 2d cylinder or spheres. r is raidus/d.

    //SpacePara(Space* space_, Vector3i bind_, string initial_diel_background, list<string>* initial_diel_list, list<double>* r_list, list<Vector2d>* center_list);
    //Build 2d cylinders with diel in the list, rest of the diel is the backgroudn diel.

    SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2);  //random rect in a region with extruded 2D geometry
    SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2);  //random rect with 3D
    SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, VectorXi* geometryPara_);
    SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2, VectorXi* geometryPara_);

    SpacePara(Space* space_, Vector3i bind_, string initial_diel, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_);
    SpacePara(Space* space_, Vector3i bind_, vector<string> initial_diel_list, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_, bool Filter_, FilterOption* Filterstats_ = NULL, string symmetry = "None", vector<double> symaxis = vector<double>{ 0.0,0.0 }); //Should be exactly the same with the previous one except of Filter

    SpacePara(Space* space_, Vector3i bind_, vector<string> initial_diel_list, vector<double> BPara_, bool Filter_ = false, FilterOption* Filterstats_ = NULL, string symmetry = "None", vector<double> symaxis = vector<double>{ 0.0,0.0 });
    SpacePara(Space* space_, Vector3i bind_, VectorXi* InputGeo, VectorXd* Inputdiel, bool Filter_ = false, FilterOption* Filterstats_ = NULL, string symmetry = "None", vector<double> symaxis = vector<double>{ 0.0,0.0 }, bool Periodic_ = false, int Lx_ = 0, int Ly_ = 0);

    void ChangeBind(Vector3i bind_);                                  //Change bind number
    void ChangeFilter();
    VectorXi cut(VectorXi* big, VectorXi* smalll);

    Space* get_space();
    VectorXi get_geometry();
    VectorXi* get_geometryPara();
    VectorXd* get_Para();
    VectorXd* get_Para_origin();
    VectorXd* get_Para_filtered();
    Vector3i* get_bind();
    VectorXi* get_Free();
    //vector<list<int>>* get_Paratogeometry();
    bool get_Filter();
    FilterOption* get_Filterstats();
    vector<vector<WeightPara>>* get_FreeWeight();
    vector<int>* get_ParaDividePos();
    vector<vector<int>>* get_Paratogeometry();
};

#endif
