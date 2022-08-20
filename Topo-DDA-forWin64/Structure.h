#ifndef TOPO_STRUCTURE_H_
#define TOPO_STRUCTURE_H_

#include <vector>

#include "Eigen/Core"

using namespace std;
using namespace Eigen;
class Structure {
private:
    VectorXi geometry;
    bool para_cond;
public:
    //--------------------------------------Dont support dependent para build-------------------------------------------------------
    Structure(VectorXi* total_space, VectorXi* geometry_, bool para_cond_ = false);            //vanilla initialization

    //Sphere
    Structure(VectorXi* total_space, double r, Vector3d center, bool para_cond_ = false);  //r: actual radius/d. center: actual center/d. In charge of Sphere 
    Structure(VectorXi* total_space, double r, double h, Vector3d center, bool para_cond_ = false);

    //Circle
    //Structure(VectorXi *total_space, string initial_diel, double r, Vector3i center, Vector3i direction, int para_); //build a circle, direction is its normalized direction in Cart. coord.

    //From file
    //Structure(VectorXi* total_space, string FileName, int para_);                                                    //Read a structure from txt file

    //--------------------------------------Support dependent para build-------------------------------------------------------------
    //Bulk
    Structure(VectorXi* total_space, Vector3d l, Vector3d center, bool para_cond_ = false);    //Ractangular
    Structure(VectorXi* total_space, Vector3d l, Vector3d center, Structure* Str, bool para_cond_ = false);
    //Duplicate
    //Structure(VectorXi *total_space, Structure *s, Vector3i direction, int times, int para_);                     //Initializa a Structure by duplicating a existing structure along a certain direction for several times. Direction is normalized and can only be alone x, y or z.
    //The original structure is not included. original structure + new structure = times * original structure. If set para=2 and original para=1, then depend on origin str as geometry_dep. If para=2 and original para=2, will copy origin geometry_dep.
    //Structure(VectorXi* total_space, Structure* s, int dep_way);     //Special copy for setting up 2-fold and 4-fold symmetry dependence in xy plane. Para auto set to 2.
    //The original str at left down corner of the xy plane(left down corner is (0,0)). dep_way=1, 2, 3 corresponds to the other blocks in clock-wise. 


    //-------------------------------------Other member functions--------------------------------------------------------------------
    VectorXi* get_geometry();
    void cut(VectorXi* big, VectorXi* smalll);
    int get_geometry_size();
    //bool sym_or_not();
    //vector<string> get_sym_condition();
    //vector<double> get_sym_axis();
    bool para_or_not();
};
#endif
