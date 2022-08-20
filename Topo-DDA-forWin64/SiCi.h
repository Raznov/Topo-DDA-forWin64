#ifndef TOPO_SICI_H_
#define TOPO_SICI_H_

#include "Eigen/Core"

using namespace Eigen;

class SiCi {
public:
    int numberSi;
    int numberCi;
    double disSi;
    double disCi;
    VectorXd Si;
    VectorXd Ci;
    SiCi();
    double get_Si(double x);
    double get_Ci(double y);

};

#endif
