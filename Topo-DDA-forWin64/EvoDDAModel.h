#ifndef TOPO_EVO_H_
#define TOPO_EVO_H_

#include "DDAModel.h"
#include "ObjDDAModel.h"

using namespace std;
using namespace Eigen;

class EvoDDAModel {
private:
    double output_time;
    CoreStructure* CStr;
    vector<DDAModel*> allModel;                    //List of DDA models sharing the same AProductCore : "Core"
    int ModelNum;                                 //number of DDA model
    string save_position;

    vector<VectorXcd> PforOrigin;
    vector<VectorXcd> PforAdjoint;
    vector<VectorXcd> PforOriginMax;
    vector<VectorXcd> PforAdjointMax;

    vector<double> objPara;
    /*list<double> MajorObjectParameters;
    list<list<double>> MinorObjectParameters;*/

    string objName;
    /*string MajorObjectFunctionName;
    double MajorObjectFunctionResult;
    list<string> MinorObjectFunctionNames;
    list<double> MinorObjectFunctionResults;*/
    vector<ObjDDAModel*> allObj;
    VectorXd Originarray;                               //Record the Obj function for partial derivative (the value before change)   
    //bool HavePenalty;
    //double PenaltyFactor;

    double MaxObj;                                //Record the historical maximum obj func
    double PreviousObj;                            //The previous obj
    int CutoffHold;
    VectorXd MaxObjarray;                         //the individual objs for each model when the average obj is maximum(not necessaily the maximum individual objs)
    double epsilon_fix;
    double epsilon_tmp;                         //The epsilon used for calculation (can be different from the fixed input epsilon)
    bool HavePathRecord;
    bool HaveOriginHeritage;
    bool HaveAdjointHeritage;
    int Stephold;

    VectorXd gradientsquare;                    //cumulative summation of gradients square. Used in Adagrad.
public:
    EvoDDAModel(string objName_, vector<double> objPara_, double epsilon_fix_, bool HavePathRecord_, bool HaveOriginHeritage_, bool HaveAdjointHeritage_, string save_position_, CoreStructure* CStr_, vector<DDAModel*> allModel_);

    //functions used to calculate partial derivatives                                 
    tuple<VectorXd, VectorXcd> devx_and_Adevxp(double epsilon, DDAModel* CurrentModel, ObjDDAModel* Obj, double origin);                       //partial derivative of obj to parameter and A to x times p
    tuple<VectorXd, VectorXcd> devx_and_Adevxp_tmp(double epsilon, DDAModel* CurrentModel, ObjDDAModel* Obj, double origin);
    VectorXcd devp(double epsilon, DDAModel* CurrentModel, ObjDDAModel* Obj, double origin);                       //partial derivative of obj to P. Size of P

    void EvoOptimization(int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method, double start_num = 0);
    void EvoOptimization(int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method, VectorXd* V_, VectorXd* S_);
    double CalTheObjForSingleStr(int MAX_ITERATION, double MAX_ERROR, int Name);                    //If you want to calculate the Obj for single DDA structure.

    //The Obj choosing function:
    ObjDDAModel* ObjFactory(string ObjectName, vector<double> ObjectParameters, DDAModel* ObjDDAModel);

    double get_output_time();
    //double L1Norm();
    VectorXd gradients_filtered(VectorXd gradients, int current_it, int Max_it);



};

#endif