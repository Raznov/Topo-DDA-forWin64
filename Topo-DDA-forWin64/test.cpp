#include "definition.h"

int main() {

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = 83; Ny = 83; Nz = 13;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    Vector3d l;
    Vector3d center;
    l << 80.0, 80.0, 12.0;
    center << 40.0, 40.0, 6.0;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);



    S = S + s1;


    double d = 25;

    double lam = 500;

    double E0 = 1.0;

    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    double epsilon = 100;

    double focus = 350;   //nm       


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 30;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = "";
    
    Vector3d n_K;
    Vector3d n_E0;

    AProductCore Core(&S, d, lam, material);

    double angle = 30 * M_PI / 180;
    n_K << 0.0, 0.0, 1.0;
    n_E0 << 1.0, 0.0, 0.0;
    DDAModel Model0(&Core, n_K, E0, n_E0);
    n_K << sin(angle), 0.0, cos(angle);
    n_E0 << cos(angle), 0.0, -sin(angle);
    DDAModel Model1(&Core, n_K, E0, n_E0);
    n_K << -sin(angle), 0.0, cos(angle);
    n_E0 << cos(angle), 0.0, sin(angle);
    DDAModel Model2(&Core, n_K, E0, n_E0);

    angle = 20 * M_PI / 180;

    n_K << sin(angle), 0.0, cos(angle);
    n_E0 << cos(angle), 0.0, -sin(angle);
    DDAModel Model3(&Core, n_K, E0, n_E0);
    n_K << -sin(angle), 0.0, cos(angle);
    n_E0 << cos(angle), 0.0, sin(angle);
    DDAModel Model4(&Core, n_K, E0, n_E0);
    //n_K << 0.0, sin(angle), cos(angle);
    //n_E0 << 1.0, 0.0, 0.0;
    //DDAModel Model3(&Core, n_K, E0, n_E0);
    //n_K << 0.0, sin(angle), -cos(angle);
    //n_E0 << 1.0, 0.0, 0.0;
    //DDAModel Model4(&Core, n_K, E0, n_E0);

    list<DDAModel*> ModelList;
    ModelList.push_back(&Model0);
    ModelList.push_back(&Model1);
    ModelList.push_back(&Model2);
    ModelList.push_back(&Model3);
    ModelList.push_back(&Model4);

    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &Core, ModelList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}





