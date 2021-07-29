#include "definition.h"
#define PI 3.14159265


int main() {

    srand((unsigned)(time(0)));

    int Nx, Ny, Nz;
    Nx = 90; Ny = 90; Nz = 20;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;

    Vector3d center;
    Vector3d l;

    d = 20;

    center << Nx / 2, Ny / 2, Nz / 2;
    l << Nx - 1, Ny - 1, Nz - 1;

    Structure s1(S.get_total_space(), l, center);



    S = S + s1;

    double lam = 700;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");

   


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;


    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;

    Vector3i bind(1, 1, 20);
    int number = 10;
    double limitx1 = 4;
    double limitx2 = 30;
    double limity1 = 4;
    double limity2 = 30;

    SpacePara spacepara(&S, bind, number, limitx1, limitx2, limity1, limity2);


    CoreStructure CStr(&spacepara, d);
    AProductCore Core(&CStr, lam, material, "LDR");
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    string save_position=".\\SiO2-rects\\";

    TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel.update_E_in_structure();
    TestModel.solve_E();
    TestModel.output_to_file(save_position, 0, 0);
    CStr.output_to_file(save_position, 0);

    int num_model = 30;

    ofstream TotalTime;
    TotalTime.open(save_position+"TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    for (int i = 0; i <= num_model-1; i++) {
        SpacePara spacepara_tmp(&S, bind, number, limitx1, limitx2, limity1, limity2, spacepara.get_geometryPara());
        CStr.UpdateStr(&spacepara_tmp);
        CStr.output_to_file(save_position, i + 1);
        TestModel.UpdateAlpha();
        TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();
        TestModel.output_to_file(save_position, 0, i+1);
    }
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();





    return 0;

}



