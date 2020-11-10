#include "definition.h"

int main() {

    int Nx, Ny, Nz;
    Nx = 83; Ny = 83; Nz = 33;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);


    Structure s1(S.get_total_space(), "Lens", 0);
    S = S + s1;

    double d = 25;
    double focus = 1217;
    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon = 100;
    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 2;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ 40 * d,40 * d,focus };
    list<list<double>*> ObjectParameters{ &ObjectParameter };



    //list<double> ObjectParameter2{r(0), r(1), r(2), l(0)*d, l(1)*d, d};

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;

    string save_position = "";

    int Name = 0;

    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material);
    TestModel.output_to_file();
    double obj = TestModel.CalTheObjForSingleStr(MAX_ITERATION_DDA, MAX_ERROR, Name);
    cout << "Objective: " << obj << endl;
    //objectives << focus << " " << d << " " << obj << endl;










    return 0;

}




