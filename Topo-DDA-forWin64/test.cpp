#include "definition.h"

int main() {

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = 83; Ny = 83; Nz = 17;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    Vector3d l;
    Vector3d center;
    l << 80.0, 80.0, 16.0;
    center << 40.0, 40.0, 8.0;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);



    S = S + s1;


    double d = 25;

    double lam = 500;

    double E0 = 1.0;

    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    double epsilon = 100;

    double focus = 450;   //nm       


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 1;

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

    CoreStructure CStr(&S, d);
    AProductCore Core(&CStr, lam, material);

    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo("AngleInfo.txt");
    ofstream nEInfo("nEInfo.txt");

    int theta_num = 2;
    //int phi_num = 36;
    VectorXd theta(theta_num);
    //VectorXd phi=VectorXd::Zero(phi_num);
    theta << 0, 10;
    //for (int i = 0; i <= phi_num - 1; i++) {
    //    phi(i) = i * 360 / phi_num;
    //}
    VectorXd phi(2);
    phi << 0, 10;
    int phi_num = 2;

    for (int i = 0; i <= theta_num - 1; i++) {
        for (int j = 0; j <= phi_num - 1; j++) {
            if (theta(i) != 0) {
                double theta_tmp = theta(i) * M_PI / 180;
                double phi_tmp = phi(j) * M_PI / 180;
                n_K << sin(theta_tmp) * cos(phi_tmp), sin(theta_tmp)* sin(phi_tmp), cos(theta_tmp);
                n_E0 = nEPerpinXZ(theta_tmp, phi_tmp);
                if (CheckPerp(n_E0, n_K) == false) {
                    cout << "----------------------------------------theta" << theta[i] << "phi" << phi[j] << "Not perpendicular---------------------------------------" << endl;
                }
                AngleInfo << theta[i] << endl;
                AngleInfo << phi[j] << endl;
                nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                DDAModel Model(&Core, n_K, E0, n_E0);
                ModelList.push_back(Model);
            }
        }
    }

    double theta_tmp = 0 * M_PI / 180;
    double phi_tmp = 0 * M_PI / 180;
    n_K << sin(theta_tmp) * cos(phi_tmp), sin(theta_tmp)* sin(phi_tmp), cos(theta_tmp);
    n_E0 = nEPerpinXZ(theta_tmp, phi_tmp);
    if (CheckPerp(n_E0, n_K) == false) {
        cout << "----------------------------------------theta" << 0 << "phi" << 0 << "Not perpendicular---------------------------------------" << endl;
    }
    AngleInfo << 0.0 << endl;
    AngleInfo << 0.0 << endl;
    nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
    DDAModel Model(&Core, n_K, E0, n_E0);
    ModelList.push_back(Model);

    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;

    list<DDAModel>::iterator it = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it));
        it++;
    }


    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &CStr, ModelpointerList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}




