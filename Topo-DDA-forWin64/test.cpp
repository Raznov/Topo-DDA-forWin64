#include "definition.h"
#define PI 3.14159265


int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();


    double ld = 375;
    double tTiO2 = 180;
    double d = 10;
    double disp = 155;

    Vector3d l;
    Vector3d center;
    l << ld / d, ld / d, tTiO2 / d;
    //l << 40.0, 40.0, 8.0;
    center << l(0) / 2, l(1) / 2, l(2) / 2;

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = round(l(0) + 1); Ny = round(l(1) + 1); Nz = round(l(2) + 1);
    cout << center << endl;
    //Nx = 23; Ny = 23; Nz = 10;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;

    //l << 20.0, 20.0, 9.0;
    //center << 10.0, 10.0, 4.5;
    Structure s1(S.get_total_space(), l, center);



    S = S + s1;

    Vector2d center1, center2;
    center1 << round(115 / d), round(115 / d);
    center2 << round(center1(0) + disp * cos(70 * PI / 180) / d), round(center1(1) + disp * sin(70 * PI / 180) / d);
    list<double> r_list{ 30 / d, 40 / d };
    list<Vector2d> center_list{ center1, center2 };
    list<string> initials_list{ "ZEROS", "ZEROS" };
    double r = 150 / d;

    Vector3i bind(1, 1, round(l(2)));
    //SpacePara spacepara(&S, bind, "ONES", "ZEROS", r);


    SpacePara spacepara(&S, bind, "ONES", &initials_list, &r_list, &center_list);


    double E0 = 1.0;


    double epsilon = 0.2;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    //double focus = (l(2) + 2) * d;
    //cout << focus << endl;

    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 500;

    list<string> ObjectFunctionNames{ "IntegratedE" };


    Vector3d n_K;
    Vector3d n_E0;

    n_K << 0.0, 0.0, -1.0;
    n_E0 << 1.0, 0.0, 0.0;


    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    double lam_min = 630;
    double lam_max = 700;
    double lam_interval = 1;
    int lam_num = round((lam_max - lam_min) / lam_interval) + 1;

    CoreStructure CStr(&spacepara, d);
    string save_position = ".\\dimer-wavelengthscan-d=10\\";
    string name = save_position + "IntegratedEwavedepend=" + to_string(d) + ".txt";

    CStr.output_to_file(save_position, 0);
    //S.show_something_about_Structures();
    ofstream fout(name);

    for (int i = 0; i <= lam_num; i++) {
        double lam = lam_min + i * lam_interval;
        Vector2cd material;
        material = Get_2_material("Air", "TiO2", lam, "nm");
        int m, n;
        double Lm, Ln;
        m = 50;
        n = 50;
        Lm = (Nx + 1) * d;
        Ln = (Ny + 1) * d;
        AProductCore Core(&CStr, lam, material, m, n, Lm, Ln, "FCD");
        DDAModel TestModel(&Core, n_K, E0, n_E0);
        TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();
        TestModel.output_to_file(save_position, 0, lam);

        VectorXcd* E_internal = TestModel.get_Einternal();
        VectorXi* R = CStr.get_R();
        N = CStr.get_N();
        double E_int = 0;
        for (int j = 0; j < N; j++) {
            if ((*R)(3 * j + 2) >= 0) {
                double E_sum_temp = 0;
                for (int k = 0; k < 3; k++) {
                    E_sum_temp += pow(abs((*E_internal)(3 * j + k)), 2);
                }

                E_int += pow(E_sum_temp, 3 / 2) * ((*(CStr.get_diel_old()))(j * 3));  
                //cout << "j: " << j << " " << E_int << endl;

            }
        }

        E_int = log(E_int * pow(d, 3));
        cout << "E_int" << E_int << endl;
        fout << lam << " " << E_int << endl;
    }
    fout.close();


    //AProductCore Core1(&CStr, lam(0), material, "LDR")






    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return 0;

}



