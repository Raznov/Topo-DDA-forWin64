

//DDAModel verification; updated after adding FCD and binding; updated on 2021-6-20
int main() {

    int Nx, Ny, Nz;
    Nx = 56; Ny = 56; Nz = 56;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;
    double r;
    Vector3d center;

    d = 2;
    r = 55 / d;
    center << Nx / 2, Ny / 2, Nz / 2;

    Structure s1(S.get_total_space(), r, center);



    S = S + s1;

    double ratioESItoG = 1 / (2.998 * pow(10, 4));
    double rationmtocm = 1 / (pow(10, 7));
    double ratioESItoGnm = (1 / (2.998 * pow(10, 4))) / (pow(10, 7));



    double lam = 525;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");

    double epsilon = 100;

    double focus = 350;   //nm       


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 100;

    //d = rationmtocm * d;
    //lam = rationmtocm * lam;
    //E0 = ratioESItoGnm * E0;
    //focus = rationmtocm * focus;


    list<string> ObjectFunctionNames{ "PointEratio" };

    double exponent = 2;
    double ratio = 4;

    //list<double> ObjectParameter1{ focus, exponent, ratio };
    list<double> ObjectParameter2{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter2 };
    string save_position = "";

    Vector3i bind(1, 1, 1);
    //SpacePara spacepara(&S, bind, "ONES", "ZEROS", r);

    SpacePara spacepara(&S, bind, "ONES");

    CoreStructure CStr(&spacepara, d);
    AProductCore Core(&CStr, lam, material, "LDR");
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel.update_E_in_structure();
    TestModel.solve_E();
    TestModel.output_to_file();

    double x, y, z;
    y = center(1) * d;
    z = center(2) * d;


    N = Core.get_N();
    VectorXcd* P = TestModel.get_P();
    VectorXi* R = Core.get_R();
    double K = 2 * M_PI / lam;
    Vector3cd E_sum = Vector3cd::Zero();
    Vector3cd E_ext = Vector3cd::Zero();

    string name = "DDAModelEfieldscand=" + to_string(d) + ".txt";
    ofstream fout(name);
    for (x = center(0) * d; x <= center(0) * d + 350; x += 2) {
        cout << "x" << x - center(0) * d << endl;
        E_ext(0) = E0 * n_E0(0) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
        E_ext(1) = E0 * n_E0(1) * (cos(K * (n_K(0) * y + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
        E_ext(2) = E0 * n_E0(2) * (cos(K * (n_K(0) * z + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
        E_sum(0) = E_ext(0);
        E_sum(1) = E_ext(1);
        E_sum(2) = E_ext(2);
        for (int i = 0; i <= N - 1; i++) {
            double rx = x - d * (*R)(3 * i);                  //R has no d in it, so needs to time d
            double ry = y - d * (*R)(3 * i + 1);
            double rz = z - d * (*R)(3 * i + 2);
            Matrix3cd A = Core.A_dic_generator(rx, ry, rz);
            E_sum(0) -= (A(0, 0) * (*P)(3 * i) + A(0, 1) * (*P)(3 * i + 1) + A(0, 2) * (*P)(3 * i + 2));
            E_sum(1) -= (A(1, 0) * (*P)(3 * i) + A(1, 1) * (*P)(3 * i + 1) + A(1, 2) * (*P)(3 * i + 2));
            E_sum(2) -= (A(2, 0) * (*P)(3 * i) + A(2, 1) * (*P)(3 * i + 1) + A(2, 2) * (*P)(3 * i + 2));
        }
        double normE = E_sum.norm();
        cout << "E" << normE << endl;
        fout << x - center(0) * d << " " << normE << endl;
    }
    fout.close();

    //EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material, AMatrixMethod);

    //TestModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}

//Template for EvoModel test
int main() {

    /*
    int Nx,Ny,Nz;
    Nx=81;Ny=81;Nz=31;
    int N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Structure s1(S.get_total_space(), "Lens", 0);
    //Structure s0(S.get_total_space(), "0", 0);
    */


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
    //l << 50.0, 50.0, 15.0;
    //center << 25.0, 25.0, 7.5;
    l << 80.0, 80.0, 12.0;
    center << 40.0, 40.0, 6.0;
    //l << 79.0, 79.0, 11.0;
    //center << 39.5, 39.5, 5.5;
    direction << 0, 0, -1;
    int times = 15;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);
    //Structure s2(S.get_total_space(), &s1, 1);
    //Structure s3(S.get_total_space(), &s1, 2);
    //Structure s4(S.get_total_space(), &s1, 3);
    //Structure s5(S.get_total_space(), &s1, direction, times, 2);
    //Structure s6(S.get_total_space(), &s2, direction, times, 2);
    //Structure s7(S.get_total_space(), &s3, direction, times, 2);
    //Structure s8(S.get_total_space(), &s4, direction, times, 2);



    S = S + s1;
    //S = S + s2;
    //S = S + s3;
    //S = S + s4;
    //S = S + s5;
    //S = S + s6;
    //S = S + s7;
    //S = S + s8;
    //S=S+s0;
    //S.show_something_about_Structures();


    double d = 25;

    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon = 100;

    double focus = 350;   //nm       

    //Vector3d r;
    //r<<center(0)*d, center(1)*d, focus;


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 200;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    //list<double> ObjectParameter1{ focus, exponent, ratio };
    list<double> ObjectParameter2{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter2 };
    string save_position = "";


    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material);

    TestModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}

//Template for EvoDDAModel test
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
    //l << 50.0, 50.0, 15.0;
    //center << 25.0, 25.0, 7.5;
    l << 80.0, 80.0, 12.0;
    center << 40.0, 40.0, 6.0;
    //l << 79.0, 79.0, 11.0;
    //center << 39.5, 39.5, 5.5;
    direction << 0, 0, -1;
    int times = 15;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);
    //Structure s2(S.get_total_space(), &s1, 1);
    //Structure s3(S.get_total_space(), &s1, 2);
    //Structure s4(S.get_total_space(), &s1, 3);
    //Structure s5(S.get_total_space(), &s1, direction, times, 2);
    //Structure s6(S.get_total_space(), &s2, direction, times, 2);
    //Structure s7(S.get_total_space(), &s3, direction, times, 2);
    //Structure s8(S.get_total_space(), &s4, direction, times, 2);



    S = S + s1;
    //S = S + s2;
    //S = S + s3;
    //S = S + s4;
    //S = S + s5;
    //S = S + s6;
    //S = S + s7;
    //S = S + s8;
    //S=S+s0;
    //S.show_something_about_Structures();


    double d = 25;

    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon = 100;

    double focus = 350;   //nm       

    //Vector3d r;
    //r<<center(0)*d, center(1)*d, focus;


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 10;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    //list<double> ObjectParameter1{ focus, exponent, ratio };
    list<double> ObjectParameter2{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter2 };
    string save_position = "";

    AProductCore Core(&S, d, lam, material);
    DDAModel Model1(&Core, n_K, E0, n_E0);

    list<DDAModel*> ModelList;
    ModelList.push_back(&Model1);

    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &Core, ModelList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}

//Template for angle sweep
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
    int MAX_ITERATION_EVO = 100;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = "";
    double angle = 10 * M_PI / 180;
    Vector3d n_K;
    Vector3d n_E0;

    AProductCore Core(&S, d, lam, material);

    n_K << 0.0, 0.0, 1.0;
    n_E0 << 1.0, 0.0, 0.0;
    DDAModel Model0(&Core, n_K, E0, n_E0);
    n_K << sin(angle), 0.0, cos(angle);
    n_E0 << cos(angle), 0.0, -sin(angle);
    DDAModel Model1(&Core, n_K, E0, n_E0);
    n_K << -sin(angle), 0.0, cos(angle);
    n_E0 << cos(angle), 0.0, sin(angle);
    DDAModel Model2(&Core, n_K, E0, n_E0);
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
    //ModelList.push_back(&Model3);
    //ModelList.push_back(&Model4);

    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &Core, ModelList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}

//Template for traditional lens

int main() {

    int Nx, Ny, Nz;
    Nx = 83; Ny = 83; Nz = 34;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);


    Structure s1(S.get_total_space(), "Lens", 0);
    S = S + s1;

    double d = 25;
    double focus = 1600;
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

//Template for list angle sweep

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
    int MAX_ITERATION_EVO = 50;

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

    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    double theta[4] = { 0,10,20,30 };
    double phi[6] = { 0,10,170,180,190,350 };
    for (int i = 0; i <= 3; i++) {
        for (int j = 0; j <= 3; j++) {
            n_K << sin(theta[i]) * cos(phi[j]), sin(theta[i])* sin(phi[j]), cos(theta[i]);
            n_E0 << cos(theta[i]) * cos(phi[j]), cos(theta[i])* sin(phi[j]), -sin(theta[i]);
            DDAModel Model(&Core, n_K, E0, n_E0);
            ModelList.push_back(Model);
        }
    }

    list<DDAModel>::iterator it = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it));
        it++;
    }


    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &Core, ModelpointerList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}

//test xz perp

int main() {
    Vector3d n_K;
    Vector3d n_E0;
    double theta[4] = { 0,10,20,30 };
    double phi[6] = { 0,10,170,180,190,350 };
    for (int i = 0; i <= 3; i++) {
        for (int j = 0; j <= 3; j++) {
            cout << "theta" << theta[i] << "phi" << phi[j] << endl;
            double theta_tmp = theta(i) * M_PI / 180;
            double phi_tmp = phi(j) * M_PI / 180;
            n_K << sin(theta_tmp) * cos(phi_tmp), sin(theta_tmp)* sin(phi_tmp), cos(theta_tmp);
            n_E0 = nEPerpinXZ(theta_tmp, phi_tmp);
            cout << "n_K" << n_K << endl;
            cout << "n_E0" << n_E0 << endl;

        }
    }
}

//theta and phi sweep

int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();



    Vector3d l;
    Vector3d center;
    l << 80.0, 80.0, 16.0;
    center << l(0) / 2, l(1) / 2, l(2) / 2;

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);
    cout << center << endl;
    //Nx = 23; Ny = 23; Nz = 10;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;

    //l << 20.0, 20.0, 9.0;
    //center << 10.0, 10.0, 4.5;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);



    S = S + s1;


    double d = 25;


    double E0 = 1.0;


    double epsilon = 100;

    double focus = (l(2) + 2) * d;   //nm       
    cout << focus << endl;

    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 50;

    list<string> ObjectFunctionNames{ "PointI" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = true;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = "";

    Vector3d n_K;
    Vector3d n_E0;



    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo("AngleInfo.txt");
    ofstream nEInfo("nEInfo.txt");

    int theta_num = 4;
    VectorXd theta(theta_num);
    theta << 0, 10, 20, 30;
    int phi_num = 10;
    VectorXd phi(phi_num);
    phi << 0, 10, 20, 160, 170, 180, 190, 200, 340, 350;
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;

    CoreStructure CStr(&S, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "SiO2", lam(0), "nm");
    AProductCore Core1(&CStr, lam(0), material);
    //material = Get_2_material("Air", "SiO2", lam(1), "nm");
    //AProductCore Core2(&CStr, lam(1), material);
    //material = Get_2_material("Air", "SiO2", lam(2), "nm");
    //AProductCore Core3(&CStr, lam(2), material);
    CorePointList.push_back(&Core1);
    //CorePointList.push_back(&Core2);
    //CorePointList.push_back(&Core3);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */
    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {
        AProductCore* Core = (*it);
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
                    if (k == 0) {
                        AngleInfo << theta[i] << endl;
                        AngleInfo << phi[j] << endl;
                        nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                    }
                    DDAModel Model(Core, n_K, E0, n_E0);
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
        if (k == 0) {
            AngleInfo << 0.0 << endl;
            AngleInfo << 0.0 << endl;
            nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
        }
        DDAModel Model(Core, n_K, E0, n_E0);
        ModelList.push_back(Model);

        it++;
    }




    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;

    list<DDAModel>::iterator it1 = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it1));
        it1++;
    }


    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, HaveOriginHeritage, HaveAdjointHeritage, PenaltyFactor, save_position, &CStr, ModelpointerList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return 0;

}

//theta and phi and lam sweep

int main() {

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = 83; Ny = 83; Nz = 17;
    //Nx = 23; Ny = 23; Nz = 10;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    Vector3d l;
    Vector3d center;
    l << 80.0, 80.0, 16.0;
    center << 40.0, 40.0, 8.0;
    //l << 20.0, 20.0, 9.0;
    //center << 10.0, 10.0, 4.5;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);



    S = S + s1;


    double d = 25;


    double E0 = 1.0;


    double epsilon = 100;

    double focus = 450;   //nm       


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



    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo("AngleInfo.txt");
    ofstream nEInfo("nEInfo.txt");

    int theta_num = 4;
    VectorXd theta(theta_num);
    theta << 0, 10, 20, 30;
    int phi_num = 10;
    VectorXd phi(phi_num);
    phi << 0, 10, 20, 160, 170, 180, 190, 200, 340, 350;
    int lam_num = 3;
    VectorXd lam(lam_num);
    lam << 450, 500, 550;

    CoreStructure CStr(&S, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "SiO2", lam(0), "nm");
    AProductCore Core1(&CStr, lam(0), material);
    material = Get_2_material("Air", "SiO2", lam(1), "nm");
    AProductCore Core2(&CStr, lam(1), material);
    material = Get_2_material("Air", "SiO2", lam(2), "nm");
    AProductCore Core3(&CStr, lam(2), material);
    CorePointList.push_back(&Core1);
    CorePointList.push_back(&Core2);
    CorePointList.push_back(&Core3);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */
    for (int k = 0; k <= lam_num - 1; k++) {
        list<AProductCore*>::iterator it = CorePointList.begin();
        AProductCore* Core = (*it);
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
                    if (k == 0) {
                        AngleInfo << theta[i] << endl;
                        AngleInfo << phi[j] << endl;
                        nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                    }
                    DDAModel Model(Core, n_K, E0, n_E0);
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
        if (k == 0) {
            AngleInfo << 0.0 << endl;
            AngleInfo << 0.0 << endl;
            nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
        }
        DDAModel Model(Core, n_K, E0, n_E0);
        ModelList.push_back(Model);

        it++;
    }




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

//with FCD
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();



    Vector3d l;
    Vector3d center;
    l << 80.0, 80.0, 16.0;
    center << l(0) / 2, l(1) / 2, l(2) / 2;

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);
    cout << center << endl;
    //Nx = 23; Ny = 23; Nz = 10;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;

    //l << 20.0, 20.0, 9.0;
    //center << 10.0, 10.0, 4.5;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);



    S = S + s1;


    double d = 25;


    double E0 = 1.0;


    double epsilon = 100;

    double focus = (l(2) + 2) * d;   //nm       
    cout << focus << endl;

    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 50;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = true;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = "";

    Vector3d n_K;
    Vector3d n_E0;



    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo("AngleInfo.txt");
    ofstream nEInfo("nEInfo.txt");

    int theta_num = 1;
    VectorXd theta(theta_num);
    theta << 0;
    int phi_num = 1;
    VectorXd phi(phi_num);
    phi << 0;
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;

    CoreStructure CStr(&S, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "SiO2", lam(0), "nm");
    AProductCore Core1(&CStr, lam(0), material, "FCD");
    //material = Get_2_material("Air", "SiO2", lam(1), "nm");
    //AProductCore Core2(&CStr, lam(1), material);
    //material = Get_2_material("Air", "SiO2", lam(2), "nm");
    //AProductCore Core3(&CStr, lam(2), material);
    CorePointList.push_back(&Core1);
    //CorePointList.push_back(&Core2);
    //CorePointList.push_back(&Core3);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */
    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {
        AProductCore* Core = (*it);
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
                    if (k == 0) {
                        AngleInfo << theta[i] << endl;
                        AngleInfo << phi[j] << endl;
                        nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                    }
                    DDAModel Model(Core, n_K, E0, n_E0);
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
        if (k == 0) {
            AngleInfo << 0.0 << endl;
            AngleInfo << 0.0 << endl;
            nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
        }
        DDAModel Model(Core, n_K, E0, n_E0);
        ModelList.push_back(Model);

        it++;
    }




    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;

    list<DDAModel>::iterator it1 = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it1));
        it1++;
    }


    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, HaveOriginHeritage, HaveAdjointHeritage, PenaltyFactor, save_position, &CStr, ModelpointerList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return 0;

}






//2021-3 with bind
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();



    Vector3d l;
    Vector3d center;
    l << 80.0, 80.0, 16.0;
    //l << 40.0, 40.0, 8.0;
    center << l(0) / 2, l(1) / 2, l(2) / 2;

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);
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

    Vector3i bind(1, 1, 1);
    SpacePara spacepara(&S, bind, "ONES");

    double d = 25;


    double E0 = 1.0;


    double epsilon = 10;

    //double focus = (l(2) + 2) * d;   //nm       
    double focus = (l(2) + 2) * d;
    cout << focus << endl;

    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 200;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = true;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = "";

    Vector3d n_K;
    Vector3d n_E0;



    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo("AngleInfo.txt");
    ofstream nEInfo("nEInfo.txt");

    int theta_num = 1;
    VectorXd theta(theta_num);
    theta << 0;
    int phi_num = 1;
    VectorXd phi(phi_num);
    phi << 0;
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;

    CoreStructure CStr(&spacepara, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "SiO2", lam(0), "nm");
    //AProductCore Core1(&CStr, lam(0), material, "LDR");

    int m, n;
    double Lm, Ln;
    m = 10;
    n = 10;
    Lm = l(0) * d * 1.25;
    Ln = l(1) * d * 1.25;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "LDR");
    //material = Get_2_material("Air", "2.5", lam(1), "nm");
    //AProductCore Core2(&CStr, lam(1), material, "LDR");
    //material = Get_2_material("Air", "2.5", lam(2), "nm");
    //AProductCore Core3(&CStr, lam(2), material, "LDR");
    CorePointList.push_back(&Core1);
    //CorePointList.push_back(&Core2);
    //CorePointList.push_back(&Core3);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {
        AProductCore* Core = (*it);
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
                    if (k == 0) {
                        AngleInfo << theta[i] << endl;
                        AngleInfo << phi[j] << endl;
                        nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                    }
                    DDAModel Model(Core, n_K, E0, n_E0);
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
        if (k == 0) {
            AngleInfo << 0.0 << endl;
            AngleInfo << 0.0 << endl;
            nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
        }
        DDAModel Model(Core, n_K, E0, n_E0);
        ModelList.push_back(Model);

        it++;
    }




    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;

    list<DDAModel>::iterator it1 = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it1));
        it1++;
    }


    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, HaveOriginHeritage, HaveAdjointHeritage, PenaltyFactor, save_position, &CStr, ModelpointerList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return 0;

}

//2021-4-8
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();



    Vector3d l;
    Vector3d center;
    l << 21.0, 21.0, 9.0;
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
    double d = 15;
    double r = 150 / d;

    Vector3i bind(1, 1, 10);
    SpacePara spacepara(&S, bind, "ONES", "ZEROS", r);

    //SpacePara spacepara(&S, bind, "RANDOM");


    double E0 = 1.0;


    double epsilon = 0.2;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    double focus = (l(2) + 2) * d;
    cout << focus << endl;

    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 88;

    list<string> ObjectFunctionNames{ "IntegratedE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ 70.0 };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 1;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = ".\\p330-lam542-beta8-TiO2-InE-circle-withinter\\";

    Vector3d n_K;
    Vector3d n_E0;



    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");

    int theta_num = 1;
    VectorXd theta(theta_num);
    theta << 0;
    int phi_num = 1;
    VectorXd phi(phi_num);
    phi << 0;
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 542;

    CoreStructure CStr(&spacepara, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "TiO2", lam(0), "nm");
    //AProductCore Core1(&CStr, lam(0), material, "LDR");

    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 22 * d;
    Ln = 22 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");
    //AProductCore Core1(&CStr, lam(0), material, "FCD");
    //material = Get_2_material("Air", "2.5", lam(1), "nm");
    //AProductCore Core2(&CStr, lam(1), material, "LDR");
    //material = Get_2_material("Air", "2.5", lam(2), "nm");
    //AProductCore Core3(&CStr, lam(2), material, "LDR");
    CorePointList.push_back(&Core1);
    //CorePointList.push_back(&Core2);
    //CorePointList.push_back(&Core3);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {
        AProductCore* Core = (*it);
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
                    if (k == 0) {
                        AngleInfo << theta[i] << endl;
                        AngleInfo << phi[j] << endl;
                        nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                    }
                    DDAModel Model(Core, n_K, E0, n_E0);
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
        if (k == 0) {
            AngleInfo << 0.0 << endl;
            AngleInfo << 0.0 << endl;
            nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
        }
        DDAModel Model(Core, n_K, E0, n_E0);
        ModelList.push_back(Model);

        it++;
    }


    cout << "nK" << n_K << endl;
    cout << "nE0" << n_E0 << endl;

    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;

    list<DDAModel>::iterator it1 = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it1));
        it1++;
    }


    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, HaveOriginHeritage, HaveAdjointHeritage, PenaltyFactor, save_position, &CStr, ModelpointerList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return 0;

}

//2021-4-8, USE this to regenerate the doubl arrow structure. Set integration objective function for z>=8 because it is 2layer.
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();



    Vector3d l;
    Vector3d center;
    l << 21.0, 21.0, 9.0;
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
    double d = 15;
    double r = 150 / d;

    Vector3i bind(1, 1, 10);
    SpacePara spacepara(&S, bind, "ONES", "ZEROS", r);

    //SpacePara spacepara(&S, bind, "RANDOM");


    double E0 = 1.0;


    double epsilon = 0.2;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    double focus = (l(2) + 2) * d;
    cout << focus << endl;

    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 88;

    list<string> ObjectFunctionNames{ "IntegratedE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ 70.0 };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 1;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = ".\\p330-lam542-beta7-TiO2-InE-2layer-circle\\";

    Vector3d n_K;
    Vector3d n_E0;



    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");

    int theta_num = 1;
    VectorXd theta(theta_num);
    theta << 0;
    int phi_num = 1;
    VectorXd phi(phi_num);
    phi << 0;
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 542;

    CoreStructure CStr(&spacepara, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "TiO2", lam(0), "nm");
    //AProductCore Core1(&CStr, lam(0), material, "LDR");

    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 22 * d;
    Ln = 22 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");
    //AProductCore Core1(&CStr, lam(0), material, "FCD");
    //material = Get_2_material("Air", "2.5", lam(1), "nm");
    //AProductCore Core2(&CStr, lam(1), material, "LDR");
    //material = Get_2_material("Air", "2.5", lam(2), "nm");
    //AProductCore Core3(&CStr, lam(2), material, "LDR");
    CorePointList.push_back(&Core1);
    //CorePointList.push_back(&Core2);
    //CorePointList.push_back(&Core3);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {
        AProductCore* Core = (*it);
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
                    if (k == 0) {
                        AngleInfo << theta[i] << endl;
                        AngleInfo << phi[j] << endl;
                        nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                    }
                    DDAModel Model(Core, n_K, E0, n_E0);
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
        if (k == 0) {
            AngleInfo << 0.0 << endl;
            AngleInfo << 0.0 << endl;
            nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
        }
        DDAModel Model(Core, n_K, E0, n_E0);
        ModelList.push_back(Model);

        it++;
    }


    cout << "nK" << n_K << endl;
    cout << "nE0" << n_E0 << endl;

    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;

    list<DDAModel>::iterator it1 = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it1));
        it1++;
    }


    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, HaveOriginHeritage, HaveAdjointHeritage, PenaltyFactor, save_position, &CStr, ModelpointerList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return 0;

}

//intermediate start from DDA
int main() {
    ofstream TotalTime;

    high_resolution_clock::time_point t_start = high_resolution_clock::now();


    string open_position = ".\\p330-lam542-beta8-TiO2-InE-circle-withinter\\";
    string name1 = open_position + "Model_output\\CoreStructure87.txt";

    ifstream fin1(name1);
    int Nx, Ny, Nz;
    int Ntmp;
    fin1 >> Nx;
    fin1 >> Ny;
    fin1 >> Nz;
    fin1 >> Ntmp;
    VectorXi geometry = VectorXi::Zero(3 * Ntmp);
    for (int i = 0; i <= Ntmp - 1; i++) {
        fin1 >> geometry(3 * i);
        fin1 >> geometry(3 * i + 1);
        fin1 >> geometry(3 * i + 2);
    }
    for (int i = 0; i <= Ntmp - 1; i++) {
        double tmp;
        fin1 >> tmp;
        fin1 >> tmp;
        fin1 >> tmp;
    }
    double d;
    fin1 >> d;
    fin1.close();


    string name2, name3, name4;
    name2 = open_position + "Model_output\\87geometryPara.txt";
    name3 = open_position + "Model_output\\87Para.txt";
    name4 = open_position + "Model_output\\87bind.txt";
    ifstream fin2(name2), fin3(name3), fin4(name4);

    VectorXi geometryPara = VectorXi::Zero(Ntmp);
    int parasize;
    fin3 >> parasize;
    VectorXd Para = VectorXd::Zero(parasize);
    Vector3i bind = Vector3i::Zero();
    for (int i = 0; i <= Ntmp - 1; i++) {
        fin2 >> geometryPara(i);
    }
    for (int i = 0; i <= parasize - 1; i++) {
        fin3 >> Para(i);
    }
    for (int i = 0; i <= 2; i++) {
        fin4 >> bind(i);
    }
    fin2.close();
    fin3.close();
    fin4.close();

    string name5, name6;
    name5 = open_position + "Model_output\\87V.txt";
    name6 = open_position + "Model_output\\87s.txt";
    ifstream fin5(name5), fin6(name6);

    VectorXd VforAdam = VectorXd::Zero(parasize);
    VectorXd SforAdam = VectorXd::Zero(parasize);
    for (int i = 0; i <= parasize - 1; i++) {
        fin5 >> VforAdam(i);
        fin6 >> SforAdam(i);
    }
    fin5.close();
    fin6.close();

    //Vector3d l;
    //Vector3d center;
    //l << 21.0, 21.0, 9.0;
    //l << 40.0, 40.0, 8.0;
    //center << l(0) / 2, l(1) / 2, l(2) / 2;

    //int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    //Nx = round(l(0) + 1); Ny = round(l(1) + 1); Nz = round(l(2) + 1);
    //cout << center << endl;
    //Nx = 23; Ny = 23; Nz = 10;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;

    //l << 20.0, 20.0, 9.0;
    //center << 10.0, 10.0, 4.5;
    //Structure s1(S.get_total_space(), l, center);
    Structure s1(S.get_total_space(), &geometry);


    S = S + s1;
    //double d = 15;
    //double r = 150 / d;

    //Vector3i bind(1, 1, 10);
    SpacePara spacepara(&S, bind, &geometryPara, &Para);
    //SpacePara spacepara(&S, bind, "ONES", "ZEROS", r);

    //SpacePara spacepara(&S, bind, "RANDOM");


    double E0 = 1.0;


    double epsilon = 0.2;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    //double focus = (l(2) + 2) * d;
    //cout << focus << endl;

    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 2;

    list<string> ObjectFunctionNames{ "IntegratedE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ 70.0 };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 1;
    list<list<double>*> ObjectParameters{ &ObjectParameter };

    string save_position = ".\\tmp1\\";




    TotalTime.open(save_position + "TotalTime.txt");

    Vector3d n_K;
    Vector3d n_E0;



    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");

    int theta_num = 1;
    VectorXd theta(theta_num);
    theta << 0;
    int phi_num = 1;
    VectorXd phi(phi_num);
    phi << 0;
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 542;

    CoreStructure CStr(&spacepara, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "TiO2", lam(0), "nm");
    //AProductCore Core1(&CStr, lam(0), material, "LDR");

    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 22 * d;
    Ln = 22 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");
    //AProductCore Core1(&CStr, lam(0), material, "FCD");
    //material = Get_2_material("Air", "2.5", lam(1), "nm");
    //AProductCore Core2(&CStr, lam(1), material, "LDR");
    //material = Get_2_material("Air", "2.5", lam(2), "nm");
    //AProductCore Core3(&CStr, lam(2), material, "LDR");
    CorePointList.push_back(&Core1);
    //CorePointList.push_back(&Core2);
    //CorePointList.push_back(&Core3);



    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {
        AProductCore* Core = (*it);
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
                    if (k == 0) {
                        AngleInfo << theta[i] << endl;
                        AngleInfo << phi[j] << endl;
                        nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                    }
                    DDAModel Model(Core, n_K, E0, n_E0);
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
        if (k == 0) {
            AngleInfo << 0.0 << endl;
            AngleInfo << 0.0 << endl;
            nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
        }
        DDAModel Model(Core, n_K, E0, n_E0);
        ModelList.push_back(Model);

        it++;
    }


    cout << "nK" << n_K << endl;
    cout << "nE0" << n_E0 << endl;

    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;

    list<DDAModel>::iterator it1 = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it1));
        it1++;
    }


    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, HaveOriginHeritage, HaveAdjointHeritage, PenaltyFactor, save_position, &CStr, ModelpointerList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam", &VforAdam, &SforAdam);





    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();


    return 0;

}

//2021-4-13
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();



    Vector3d l;
    Vector3d center;
    l << 21.0, 21.0, 9.0;
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
    double d = 15;
    double r = 150 / d;

    Vector3i bind(1, 1, 10);
    //SpacePara spacepara(&S, bind, "ONES", "ZEROS", r);

    SpacePara spacepara(&S, bind, "RANDOM");


    double E0 = 1.0;


    double epsilon = 0.2;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    double focus = (l(2) + 2) * d;
    cout << focus << endl;

    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 500;

    list<string> ObjectFunctionNames{ "IntegratedE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ 70.0 };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 1;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = ".\\p330-lam542-beta7-TiO2-InE-rand3\\";

    Vector3d n_K;
    Vector3d n_E0;



    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");

    int theta_num = 1;
    VectorXd theta(theta_num);
    theta << 0;
    int phi_num = 1;
    VectorXd phi(phi_num);
    phi << 0;
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 542;

    CoreStructure CStr(&spacepara, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "TiO2", lam(0), "nm");
    //AProductCore Core1(&CStr, lam(0), material, "LDR");

    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 22 * d;
    Ln = 22 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");
    //AProductCore Core1(&CStr, lam(0), material, "FCD");
    //material = Get_2_material("Air", "2.5", lam(1), "nm");
    //AProductCore Core2(&CStr, lam(1), material, "LDR");
    //material = Get_2_material("Air", "2.5", lam(2), "nm");
    //AProductCore Core3(&CStr, lam(2), material, "LDR");
    CorePointList.push_back(&Core1);
    //CorePointList.push_back(&Core2);
    //CorePointList.push_back(&Core3);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {
        AProductCore* Core = (*it);
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
                    if (k == 0) {
                        AngleInfo << theta[i] << endl;
                        AngleInfo << phi[j] << endl;
                        nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                    }
                    DDAModel Model(Core, n_K, E0, n_E0);
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
        if (k == 0) {
            AngleInfo << 0.0 << endl;
            AngleInfo << 0.0 << endl;
            nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
        }
        DDAModel Model(Core, n_K, E0, n_E0);
        ModelList.push_back(Model);

        it++;
    }


    cout << "nK" << n_K << endl;
    cout << "nE0" << n_E0 << endl;

    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;

    list<DDAModel>::iterator it1 = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it1));
        it1++;
    }


    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, HaveOriginHeritage, HaveAdjointHeritage, PenaltyFactor, save_position, &CStr, ModelpointerList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return 0;

}

//2021-4-19 Find integratedE for 450 period
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();



    Vector3d l;
    Vector3d center;
    l << 31.0, 31.0, 9.0;
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
    double d = 15;
    double r = 150 / d;

    Vector3i bind(1, 1, 10);
    //SpacePara spacepara(&S, bind, "ONES", "ZEROS", r);

    SpacePara spacepara(&S, bind, "RANDOM");


    double E0 = 1.0;


    double epsilon = 0.2;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    double focus = (l(2) + 2) * d;
    cout << focus << endl;

    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 500;

    list<string> ObjectFunctionNames{ "IntegratedE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ 70.0 };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 1;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = ".\\p480-lam542-beta7-TiO2-InE-rand1\\";

    Vector3d n_K;
    Vector3d n_E0;



    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");

    int theta_num = 1;
    VectorXd theta(theta_num);
    theta << 0;
    int phi_num = 1;
    VectorXd phi(phi_num);
    phi << 0;
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 542;

    CoreStructure CStr(&spacepara, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "TiO2", lam(0), "nm");
    //AProductCore Core1(&CStr, lam(0), material, "LDR");

    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 32 * d;
    Ln = 32 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");
    //AProductCore Core1(&CStr, lam(0), material, "FCD");
    //material = Get_2_material("Air", "2.5", lam(1), "nm");
    //AProductCore Core2(&CStr, lam(1), material, "LDR");
    //material = Get_2_material("Air", "2.5", lam(2), "nm");
    //AProductCore Core3(&CStr, lam(2), material, "LDR");
    CorePointList.push_back(&Core1);
    //CorePointList.push_back(&Core2);
    //CorePointList.push_back(&Core3);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {
        AProductCore* Core = (*it);
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
                    if (k == 0) {
                        AngleInfo << theta[i] << endl;
                        AngleInfo << phi[j] << endl;
                        nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                    }
                    DDAModel Model(Core, n_K, E0, n_E0);
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
        if (k == 0) {
            AngleInfo << 0.0 << endl;
            AngleInfo << 0.0 << endl;
            nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
        }
        DDAModel Model(Core, n_K, E0, n_E0);
        ModelList.push_back(Model);

        it++;
    }


    cout << "nK" << n_K << endl;
    cout << "nE0" << n_E0 << endl;

    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;

    list<DDAModel>::iterator it1 = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it1));
        it1++;
    }


    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, HaveOriginHeritage, HaveAdjointHeritage, PenaltyFactor, save_position, &CStr, ModelpointerList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return 0;

}

//2021-4-20 for 4500 ObjE
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();



    Vector3d l;
    Vector3d center;
    l << 180.0, 180.0, 4.0;
    //l << 40.0, 40.0, 8.0;
    center << l(0) / 2, l(1) / 2, l(2) / 2;

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);
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

    Vector3i bind(1, 1, 1);
    SpacePara spacepara(&S, bind, "ONES");

    double d = 25;


    double E0 = 1.0;


    double epsilon = 10;

    //double focus = (l(2) + 2) * d;   //nm       
    double focus = (l(2) + 2) * d;
    cout << focus << endl;

    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 400;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = true;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = ".\\thick100-phi0theta0-lam500-size4500-300it-Adam\\";

    Vector3d n_K;
    Vector3d n_E0;



    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo("AngleInfo.txt");
    ofstream nEInfo("nEInfo.txt");

    int theta_num = 1;
    VectorXd theta(theta_num);
    theta << 0;
    int phi_num = 1;
    VectorXd phi(phi_num);
    phi << 0;
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;

    CoreStructure CStr(&spacepara, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "SiO2", lam(0), "nm");
    AProductCore Core1(&CStr, lam(0), material, "LDR");

    int m, n;
    double Lm, Ln;
    m = 10;
    n = 10;
    Lm = l(0) * d * 1.25;
    Ln = l(1) * d * 1.25;
    //AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "LDR");
    //material = Get_2_material("Air", "2.5", lam(1), "nm");
    //AProductCore Core2(&CStr, lam(1), material, "LDR");
    //material = Get_2_material("Air", "2.5", lam(2), "nm");
    //AProductCore Core3(&CStr, lam(2), material, "LDR");
    CorePointList.push_back(&Core1);
    //CorePointList.push_back(&Core2);
    //CorePointList.push_back(&Core3);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {
        AProductCore* Core = (*it);
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
                    if (k == 0) {
                        AngleInfo << theta[i] << endl;
                        AngleInfo << phi[j] << endl;
                        nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                    }
                    DDAModel Model(Core, n_K, E0, n_E0);
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
        if (k == 0) {
            AngleInfo << 0.0 << endl;
            AngleInfo << 0.0 << endl;
            nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
        }
        DDAModel Model(Core, n_K, E0, n_E0);
        ModelList.push_back(Model);

        it++;
    }




    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;

    list<DDAModel>::iterator it1 = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it1));
        it1++;
    }


    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, HaveOriginHeritage, HaveAdjointHeritage, PenaltyFactor, save_position, &CStr, ModelpointerList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return 0;

}

//2021-6-20 feeding heterodimer back d=15
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();


    double ld = 375;
    double tTiO2 = 180;
    double d = 15;
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
    center1 << 100 / d, 100 / d;
    center2 << (100 + disp * cos(70 * PI / 180)) / d, (100 + disp * sin(70 * PI / 180)) / d;
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
    double lam_max = 730;
    double lam_interval = 1;
    int lam_num = round((lam_max - lam_min) / lam_interval) + 1;

    CoreStructure CStr(&spacepara, d);
    string save_position = ".\\dimer-wavelengthscan-d=15\\";
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
        double E_int = 0.0;
        for (int j = 0; j < N; j++) {
            if ((*R)(3 * j + 2) >= 0) {
                double E_sum_temp = 0;
                for (int k = 0; k < 3; k++) {
                    E_sum_temp += pow(abs((*E_internal)(3 * j + k)), 2);
                }

                E_int += pow(E_sum_temp, 3 / 2) * ((*(CStr.get_diel_old()))(j * 3));  //prevent nan result for devp calculation.


            }
        }

        E_int = log(E_int * pow(d, 3));   //unit: nm^3
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

//2021-6-20 feeding heterodimer back d=10
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
    double lam_max = 660;
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

                E_int += pow(E_sum_temp, 2) * ((*(CStr.get_diel_old()))(j * 3) + 0.0001) / 4.0;  //prevent nan result for devp calculation.
                //cout << "j: " << j << " " << E_int << endl;

            }
        }

        E_int = log(E_int);
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

//2021-6-21 import structure from binary files converted from COMSOL and do DDA calculation
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();


    string open_position = ".\\170-air-d=10\\";
    string name1 = open_position + "CoreStructre-fromComsol-170air.txt";

    ifstream fin1(name1);
    int Nx, Ny, Nz;
    int Ntmp;
    fin1 >> Nx;
    fin1 >> Ny;
    fin1 >> Nz;
    fin1 >> Ntmp;
    VectorXi geometry = VectorXi::Zero(3 * Ntmp);
    for (int i = 0; i <= Ntmp - 1; i++) {
        fin1 >> geometry(3 * i);
        fin1 >> geometry(3 * i + 1);
        fin1 >> geometry(3 * i + 2);
    }
    double d;
    fin1 >> d;
    fin1.close();


    string name2, name3, name4;
    name2 = open_position + "geometryPara-fromComsol-170air.txt";
    name3 = open_position + "Para-fromComsol-170air.txt";
    name4 = open_position + "bind-fromComsol-170air.txt";
    ifstream fin2(name2), fin3(name3), fin4(name4);

    VectorXi geometryPara = VectorXi::Zero(Ntmp);
    int parasize;
    fin3 >> parasize;
    VectorXd Para = VectorXd::Zero(parasize);
    Vector3i bind = Vector3i::Zero();
    for (int i = 0; i <= Ntmp - 1; i++) {
        fin2 >> geometryPara(i);
    }
    for (int i = 0; i <= parasize - 1; i++) {
        fin3 >> Para(i);
    }
    for (int i = 0; i <= 2; i++) {
        fin4 >> bind(i);
    }
    fin2.close();
    fin3.close();
    fin4.close();


    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Structure s1(S.get_total_space(), &geometry);

    S = S + s1;


    SpacePara spacepara(&S, bind, &geometryPara, &Para);


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

    double lam_min = 600;
    double lam_max = 680;
    double lam_interval = 1;
    int lam_num = round((lam_max - lam_min) / lam_interval) + 1;

    CoreStructure CStr(&spacepara, d);
    string save_position = ".\\170-air-d=10\\";
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
        double E_int = 0.0;
        for (int j = 0; j < N; j++) {
            if ((*R)(3 * j + 2) >= 0) {
                double E_sum_temp = 0;
                for (int k = 0; k < 3; k++) {
                    E_sum_temp += pow(abs((*E_internal)(3 * j + k)), 2);
                }

                E_int += pow(E_sum_temp, 3 / 2) * ((*(CStr.get_diel_old()))(j * 3));  //prevent nan result for devp calculation.


            }
        }

        E_int = log(E_int * pow(d, 3));   //unit: nm^3
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

//random rect generation 2021-7-31
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
    int number = 6;
    double limitx1 = 4;
    double limitx2 = 25;
    double limity1 = 4;
    double limity2 = 25;

    SpacePara spacepara(&S, bind, number, limitx1, limitx2, limity1, limity2);


    CoreStructure CStr(&spacepara, d);
    AProductCore Core(&CStr, lam, material, "LDR");
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    string save_position = ".\\SiO2-rects\\";

    TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel.update_E_in_structure();
    TestModel.solve_E();
    //TestModel.output_to_file(save_position, 0, 0);
    //CStr.output_to_file(save_position, 0);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    int num_model = 15000;
    int start_num = 0;

    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    for (int i = 0; i <= num_model - 1; i++) {
        SpacePara spacepara_tmp(&S, bind, number, limitx1, limitx2, limity1, limity2, spacepara.get_geometryPara());
        CStr.UpdateStr(&spacepara_tmp);
        CStr.output_to_file(save_position, start_num + i + 1, "Simple");
        TestModel.UpdateAlpha();
        TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();
        TestModel.output_to_file(save_position, start_num + i + 1);
    }
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();





    return 0;

}

