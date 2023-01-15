﻿#define PI 3.14159265
#define _USE_MATH_DEFINES

#include <chrono>
#include <fstream>
#include <iostream>

#include "EvoDDAModel.h"
#include "filterReader.h"
#include "INIReader.h"
#include "symReader.h"
#include "ObjReader.h"
#include "Tools.h"


//---Tmp including head files
#include "Space.h"
#include "SpacePara.h"


using namespace std::chrono;

ObjDDAModel* ObjFactoryOut(string ObjectName, vector<double> ObjectParameters, DDAModel* ObjDDAModel) {
    /*if (HavePenalty) {
        cout << "Using L1 Penalty with Penalty Factor " << PenaltyFactor << endl;
    }*/
    if (ObjectName == "PointE") {
        return new ObjPointEDDAModel(ObjectParameters, ObjDDAModel);
    }
    if (ObjectName == "IntegratedE") {
        return new ObjIntegratedEDDAModel(ObjectParameters, ObjDDAModel);
    }
    cout << "NOT A LEGIT OBJECTIVE NAME!" << endl;
    return new ObjPointEDDAModel(ObjectParameters, ObjDDAModel);
}

void task() {
    INIReader reader1("task.ini");

    if (reader1.ParseError() != 0) {
        std::cout << "Can't load 'task.ini'\n";
        return;
    }

    map<string, string> findtask;
    findtask.insert(pair<string, string>("DDA verify", "DDA_verify_path"));
    findtask.insert(pair<string, string>("DDA input", "DDA_input_path"));
    findtask.insert(pair<string, string>("EvoOpt", "EvoOpt_path"));
    findtask.insert(pair<string, string>("EvoOpt 2D input", "EvoOpt_2D_input_path"));
    findtask.insert(pair<string, string>("EvoOpt 2D input periodic", "EvoOpt_2D_input_periodic_path"));
    findtask.insert(pair<string, string>("NN data generate", "NN_path"));

    string tasktype = reader1.Get("Model Name", "name", "UNKNOWN");
    string tasktypepath = findtask[tasktype];
    string path = reader1.Get("Path", tasktypepath, "UNKNOWN");

    cout << path << endl;
    INIReader reader2(path);
    if (tasktype == "DDA verify") {
        high_resolution_clock::time_point t_start = high_resolution_clock::now();
        int Nx, Ny, Nz;
        Nx = reader2.GetInteger("Geometry", "Nx", -1); Ny = reader2.GetInteger("Geometry", "Ny", -1); Nz = reader2.GetInteger("Geometry", "Nz", -1);
        int N = 0;
        VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
        vector<Structure> ln;
        Space S(&total_space, Nx, Ny, Nz, N, &ln);

        double d;
        double r;
        Vector3d center;

        d = reader2.GetInteger("Grid", "d", -1);
        r = reader2.GetInteger("Geometry", "r", -1) / d;
        center(0) = Nx / 2;
        center(1) = Ny / 2;
        center(2) = Nz / 2;

        Structure s1(S.get_total_space(), r, center);
        S = S + s1;

        double lam = reader2.GetInteger("Input field", "lam", -1);
        Vector3d n_K;
        n_K(0) = 0.0;
        n_K(1) = 0.0;
        n_K(2) = 1.0;
        double E0 = 1.0;
        Vector3d n_E0;
        n_E0(0) = 1.0;
        n_E0(1) = 0.0;
        n_E0(2) = 0.0;
        list<string> mat_l{ reader2.Get("Material", "material1", "UNKNOWN"), reader2.Get("Material", "material2", "UNKNOWN") };
        VectorXcd material = Get_X_material(mat_l, lam, "nm");



        int MAX_ITERATION_DDA = reader2.GetInteger("DDA iteration", "MAX_ITERATION_DDA", -1);
        double MAX_ERROR = reader2.GetFloat("DDA iteration", "MAX_ERROR", -1);

        Vector3i bind(1, 1, 1);


        VectorXi* s1geometry = s1.get_geometry();
        vector<string> FParaInit{};
        vector<double> BParal{ 1.0 };
        SpacePara spacepara(&S, bind, FParaInit, BParal);

        CoreStructure CStr(&spacepara, d);
        double nback = sqrt(real(material(0)));
        AProductCore Core(&CStr, lam, material, nback, "LDR");
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

        high_resolution_clock::time_point t_end = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(t_end - t_start).count();
        cout << "total time: " << duration << "ms" << endl;


        return;
    }
    else if (tasktype == "DDA input") {


        VectorXi inputGeo;
        VectorXd inputDiel;
        string pathCommonData = reader2.Get("Geometry", "pathCommonData", "UNKNOWN");
        string pathPara = reader2.Get("Geometry", "pathPara", "UNKNOWN");
        tie(inputGeo, inputDiel) = getInputStr(pathCommonData, pathPara);

        string save_position = reader2.Get("Output", "saveDir", "UNKNOWN");       //output file

        /*MatrixXi scope = find_scope_3_dim(&inputGeo);*/
        int Nx, Ny, Nz;
        /*Nx = scope(0, 1) - scope(0, 0) + 1;
        Ny = scope(1, 1) - scope(1, 0) + 1;
        Nz = scope(2, 1) - scope(2, 0) + 1;*/
        tie(Nx, Ny, Nz) = getInputNs(pathCommonData);
        int N = 0;
        VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
        vector<Structure> ln;
        Space S(&total_space, Nx, Ny, Nz, N, &ln);

        double d;
        Vector3d center;
        Vector3d l;
        d = stod(reader2.Get("Grid", "d", "UNKNOWN"));
        center(0) = Nx / 2;
        center(1) = Ny / 2;
        center(2) = Nz / 2;
        /*l << Nx - 1, Ny - 1, Nz - 1;
        Structure s1(S.get_total_space(), l, center);*/
        Structure s1(S.get_total_space(), &inputGeo);

        S = S + s1;

        vector<double> lamSweep = ReadLam(reader2.Get("Input field", "lam", "UNKNOWN"));


        Vector3d n_K;
        n_K(0) = reader2.GetFloat("Input field", "nKx", 0.0);
        n_K(1) = reader2.GetFloat("Input field", "nKy", 0.0);
        n_K(2) = reader2.GetFloat("Input field", "nKz", 0.0);
        double E0 = 1.0;
        Vector3d n_E0;
        n_E0(0) = reader2.GetFloat("Input field", "nEx", 0.0);
        n_E0(1) = reader2.GetFloat("Input field", "nEy", 0.0);
        n_E0(2) = reader2.GetFloat("Input field", "nEz ", 0.0);

        list<string> mat_l = ReadMat(reader2.Get("Material", "materialnames", "UNKNOWN"));
        int MAX_ITERATION_DDA = reader2.GetInteger("DDA iteration", "MAX_ITERATION_DDA", -1);
        double MAX_ERROR = reader2.GetFloat("DDA iteration", "MAX_ERROR", -1);
        Vector3i bind(1, 1, 1);
        /*vector<string> FParaInit{};
        vector<double> BParal{ 1.0 };*/
        SpacePara spacepara(&S, bind, &inputGeo, &inputDiel);
        CoreStructure CStr(&spacepara, d);

        string objPath = save_position + "Obj200-150-80.txt";
        ofstream objOut;
        objOut.open(objPath);

        string objPath1 = save_position + "Obj200-150-200.txt";
        ofstream objOut1;
        objOut1.open(objPath1);

        string objPath2 = save_position + "Obj200-150-300.txt";
        ofstream objOut2;
        objOut2.open(objPath2);

        string objPath3 = save_position + "Obj200-150-400.txt";
        ofstream objOut3;
        objOut3.open(objPath3);

        string timePath = save_position + "Time.txt";
        ofstream timeCount;
        timeCount.open(timePath);
        for (int i = 0; i < lamSweep.size(); i++) {
            high_resolution_clock::time_point t_start = high_resolution_clock::now();
            double lam = lamSweep[i];
            VectorXcd material = Get_X_material(mat_l, lam, "nm");
            double nback = sqrt(real(material(0)));
            AProductCore Core(&CStr, lam, material, nback, "FCD");
            DDAModel TestModel(&Core, n_K, E0, n_E0);

            TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
            TestModel.update_E_in_structure();
            TestModel.solve_E();
            string savePosition = reader2.Get("Output", "saveDir", "UNKNOWN");
            TestModel.output_to_file(savePosition + "Model_output\\", i, 0);
            ObjReader objReader(reader2);
            string objName = objReader.GetObjName();
            vector<double> objPara = objReader.GetObjPara();

            ObjDDAModel* objEvas = ObjFactoryOut(objName, objPara, &TestModel);
            objEvas->Reset();
            objOut << objEvas->GetVal() << endl;

            ObjDDAModel* objEvas1 = ObjFactoryOut("PointE", vector<double>{200, 150, 200}, & TestModel);
            objEvas1->Reset();
            objOut1 << objEvas1->GetVal() << endl;

            ObjDDAModel* objEvas2 = ObjFactoryOut("PointE", vector<double>{200, 150, 300}, &TestModel);
            objEvas2->Reset();
            objOut2 << objEvas2->GetVal() << endl;

            ObjDDAModel* objEvas3 = ObjFactoryOut("PointE", vector<double>{200, 150, 400}, &TestModel);
            objEvas3->Reset();
            objOut3 << objEvas3->GetVal() << endl;

            CStr.output_to_file(save_position + "CoreStructure\\", i, "simple");
            high_resolution_clock::time_point t_end = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(t_end - t_start).count();
            timeCount << duration / 1000 << endl;
            
        }

        
        string nameCommonData = save_position + "commondata.txt";
        ofstream Common;
        Common.open(nameCommonData);
        Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
        Common << (spacepara.get_geometry()) << endl;
        Common << d << endl;
        Common << n_E0 << endl;
        Common << n_K << endl;
        return;

    }
    else if (tasktype == "EvoOpt 2D input") {
        VectorXi inputGeo;
        VectorXd inputDiel;
        string pathCommonData = reader2.Get("Geometry", "pathCommonData", "UNKNOWN");
        string pathPara = reader2.Get("Geometry", "pathPara", "UNKNOWN");
        tie(inputGeo, inputDiel) = getInputStr(pathCommonData, pathPara);

        string save_position = reader2.Get("Output", "saveDir", "UNKNOWN");       //output file

        int Nx, Ny, Nz;
        tie(Nx, Ny, Nz) = getInputNs(pathCommonData);
        int N = 0;
        VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
        vector<Structure> ln;
        Space S(&total_space, Nx, Ny, Nz, N, &ln);

        double d;
        Vector3d center;
        Vector3d l;
        d = stod(reader2.Get("Grid", "d", "UNKNOWN"));
        center(0) = Nx / 2;
        center(1) = Ny / 2;
        center(2) = Nz / 2;
        /*l << Nx - 1, Ny - 1, Nz - 1;
        Structure s1(S.get_total_space(), l, center);*/
        Structure s1(S.get_total_space(), &inputGeo, true);

        S = S + s1;

        double lam = reader2.GetInteger("Input field", "lam", -1);
        Vector3d n_K;
        n_K(0) = reader2.GetFloat("Input field", "nKx", 0.0);
        n_K(1) = reader2.GetFloat("Input field", "nKy", 0.0);
        n_K(2) = reader2.GetFloat("Input field", "nKz", 0.0);
        double E0 = 1.0;
        Vector3d n_E0;
        n_E0(0) = reader2.GetFloat("Input field", "nEx", 0.0);
        n_E0(1) = reader2.GetFloat("Input field", "nEy", 0.0);
        n_E0(2) = reader2.GetFloat("Input field", "nEz ", 0.0);

        list<string> mat_l{ reader2.Get("Material", "material1", "UNKNOWN"), reader2.Get("Material", "material2", "UNKNOWN") };
        VectorXcd material = Get_X_material(mat_l, lam, "nm");

        int MAX_ITERATION_DDA = reader2.GetInteger("DDA iteration", "MAX_ITERATION_DDA", -1);
        int MAX_ITERATION_EVO = reader2.GetInteger("Evo Option", "MAX_ITERATION_EVO", -1);
        double MAX_ERROR = reader2.GetFloat("DDA iteration", "MAX_ERROR", -1);

        MatrixXi scope = find_scope_3_dim(&inputGeo);
        int bindz = scope(2, 1) - scope(2, 0) + 1;
        Vector3i bind(1, 1, bindz);

        filterReader readTheFilter(reader2);
        bool filter = readTheFilter.getFilter();
        vector<filterinfo> filterList = readTheFilter.getFilterList();
        
        FilterOption filterOpt(readTheFilter.getBetaMin(), readTheFilter.getBetaMax(), readTheFilter.getIta(), readTheFilter.getBetaType(), filterList);

        symReader readSymmetry(reader2);
        string symmetry = readSymmetry.getSymmetry();
        vector<double> symAxis = readSymmetry.getSymAxis();
        SpacePara spacepara(&S, bind, &inputGeo, &inputDiel, filter, &filterOpt, symmetry, symAxis);
        CoreStructure CStr(&spacepara, d);
        double nback = sqrt(real(material(0)));
        AProductCore Core(&CStr, lam, material, nback, "LDR");
        DDAModel TestModel(&Core, n_K, E0, n_E0);

        ObjReader objReader(reader2);
        string objName = objReader.GetObjName();
        vector<double> objPara=objReader.GetObjPara();  //Focal spot postition.
        //vector<DDAModel> models;
        vector<DDAModel*> modelPtrs;
        modelPtrs.push_back(&TestModel);

        bool HavePathRecord = false;
        bool HaveOriginHeritage = false;
        bool HaveAdjointHeritage = false;
        double epsilon = reader2.GetInteger("Evo Option", "epsilon", 1);
        EvoDDAModel evoModel(objName, objPara, epsilon, HavePathRecord, HaveOriginHeritage, HaveAdjointHeritage, save_position, &CStr, modelPtrs);
        evoModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");
        /*TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();

        string savePosition = reader2.Get("Output", "saveDir", "UNKNOWN");
        TestModel.output_to_file(savePosition + "Model_output\\", 0);*/

        string nameCommonData = save_position + "commondata.txt";
        ofstream Common;
        Common.open(save_position + "commondata.txt");
        Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
        Common << (spacepara.get_geometry()) << endl;
        Common << d << endl;
        Common << n_E0 << endl;
        Common << n_K << endl;
        return;


    }
    else if (tasktype == "EvoOpt 2D input periodic") {
    VectorXi inputGeo;
    VectorXd inputDiel;
    string pathCommonData = reader2.Get("Geometry", "pathCommonData", "UNKNOWN");
    string pathPara = reader2.Get("Geometry", "pathPara", "UNKNOWN");
    tie(inputGeo, inputDiel) = getInputStr(pathCommonData, pathPara);

    string save_position = reader2.Get("Output", "saveDir", "UNKNOWN");       //output file

    int Nx, Ny, Nz;
    tie(Nx, Ny, Nz) = getInputNs(pathCommonData);
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    vector<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;
    Vector3d center;
    Vector3d l;
    d = stod(reader2.Get("Grid", "d", "UNKNOWN"));
    center(0) = Nx / 2;
    center(1) = Ny / 2;
    center(2) = Nz / 2;
    /*l << Nx - 1, Ny - 1, Nz - 1;
    Structure s1(S.get_total_space(), l, center);*/
    Structure s1(S.get_total_space(), &inputGeo, true);

    S = S + s1;

    double lam = reader2.GetInteger("Input field", "lam", -1);
    Vector3d n_K;
    n_K(0) = reader2.GetFloat("Input field", "nKx", 0.0);
    n_K(1) = reader2.GetFloat("Input field", "nKy", 0.0);
    n_K(2) = reader2.GetFloat("Input field", "nKz", 0.0);
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0(0) = reader2.GetFloat("Input field", "nEx", 0.0);
    n_E0(1) = reader2.GetFloat("Input field", "nEy", 0.0);
    n_E0(2) = reader2.GetFloat("Input field", "nEz ", 0.0);

    list<string> mat_l{ reader2.Get("Material", "material1", "UNKNOWN"), reader2.Get("Material", "material2", "UNKNOWN") };
    VectorXcd material = Get_X_material(mat_l, lam, "nm");

    int MAX_ITERATION_DDA = reader2.GetInteger("DDA iteration", "MAX_ITERATION_DDA", -1);
    int MAX_ITERATION_EVO = reader2.GetInteger("Evo Option", "MAX_ITERATION_EVO", -1);
    double MAX_ERROR = reader2.GetFloat("DDA iteration", "MAX_ERROR", -1);

    MatrixXi scope = find_scope_3_dim(&inputGeo);
    int bindz = scope(2, 1) - scope(2, 0) + 1;
    Vector3i bind(1, 1, bindz);

    filterReader readTheFilter(reader2);
    bool filter = readTheFilter.getFilter();
    vector<filterinfo> filterList = readTheFilter.getFilterList();
    int m, n;
    int Lm, Ln;
    m = reader2.GetInteger("Periodicity Option", "m", -1);
    n = reader2.GetInteger("Periodicity Option", "n", -1);
    Lm = reader2.GetInteger("Periodicity Option", "Lx", -1);
    Ln = reader2.GetInteger("Periodicity Option", "Ly", -1);
    FilterOption filterOpt(readTheFilter.getBetaMin(), readTheFilter.getBetaMax(), readTheFilter.getIta(), readTheFilter.getBetaType(), filterList);

    symReader readSymmetry(reader2);
    string symmetry = readSymmetry.getSymmetry();
    vector<double> symAxis = readSymmetry.getSymAxis();
    SpacePara spacepara(&S, bind, &inputGeo, &inputDiel, filter, &filterOpt, symmetry, symAxis, readTheFilter.getPeriodic(), Lm, Ln);
    CoreStructure CStr(&spacepara, d);
    double nback = sqrt(real(material(0)));
    
    AProductCore Core(&CStr, lam, material, nback, m, n, Lm* d, Ln* d, "FCD");
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    ObjReader objReader(reader2);
    string objName = objReader.GetObjName();
    vector<double> objPara = objReader.GetObjPara();  //Focal spot postition.
    //vector<DDAModel> models;
    vector<DDAModel*> modelPtrs;
    modelPtrs.push_back(&TestModel);

    bool HavePathRecord = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double epsilon = reader2.GetInteger("Evo Option", "epsilon", 1);
    EvoDDAModel evoModel(objName, objPara, epsilon, HavePathRecord, HaveOriginHeritage, HaveAdjointHeritage, save_position, &CStr, modelPtrs);
    evoModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");
    /*TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel.update_E_in_structure();
    TestModel.solve_E();

    string savePosition = reader2.Get("Output", "saveDir", "UNKNOWN");
    TestModel.output_to_file(savePosition + "Model_output\\", 0);*/

    string nameCommonData = save_position + "commondata.txt";
    ofstream Common;
    Common.open(save_position + "commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;
    return;


    }
    else if (tasktype == "NN data generate") {
    srand((unsigned)(time(0)));

    int Nx, Ny, Nz;
    Nx = 32; Ny = 32; Nz = 8;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    vector<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;

    Vector3d center;
    Vector3d l;

    d = 20;

    center << double(Nx - 1) / 2, double(Ny - 1) / 2, double(Nz - 1) / 2;
    l << Nx - 1, Ny - 1, Nz - 1;

    Structure s1(S.get_total_space(), l, center);



    S = S + s1;

    double lam = 650;
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

    Vector3i bind(1, 1, 1);
    int maxnumber = 15;
    int minnumber = 6;
    int number = 0;
    double limitx1 = 2;
    double limitx2 = 5;
    double limity1 = 2;
    double limity2 = 5;
    double limitz1 = 2;
    double limitz2 = 5;

    SpacePara spacepara(&S, bind, maxnumber, limitx1, limitx2, limity1, limity2, limitz1, limitz2);

    double nback = sqrt(real(material(0)));
    CoreStructure CStr(&spacepara, d);
    AProductCore Core(&CStr, lam, material, nback, "LDR");

    DDAModel TestModel(&Core, n_K, E0, n_E0, false);


    string save_position = reader2.Get("Output", "saveDir", "UNKNOWN");

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

    int num_model = 10;

    int start_num = 0;

    double paralow = 0.05;
    double parahigh = 1.0;

    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    for (int i = 0; i <= num_model - 1; i++) {
        
        number = round(((double)rand() / RAND_MAX) * (maxnumber - minnumber + 1)) + minnumber;
        SpacePara spacepara_tmp(&S, bind, number, limitx1, limitx2, limity1, limity2, limitz1, limitz2, spacepara.get_geometryPara(), paralow, parahigh);
        CStr.UpdateStr(&spacepara_tmp);
        CStr.output_to_file(save_position + "CoreStructure\\", start_num + i + 1, "Simple");
        TestModel.UpdateAlpha();
        TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();
        TestModel.output_to_file(save_position + "Model_output\\", start_num + i + 1);
        cout << start_num + i + 1 << endl;
    }
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();





    return;
    }
    else {

    }



}

int main() {

    task();

}








