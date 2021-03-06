#include "definition.h"

CoreStructure::CoreStructure(SpacePara* spacepara_, double d_) {

    spacepara = spacepara_;
    space = (*spacepara).get_space();
    d = d_;

    cout << "(d=" << d << ") " << endl;

    tie(Nx, Ny, Nz, N) = (*space).get_Ns();

    //-----------------------------------------------------------------Input strs-------------------------------------------------------------
    R = (*spacepara).get_geometry();
    VectorXi* geometryPara = (*spacepara).get_geometryPara();
    VectorXd* Para = (*spacepara).get_Para();
    //---------------------------------------------------initial diel------------------------------------
    diel_old = VectorXd::Zero(3 * N);
    diel_old_max = diel_old;
    for (int i = 0; i <= N - 1; i++) {
        double dieltmp = (*Para)((*geometryPara)(i));
        diel_old(3 * i) = dieltmp;
        diel_old(3 * i + 1) = dieltmp;
        diel_old(3 * i + 2) = dieltmp;
    }
}

void CoreStructure::UpdateStr(VectorXd step) {

    cout << "step in UpdateStr" << step.mean() << endl;
    VectorXi* geometryPara = (*spacepara).get_geometryPara();
    VectorXd* Para = (*spacepara).get_Para();
    int Parasize = (*Para).size();
    if (Parasize != step.size()) {
        cout << "ERROR: In CoreStructure::UpdateStr(VectorXd step), step.size!=Para.size";
    }
    for (int i = 0; i <= Parasize - 1; i++) {
        (*Para)(i) += step(i);
        if ((*Para)(i) >= 1) {
            (*Para)(i) = 1;
        }
        if ((*Para)(i) <= 0) {
            (*Para)(i) = 0;
        }
    }
    for (int i = 0; i <= N - 1; i++) {
        int position = (*geometryPara)(i);
        double value = (*Para)(position);
        diel_old(3 * i) = value;
        diel_old(3 * i + 1) = value;
        diel_old(3 * i + 2) = value;
    }
}

void CoreStructure::output_to_file() {

    ofstream fout("CoreStructure.txt");
    fout << Nx << endl << Ny << endl << Nz << endl << N << endl;
    fout << R << endl;
    fout << diel_old << endl;
    fout << d << endl;
    fout.close();
}

void CoreStructure::output_to_file(string save_position, int iteration) {

    string name;
    //name=save_position+"AProductCoreit" + to_string(iteration) + ".txt";
    name = save_position + "Model_output\\CoreStructure" + to_string(iteration) + ".txt";
    ofstream fout(name);
    fout << Nx << endl << Ny << endl << Nz << endl << N << endl;
    fout << R << endl;
    fout << diel_old << endl;
    fout << d << endl;
    fout.close();

    //-------------------------------Following for intermediate restart (output the structure and restart as initial in another calculation)--------------------
    /*
    if ((iteration == 87) || (iteration == 88)) {
        string name1, name2, name3;
        name1 = save_position + to_string(iteration) + "geometryPara.txt";
        name2 = save_position + to_string(iteration) + "Para.txt";
        name3 = save_position + to_string(iteration) + "bind.txt";
        ofstream fout1(name1), fout2(name2), fout3(name3);
        fout1 << *((*spacepara).get_geometryPara()) << endl;
        fout2 << (*((*spacepara).get_Para())).size() << endl;
        fout2 << *((*spacepara).get_Para()) << endl;
        fout3 << *((*spacepara).get_bind()) << endl;
        fout1.close();
        fout2.close();
        fout3.close();
    }
    */
}

int CoreStructure::get_N() {
    return N;
}
int CoreStructure::get_Nx() {
    return Nx;
}
int CoreStructure::get_Ny() {
    return Ny;
}
int CoreStructure::get_Nz() {
    return Nz;
}
VectorXi* CoreStructure::get_R() {
    return &R;
}
double CoreStructure::get_d() {
    return d;
}
SpacePara* CoreStructure::get_spacepara() {
    return spacepara;
}
VectorXd* CoreStructure::get_diel_old() {
    return &diel_old;
}
VectorXd* CoreStructure::get_diel_old_max() {
    return &diel_old_max;
}

/*
tuple<list<int>, list<int>, list<int>, list<int>> CoreStructure::get_para_info() {
    return make_tuple(para_nums, para_starts, para_dep_nums, para_dep_starts);
}

list<list<int>>* CoreStructure::get_PositionDep() {
    return &PositionDep;
}
VectorXi* CoreStructure::get_PositionPara() {
    return &PositionPara;
}
list<int>* CoreStructure::get_para_nums() {
    return &para_nums;
}
list<int>* CoreStructure::get_para_starts() {
    return &para_starts;
}
list<int>* CoreStructure::get_para_dep_nums() {
    return &para_dep_nums;
}
list<int>* CoreStructure::get_para_dep_starts() {
    return &para_dep_starts;
}
*/