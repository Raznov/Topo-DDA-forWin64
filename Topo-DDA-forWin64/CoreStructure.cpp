#include "definition.h"

CoreStructure::CoreStructure(SpacePara* spacepara_, double d_) {

    spacepara = spacepara_;
    space = (*spacepara).get_space();
    d = d_;

    cout << "(d=" << d << ") " << endl;

    tie(Nx, Ny, Nz, N) = (*space).get_Ns();

    //-----------------------------------------------------------------Input strs-------------------------------------------------------------
    R = (*spacepara).get_geometry();
    Vectori* geometryPara = (*spacepara).get_geometryPara();
    Vectord* Para = (*spacepara).get_Para();
    //---------------------------------------------------initial diel------------------------------------
    diel_old = Vectord(3 * N);
    diel_old_max = diel_old;
    for (int i = 0; i <= N - 1; i++) {
        double dieltmp = (*Para)((*geometryPara)(i));
        diel_old(3 * i) = dieltmp;
        diel_old(3 * i + 1) = dieltmp;
        diel_old(3 * i + 2) = dieltmp;
    }
}

void CoreStructure::UpdateStr(Vectord step, int current_it, int Max_it) {
    cout << "step in UpdateStr" << vecmean(step) << endl;
    Vectori* geometryPara = (*spacepara).get_geometryPara();
    Vectord* Para = (*spacepara).get_Para();
    Vectori* Free = (*spacepara).get_Free();

    int Parasize = (*Free).size();
    if (Parasize != step.size()) {
        cout << "ERROR: In CoreStructure::UpdateStr(Vectord step), step.size!=FreePara.size";
        throw 1;
    }
    

    if ((*spacepara).get_Filter()) {
        //When there is filter
        Vectord* Para_origin = (*spacepara).get_Para_origin();
        Vectord* Para_filtered = (*spacepara).get_Para_filtered();
        FilterOption* Filterstats = (*spacepara).get_Filterstats();
        (*Filterstats).update_beta(current_it, Max_it);                  //Update beta value according to current iteration
        
        for (int i = 0; i <= Parasize - 1; i++) {
            (*Para_origin)((*Free)(i)) += step(i);
            if ((*Para_origin)((*Free)(i)) >= 1) {
                (*Para_origin)((*Free)(i)) = 1;
            }
            if ((*Para_origin)((*Free)(i)) <= 0) {
                (*Para_origin)((*Free)(i)) = 0;
            }
        }
        
        cout << "Beta at iteration " << current_it << " is " << (*Filterstats).get_beta() << endl;

        vector<vector<WeightPara>>* FreeWeight = (*spacepara).get_FreeWeight();
        for (int i = 0; i <= Parasize - 1; i++) {
            int weightnum = ((*FreeWeight)[i]).size();
            double numerator = 0.0;
            double denominator = 0.0;
            for (int j = 0; j <= weightnum - 1; j++) {
                numerator += ((*FreeWeight)[i][j].weight) * (*Para_origin)((*FreeWeight)[i][j].position);
                denominator += ((*FreeWeight)[i][j].weight);
            }
            (*Para_filtered)((*Free)(i)) = numerator / denominator;

            double Para_physical = (*Filterstats).SmoothDensity((*Para_filtered)((*Free)(i)));     
            (*Para)((*Free)(i)) = Para_physical;

        }
    }
    else {//When there is no filter
        for (int i = 0; i <= Parasize - 1; i++) {
            (*Para)((*Free)(i)) += step(i);
            if ((*Para)((*Free)(i)) >= 1) {
                (*Para)((*Free)(i)) = 1;
            }
            if ((*Para)((*Free)(i)) <= 0) {
                (*Para)((*Free)(i)) = 0;
            }
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

void CoreStructure::UpdateStr(SpacePara* spacepara_){
  
    Vectori* geometryPara = (*spacepara).get_geometryPara();
    Vectord* Para = (*spacepara).get_Para();
    Vectori* Free = (*spacepara).get_Free();
    Vectord* Parareplace = (*spacepara_).get_Para();
    Vectori* Freereplace = (*spacepara_).get_Free();

    
    for (int i = 0; i <= (*Para).size() - 1; i++) {
        (*Para)(i) = (*Parareplace)(i);
    }
    for (int i = 0; i <= (*Free).size() - 1; i++) {
        (*Free)(i) = (*Freereplace)(i);
    }
    for (int i = 0; i <= N - 1; i++) {
        int position = (*geometryPara)(i);
        double value = (*Para)(position);
        diel_old(3 * i) = value;
        diel_old(3 * i + 1) = value;
        diel_old(3 * i + 2) = value;
    }
}

void CoreStructure::UpdateStrSingle(int idx, double value) {
    //vector<list<int>>* Paratogeometry = (*spacepara).get_Paratogeometry();

    /*
    list<int>::iterator it = (*Paratogeometry)[idx].begin();
    for (int j = 0; j <= (*Paratogeometry)[idx].size() - 1; j++) {
        int position = *it;
        diel_old(3 * position) = value;
        diel_old(3 * position + 1) = value;
        diel_old(3 * position + 2) = value;
        it++;
    }
    */
    diel_old(3 * idx) = value;
    diel_old(3 * idx + 1) = value;
    diel_old(3 * idx + 2) = value;

}

void CoreStructure::output_to_file() {

    ofstream fout("CoreStructure.txt");
    fout << Nx << endl << Ny << endl << Nz << endl << N << endl;
    vectofile(fout, R);
    vectofile(fout, diel_old);
    fout << d << endl;
    fout.close();
}

void CoreStructure::output_to_file(string save_position, int iteration, string mode) {

    if (mode == "normal") {
        string name;
        //name=save_position+"AProductCoreit" + to_string(iteration) + ".txt";
        name = save_position + "CoreStructure" + to_string(iteration) + ".txt";
        //name = save_position + "CoreStructure_verify\\CoreStructure" + to_string(iteration) + ".txt";
        ofstream fout(name);
        fout << Nx << endl << Ny << endl << Nz << endl << N << endl;
        vectofile(fout, R);
        vectofile(fout, diel_old);
        fout << d << endl;
        fout.close();
    }
    else {
        string name;
        //name=save_position+"AProductCoreit" + to_string(iteration) + ".txt";
        name = save_position + "CoreStructure" + to_string(iteration) + ".txt";
        //name = save_position + "CoreStructure_verify\\CoreStructure" + to_string(iteration) + ".txt";
        ofstream fout(name);
        vectofile(fout, diel_old);
        fout.close();
    }
    

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
Vectori* CoreStructure::get_R() {
    return &R;
}
double CoreStructure::get_d() {
    return d;
}
SpacePara* CoreStructure::get_spacepara() {
    return spacepara;
}
Vectord* CoreStructure::get_diel_old() {
    return &diel_old;
}
Vectord* CoreStructure::get_diel_old_max() {
    return &diel_old_max;
}

/*
tuple<list<int>, list<int>, list<int>, list<int>> CoreStructure::get_para_info() {
    return make_tuple(para_nums, para_starts, para_dep_nums, para_dep_starts);
}

list<list<int>>* CoreStructure::get_PositionDep() {
    return &PositionDep;
}
Vectori* CoreStructure::get_PositionPara() {
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