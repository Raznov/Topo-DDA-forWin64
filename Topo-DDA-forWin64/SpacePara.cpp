#include "definition.h"

VectorXi SpacePara::cut(VectorXi* big, VectorXi* smalll) {

    int number_origin = round((*smalll).size() / 3);
    MatrixXi big_scope = find_scope_3_dim(big);
    //cout<<"big_scope "<<big_scope<<endl;
    list<int> positions_in;
    int number_out = 0;
    //cout<<"small_scope "<<find_scope_3_dim(small)<<endl;
    for (int i = 0; i <= number_origin - 1; i++) {
        if (((*smalll)(3 * i) < big_scope(0, 0)) || ((*smalll)(3 * i) > big_scope(0, 1)) ||
            ((*smalll)(3 * i + 1) < big_scope(1, 0)) || ((*smalll)(3 * i + 1) > big_scope(1, 1)) ||
            ((*smalll)(3 * i + 2) < big_scope(2, 0)) || ((*smalll)(3 * i + 2) > big_scope(2, 1))) {
            number_out += 1;
        }
        else {
            positions_in.push_back(i);
        }
    }
    int number_in = positions_in.size();
    VectorXi geometry = VectorXi::Zero(3 * number_in);
    for (int i = 0; i <= number_in - 1; i++) {
        int j = positions_in.front();
        positions_in.pop_front();
        geometry(3 * i) = (*smalll)(3 * j);
        geometry(3 * i + 1) = (*smalll)(3 * j + 1);
        geometry(3 * i + 2) = (*smalll)(3 * j + 2);
    }
    if (number_out == 0) {
        cout << "The geometry you built is entirely in the space." << endl;
        cout << "number_origin " << number_origin << endl;
        cout << "number_real " << number_in << endl;
    }
    else {
        cout << "The geometry you built is at least partially outside of space." << endl;
        cout << "number_origin " << number_origin << endl;
        cout << "number_out " << number_out << endl;
        cout << "number_real " << number_in << endl;
    }

    return geometry;
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, VectorXi* geometryPara_, VectorXd* Para_) {
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = *geometryPara_;
    Para = *Para_;
    
    

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel) {
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1+j)= (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = VectorXi::Zero(N);
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    //cout << "scope" << endl << scope << endl;
    //cout << "(Npara, Nparay, Nparaz) " << "(" << Nparax << ", " << Nparay << ", " << Nparaz << ")" << endl;
    //cout << "Npara: " << Npara << endl;
    Para = initial_diel_func(initial_diel, Npara);
    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos= paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
    }

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel_center, string initial_diel_ring, double r) {
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = VectorXi::Zero(N);
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;

    double xcenter, ycenter;
    xcenter = double(scope(0, 1) + scope(0, 0)) / 2;
    ycenter = double(scope(1, 1) + scope(1, 0)) / 2;
    
    cout << "xcenter"<<xcenter << endl;
    cout << "ycenter"<<ycenter << endl;
    Para = VectorXd::Zero(Npara);
    for (int i = 0; i <= Nparax - 1; i++) {
        for (int j = 0; j <= Nparay - 1; j++) {
            for (int k = 0; k <= Nparaz - 1; k++) {
                int pos = k + Nparaz * (j + Nparay * i);
                double x, y;      //center of one 2D para region
                x = bind(0) * (2 * i + 1) / 2;
                y = bind(1) * (2 * j + 1) / 2;
                double rx, ry;
                rx = x - xcenter;
                ry = y - ycenter;
                if (rx * rx + ry * ry < r * r) {
                    Para(pos) = initial_diel_func(initial_diel_center);
                }
                else {
                    Para(pos) = initial_diel_func(initial_diel_ring);
                }

            }
        }
    }


    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
    }

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel_background, list<string>* initial_diel_list, list<double>* r_list, list<Vector2d>* center_list){
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = VectorXi::Zero(N);
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    Para = VectorXd::Zero(Npara);

    for (int i = 0; i <= Npara - 1; i++) {                          //First inital all to background
        Para(i) = initial_diel_func(initial_diel_background);
    }

    list<string>::iterator it_diel = (*initial_diel_list).begin();
    list<double>::iterator it_r = (*r_list).begin();
    list<Vector2d>::iterator it_center = (*center_list).begin();

    for (int itnum = 0; itnum <= (*initial_diel_list).size() - 1; itnum++) {
        string initial_diel = (*it_diel);
        double r = (*it_r);
        Vector2d center = (*it_center);
        double xcenter, ycenter;
        xcenter = center(0);
        ycenter = center(1);
        cout << "xcenter" << xcenter << endl;
        cout << "ycenter" << ycenter << endl;
 

        for (int i = 0; i <= Nparax - 1; i++) {
            for (int j = 0; j <= Nparay - 1; j++) {
                for (int k = 0; k <= Nparaz - 1; k++) {
                    int pos = k + Nparaz * (j + Nparay * i);
                    double x, y;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;
                    double rx, ry;
                    rx = x - xcenter;
                    ry = y - ycenter;
                    if (rx * rx + ry * ry < (r+0.1) * (r+0.1)) {
                        Para(pos) = initial_diel_func(initial_diel);
                    }
                }
            }
        }

        it_diel++;
        it_r++;
        it_center++;


    }
    

    
    
    


    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
    }

}

void SpacePara::ChangeBind(Vector3i bind_) {
    bind = bind_;
    VectorXi geometryPara_before = geometryPara;
    VectorXd Para_before = Para;
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil((scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil((scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil((scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;

    vector<vector<int>> Paratogeometry(Npara, vector<int>(2, 0));
    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
        Paratogeometry[pos][0] += 1;
        Paratogeometry[pos][1] += Para_before(geometryPara_before(i));
    }

    for (int i = 0; i <= Npara - 1; i++) {
        Para(i) = Paratogeometry[i][1] / Paratogeometry[i][0];
    }

    


}

Space* SpacePara::get_space() {
    return space;
}

VectorXi SpacePara::get_geometry() {
    return geometry;
}

VectorXi* SpacePara::get_geometryPara() {
    return &geometryPara;
}

VectorXd* SpacePara::get_Para() {
    return &Para;
}

Vector3i* SpacePara::get_bind() {
    return &bind;
}