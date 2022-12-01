#include <iostream>

#include "ObjDDAModel.h"

double SmoothDensity(double input, double ita, double beta) {
    double result = 0.0;
    if (input <= ita && input >= 0.0) {
        return ita * (exp(-beta * (1 - input / ita)) - (1 - input / ita) * exp(-beta));
    }
    else if (input > ita && input <= 1.0) {
        return (1 - ita) * (1 - exp(-beta * (input - ita) / (1 - ita)) + (input - ita) / (1 - ita) * exp(-beta)) + ita;
    }
    else {
        cout << "ERROR: FilterOption::SmoothDensity(double input)--input out of range" << endl;
        throw 1;
    }
}


ObjPointEDDAModel::ObjPointEDDAModel(vector<double> parameters, DDAModel *model_){
    VectorXd PointEParameters = VectorXd::Zero((parameters).size());
    //auto it=(parameters).begin();
    for(int i=0;i<=int((parameters).size()-1);i++){
        PointEParameters(i) = parameters[i];
    }
    //Have_Penalty = HavePenalty_;
    x=PointEParameters(0);
    y=PointEParameters(1);
    z=PointEParameters(2);
    Have_Devx = false;
    model = model_;
    //evomodel = evomodel_;
    AProductCore* Core = (*model).get_Core();
    d = (*Core).get_d();
    N = (*Core).get_N();
    P = (*model).get_P();
    R = (*Core).get_R();
    Vector3d n_E0 = (*model).get_nE0();
    Vector3d n_K = (*model).get_nK();
    double E0 = (*model).get_E0();
    double lam = (*Core).get_lam();
    cout << "lam" << lam << endl;
    double K = (*Core).get_K();
    E_sum = Vector3cd::Zero();                                                                             //是不是E_sum忘了加E_ext了？ It is actually in Rest.
    E_ext = Vector3cd::Zero();
    E_ext(0) = E0*n_E0(0)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
    E_ext(1) = E0*n_E0(1)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
    E_ext(2) = E0*n_E0(2)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);   
    // cout << R(3*5444+2) << "this" << endl;
}

void ObjPointEDDAModel::SingleResponse(int idx, bool deduction){
    //VectorXcd P = model->get_P();
    //VectorXi R = model->get_R();
    double rx=x-d*(*R)(3*idx);                  //R has no d in it, so needs to time d
    double ry=y-d*(*R)(3*idx+1);
    double rz=z-d*(*R)(3*idx+2);
    //cout << rx << "," << ry << "," << rz << idx << endl;
    AProductCore* Core = (*model).get_Core();
    Matrix3cd A=(*Core).A_dic_generator(rx,ry,rz);
    if (deduction == false){
        E_sum(0)-=(A(0,0)*(*P)(3*idx)+A(0,1)*(*P)(3*idx+1)+A(0,2)*(*P)(3*idx+2));
        E_sum(1)-=(A(1,0)*(*P)(3*idx)+A(1,1)*(*P)(3*idx+1)+A(1,2)*(*P)(3*idx+2));
        E_sum(2)-=(A(2,0)*(*P)(3*idx)+A(2,1)*(*P)(3*idx+1)+A(2,2)*(*P)(3*idx+2));
    }
    else{
        E_sum(0)+=(A(0,0)*(*P)(3*idx)+A(0,1)*(*P)(3*idx+1)+A(0,2)*(*P)(3*idx+2));
        E_sum(1)+=(A(1,0)*(*P)(3*idx)+A(1,1)*(*P)(3*idx+1)+A(1,2)*(*P)(3*idx+2));
        E_sum(2)+=(A(2,0)*(*P)(3*idx)+A(2,1)*(*P)(3*idx+1)+A(2,2)*(*P)(3*idx+2));    
    }
}

double ObjPointEDDAModel::GroupResponse(){
    /*if (Have_Penalty){
        return (E_sum).norm()-(*evomodel).L1Norm();
    }*/
    /*else{*/
    return (E_sum).norm();
    /*}*/
    
}

double ObjPointEDDAModel::GetVal(){
    Reset();
    for (int idx=0;idx<N;idx++){
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void ObjPointEDDAModel::Reset(){
    E_sum(0) = E_ext(0);
    E_sum(1) = E_ext(1);
    E_sum(2) = E_ext(2);
}


ObjIntegratedEDDAModel::ObjIntegratedEDDAModel(vector<double> parameters, DDAModel* model_) {
    VectorXd PointEParameters = VectorXd::Zero((parameters).size());
    //auto it=(parameters).begin();
    for (int i = 0; i <= int((parameters).size() - 1); i++) {
        PointEParameters(i) = parameters[i];
    }
    //Have_Penalty = HavePenalty_;
    powNum = int(round(PointEParameters(0)));
    xMin = int(round(PointEParameters(1)));
    xMax = int(round(PointEParameters(2)));
    yMin = int(round(PointEParameters(3)));
    yMax = int(round(PointEParameters(4)));
    zMin = int(round(PointEParameters(5)));
    zMax = int(round(PointEParameters(6)));
    ita = PointEParameters(7);
    beta = PointEParameters(8);
    Have_Devx = true;
    model = model_;
    //evomodel = evomodel_;
    AProductCore* Core = (*model).get_Core();
    d = (*Core).get_d();
    N = (*Core).get_N();
    P = (*model).get_P();
    R = (*Core).get_R();
    al = (*model).get_al();
    diel_old = (*Core).get_diel_old();
    Vector3d n_E0 = (*model).get_nE0();
    Vector3d n_K = (*model).get_nK();
    double E0 = (*model).get_E0();
    double lam = (*Core).get_lam();
    cout << "lam" << lam << endl;
    double K = (*Core).get_K();
    E_int = 0.0;
    // cout << R(3*5444+2) << "this" << endl;
}

void ObjIntegratedEDDAModel::SingleResponse(int idx, bool deduction) {
    //VectorXcd P = model->get_P();
    //VectorXi R = model->get_R();
    if ((xMin <= (*R)(3 * idx) <= xMax)&&(yMin <= (*R)(3 * idx + 1) <= yMax)&&(zMin <= (*R)(3 * idx + 2) <= zMax)) {
        double factor = SmoothDensity((*diel_old)(3 * idx), ita, beta);
        if (deduction == false) {
            E_int += factor * pow(abs((*al)(3 * idx) * (*P)(3 * idx)), powNum);
            E_int += factor * pow(abs((*al)(3 * idx + 1) * (*P)(3 * idx + 1)), powNum);
            E_int += factor * pow(abs((*al)(3 * idx + 2) * (*P)(3 * idx + 2)), powNum);
        }
        else {
            E_int -= factor * pow(abs((*al)(3 * idx) * (*P)(3 * idx)), powNum);
            E_int -= factor * pow(abs((*al)(3 * idx + 1) * (*P)(3 * idx + 1)), powNum);
            E_int -= factor * pow(abs((*al)(3 * idx + 2) * (*P)(3 * idx + 2)), powNum);
        }
    }
    
    
    
}

double ObjIntegratedEDDAModel::GroupResponse() {
    /*if (Have_Penalty){
        return (E_sum).norm()-(*evomodel).L1Norm();
    }*/
    /*else{*/
    return log(E_int);
    /*}*/

}

double ObjIntegratedEDDAModel::GetVal() {
    Reset();
    for (int idx = 0; idx < N; idx++) {
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void ObjIntegratedEDDAModel::Reset() {
    E_int = 0.0;
}





//ObjIntegratedEDDAModel::ObjIntegratedEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    Have_Penalty = HavePenalty_;
//    Have_Devx = false;
//    model = model_;
//    evomodel = evomodel_;
//    AProductCore* Core = (*model).get_Core();
//    N = (*Core).get_N();
//    Nx = (*Core).get_Nx();
//    Ny = (*Core).get_Ny();
//    Nz = (*Core).get_Nz();
//    d = (*Core).get_d();
//    al = (*model).get_al();
//    P = (*model).get_P();
//    E = VectorXcd::Zero(N * 3);
//    R = (*Core).get_R();
//    E_int = 0;
//}
//
//void ObjIntegratedEDDAModel::SingleResponse(int idx, bool deduction) {
//    return;
//}
//
//double ObjIntegratedEDDAModel::GroupResponse() {
//    for (int i = 0; i < N; i++) {
//        E(i * 3) = (*al)(i * 3) * (*P)(i * 3);
//        E(i * 3 + 1) = (*al)(i * 3 + 1) * (*P)(i * 3 + 1);
//        E(i * 3 + 2) = (*al)(i * 3 + 2) * (*P)(i * 3 + 2);
//    }
//    E_int = 0;
//    double diel_sum = 0;
//    
//    for (int i = 0; i < N; i++) {
//        if ((*R)(3*i+2) >= 27) {
//            double E_sum_temp = 0;
//            for (int j = 0; j < 3; j++) {
//                E_sum_temp += pow(abs(E(3 * i + j)), 2);
//            }
//            diel_sum += (*((*model).get_Core()->get_diel_old()))(i * 3);
//            E_int += pow(E_sum_temp, 2) * ((*((*model).get_Core()->get_diel_old()))(i * 3) + 0.0001) / 4.0; //prevent nan result for devp calculation.
//        }
//    }
//    
//    E_int = log(E_int);
//
//
//    // E_int /= diel_sum;
//    if (Have_Penalty) {
//        return E_int - (*evomodel).L1Norm();
//    }
//    else {
//        return E_int;
//    }
//
//}
//
//double ObjIntegratedEDDAModel::GetVal() {
//    Reset();
//    for (int i = 0; i < N; i++) {
//        E(i * 3) = (*al)(i * 3) * (*P)(i * 3);
//        E(i * 3 + 1) = (*al)(i * 3 + 1) * (*P)(i * 3 + 1);
//        E(i * 3 + 2) = (*al)(i * 3 + 2) * (*P)(i * 3 + 2);
//    }
//    return GroupResponse();
//}
//
//void ObjIntegratedEDDAModel::Reset() {
//    E_int = 0;
//    E = VectorXcd::Zero(N * 3);
//}



//ObjPointListEDDAModel::ObjPointListEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    VectorXd PointEParameters = VectorXd::Zero((parameters).size());
//    list<double>::iterator it = (parameters).begin();
//    if ((parameters.size() % 3) != 0) {
//        cout << "ERROR: (parameters.size() % 3) != 0" << endl;
//    }
//    PNum = int(round(parameters.size() / 3));
//    x = VectorXd::Zero(PNum);
//    y = VectorXd::Zero(PNum);
//    z = VectorXd::Zero(PNum);
//    E_sum.resize(PNum, 3);
//    E_ext.resize(PNum, 3);
//    for (int i = 0; i <= PNum - 1; i++) {
//        x(i) = (*it);
//        it++;
//        y(i) = (*it);
//        it++;
//        z(i) = (*it);
//        it++;
//    }
//    Have_Penalty = HavePenalty_;
//    
//    Have_Devx = false;
//    model = model_;
//    evomodel = evomodel_;
//    AProductCore* Core = (*model).get_Core();
//    d = (*Core).get_d();
//    N = (*Core).get_N();
//    P = (*model).get_P();
//    R = (*Core).get_R();
//    Vector3d n_E0 = (*model).get_nE0();
//    Vector3d n_K = (*model).get_nK();
//    double E0 = (*model).get_E0();
//    double lam = (*Core).get_lam();
//    cout << "lam" << lam << endl;
//    double K = (*Core).get_K();
//
//    for (int i = 0; i <= PNum - 1; i++) {
//        E_ext(i, 0) = E0 * n_E0(0) * (cos(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) + sin(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) * 1i);
//        E_ext(i, 1) = E0 * n_E0(1) * (cos(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) + sin(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) * 1i);
//        E_ext(i, 2) = E0 * n_E0(2) * (cos(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) + sin(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) * 1i);
//    }
//    
//    // cout << R(3*5444+2) << "this" << endl;
//}
//
//void ObjPointListEDDAModel::SingleResponse(int idx, bool deduction) {
//    //VectorXcd P = model->get_P();
//    //VectorXi R = model->get_R();
//    AProductCore* Core = (*model).get_Core();
//    for (int i = 0; i <= PNum - 1; i++) {
//        double rx = x(i) - d * (*R)(3 * idx);                  //R has no d in it, so needs to time d
//        double ry = y(i) - d * (*R)(3 * idx + 1);
//        double rz = z(i) - d * (*R)(3 * idx + 2);
//        Matrix3cd A = (*Core).A_dic_generator(rx, ry, rz);
//        if (deduction == false) {
//            E_sum(i, 0) -= (A(0, 0) * (*P)(3 * idx) + A(0, 1) * (*P)(3 * idx + 1) + A(0, 2) * (*P)(3 * idx + 2));
//            E_sum(i, 1) -= (A(1, 0) * (*P)(3 * idx) + A(1, 1) * (*P)(3 * idx + 1) + A(1, 2) * (*P)(3 * idx + 2));
//            E_sum(i, 2) -= (A(2, 0) * (*P)(3 * idx) + A(2, 1) * (*P)(3 * idx + 1) + A(2, 2) * (*P)(3 * idx + 2));
//        }
//        else {
//            E_sum(i, 0) += (A(0, 0) * (*P)(3 * idx) + A(0, 1) * (*P)(3 * idx + 1) + A(0, 2) * (*P)(3 * idx + 2));
//            E_sum(i, 1) += (A(1, 0) * (*P)(3 * idx) + A(1, 1) * (*P)(3 * idx + 1) + A(1, 2) * (*P)(3 * idx + 2));
//            E_sum(i, 2) += (A(2, 0) * (*P)(3 * idx) + A(2, 1) * (*P)(3 * idx + 1) + A(2, 2) * (*P)(3 * idx + 2));
//        }
//    }
//
//    
//    //cout << rx << "," << ry << "," << rz << idx << endl;
//    
//    
//    
//}
//
//double ObjPointListEDDAModel::GroupResponse() {
//    double result = 0;
//    if (Have_Penalty) {
//        for (int i = 0; i <= PNum - 1; i++) {
//            Vector3cd tmp;
//            tmp(0) = E_sum(i, 0);
//            tmp(1) = E_sum(i, 1);
//            tmp(2) = E_sum(i, 2);
//            //result = result+ sqrt(pow(tmp(0).real(), 2) + pow(tmp(0).imag(), 2)+ pow(tmp(1).real(), 2) + pow(tmp(1).imag(), 2)+ pow(tmp(2).real(), 2) + pow(tmp(2).imag(), 2));
//            result += tmp.norm();
//        }
//        
//        return result/PNum - (*evomodel).L1Norm();
//    }
//    else {
//        for (int i = 0; i <= PNum - 1; i++) {
//            Vector3cd tmp;
//            tmp(0) = E_sum(i, 0);
//            tmp(1) = E_sum(i, 1);
//            tmp(2) = E_sum(i, 2);
//            //result = result + sqrt(pow(tmp(0).real(), 2) + pow(tmp(0).imag(), 2) + pow(tmp(1).real(), 2) + pow(tmp(1).imag(), 2) + pow(tmp(2).real(), 2) + pow(tmp(2).imag(), 2));
//            result += tmp.norm();
//        }
//        return result / PNum;
//    }
//
//}
//
//double ObjPointListEDDAModel::GetVal() {
//    Reset();
//    for (int idx = 0; idx < N; idx++) {
//        SingleResponse(idx, false);
//        // cout << E_sum(0) << endl;
//    }
//    return GroupResponse();
//}
//
//void ObjPointListEDDAModel::Reset() {
//    for (int i = 0; i <= PNum - 1; i++) {
//        E_sum(i, 0) = E_ext(i, 0);
//        E_sum(i, 1) = E_ext(i, 1);
//        E_sum(i, 2) = E_ext(i, 2);
//    }
//    
//}
//
//
//
//
//ObjPointIDDAModel::ObjPointIDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    VectorXd PointEParameters = VectorXd::Zero((parameters).size());
//    list<double>::iterator it = (parameters).begin();
//    for (int i = 0; i <= int((parameters).size() - 1); i++) {
//        PointEParameters(i) = (*it);
//        it++;
//    }
//    Have_Penalty = HavePenalty_;
//    x = PointEParameters(0);
//    y = PointEParameters(1);
//    z = PointEParameters(2);
//    Have_Devx = false;
//    model = model_;
//    evomodel = evomodel_;
//    AProductCore* Core = (*model).get_Core();
//    d = (*Core).get_d();
//    N = (*Core).get_N();
//    P = (*model).get_P();
//    R = (*Core).get_R();
//    Vector3d n_E0 = (*model).get_nE0();
//    Vector3d n_K = (*model).get_nK();
//    double E0 = (*model).get_E0();
//    double lam = (*Core).get_lam();
//    cout << "lam" << lam << endl;
//    double K = (*Core).get_K();
//    E_sum = Vector3cd::Zero();
//    E_ext = Vector3cd::Zero();
//    E_ext(0) = E0 * n_E0(0) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
//    E_ext(1) = E0 * n_E0(1) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
//    E_ext(2) = E0 * n_E0(2) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
//    // cout << R(3*5444+2) << "this" << endl;
//}
//
//void ObjPointIDDAModel::SingleResponse(int idx, bool deduction) {
//    //VectorXcd P = model->get_P();
//    //VectorXi R = model->get_R();
//    double rx = x - d * (*R)(3 * idx);                  //R has no d in it, so needs to time d
//    double ry = y - d * (*R)(3 * idx + 1);
//    double rz = z - d * (*R)(3 * idx + 2);
//    //cout << rx << "," << ry << "," << rz << idx << endl;
//    AProductCore* Core = (*model).get_Core();
//    Matrix3cd A = (*Core).A_dic_generator(rx, ry, rz);
//    if (deduction == false) {
//        E_sum(0) -= (A(0, 0) * (*P)(3 * idx) + A(0, 1) * (*P)(3 * idx + 1) + A(0, 2) * (*P)(3 * idx + 2));
//        E_sum(1) -= (A(1, 0) * (*P)(3 * idx) + A(1, 1) * (*P)(3 * idx + 1) + A(1, 2) * (*P)(3 * idx + 2));
//        E_sum(2) -= (A(2, 0) * (*P)(3 * idx) + A(2, 1) * (*P)(3 * idx + 1) + A(2, 2) * (*P)(3 * idx + 2));
//    }
//    else {
//        E_sum(0) += (A(0, 0) * (*P)(3 * idx) + A(0, 1) * (*P)(3 * idx + 1) + A(0, 2) * (*P)(3 * idx + 2));
//        E_sum(1) += (A(1, 0) * (*P)(3 * idx) + A(1, 1) * (*P)(3 * idx + 1) + A(1, 2) * (*P)(3 * idx + 2));
//        E_sum(2) += (A(2, 0) * (*P)(3 * idx) + A(2, 1) * (*P)(3 * idx + 1) + A(2, 2) * (*P)(3 * idx + 2));
//    }
//}
//
//double ObjPointIDDAModel::GroupResponse() {
//    if (Have_Penalty) {
//        return pow((E_sum).norm(),2) - (*evomodel).L1Norm();
//    }
//    else {
//        return pow((E_sum).norm(),2);
//    }
//
//}
//
//double ObjPointIDDAModel::GetVal() {
//    Reset();
//    for (int idx = 0; idx < N; idx++) {
//        SingleResponse(idx, false);
//        // cout << E_sum(0) << endl;
//    }
//    return GroupResponse();
//}
//
//void ObjPointIDDAModel::Reset() {
//    E_sum(0) = E_ext(0);
//    E_sum(1) = E_ext(1);
//    E_sum(2) = E_ext(2);
//}
//
//
//

//
//
//
//ObjMidAvgEDDAModel::ObjMidAvgEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    VectorXd PointEParameters = VectorXd::Zero((parameters).size());
//    list<double>::iterator it = (parameters).begin();
//    for (int i = 0; i <= int((parameters).size() - 1); i++) {
//        PointEParameters(i) = (*it);
//        it++;
//    }
//    r = PointEParameters(0);
//    
//    Have_Penalty = HavePenalty_;
//    Have_Devx = false;
//    model = model_;
//    evomodel = evomodel_;
//    AProductCore* Core = (*model).get_Core();
//    N = (*Core).get_N();
//    Nx = (*Core).get_Nx();
//    Ny = (*Core).get_Ny();
//    Nz = (*Core).get_Nz();
//    d = (*Core).get_d();
//    al = (*model).get_al();
//    P = (*model).get_P();
//    E = VectorXcd::Zero(N * 3);
//    R = (*Core).get_R();
//
//}
//
//void ObjMidAvgEDDAModel::SingleResponse(int idx, bool deduction) {
//    return;
//}
//
//double ObjMidAvgEDDAModel::GroupResponse() {
//    for (int i = 0; i < N; i++) {
//        E(i * 3) = (*al)(i * 3) * (*P)(i * 3);
//        E(i * 3 + 1) = (*al)(i * 3 + 1) * (*P)(i * 3 + 1);
//        E(i * 3 + 2) = (*al)(i * 3 + 2) * (*P)(i * 3 + 2);
//    }
//
//    double avg_num = 0;
//
//    for (int i = 0; i < N; i++) {
//        if ((48<=(*R)(3 * i) * d<=118 && 48 <= (*R)(3 * i+1) * d <= 118) || ((212 <= (*R)(3 * i) * d <= 282 && 212 <= (*R)(3 * i + 1) * d <= 282))) {
//            double E_sum_temp = 0;
//            
//            for (int j = 0; j < 3; j++) {
//                E_sum_temp += pow(abs(E(3 * i + j)), 2);
//            }
//            E_sum_temp = sqrt(E_sum_temp);
//
//            E_avg += E_sum_temp * (pow((*((*model).get_Core()->get_diel_old()))(i * 3), 1) + 0.0001);
//            avg_num += 1;
//        }
//    }
//
//    E_avg = E_avg;
//
//
//    // E_int /= diel_sum;
//    if (Have_Penalty) {
//        return E_avg - (*evomodel).L1Norm();
//    }
//    else {
//        return E_avg;
//    }
//
//}
//
//double ObjMidAvgEDDAModel::GetVal() {
//    Reset();
//    for (int i = 0; i < N; i++) {
//        E(i * 3) = (*al)(i * 3) * (*P)(i * 3);
//        E(i * 3 + 1) = (*al)(i * 3 + 1) * (*P)(i * 3 + 1);
//        E(i * 3 + 2) = (*al)(i * 3 + 2) * (*P)(i * 3 + 2);
//    }
//    return GroupResponse();
//}
//
//void ObjMidAvgEDDAModel::Reset() {
//    E_avg = 0;
//    E = VectorXcd::Zero(N * 3);
//}
//
//
//Objscattering0D::Objscattering0D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    Paralength = (parameters).size();
//    VectorXd FOMParameters = VectorXd::Zero(Paralength);
//    list<double>::iterator it = (parameters).begin();
//    for (int i = 0; i <= int(Paralength - 1); i++) {
//        FOMParameters(i) = (*it);
//        it++;
//    }
//    if (Paralength % 3 != 0) {
//        cout << "FOMscattering2D ERROR: parameter must be times of 3." << endl;
//        throw 1;
//    }
//
//    Have_Penalty = HavePenalty_;
//    Have_Devx = false;
//    model = model_;
//    evomodel = evomodel_;
//
//    AProductCore* Core = (*model).get_Core();
//    d = (*Core).get_d();
//    N = (*Core).get_N();                   //Number of dipoles
//    P = (*model).get_P();
//    R = (*Core).get_R();
//    Vector3d n_E0 = (*model).get_nE0();
//    Vector3d n_K = (*model).get_nK();
//    E0 = (*model).get_E0();
//    double lam = (*Core).get_lam();
//    cout << "lam" << lam << endl;
//    K = (*Core).get_K();
//    
//
//    for (int i = 0; i <= int(round(Paralength / 3) - 1); i++) {
//        Vector3d n_K_tmp;
//        n_K_tmp(0) = FOMParameters(3 * i);
//        n_K_tmp(1) = FOMParameters(3 * i + 1);
//        n_K_tmp(2) = FOMParameters(3 * i + 2);
//        n_K_s_l.push_back(n_K_tmp);
//
//        Matrix3d FconstM;
//        double nkx = n_K_tmp(0);
//        double nky = n_K_tmp(1);
//        double nkz = n_K_tmp(2);
//        double K3 = pow(K, 3);
//        FconstM(0, 0) = K3 * (1 - nkx * nkx);
//        FconstM(0, 1) = -K3 * nkx * nky;
//        FconstM(0, 2) = -K3 * nkx * nkz;
//        FconstM(1, 1) = K3 * (1 - nky * nky);
//        FconstM(1, 2) = -K3 * nky * nkz;
//        FconstM(2, 2) = K3 * (1 - nkz * nkz);
//        FconstM(1, 0) = FconstM(0, 1);
//        FconstM(2, 0) = FconstM(0, 2);
//        FconstM(2, 1) = FconstM(1, 2);
//        FconstM_l.push_back(FconstM);
//
//        Vector3cd PSum_tmp;
//        PSum_tmp = Vector3cd::Zero();
//        PSum_l.push_back(PSum_tmp);
//        
//    }
//
//    
//
//
//
//}
//
//void Objscattering0D::SingleResponse(int idx, bool deduction) {
//    //list<Matrix3d>::iterator it1 = (FconstM_l).begin();
//    list<Vector3d>::iterator it2 = (n_K_s_l).begin();
//    list<Vector3cd>::iterator it3 = (PSum_l).begin();
//    for (int i = 0; i <= int(round(Paralength / 3) - 1); i++) {
//        //Matrix3d FconstM = (*it1);
//        Vector3d n_K_s = (*it2);
//        double nkx = n_K_s(0);
//        double nky = n_K_s(1);
//        double nkz = n_K_s(2);
//        double phaseterm = -d * K * (nkx * ((*R)(3 * idx) - 32.5) + nky * ((*R)(3 * idx + 1) - 32.5) + nkz * ((*R)(3 * idx + 2) - 32.5));  //From equation 17. Time term will be eliminated
//        complex<double> phase = cos(phaseterm) + sin(phaseterm) * 1i;
//        if (deduction == false) {
//            (*it3)(0) += (*P)(3 * idx) * phase;
//            (*it3)(1) += (*P)(3 * idx + 1) * phase;
//            (*it3)(2) += (*P)(3 * idx + 2) * phase;
//
//        }
//        else {
//            (*it3)(0) -= (*P)(3 * idx) * phase;
//            (*it3)(1) -= (*P)(3 * idx + 1) * phase;
//            (*it3)(2) -= (*P)(3 * idx + 2) * phase;
//        }
//        it2++;
//        it3++;
//    }
//}
//
//double Objscattering0D::GroupResponse() {
//    if (Have_Penalty) {
//        return log10(this->FTUCnsquare() / (pow(K, 2) * pow(E0, 2))) - (*evomodel).L1Norm();
//    }
//    else {
//        return log10(this->FTUCnsquare() / (pow(K, 2) * pow(E0, 2)));                           //K does not depend on scattering angle, so it is fine to divide it here.
//    }
//
//}
//
//double Objscattering0D::GetVal() {
//    Reset();
//    for (int idx = 0; idx < N; idx++) {
//        SingleResponse(idx, false);
//    }
//    return GroupResponse();
//}
//
//void Objscattering0D::Reset() {
//    for (list<Vector3cd>::iterator it = PSum_l.begin(); it != PSum_l.end(); it++) {
//        (*it) = Vector3cd::Zero();
//    }
//}
//
//double Objscattering0D::FTUCnsquare() {
//
//    list<Matrix3d>::iterator it1 = (FconstM_l).begin();
//    list<Vector3cd>::iterator it2 = (PSum_l).begin();
//    double result = 0.0;
//    for (int i = 0; i <= int(round(Paralength / 3) - 1); i++) {
//        Vector3cd FTUC;
//        FTUC(0) = (*it1)(0, 0) * (*it2)(0) + (*it1)(0, 1) * (*it2)(1) + (*it1)(0, 2) * (*it2)(2);
//        FTUC(1) = (*it1)(1, 0) * (*it2)(0) + (*it1)(1, 1) * (*it2)(1) + (*it1)(1, 2) * (*it2)(2);
//        FTUC(2) = (*it1)(2, 0) * (*it2)(0) + (*it1)(2, 1) * (*it2)(1) + (*it1)(2, 2) * (*it2)(2);
//
//
//        it1++;
//        it2++;
//        result += norm(FTUC(0)) + norm(FTUC(1)) + norm(FTUC(2));                                 //In C++, norm is the square of magnitude.
//        
//
//    }
//
//    return result/double(round(Paralength / 3));
//    
//}
//
//
//
//
//
//
//
//
//
//
//Objscattering2D::Objscattering2D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    Paralength = (parameters).size();
//    VectorXd FOMParameters = VectorXd::Zero(Paralength);
//    list<double>::iterator it = (parameters).begin();
//    for (int i = 0; i <= int(Paralength - 1); i++) {
//        FOMParameters(i) = (*it);
//        it++;
//    }
//    if (Paralength % 2 != 0) {
//        cout << "Objscattering2D ERROR: parameter must be times of 2." << endl;
//        throw 1;
//    }
//
//    Have_Penalty = HavePenalty_;
//    Have_Devx = false;
//    model = model_;
//    evomodel = evomodel_;
//
//    AProductCore* Core = (*model).get_Core();
//    d = (*Core).get_d();
//    N = (*Core).get_N();                   //Number of dipoles
//    P = (*model).get_P();
//    R = (*Core).get_R();
//    Vector3d n_E0 = (*model).get_nE0();
//    Vector3d n_K = (*model).get_nK();
//    E0 = (*model).get_E0();
//    double lam = (*Core).get_lam();
//    cout << "lam" << lam << endl;
//    K = (*Core).get_K();
//
//    double Lm = (*Core).get_Lm();
//    double Ln = (*Core).get_Ln();
//    ATUC = Lm * Ln;
//
//
//    for (int i = 0; i <= int(round(Paralength / 2) - 1); i++) {
//        Vector3d n_K_tmp;
//        n_K_tmp(0) = 2 * M_PI * FOMParameters(2 * i) / (Lm * K) + n_K(0);
//        n_K_tmp(1) = 2 * M_PI * FOMParameters(2 * i + 1) / (Ln * K) + n_K(1);
//        n_K_tmp(2) = sqrt(1 - pow(n_K_tmp(0), 2) - pow(n_K_tmp(1), 2));
//        n_K_s_l.push_back(n_K_tmp);
//        cout << n_K_tmp << endl;
//        Matrix3d FconstM;
//        double nkx = n_K_tmp(0);
//        double nky = n_K_tmp(1);
//        double nkz = n_K_tmp(2);
//        double K3 = pow(K, 3);
//        FconstM(0, 0) = K3 * (1 - nkx * nkx);
//        FconstM(0, 1) = -K3 * nkx * nky;
//        FconstM(0, 2) = -K3 * nkx * nkz;
//        FconstM(1, 1) = K3 * (1 - nky * nky);
//        FconstM(1, 2) = -K3 * nky * nkz;
//        FconstM(2, 2) = K3 * (1 - nkz * nkz);
//        FconstM(1, 0) = FconstM(0, 1);
//        FconstM(2, 0) = FconstM(0, 2);
//        FconstM(2, 1) = FconstM(1, 2);
//        FconstM_l.push_back(FconstM);
//
//        Vector3cd PSum_tmp;
//        PSum_tmp = Vector3cd::Zero();
//        PSum_l.push_back(PSum_tmp);
//
//    }
//
//
//
//
//
//}
//
//void Objscattering2D::SingleResponse(int idx, bool deduction) {
//    //list<Matrix3d>::iterator it1 = (FconstM_l).begin();
//    list<Vector3d>::iterator it2 = (n_K_s_l).begin();
//    list<Vector3cd>::iterator it3 = (PSum_l).begin();
//    for (int i = 0; i <= int(round(Paralength / 2) - 1); i++) {
//        //Matrix3d FconstM = (*it1);
//        Vector3d n_K_s = (*it2);
//        double nkx = n_K_s(0);
//        double nky = n_K_s(1);
//        double nkz = n_K_s(2);
//        double phaseterm = -d * K * (nkx * ((*R)(3 * idx)) + nky * ((*R)(3 * idx + 1)) + nkz * ((*R)(3 * idx + 2)));  //From equation 17. Time term will be eliminated
//        complex<double> phase = cos(phaseterm) + sin(phaseterm) * 1i;
//        if (deduction == false) {
//            (*it3)(0) += (*P)(3 * idx) * phase;
//            (*it3)(1) += (*P)(3 * idx + 1) * phase;
//            (*it3)(2) += (*P)(3 * idx + 2) * phase;
//
//        }
//        else {
//            (*it3)(0) -= (*P)(3 * idx) * phase;
//            (*it3)(1) -= (*P)(3 * idx + 1) * phase;
//            (*it3)(2) -= (*P)(3 * idx + 2) * phase;
//        }
//        it2++;
//        it3++;
//    }
//}
//
//double Objscattering2D::GroupResponse() {
//    if (Have_Penalty) {
//        return -(this->FTUCnsquareoversinal() * (pow(2 * M_PI, 2)) / (pow(K, 4) * pow(E0, 2) * pow(ATUC, 2)) - (*evomodel).L1Norm());
//    }
//    else {
//        return -this->FTUCnsquareoversinal() * (pow(2 * M_PI, 2)) / (pow(K, 4) * pow(E0, 2) * pow(ATUC, 2));                           //K does not depend on scattering angle, so it is fine to divide it here.
//    }
//
//}
//
//double Objscattering2D::GetVal() {
//    Reset();
//    for (int idx = 0; idx < N; idx++) {
//        SingleResponse(idx, false);
//    }
//    return GroupResponse();
//}
//
//void Objscattering2D::Reset() {
//    for (list<Vector3cd>::iterator it = PSum_l.begin(); it != PSum_l.end(); it++) {
//        (*it) = Vector3cd::Zero();
//    }
//}
//
//double Objscattering2D::FTUCnsquareoversinal() {
//
//    list<Matrix3d>::iterator it1 = (FconstM_l).begin();
//    list<Vector3cd>::iterator it2 = (PSum_l).begin();
//    list<Vector3d>::iterator it3 = n_K_s_l.begin();
//    double result = 0.0;
//    for (int i = 0; i <= int(round(Paralength / 2) - 1); i++) {
//        double ksz = (*it3)(2);
//        Vector3cd FTUC;
//        FTUC(0) = (*it1)(0, 0) * (*it2)(0) + (*it1)(0, 1) * (*it2)(1) + (*it1)(0, 2) * (*it2)(2);
//        FTUC(1) = (*it1)(1, 0) * (*it2)(0) + (*it1)(1, 1) * (*it2)(1) + (*it1)(1, 2) * (*it2)(2);
//        FTUC(2) = (*it1)(2, 0) * (*it2)(0) + (*it1)(2, 1) * (*it2)(1) + (*it1)(2, 2) * (*it2)(2);
//
//
//        it1++;
//        it2++;
//        it3++;
//        result += (norm(FTUC(0)) + norm(FTUC(1)) + norm(FTUC(2))) / (ksz * ksz);                                 //In C++, norm is the square of magnitude.
//
//
//    }
//
//    return result / double(round(Paralength / 2));
//
//}
//
//
//
//
//Objreflect2D::Objreflect2D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    Paralength = (parameters).size();
//    VectorXd FOMParameters = VectorXd::Zero(Paralength);
//    list<double>::iterator it = (parameters).begin();
//    for (int i = 0; i <= int(Paralength - 1); i++) {
//        FOMParameters(i) = (*it);
//        it++;
//    }
//    if (Paralength % 2 != 0) {
//        cout << "Objreflect2D ERROR: parameter must be times of 2." << endl;
//        throw 1;
//    }
//
//    Have_Penalty = HavePenalty_;
//    Have_Devx = false;
//    model = model_;
//    evomodel = evomodel_;
//
//    AProductCore* Core = (*model).get_Core();
//    d = (*Core).get_d();
//    N = (*Core).get_N();                   //Number of dipoles
//    P = (*model).get_P();
//    R = (*Core).get_R();
//    Vector3d n_E0 = (*model).get_nE0();
//    Vector3d n_K = (*model).get_nK();
//    E0 = (*model).get_E0();
//    double lam = (*Core).get_lam();
//    cout << "lam" << lam << endl;
//    K = (*Core).get_K();
//
//    double Lm = (*Core).get_Lm();
//    double Ln = (*Core).get_Ln();
//    ATUC = Lm * Ln;
//
//
//    for (int i = 0; i <= int(round(Paralength / 2) - 1); i++) {
//        Vector3d n_K_tmp;
//        n_K_tmp(0) = 2 * M_PI * FOMParameters(2 * i) / (Lm * K) + n_K(0);
//        n_K_tmp(1) = 2 * M_PI * FOMParameters(2 * i + 1) / (Ln * K) + n_K(1);
//        n_K_tmp(2) = sqrt(1 - pow(n_K_tmp(0), 2) - pow(n_K_tmp(1), 2));
//        cout << n_K_tmp << endl;
//        n_K_s_l.push_back(n_K_tmp);
//
//        Matrix3d FconstM;
//        double nkx = n_K_tmp(0);
//        double nky = n_K_tmp(1);
//        double nkz = n_K_tmp(2);
//        double K3 = pow(K, 3);
//        FconstM(0, 0) = K3 * (1 - nkx * nkx);
//        FconstM(0, 1) = -K3 * nkx * nky;
//        FconstM(0, 2) = -K3 * nkx * nkz;
//        FconstM(1, 1) = K3 * (1 - nky * nky);
//        FconstM(1, 2) = -K3 * nky * nkz;
//        FconstM(2, 2) = K3 * (1 - nkz * nkz);
//        FconstM(1, 0) = FconstM(0, 1);
//        FconstM(2, 0) = FconstM(0, 2);
//        FconstM(2, 1) = FconstM(1, 2);
//        FconstM_l.push_back(FconstM);
//
//        Vector3cd PSum_tmp;
//        PSum_tmp = Vector3cd::Zero();
//        PSum_l.push_back(PSum_tmp);
//
//    }
//
//
//
//
//
//}
//
//void Objreflect2D::SingleResponse(int idx, bool deduction) {
//    //list<Matrix3d>::iterator it1 = (FconstM_l).begin();
//    list<Vector3d>::iterator it2 = (n_K_s_l).begin();
//    list<Vector3cd>::iterator it3 = (PSum_l).begin();
//    for (int i = 0; i <= int(round(Paralength / 2) - 1); i++) {
//        //Matrix3d FconstM = (*it1);
//        Vector3d n_K_s = (*it2);
//        double nkx = n_K_s(0);
//        double nky = n_K_s(1);
//        double nkz = n_K_s(2);
//        double phaseterm = -d * K * (nkx * ((*R)(3 * idx)) + nky * ((*R)(3 * idx + 1)) + nkz * ((*R)(3 * idx + 2)));  //From equation 17. Time term will be eliminated
//        complex<double> phase = cos(phaseterm) + sin(phaseterm) * 1i;
//        if (deduction == false) {
//            (*it3)(0) += (*P)(3 * idx) * phase;
//            (*it3)(1) += (*P)(3 * idx + 1) * phase;
//            (*it3)(2) += (*P)(3 * idx + 2) * phase;
//
//        }
//        else {
//            (*it3)(0) -= (*P)(3 * idx) * phase;
//            (*it3)(1) -= (*P)(3 * idx + 1) * phase;
//            (*it3)(2) -= (*P)(3 * idx + 2) * phase;
//        }
//        it2++;
//        it3++;
//    }
//}
//
//double Objreflect2D::GroupResponse() {
//    if (Have_Penalty) {
//        return log10(this->FTUCnsquareoversinal() * (pow(2 * M_PI, 2)) / (pow(K, 4) * pow(E0, 2) * pow(ATUC, 2))) - (*evomodel).L1Norm();
//    }
//    else {
//        return log10(this->FTUCnsquareoversinal() * (pow(2 * M_PI, 2)) / (pow(K, 4) * pow(E0, 2) * pow(ATUC, 2)));                           //K does not depend on scattering angle, so it is fine to divide it here.
//    }
//
//}
//
//double Objreflect2D::GetVal() {
//    Reset();
//    for (int idx = 0; idx < N; idx++) {
//        SingleResponse(idx, false);
//    }
//    return GroupResponse();
//}
//
//void Objreflect2D::Reset() {
//    for (list<Vector3cd>::iterator it = PSum_l.begin(); it != PSum_l.end(); it++) {
//        (*it) = Vector3cd::Zero();
//    }
//}
//
//double Objreflect2D::FTUCnsquareoversinal() {
//
//    list<Matrix3d>::iterator it1 = (FconstM_l).begin();
//    list<Vector3cd>::iterator it2 = (PSum_l).begin();
//    list<Vector3d>::iterator it3 = n_K_s_l.begin();
//    double result = 0.0;
//    for (int i = 0; i <= int(round(Paralength / 2) - 1); i++) {
//        double ksz = (*it3)(2);
//        Vector3cd FTUC;
//        FTUC(0) = (*it1)(0, 0) * (*it2)(0) + (*it1)(0, 1) * (*it2)(1) + (*it1)(0, 2) * (*it2)(2);
//        FTUC(1) = (*it1)(1, 0) * (*it2)(0) + (*it1)(1, 1) * (*it2)(1) + (*it1)(1, 2) * (*it2)(2);
//        FTUC(2) = (*it1)(2, 0) * (*it2)(0) + (*it1)(2, 1) * (*it2)(1) + (*it1)(2, 2) * (*it2)(2);
//
//
//        it1++;
//        it2++;
//        it3++;
//        result += (norm(FTUC(0)) + norm(FTUC(1)) + norm(FTUC(2))) / (ksz * ksz);                                 //In C++, norm is the square of magnitude.
//
//
//    }
//
//    return result / double(round(Paralength / 2));
//
//}
//
//
//
//ObjAbsbyfar::ObjAbsbyfar(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    Paralength = (parameters).size();
//    VectorXd FOMParameters = VectorXd::Zero(Paralength);
//    list<double>::iterator it = (parameters).begin();
//    for (int i = 0; i <= int(Paralength - 1); i++) {
//        FOMParameters(i) = (*it);
//        it++;
//    }
//    if (Paralength % 2 != 0) {
//        cout << "Objscattering2D ERROR: parameter must be times of 2." << endl;
//        throw 1;
//    }
//
//    Have_Penalty = HavePenalty_;
//    Have_Devx = false;
//    model = model_;
//    evomodel = evomodel_;
//
//    AProductCore* Core = (*model).get_Core();
//    d = (*Core).get_d();
//    N = (*Core).get_N();                   //Number of dipoles
//    P = (*model).get_P();
//    R = (*Core).get_R();
//    Vector3d n_E0 = (*model).get_nE0();
//    Vector3d n_K = (*model).get_nK();
//    E0 = (*model).get_E0();
//    double lam = (*Core).get_lam();
//    cout << "lam" << lam << endl;
//    K = (*Core).get_K();
//
//    double Lm = (*Core).get_Lm();
//    double Ln = (*Core).get_Ln();
//    ATUC = Lm * Ln;
//
//
//    for (int i = 0; i <= int(round(Paralength / 2) - 1); i++) {
//        Vector3d n_K_tmp;
//        n_K_tmp(0) = 2 * M_PI * FOMParameters(2 * i) / (Lm * K) + n_K(0);
//        n_K_tmp(1) = 2 * M_PI * FOMParameters(2 * i + 1) / (Ln * K) + n_K(1);
//        n_K_tmp(2) = sqrt(1 - pow(n_K_tmp(0), 2) - pow(n_K_tmp(1), 2));
//        n_K_s_l.push_back(n_K_tmp);
//        cout << n_K_tmp << endl;
//        Matrix3d FconstM;
//        double nkx = n_K_tmp(0);
//        double nky = n_K_tmp(1);
//        double nkz = n_K_tmp(2);
//        double K3 = pow(K, 3);
//        FconstM(0, 0) = K3 * (1 - nkx * nkx);
//        FconstM(0, 1) = -K3 * nkx * nky;
//        FconstM(0, 2) = -K3 * nkx * nkz;
//        FconstM(1, 1) = K3 * (1 - nky * nky);
//        FconstM(1, 2) = -K3 * nky * nkz;
//        FconstM(2, 2) = K3 * (1 - nkz * nkz);
//        FconstM(1, 0) = FconstM(0, 1);
//        FconstM(2, 0) = FconstM(0, 2);
//        FconstM(2, 1) = FconstM(1, 2);
//        FconstM_l.push_back(FconstM);
//
//        Vector3cd PSum_tmp;
//        PSum_tmp = Vector3cd::Zero();
//        PSum_l.push_back(PSum_tmp);
//
//    }
//
//
//
//
//
//}
//
//void ObjAbsbyfar::SingleResponse(int idx, bool deduction) {
//    //list<Matrix3d>::iterator it1 = (FconstM_l).begin();
//    list<Vector3d>::iterator it2 = (n_K_s_l).begin();
//    list<Vector3cd>::iterator it3 = (PSum_l).begin();
//    for (int i = 0; i <= int(round(Paralength / 2) - 1); i++) {
//        //Matrix3d FconstM = (*it1);
//        Vector3d n_K_s = (*it2);
//        double nkx = n_K_s(0);
//        double nky = n_K_s(1);
//        double nkz = n_K_s(2);
//        double phaseterm = -d * K * (nkx * ((*R)(3 * idx)) + nky * ((*R)(3 * idx + 1)) + nkz * ((*R)(3 * idx + 2)));  //From equation 17. Time term will be eliminated
//        complex<double> phase = cos(phaseterm) + sin(phaseterm) * 1i;
//        if (deduction == false) {
//            (*it3)(0) += (*P)(3 * idx) * phase;
//            (*it3)(1) += (*P)(3 * idx + 1) * phase;
//            (*it3)(2) += (*P)(3 * idx + 2) * phase;
//
//        }
//        else {
//            (*it3)(0) -= (*P)(3 * idx) * phase;
//            (*it3)(1) -= (*P)(3 * idx + 1) * phase;
//            (*it3)(2) -= (*P)(3 * idx + 2) * phase;
//        }
//        it2++;
//        it3++;
//    }
//}
//
//double ObjAbsbyfar::GroupResponse() {
//    if (Have_Penalty) {
//        return 1.0 - (this->FTUCnsquareoversinal() * (pow(2 * M_PI, 2)) / (pow(K, 4) * pow(E0, 2) * pow(ATUC, 2)) - (*evomodel).L1Norm());
//    }
//    else {
//        return 1.0 - this->FTUCnsquareoversinal() * (pow(2 * M_PI, 2)) / (pow(K, 4) * pow(E0, 2) * pow(ATUC, 2));                           //K does not depend on scattering angle, so it is fine to divide it here.
//    }
//
//}
//
//double ObjAbsbyfar::GetVal() {
//    Reset();
//    for (int idx = 0; idx < N; idx++) {
//        SingleResponse(idx, false);
//    }
//    return GroupResponse();
//}
//
//void ObjAbsbyfar::Reset() {
//    for (list<Vector3cd>::iterator it = PSum_l.begin(); it != PSum_l.end(); it++) {
//        (*it) = Vector3cd::Zero();
//    }
//}
//
//double ObjAbsbyfar::FTUCnsquareoversinal() {
//
//    list<Matrix3d>::iterator it1 = (FconstM_l).begin();
//    list<Vector3cd>::iterator it2 = (PSum_l).begin();
//    list<Vector3d>::iterator it3 = n_K_s_l.begin();
//    double result = 0.0;
//    for (int i = 0; i <= int(round(Paralength / 2) - 1); i++) {
//        double ksz = (*it3)(2);
//        Vector3cd FTUC;
//        FTUC(0) = (*it1)(0, 0) * (*it2)(0) + (*it1)(0, 1) * (*it2)(1) + (*it1)(0, 2) * (*it2)(2);
//        FTUC(1) = (*it1)(1, 0) * (*it2)(0) + (*it1)(1, 1) * (*it2)(1) + (*it1)(1, 2) * (*it2)(2);
//        FTUC(2) = (*it1)(2, 0) * (*it2)(0) + (*it1)(2, 1) * (*it2)(1) + (*it1)(2, 2) * (*it2)(2);
//
//
//        it1++;
//        it2++;
//        it3++;
//        result += (norm(FTUC(0)) + norm(FTUC(1)) + norm(FTUC(2))) / (ksz * ksz);                                 //In C++, norm is the square of magnitude.
//
//
//    }
//
//    return result / double(round(Paralength / 2));
//
//}
//
//
//
//
//
//
//ObjAbs::ObjAbs(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    Have_Penalty = HavePenalty_;
//    Have_Devx = true;
//    model = model_;
//    evomodel = evomodel_;
//    AProductCore* Core = (*model).get_Core();
//    N = (*Core).get_N();
//    Nx = (*Core).get_Nx();
//    Ny = (*Core).get_Ny();
//    Nz = (*Core).get_Nz();
//    d = (*Core).get_d();
//    al = (*model).get_al();
//    P = (*model).get_P();
//    E = VectorXcd::Zero(N * 3);
//    R = (*Core).get_R();
//    Cabs = 0.0;
//    E0 = (*model).get_E0();
//    double lam = (*Core).get_lam();
//    cout << "lam" << lam << endl;
//    K = (*Core).get_K();
//    K3 = pow(K, 3);
//}
//
//void ObjAbs::SingleResponse(int idx, bool deduction) {
//    complex<double> al_tmp = (*al)(3 * idx);
//    Vector3cd P_tmp((*P)(3 * idx), (*P)(3 * idx + 1), (*P)(3 * idx + 2));
//
//    if (deduction == false) {
//        complex<double> tmp = (al_tmp * P_tmp).dot(P_tmp);                    
//        Cabs += (tmp.imag() - (2.0 / 3.0) * K3 * (P_tmp.dot(P_tmp)).real());
//    }
//    else {
//        complex<double> tmp = (al_tmp * P_tmp).dot(P_tmp);                    //Eigen dot product is conjugate linear in the first variable, and linear in the second one.
//        Cabs -= (tmp.imag() - (2.0 / 3.0) * K3 * (P_tmp.dot(P_tmp)).real());  //The Eigen Vector norm turns out to be the square-root term, instead of the one with square like C++ complex number norm.
//    }
//    return;
//}
//
//double ObjAbs::GroupResponse() {
//
//    if (Have_Penalty) {
//        return log10(Cabs * (4 * M_PI * K / pow(E0, 2))) - (*evomodel).L1Norm();
//    }
//    else {
//        return log10(Cabs * (4 * M_PI * K / pow(E0, 2)));
//    }
//
//}
//
//double ObjAbs::GetVal() {
//    Reset();
//    for (int idx = 0; idx < N; idx++) {
//        SingleResponse(idx, false);
//    }
//    return GroupResponse();
//}
//
//void ObjAbs::Reset() {
//    Cabs = 0.0;
//}
//
//
//
//ObjAbsPartial::ObjAbsPartial(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    Have_Penalty = HavePenalty_;
//    Have_Devx = true;
//    model = model_;
//    evomodel = evomodel_;
//    AProductCore* Core = (*model).get_Core();
//    N = (*Core).get_N();
//    Nx = (*Core).get_Nx();
//    Ny = (*Core).get_Ny();
//    Nz = (*Core).get_Nz();
//    d = (*Core).get_d();
//    al = (*model).get_al();
//    P = (*model).get_P();
//    E = VectorXcd::Zero(N * 3);
//    R = (*Core).get_R();
//    Cabs = 0.0;
//    E0 = (*model).get_E0();
//    double lam = (*Core).get_lam();
//    cout << "lam" << lam << endl;
//    K = (*Core).get_K();
//    K3 = pow(K, 3);
//    SpacePara* spaceparatmp = (*((*Core).get_CStr())).get_spacepara();
//    vector<int>* ParaDividePos = (*spaceparatmp).get_ParaDividePos();
//    int NPara = (*((*spaceparatmp).get_Para())).size();
//    int numPara = (*ParaDividePos).size();
//    list<double>::iterator it = parameters.begin();
//    vector<int> startpos;
//    vector<int> endpos;
//    while (it != parameters.end()) {
//        int geometrylabel = int(round(*it));
//        if (geometrylabel >= numPara || geometrylabel < 0) {
//            cout << "ERROR:ObjAbsPartial--geometrylabel out of range" << endl;
//            throw 1;
//        }
//        startpos.push_back((*ParaDividePos)[geometrylabel]);
//        cout << "startpos.back(): " << startpos.back() << endl;
//        if (geometrylabel == numPara - 1) {
//            endpos.push_back(NPara - 1);
//        }
//        else {
//            endpos.push_back((*ParaDividePos)[geometrylabel + 1] - 1);
//        }
//        cout << "endpos.back(): " << endpos.back() << endl;
//        it++;
//    }
//    vector<vector<int>>* Paratogeometry = (*spaceparatmp).get_Paratogeometry();
//    
//    for (int i = 0; i < startpos.size(); i++) {
//        for (int j = startpos[i]; j <= endpos[i]; j++) {
//            for (int k = 0; k < (*Paratogeometry)[j].size(); k++) {
//                integralpos.insert((*Paratogeometry)[j][k]);
//                //cout << (*Paratogeometry)[j][k] << endl;
//            }
//        }
//        
//    }
//
//
//}
//
//void ObjAbsPartial::SingleResponse(int idx, bool deduction) {
//    if (integralpos.count(idx)) {
//        complex<double> al_tmp = (*al)(3 * idx);
//        Vector3cd P_tmp((*P)(3 * idx), (*P)(3 * idx + 1), (*P)(3 * idx + 2));
//
//        if (deduction == false) {
//            complex<double> tmp = (al_tmp * P_tmp).dot(P_tmp);
//            Cabs += (tmp.imag() - (2.0 / 3.0) * K3 * (P_tmp.dot(P_tmp)).real());
//        }
//        else {
//            complex<double> tmp = (al_tmp * P_tmp).dot(P_tmp);                    //Eigen dot product is conjugate linear in the first variable, and linear in the second one.
//            Cabs -= (tmp.imag() - (2.0 / 3.0) * K3 * (P_tmp.dot(P_tmp)).real());  //The Eigen Vector norm turns out to be the square-root term, instead of the one with square like C++ complex number norm.
//        }
//        return;
//    }
//
//    
//}
//
//double ObjAbsPartial::GroupResponse() {
//
//    if (Have_Penalty) {
//        return log10(Cabs * (4 * M_PI * K / pow(E0, 2))) - (*evomodel).L1Norm();
//    }
//    else {
//        return log10(Cabs * (4 * M_PI * K / pow(E0, 2)));
//    }
//
//}
//
//double ObjAbsPartial::GetVal() {
//    Reset();
//    for (int idx = 0; idx < N; idx++) {
//        SingleResponse(idx, false);
//    }
//    return GroupResponse();
//}
//
//void ObjAbsPartial::Reset() {
//    Cabs = 0.0;
//}
///*
//Objscattering2D::Objscattering2D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    VectorXd PointEParameters = VectorXd::Zero((parameters).size());
//    list<double>::iterator it = (parameters).begin();
//    for (int i = 0; i <= int((parameters).size() - 1); i++) {
//        PointEParameters(i) = (*it);
//        it++;
//    }
//    Have_Penalty = HavePenalty_;
//    x = PointEParameters(0);
//    y = PointEParameters(1);
//    z = PointEParameters(2);
//    Have_Devx = false;
//    model = model_;
//    evomodel = evomodel_;
//    AProductCore* Core = (*model).get_Core();
//    d = (*Core).get_d();
//    N = (*Core).get_N();
//    P = (*model).get_P();
//    R = (*Core).get_R();
//    Vector3d n_E0 = (*model).get_nE0();
//    Vector3d n_K = (*model).get_nK();
//    double E0 = (*model).get_E0();
//    double lam = (*Core).get_lam();
//    cout << "lam" << lam << endl;
//    double K = 2 * M_PI / lam;
//    E_sum = Vector3cd::Zero();
//    E_ext = Vector3cd::Zero();
//    E_ext(0) = E0 * n_E0(0) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
//    E_ext(1) = E0 * n_E0(1) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
//    E_ext(2) = E0 * n_E0(2) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
//    // cout << R(3*5444+2) << "this" << endl;
//}
//
//void Objscattering2D::SingleResponse(int idx, bool deduction) {
//}
//
//double Objscattering2D::GroupResponse() {
//}
//
//double Objscattering2D::GetVal() {
//    Reset();
//    for (int idx = 0; idx < N; idx++) {
//        SingleResponse(idx, false);
//        // cout << E_sum(0) << endl;
//    }
//    return GroupResponse();
//}
//
//void Objscattering2D::Reset() {
//    E_sum(0) = E_ext(0);
//    E_sum(1) = E_ext(1);
//    E_sum(2) = E_ext(2);
//}
//
//double Objscattering2D::FTUCnsquare() {
//    double 
//}
//*/
//
//
//
//ObjAbsPartialzslice::ObjAbsPartialzslice(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    Have_Penalty = HavePenalty_;
//    Have_Devx = true;
//    model = model_;
//    evomodel = evomodel_;
//    AProductCore* Core = (*model).get_Core();
//    N = (*Core).get_N();
//    Nx = (*Core).get_Nx();
//    Ny = (*Core).get_Ny();
//    Nz = (*Core).get_Nz();
//    d = (*Core).get_d();
//    al = (*model).get_al();
//    P = (*model).get_P();
//    E = VectorXcd::Zero(N * 3);
//    R = (*Core).get_R();
//    Cabs = 0.0;
//    E0 = (*model).get_E0();
//    double lam = (*Core).get_lam();
//    cout << "lam" << lam << endl;
//    K = (*Core).get_K();
//    K3 = pow(K, 3);
//    SpacePara* spaceparatmp = (*((*Core).get_CStr())).get_spacepara();
//    vector<int>* ParaDividePos = (*spaceparatmp).get_ParaDividePos();
//    int NPara = (*((*spaceparatmp).get_Para())).size();
//    int numPara = (*ParaDividePos).size();
//    list<double>::iterator it = parameters.begin();
//    vector<int> startpos;
//    vector<int> endpos;
//    while (it != parameters.end()) {
//        int geometrylabel = int(round(*it));
//        if (geometrylabel >= numPara || geometrylabel < 0) {
//            cout << "ERROR:ObjAbsPartial--geometrylabel out of range" << endl;
//            throw 1;
//        }
//        startpos.push_back((*ParaDividePos)[geometrylabel]);
//        cout << "startpos.back(): " << startpos.back() << endl;
//        if (geometrylabel == numPara - 1) {
//            endpos.push_back(NPara - 1);
//        }
//        else {
//            endpos.push_back((*ParaDividePos)[geometrylabel + 1] - 1);
//        }
//        cout << "endpos.back(): " << endpos.back() << endl;
//        it++;
//    }
//    vector<vector<int>>* Paratogeometry = (*spaceparatmp).get_Paratogeometry();
//
//    for (int i = 0; i < startpos.size(); i++) {
//        for (int j = startpos[i]; j <= endpos[i]; j++) {
//            for (int k = 0; k < (*Paratogeometry)[j].size(); k++) {
//                integralpos.insert((*Paratogeometry)[j][k]);
//                //cout << (*Paratogeometry)[j][k] << endl;
//            }
//        }
//
//    }
//
//
//
//}
//
//void ObjAbsPartialzslice::SingleResponse(int idx, bool deduction) {
//    if (integralpos.count(idx)) {
//        complex<double> al_tmp = (*al)(3 * idx);
//        Vector3cd P_tmp((*P)(3 * idx), (*P)(3 * idx + 1), (*P)(3 * idx + 2));
//
//        if (zslices.count((*R)(3 * idx + 2))) {
//            if (deduction == false) {
//                complex<double> tmp = (al_tmp * P_tmp).dot(P_tmp);
//                Cabs += (tmp.imag() - (2.0 / 3.0) * K3 * (P_tmp.dot(P_tmp)).real());
//            }
//            else {
//                complex<double> tmp = (al_tmp * P_tmp).dot(P_tmp);                    //Eigen dot product is conjugate linear in the first variable, and linear in the second one.
//                Cabs -= (tmp.imag() - (2.0 / 3.0) * K3 * (P_tmp.dot(P_tmp)).real());  //The Eigen Vector norm turns out to be the square-root term, instead of the one with square like C++ complex number norm.
//            }
//        }
//
//        
//        return;
//    }
//
//
//}
//
//double ObjAbsPartialzslice::GroupResponse() {
//
//    if (Have_Penalty) {
//        return log10(Cabs * (4 * M_PI * K / pow(E0, 2))) - (*evomodel).L1Norm();
//    }
//    else {
//        return log10(Cabs * (4 * M_PI * K / pow(E0, 2)));
//    }
//
//}
//
//double ObjAbsPartialzslice::GetVal() {
//    Reset();
//    for (int idx = 0; idx < N; idx++) {
//        SingleResponse(idx, false);
//    }
//    return GroupResponse();
//}
//
//void ObjAbsPartialzslice::Reset() {
//    Cabs = 0.0;
//}
//
//
//
//
//ObjIntegrateEPartial::ObjIntegrateEPartial(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
//    Have_Penalty = HavePenalty_;
//    Have_Devx = true;
//    model = model_;
//    evomodel = evomodel_;
//    AProductCore* Core = (*model).get_Core();
//    N = (*Core).get_N();
//    Nx = (*Core).get_Nx();
//    Ny = (*Core).get_Ny();
//    Nz = (*Core).get_Nz();
//    d = (*Core).get_d();
//    al = (*model).get_al();
//    P = (*model).get_P();
//    E = VectorXcd::Zero(N * 3);
//    R = (*Core).get_R();
//    Total = 0.0;
//    E0 = (*model).get_E0();
//    double lam = (*Core).get_lam();
//    cout << "lam" << lam << endl;
//    K = (*Core).get_K();
//    K3 = pow(K, 3);
//    SpacePara* spaceparatmp = (*((*Core).get_CStr())).get_spacepara();
//    vector<int>* ParaDividePos = (*spaceparatmp).get_ParaDividePos();
//    int NPara = (*((*spaceparatmp).get_Para())).size();
//    int numPara = (*ParaDividePos).size();
//    list<double>::iterator it = parameters.begin();
//    vector<int> startpos;
//    vector<int> endpos;
//    while (it != parameters.end()) {
//        int geometrylabel = int(round(*it));
//        if (geometrylabel >= numPara || geometrylabel < 0) {
//            cout << "ERROR:ObjAbsPartial--geometrylabel out of range" << endl;
//            throw 1;
//        }
//        startpos.push_back((*ParaDividePos)[geometrylabel]);
//        if (geometrylabel == numPara - 1) {
//            endpos.push_back(NPara - 1);
//        }
//        else {
//            endpos.push_back((*ParaDividePos)[geometrylabel + 1] - 1);
//        }
//        it++;
//    }
//    vector<vector<int>>* Paratogeometry = (*spaceparatmp).get_Paratogeometry();
//
//    for (int i = 0; i < startpos.size(); i++) {
//        for (int j = startpos[i]; j <= endpos[i]; j++) {
//            for (int k = 0; k < (*Paratogeometry)[j].size(); k++) {
//                integralpos.insert((*Paratogeometry)[j][k]);
//                //cout << (*Paratogeometry)[j][k] << endl;
//            }
//        }
//
//    }
//
//
//}
//
//void ObjIntegrateEPartial::SingleResponse(int idx, bool deduction) {
//    if (integralpos.count(idx)) {
//        complex<double> al_tmp = (*al)(3 * idx);
//        Vector3cd P_tmp((*P)(3 * idx), (*P)(3 * idx + 1), (*P)(3 * idx + 2));
//        Vector3cd E_tmp = al_tmp * P_tmp;
//
//        if (deduction == false) {
//            double E_tmp_square = 0.0;
//            for (int i = 0; i < 3; i++) {
//                E_tmp_square += norm(E_tmp(i));
//            }
//            Total += pow(E_tmp_square, 2);
//        }
//        else {
//            double E_tmp_square = 0.0;
//            for (int i = 0; i < 3; i++) {
//                E_tmp_square += norm(E_tmp(i));
//            }
//            Total -= pow(E_tmp_square, 2);
//        }
//        return;
//    }
//
//
//}
//
//double ObjIntegrateEPartial::GroupResponse() {
//
//    if (Have_Penalty) {
//        return log10(Total) - (*evomodel).L1Norm();
//    }
//    else {
//        return log10(Total);
//    }
//
//}
//
//double ObjIntegrateEPartial::GetVal() {
//    Reset();
//    for (int idx = 0; idx < N; idx++) {
//        SingleResponse(idx, false);
//    }
//    return GroupResponse();
//}
//
//void ObjIntegrateEPartial::Reset() {
//    Total = 0.0;
//}





