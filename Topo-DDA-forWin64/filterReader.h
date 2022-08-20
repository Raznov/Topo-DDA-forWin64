#ifndef TOPO_FILTERREAD_H_
#define TOPO_FILTERREAD_H_
#include <vector>

#include "INIReader.h"

#include "filterOption.h"

using namespace std;

class filterReader {
private:
    bool filter;
    vector<filterinfo> filterList;
    string betaType;
    double betaMin;
    double betaMax;
    double ita;
public:
    filterReader(INIReader reader);
    vector<filterinfo> readFilterList(string input);
    bool getFilter();
    vector<filterinfo> getFilterList();
    string getBetaType();
    double getBetaMin();
    double getBetaMax();
    double getIta();
};
#endif
