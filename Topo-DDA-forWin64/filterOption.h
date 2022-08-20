#ifndef TOPO_FILTEROPT_H_
#define TOPO_FILTEROPT_H_
#include <string>
#include <vector>
using namespace std;

struct filterinfo {
    int iteration;
    double rfilter;
};

class FilterOption {
private:
    double beta;
    double ita;
    string beta_type;
    double beta_min;
    double beta_max;
    double rfilter;
    vector<filterinfo> rfilterlist;
    bool fixit;
    int MAX_ITERATION_FIXED;

public:
    FilterOption(double beta_min_, double beta_max_, double ita_, string beta_type_, vector<filterinfo> rfilterlist_, bool fixit_ = false, int MAX_ITERATION_FIXED_ = 100);
    double get_beta();
    double get_ita();
    double get_rfilter();
    void update_beta(const int iteration, const int Max_iteration);
    double SmoothDensity(double input);
    bool filterchange(int iteration);
};
#endif
