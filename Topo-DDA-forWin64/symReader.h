#ifndef TOPO_SYMREADER_H_
#define TOPO_SYMREADER_H_
#include <string>
#include <vector>

#include "INIReader.h"

using namespace std;

class symReader {
private:
	bool hasSym;
	string symmetry;
	vector<double> symAxis;
public:
	symReader(INIReader reader);
	vector<double> readSymAxis(string input);
	bool getHasSym();
	string getSymmetry();
	vector<double> getSymAxis();
};

#endif
