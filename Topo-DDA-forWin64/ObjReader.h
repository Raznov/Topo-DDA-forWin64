#ifndef TOPO_OBJREADER_H_
#define TOPO_OBJREADER_H_
#include <string>
#include <vector>

#include "INIReader.h"

using namespace std;

class ObjReader {
private:
	string objName;
	vector<double> objPara;
public:
	ObjReader(INIReader reader);
	vector<double> ReadObjPara(string input);
	string GetObjName();
	vector<double> GetObjPara();
};

#endif