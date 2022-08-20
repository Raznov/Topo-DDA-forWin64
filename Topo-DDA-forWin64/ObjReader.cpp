#include "ObjReader.h"
#include "Tools.h"
ObjReader::ObjReader(INIReader reader) {
	objName = reader.Get("Obj Option", "objName", "UNKNOWN");
	objPara = ReadObjPara(reader.Get("Obj Option", "objPara", "UNKNOWN"));
	return;
}

vector<double> ObjReader::ReadObjPara(string input) {
	vector<string> split1 = splitInputStr(input, "/");
	vector<double> result;
	for (int i = 0; i < split1.size(); i++) {
		result.push_back(stod(split1[i]));
	}
	return result;
}

string ObjReader::GetObjName() {
	return objName;
}

vector<double> ObjReader::GetObjPara() {
	return objPara;
}
