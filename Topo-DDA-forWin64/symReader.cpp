#include "symReader.h"
#include "Tools.h"

symReader::symReader(INIReader reader) {
	string tmp = reader.Get("Symmetry Option", "hasSym", "UNKNOWN");
	if (tmp != "True") {
		hasSym = false;
		symmetry = "None";
	}
	else {
		hasSym = true;
		symmetry = reader.Get("Symmetry Option", "symmetry", "UNKNOWN");
		symAxis = readSymAxis(reader.Get("Symmetry Option", "symAxis", "UNKNOWN"));
	}


	return;
}

vector<double> symReader::readSymAxis(string input) {
	vector<string> split1 = splitInputStr(input, "/");
	vector<double> result;
	for (int i = 0; i < split1.size(); i++) {
		result.push_back(stod(split1[i]));
	}
	return result;
}

bool symReader::getHasSym() {
	return hasSym;
}

string symReader::getSymmetry() {
	return symmetry;
}

vector<double> symReader::getSymAxis() {
	return symAxis;
}
