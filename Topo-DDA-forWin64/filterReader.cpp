#include <iostream>

#include "filterReader.h"
#include "Tools.h"

filterReader::filterReader(INIReader reader) {
	string tmp = reader.Get("Filter Option", "filter", "UNKNOWN");
	if (tmp != "True") {
		filter = false;
	}
	else {
		filter = true;
        /*cout << reader.Get("Filter Option", "filterInfo", "UNKNOWN") << endl;*/
        filterList = readFilterList(reader.Get("Filter Option", "filterInfo", "UNKNOWN"));
        betaType = reader.Get("Filter Option", "betaType", "UNKNOWN");
        betaMin = reader.GetReal("Filter Option", "betaMin", 0.0);
        betaMax = reader.GetReal("Filter Option", "betaMax", 50.0);
        ita = reader.GetReal("Filter Option", "ita", 0.5);
	}
    

    return;
}

vector<filterinfo> filterReader::readFilterList(string input) {
    vector<string> split1 = splitInputStr(input, "/");
    vector<filterinfo> result;
    for (int i = 0; i < split1.size(); i++) {
        vector<string> split2 = splitInputStr(split1[i], ",");
        if (split2.size() != 2) {
            cout << "ERROR:vector<filterinfo> getFilterList(string input)--Input filter parameter size wrong" << endl;
            throw 1;
        }
        result.push_back(filterinfo{ stoi(split2[0]), stod(split2[1]) });
        /*cout << stoi(split2[0]) << endl;
        cout << stod(split2[1]) << endl;*/
    }
    return result;
}

bool filterReader::getFilter() {
    return filter;
}

vector<filterinfo> filterReader::getFilterList() {
    return filterList;
}

string filterReader::getBetaType() {
    return betaType;
}

double filterReader::getBetaMin() {
    return betaMin;
}

double filterReader::getBetaMax() {
    return betaMax;
}

double filterReader::getIta() {
    return ita;
}