/*
    for some common define
    and common method
*/
#ifndef DEFINE_H
#define DEFINE_H

#include <vector>
#include <string>
#include <sys/time.h>
#include "cudd.h"
#include "meddly.h"

#define DEBUG 0
#define GENERAL 1
#define BSTORE 1

#define DEFAULT_TEST_DIR "./data/testData/"
#define DEFAULT_TEST_DAG_FILE "testDag"
#define DEFAULT_TEST_CHK_FILE "testCHK"
#define DEFAULT_TEST_DIS_FILE "testDIS"

#define DEFAULT_EXPERIMENT_DIR "./data/experimentData/"
#define KNAP_EXPERIMENT_DAG_FILE "expKnapDag"
#define KNAP_EXPERIMENT_CHK1_3_1_FILE "expKnapCHK1"
#define KNAP_EXPERIMENT_CHK2_5_5_FILE "expKnapCHK2"
#define KNAP_EXPERIMENT_CHK3_80_20_FILE "expKnapCHK3"
#define KNAP_EXPERIMENT_DIS1_50_1_FILE "expKnapDIS1"
#define KNAP_EXPERIMENT_DIS2_60_70_FILE "expKnapDIS2"
#define KNAP_EXPERIMENT_DIS3_100_80_FILE "expKnapDIS3"

#define CIT_EXPERIMENT_DAG_FILE "expCitDag"
#define CIT_EXPERIMENT_CHK1_3_1_FILE "expCitCHK1"
#define CIT_EXPERIMENT_CHK2_5_5_FILE "expCitCHK2"
#define CIT_EXPERIMENT_CHK3_8_20_FILE "expCitCHK3"
#define CIT_EXPERIMENT_DIS1_50_1_FILE "expCitDIS1"
#define CIT_EXPERIMENT_DIS2_60_70_FILE "expCitDIS2"
#define CIT_EXPERIMENT_DIS3_100_80_FILE "expCitDIS3"

#define ED_EXPERIMENT_DAG_FILE "expEdDag"
#define ED_EXPERIMENT_CHK1_3_1_FILE "expEdCHK1"
#define ED_EXPERIMENT_CHK2_5_5_FILE "expEdCHK2"
#define ED_EXPERIMENT_CHK3_8_20_FILE "expEdCHK3"
#define ED_EXPERIMENT_DIS1_50_1_FILE "expEdDIS1"
#define ED_EXPERIMENT_DIS2_60_70_FILE "expEdDIS2"
#define ED_EXPERIMENT_DIS3_100_80_FILE "expEdDIS3"

#define POS_EXPERIMENT_DAG_FILE "expPosDag"
#define POS_EXPERIMENT_CHK1_3_1_FILE "expPosCHK1"
#define POS_EXPERIMENT_CHK2_5_5_FILE "expPosCHK2"
#define POS_EXPERIMENT_CHK3_8_20_FILE "expPosCHK3"
#define POS_EXPERIMENT_DIS1_50_1_FILE "expPosDIS1"
#define POS_EXPERIMENT_DIS2_60_70_FILE "expPosDIS2"
#define POS_EXPERIMENT_DIS3_100_80_FILE "expPosDIS3"

#define DCKP_RANDOM_EXPERIMENT_DAG_FILE "expDckpRandom"
#define DCKP_RANDOM_EXPERIMENT_DIS1_1_FILE "expDckpRandomDis1"
#define DCKP_RANDOM_EXPERIMENT_DIS2_2_FILE "expDckpRandomDis2"
#define DCKP_RANDOM_EXPERIMENT_DIS3_5_FILE "expDckpRandomDis3"

#define DCKP_CORREL_EXPERIMENT_DAG_FILE "expDckpCorrel"
#define DCKP_CORREL_EXPERIMENT_DIS1_1_FILE "expDckpCorrelDis1"
#define DCKP_CORREL_EXPERIMENT_DIS2_2_FILE "expDckpCorrelDis2"
#define DCKP_CORREL_EXPERIMENT_DIS3_5_FILE "expDckpCorrelDis3"

#define DEFAULT_PICTURE_DIR "./picture/"
#define DEFAULT_PICTURE_NAME "testDDPicture"
#define DEFAULT_MDD_PICTURE_NAME "testMDDPicture"

#define DEFAULT_ORIGNAL_DATA_DIR "./orignalData/"
#define DEFAULT_CIT_DAG_ORIGNAL_DATA_FILE "originalScCitDAG"


using namespace MEDDLY;
using namespace std;

enum dagProblem {
    knap,
    other,
    unKnow
};

struct sortBddNode {  
    bool operator()(DdNode* a, DdNode* b) const {
        return Cudd_NodeReadIndex(a) > Cudd_NodeReadIndex(b);
    }  
};

class commonTool {
    public:

    static inline vector<int> sepnum(string line);

    static inline bool pairSecondAsendSort(pair<int, int> p1, pair<int, int> p2);

    /*
        find the edge/group position in var2eg
        if not exist, return -1
    */
    static inline int findExistEG(int eg, vector<int>* var2gro);

    static inline bool comparePairs(const std::pair<int, int>& a, const std::pair<int, int>& b);

    static inline void compareTtable(vector<unordered_map<DdNode*, int>>* bddT, vector<unordered_map<node_handle, int>>* mddT, int numOfDiffT = -1);

    static inline void compareEgPath(vector<int>* bcsEgPath, vector<int>* mcsEgpath);

};

/*
    inline method
*/

inline vector<int> commonTool::sepnum(string line) {
    vector<int> num;
	string nstring;
	for (char c : line) {
		if (c == ' ') {
			num.push_back(atoi(nstring.c_str()));
			nstring = "";
			continue;
		}
		nstring.push_back(c);
	}
	num.push_back(atoi(nstring.c_str()));
	return num;
}

inline bool commonTool::pairSecondAsendSort(pair<int, int> p1, pair<int, int> p2) {
    return p1.second < p2.second;
}

inline int commonTool::findExistEG(int eg, vector<int>* var2eg) {
    eg = abs(eg);
    for(int i = 0; i < var2eg->size(); i++) {
        if((*var2eg)[i] == eg)
            return i;
    }
    return -1;
}

inline bool commonTool::comparePairs(const std::pair<int, int>& a, const std::pair<int, int>& b) {
    return a.first > b.first;
}

inline void commonTool::compareTtable(vector<unordered_map<DdNode*, int>>* bddT, vector<unordered_map<node_handle, int>>* mddT, int numOfDiffT) {
    int n = 0;
    cout << "[COMPARE T TABLES]:\n";
    if(bddT->size() != mddT->size()) {
        cerr << "[ERROR]T tables size is not same.";
        exit(EXIT_FAILURE);
    }
    bool isDiff = false;
    for(int i = 0; i < bddT->size(); i++) {
        if(bddT->at(i).size() == mddT->at(i).size()) 
            continue;
        isDiff = true;
        n++;
        cout << "  T[" << i << "]" << "is different.\n";
        cout << "\tbdd T size is " << bddT->at(i).size() << ".\n";
        for(auto tp : bddT->at(i)){
			cout << "\t\tT["<< i << "][" << tp.first << "] = " << tp.second << "\n";
		}
        cout << "\tmdd T size is " << mddT->at(i).size() << ".\n";
        for(auto tp : mddT->at(i)){
			cout << "\t\tT["<< i << "][" << tp.first << "] = " << tp.second << "\n";
		}
        // decide whether only print first difference.
        if(numOfDiffT != -1) {
            if(numOfDiffT == n)
                break;
        }
    }
    if(!isDiff)
        cout << "\tT table is same." << endl;
}

inline void commonTool::compareEgPath(vector<int>* p1, vector<int>* p2) {
    int size = min(p1->size(), p2->size());
    bool isDiff = true;
    cout << "[COMPARE EG PATHS]:\n";
    if(size == 0) {
        cerr << "\t[ERROR]Using the method under debug mode." << endl;
        return;
    }
    for(int i = 0; i < size; i++) {
        if(p1->at(i) == p2->at(i)) continue;
        isDiff = false;
        cout << "\tDiff from " << i << " th path:\n";
        cout << "\t\t bcs path eg = " << p1->at(i) << " mcs path eg = " << p2->at(i) << "\n";
    }
    if(isDiff)
        cout << "\tPaths are same." << endl;
}

#endif