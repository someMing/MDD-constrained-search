#ifndef TESTDATA_H
#define TESTDATA_H

#include <string>
#include <stack>
#include <limits.h>
#include <algorithm>
#include <fstream>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include "dag.h"

using namespace std;

// only use for generate experiment data
struct sampleNode {
    int id;
    // target id, weight
    vector<pair<sampleNode*, int>> neighbor;
};

class TestData {
private:
    // safe model deny any generate to avoid cover
    bool safeModel;
    string dirName;
    string dagFileName;
    string chkFileName;
    string disFileName;
    dagProblem problemType;
    vector<pair<int,int>> bestGroup;
    int dagGroNum;

public:
    TestData(string dirName = DEFAULT_TEST_DIR, 
        string dagFileName = DEFAULT_TEST_DAG_FILE, 
        string chkFileName = DEFAULT_TEST_CHK_FILE,
        string disFileName = DEFAULT_TEST_DIS_FILE);

    void setChkName(string newName);

    void setDisName(string newName);

    void setDirName(string newName);

    void setDagName(string newName);

    string getDagFileName();

    string getChkFileName();

    string getDisFileName();

    void allowGenerate();

    template <typename dagType>
    DAG<dagType>* readDag(bool isDetial = false);

    template <typename dagType>
    DAG<dagType>* readPosDag(bool isDetial = false);

    static vector<sampleNode*> topologyOrder(vector<sampleNode*>& orignalList, vector<int> in);

    /*
        get a shortest/longest path of the DAG
        NOTICE: if need detial print, please make sure vertex type can 
        be print in cout
    */
    template <typename dagType>
    static int DP_path(DAG<dagType>* dag, bool isShort = true, bool isDetail = false);
};

/*
    template method
*/
template <typename dagType>
int TestData::DP_path(DAG<dagType>* dag, bool isShort, bool isDetail) {
    int numVex = dag->getVertexNum();
    int limit = isShort ? INT_MAX : INT_MIN;
	vector<int> T(numVex, limit);
	vector<Edge*> B(numVex);
    int step = 0;

    // DP algorithm
	T[0] = 0;
	for (int s = 0; s < numVex; s++) {
        vector<Edge>* vexEdge = dag->getEdgeList(s);
        Edge* curEdge;
        for(int i = 0; i < vexEdge->size(); i++) {
            step++;
            curEdge = &((*vexEdge)[i]);
            int vt = curEdge->target();
            int wt = curEdge->weight();
            if((isShort && T[vt] > T[s] + wt) || (!isShort && T[vt] < T[s] + wt)) {
                // replace
                T[vt] = T[s] + wt;
                B[vt] = curEdge;
            }
        }
	}
    string pathType = isShort ? "Short" : "Long";
    cout << "[DP " << pathType << "est" << " Path's Length is:" 
        << T[numVex-1] << "]\n";
    cout << "\tstep number is:" << step << endl;
    if(!isDetail) return T[numVex-1];

    // detail print
    stack<int> st;
	int v = numVex - 1;
	while (v > 0) {
        st.push(B[v]->target());
		v = B[v]->source();
	}
    cout << "\t" << 0;
    while (!st.empty()){
		cout << " --> " << st.top();
		st.pop();
	}
    cout << endl;
    return T[numVex-1];
}

template <typename dagType>
DAG<dagType>* TestData::readDag(bool isDetail) {
    problemType = knap;
    string fileName = dirName + dagFileName;
    bestGroup.clear();
    ifstream in(fileName, ios::in);
	string line;
    vector<int> intV;
    cout << "[READ KNAP FILE:" << fileName << "]\n";
    while (getline(in, line)) {
		if (line == "#numVex") {
            if(isDetail) 
                cout << "\tis reading number of vex ...\n";
            break;
        }
	}
    getline(in, line);
	int numVex = commonTool::sepnum(line)[0];
	cout << "\tvertex num is:" << numVex << "\n";

    while (getline(in, line)) {
		if (line == "#bestGroup") {
            if(isDetail) 
                cout << "\tis reading best group ...\n";
            break;
        }
	}
    while (getline(in, line)){
        if(line == "") continue;
		intV = commonTool::sepnum(line);
		bestGroup.push_back({intV[0],intV[1]});
	}
    dagGroNum = bestGroup.size();
    cout << "\tthe number of edge group is:" << dagGroNum << endl;
    DAG<string>* dag = new DAG<string>(numVex, dagGroNum);

    in.clear();
    in.seekg(0, ios::beg);
	getline(in, line);
    if(line == "#edge") {
        if(isDetail)
            cout << "\tis reading edge ...\n";
    }
	// reading edge
    int nOfe = 0;
	while (getline(in, line)) {
		if (line == "#numVex") {
			getline(in, line);
			getline(in, line);
			break;
		}
		// remain edge
		intV = commonTool::sepnum(line);
		dag->setEdge(intV[0], intV[1], intV[2], intV[3]);
		nOfe++;
	}
	cout << "\tedge num is:" << nOfe << "\n";
	// reading item
    int nOfit = 0;
	int befItem = -1;
	if (isDetail && line == "#itemSep") cout << "\tis reading item ...\n";
	while (getline(in, line)) {   
		// reading best
		if (line == "#bestGroup") {
			if(isDetail)
                cout << "\tread finished\n";
			break;
		}
		// reading res
		if (line == "#maxValue") {
            cout << "\tthe number of node group is:" << nOfit << "\n";
            if(isDetail)
                cout << "\tis reading max value ...\n";
			getline(in, line);
			cout << "\tthe max value result by DP is:" << commonTool::sepnum(line)[0] << "\n";
			continue;
		}
		// set node
		int aftItem = commonTool::sepnum(line)[0];
		for (int i = befItem + 1; i <= aftItem; i++) {
            // for knap data set
            // final item means result, not a true item
            // for example n = 200, item[0,199] means real item and item200 means result layer
            // so node group number = item group num + 2
            // if no need node gro, just make this [0, numVex-1] and the node gro will be [1, numVex];
            dag->setVertex("item" + to_string(nOfit), i);
            dag->setNodeGro(i, nOfit+1);
		}
		befItem = aftItem;
		nOfit++;
	}
	in.close();
	return dag;
}

template <typename dagType>
DAG<dagType>* TestData::readPosDag(bool isDetail) {
    problemType = knap;
    string fileName = dirName + dagFileName;
    bestGroup.clear();
    ifstream in(fileName, ios::in);
	string line;
    vector<int> intV;
    unordered_set<int> nodeGro;
    cout << "[READ POS FILE:" << fileName << "]\n";
    getline(in, line);
    if (line.find("#numVex") != string::npos) {
        if(isDetail) cout << "\tis reading number of vex ...\n";
    }
	getline(in, line);
	int numVex = commonTool::sepnum(line)[0];
	cout << "\tvertex num is:" << numVex << "\n";

    getline(in, line);
    if (line.find("#numEdgeGroup") != string::npos) {
        if(isDetail) cout << "\tis reading number of edge groups ...\n";
    }
	getline(in, line);
	dagGroNum= commonTool::sepnum(line)[0];
	cout << "\tthe number of edge group is:" << dagGroNum << endl;
    DAG<string>* dag = new DAG<string>(numVex, dagGroNum);

    getline(in, line);
    if (line.find("#node") != string::npos) {
        if(isDetail) cout << "\tis reading info of node ...\n";
    }
    while(getline(in, line)) {
        if(line.find("#edge") != string::npos) {
            if(isDetail) cout << "\tis reading edge ...\n";
            break;
        }
        intV = commonTool::sepnum(line);
		dag->setNodeGro(intV[0], intV[1]);
        nodeGro.emplace(intV[1]);
    }
    cout << "\tthe number of node group is:" << nodeGro.size() << "\n";

    int nOfe = 0;
    while(getline(in, line)) {
        if(line.find("#bestGroup") != string::npos) {
            if(isDetail) cout << "\tis reading best groups ...\n";
            break;
        }
        intV = commonTool::sepnum(line);
		dag->setEdge(intV[0], intV[1], intV[2], intV[3]);
		nOfe++;
    }
    cout << "\tedge num is:" << nOfe << "\n";

    while(getline(in, line)) {
        if(line == "") {
            cout << "\tread finished\n";
            break;
        }
        intV = commonTool::sepnum(line);
		bestGroup.push_back({intV[0],intV[1]});
    }
    
    in.close();
    return dag;
}


#endif