#ifndef BCSALGORITHM_H
#define BCSALGORITHM_H
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include "dag.h"
#include "cudd.h"

template <typename dagType>
class BcsAlgorithm {
private:
    bool isShort;
    DAG<dagType>* dag;
    vector<int> var2eg;
    dagProblem problemType;
    string consFileName;
    DdManager* gbm;
    DdNode* consRoot;
	vector<unordered_map<DdNode*, int>> T;
	vector<int> egPath;
	bool isOutput;
    int numOfvars;
    int numOfDDnodes;
    int DDwidth;
    int tableLsize;
    int tableBsize;
    int numOfsteps; 
    long slideCount;
    double ddTime;
    double searchTime;

public:
    BcsAlgorithm(DAG<dagType>* dag, 
        string consFileName, 
        bool isShort = true,
		bool isOut = true);

	void setIsShort(bool newBool);

    int BCS(bool isDetail = false);

    double setCons(string consNewName);

	void printCons(int pr = 1, int n = 2);

	DdManager* getDdManager();
	
	DdNode* getConsRootNode();

	void writePictureConsBDD(string bddName = DEFAULT_PICTURE_NAME);

	int getConsTreeWidth();

	vector<unordered_map<DdNode*, int>>* getTtable();

	vector<int>* getEgPath();	

	void updateExpInfo(int &var, int &nodes, int &width, int& states, int &steps, long &traArcs, double &ddtime, double &searchtime, int &ddArcs);

private:
    inline DdNode* followBDD(DdNode* n, int eg);
};

/*
    public method
*/
template <typename dagType>
BcsAlgorithm<dagType>::BcsAlgorithm(DAG<dagType>* g, string consFN, bool isS, bool isO) {
    dag = g;
    consFileName = consFN;
    isShort = isS;
	isOutput = isO;
    setCons(consFileName);
}

template <typename dagType>
void BcsAlgorithm<dagType>::setIsShort(bool nb) {
    isShort = nb;
}

template <typename dagType>
vector<unordered_map<DdNode*, int>>* BcsAlgorithm<dagType>::getTtable() {
    return &T;
}

template <typename dagType>
int BcsAlgorithm<dagType>::BCS(bool isDetail) {
	struct timeval t1,t2;
		gettimeofday(&t1,NULL);
    int numOfvex = dag->getVertexNum();
	int step  = 0;
	int vt, wt, gro;
	egPath.clear();
	T.clear();
	T = vector<unordered_map<DdNode*, int>>(numOfvex);
#if BSTORE
	vector<unordered_map<DdNode*, pair<Edge*, DdNode*>>> B(numOfvex);
#endif
	// now BDD node
	DdNode* n = consRoot;
	DdNode* np;
	slideCount = 0;
	T[0][n] = 0;
	// BCS for edge [1...|E|] for (n,l) stored in T[vs]
	for (int s = 0; s < numOfvex; s++) {
		// use link table to travel edge
		vector<Edge>* vexEdge = dag->getEdgeList(s);
        Edge* curEdge;
		for (int i = 0; i < vexEdge->size(); i++) {
			curEdge = &((*vexEdge)[i]);
			vt = curEdge->target();
			wt = curEdge->weight();
			gro = curEdge->group();
			// continue
			for (auto nl : T[s]) {
				step++;
				DdNode* n = nl.first;
				int l = nl.second;
				np = followBDD(n, gro);
				// if n = nonT continue
				if (Cudd_IsConstant(np)) {
					if (Cudd_V(np) == 0) 
						continue;
				}else {
					// if label(n') = ei
					if (var2eg[Cudd_NodeReadIndex(np)] == gro) {
						slideCount++;
						np = Cudd_T(np);
					}
				}
				// if T[vt][n'] < l + wt
				if (T[vt].count(np) != 0) {
					// find other path
					if ((isShort && T[vt][np] > l + wt)||
							(!isShort && T[vt][np] < l + wt)) {
						T[vt][np] = l + wt;
						B[vt][np] = make_pair(curEdge, n);
					}
				}else {
					// haven't be traveled
					T[vt][np] = l + wt;
					B[vt][np] = make_pair(curEdge, n);
				}
			}
		}
	}
	DdNode* Tnode;
	int res = -1;
	bool reach = false;
	for(auto nl : T[numOfvex-1]) {
		Tnode = nl.first;
		if(Cudd_IsConstant(Tnode)) {
			if(Cudd_V(Tnode) == 1) {
				reach = true;
				break;
			}
		}
	}
	if(reach) {
		res =  T[numOfvex-1][Tnode];
		stack<int> egStack;
		auto p = B[numOfvex - 1][Tnode];
		while (p.first->source() != 0) {
            egStack.push(p.first->group());
            p = B[p.first->source()][p.second];
        }
		egStack.push(p.first->group());
	}
		gettimeofday(&t2,NULL);
	double timeUse = (t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0)*1000;
	if(isOutput) {
		string pathType = isShort ? "Short" : "Long";
		cout << "[BCS " << pathType << "est" << " Path's Length is:" ;
		if(reach) {
			cout << res;
		} else {
			cout << "UNREACHED";
		}
		cout << "]\n";
		cout << "\tstep number is:" << step << "\n";
		cout << "\tbcs time is:" << timeUse << "ms" << endl;
		if(slideCount != 0)
			cout << "\tbdd slide number is:" << slideCount << "\n";
	}
	// update info
    numOfsteps = step;
    searchTime = timeUse;
	numOfDDnodes = Cudd_DagSize(consRoot);
    tableLsize = 0;
    tableBsize = 0;
    for(int i = 0; i < T.size(); i++) {
        tableLsize += T[i].size();
    }
    for(int i = 0; i < B.size(); i++) {
        tableBsize += B[i].size();
    }
	// print path of node
#if DEBUG
	if(reach) {
        stack<int> pathStack;
        stack<int> egStack;
		stack<int> ngStack;
        auto p = B[numOfvex - 1][Tnode];
        cout << "\tprint path:";
        pathStack.push(p.first->target());
		ngStack.push(dag->getGroupFromNode(p.first->target()));
        //cout << " " << p.first->target() << " ";
        while (p.first->source() != 0) {
            pathStack.push(p.first->source());
            egStack.push(p.first->group());
			ngStack.push(dag->getGroupFromNode(p.first->source()));
            //cout << p.first->source() << " ";
            p = B[p.first->source()][p.second];
        }
        pathStack.push(0);
        egStack.push(p.first->group());
		ngStack.push(dag->getGroupFromNode(p.first->source()));
        while (!pathStack.empty()) {
            if(isDetail) cout << pathStack.top() << " ";
			//egPath.push_back(pathStack.top());
            pathStack.pop();
        }
        cout << "\n\tprint edge group:";
        while (!egStack.empty()) {
            if(isDetail) cout << egStack.top() << " ";
			egPath.push_back(egStack.top());
            egStack.pop();
        }
		cout << "\n\tprint node group:";
        while (!ngStack.empty()) {
            if(isDetail) cout << ngStack.top() << " ";
            ngStack.pop();
        }
        cout << "\n";
    } 
	/*
	if(reach) {
		auto p = B[numOfvex - 1][Tnode];
		cout << "\tprint path:";
		cout << "\t" << dag->getVertex(p.first->target()) << ":" << p.first->target() << " <- ";
		while (p.first->source() != 0) {
			cout << dag->getVertex(p.first->source()) << ":" << p.first->source() << " <- ";
			p = B[p.first->source()][p.second];
		}
		cout << "source:0" <<"\n";
	}
	*/
#endif
	// print the table of T[v][n]
	if(!isDetail) return res;
	cout << "\tprint the T table:\n";
	for (int s = 0; s < numOfvex; s++){
		cout << "\t#\n";
		for(auto tp : T[s]){
			cout << "\tT["<< s << "][" << tp.first << "] = " << tp.second << "\n";
		}
	}
	cout << endl;
	return res;
}

template <typename dagType>
double BcsAlgorithm<dagType>::setCons(string consNewName) {
	struct timeval t1,t2;
		gettimeofday(&t1,NULL);
    gbm = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    consFileName = consNewName;
    var2eg.clear();
    // follow the order of group
	string line;
	ifstream in(consFileName, ios::in);
	bool isfirstn = true;
	bool isfirstr = true;
	bool isCnf;
	getline(in, line);
	if(line.find("cnf") != string::npos) {
		isCnf = true;
	}else {
		if(line.find("dnf") != string::npos) {
			isCnf = false;
		}else {
			cerr << "[ERROR]Input error format constrain from file:" 
				<< consFileName << "]" << endl;
			return 0;
		}
	}
    vector<int> clause = commonTool::sepnum(line);
	int numofcla = clause[3];
	int countOfindex = 0;
	int posOfv2eg;
	DdNode* tmp1;
	DdNode* tmp2;
	DdNode* root;
    // read all group index
	for(int j = 0; j < numofcla; j++){
		getline(in, line);
		clause = commonTool::sepnum(line);
		for(int i : clause){
			if(i==0) break;
			if(commonTool::findExistEG(i, &var2eg) != -1){
				continue;
			}else{
				var2eg.push_back(abs(i));
			}
		}
	}
    // ascend order to assign var index
	sort(var2eg.begin(),var2eg.end(),less<int>());
	for(int i = 0; i < var2eg.size(); i++){
		DdNode* nnode = Cudd_bddIthVar(gbm, i);
	}
    // struct BDD
	in.seekg(in.beg);
	getline(in, line);
	for(int j = 0; j < numofcla; j++){
		getline(in, line);
		clause = commonTool::sepnum(line);
		for(int i : clause){
			//cerr << "Now i is : " << i << "\n";
			if(i == 0) break;
			DdNode* nnode;
			if((posOfv2eg=commonTool::findExistEG(i, &var2eg)) == -1){
				// gro haven't exist
				cerr << "[ERROR] Unmarked var while read bdd\n";
                return 0;
			}else{
				// gro have existed
				//cerr << posOfn2g << "\n";
				nnode = Cudd_bddIthVar(gbm, posOfv2eg);
			}
			// true or false
			if(i < 0){
				nnode = Cudd_Not(nnode);
				Cudd_Ref(nnode);
			}
			//nnode = Cudd_bddIthVar(gbm, countOfindex++);
			//index2group.push_back(i);
			if(isfirstn){
				tmp1 = nnode;
				tmp2 = nnode;
				isfirstn = false;
			}else{
				tmp1 = isCnf ? Cudd_bddOr(gbm, nnode, tmp1) : Cudd_bddAnd(gbm, nnode, tmp1);
				Cudd_Ref(tmp1);
				Cudd_RecursiveDeref(gbm, tmp2);
				tmp2 = tmp1;
			}
		}
		isfirstn = true;
		if(isfirstr){
			root = tmp1;
			isfirstr = false;
		}else{
			tmp2 = root;
			root = isCnf ? Cudd_bddAnd(gbm, root, tmp1) : Cudd_bddOr(gbm, root, tmp1);
			Cudd_Ref(root);
			Cudd_RecursiveDeref(gbm, tmp1);
			Cudd_RecursiveDeref(gbm, tmp2);
		}
	}
	consRoot = Cudd_BddToAdd(gbm, root);
		gettimeofday(&t2,NULL);
	double timeUse = (t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0)*1000;
	if(isOutput) {
		cout << "[SET NEW CONS FROM FILE:" << consFileName << "]\n";
		cout << "\tbdd time:" << timeUse << "ms" << endl;
	}
	numOfvars = var2eg.size();
    ddTime = timeUse;
    DDwidth = getConsTreeWidth();

	return timeUse;
}

template <typename dagType>
void BcsAlgorithm<dagType>::printCons(int pr, int n) {
	cout << "[PRINT CONS BDD IN:" << consFileName << "]\n";
	cout << "\tDdManager nodes number:" << Cudd_ReadNodeCount(gbm) << "\n";
	cout << "\tDdManager var number:" << Cudd_ReadSize(gbm) << "\n";
	cout << "\tDdManager reorderings:" << Cudd_ReadReorderings(gbm) << "\n";
	cout << "\tDdManager memory:" << Cudd_ReadMemoryInUse(gbm) << "\n";
	cout << "\tcons bdd detail";
    Cudd_PrintDebug(gbm, consRoot, n, pr);
}

#if DEBUG
template <typename dagType>
DdManager* BcsAlgorithm<dagType>::getDdManager() {
	return gbm;
}

template <typename dagType>
DdNode* BcsAlgorithm<dagType>::getConsRootNode() {
	return consRoot;
}

#endif

/*
    private method
*/

template <typename dagType>
inline DdNode* BcsAlgorithm<dagType>::followBDD(DdNode* n, int eg) {
    if (Cudd_IsConstant(n)) return n;
	while (var2eg[Cudd_NodeReadIndex(n)] < eg) {
		slideCount++;
		n = Cudd_E(n);
		if (Cudd_IsConstant(n)) break;
	}
	return n;
}


template <typename dagType>
void BcsAlgorithm<dagType>::writePictureConsBDD(string bddName) {
	string fileName = DEFAULT_PICTURE_DIR + bddName + ".dot";
	char** inputNames;
    char** outputName = new char* [1];
    outputName[0] = const_cast<char*>(bddName.c_str());
    DdNode **outputs = new DdNode* [1];
	outputs[0] = consRoot;
    FILE *f = fopen(fileName.c_str(), "w");
	Cudd_DumpDot(gbm, 1, outputs, inputNames, outputName, f);
	cout << "[WRITE PICTURE BDD IN FILE:" << fileName << "]" << endl;
}

template <typename dagType>
int BcsAlgorithm<dagType>::getConsTreeWidth() {
	int varNum = var2eg.size();
	int width = 0;
	priority_queue<DdNode*, vector<DdNode*>, sortBddNode> nodeQue;
	unordered_set<DdNode*> nodeSet;
	nodeQue.push(consRoot);
	nodeSet.emplace(consRoot);  

	for(int i = 0; i < varNum; i++) {
		int levelWidth;
		/*
		priority_queue<DdNode*, vector<DdNode*>, sortBddNode> nodeQue_copy = nodeQue;
		while(!nodeQue_copy.empty()) {
			cout << Cudd_NodeReadIndex(nodeQue_copy.top()) << " ";
			nodeQue_copy.pop();
		}
		cout << endl;
		*/
		while(Cudd_NodeReadIndex(nodeQue.top()) == i) {
			//cerr << "\ttop index = " << Cudd_NodeReadIndex(nodeQue.top()) << endl;
			DdNode* Pnode = nodeQue.top();
			DdNode* Tnode = Cudd_T(Pnode);
			DdNode* Enode = Cudd_E(Pnode);
			nodeQue.pop();
			if(nodeSet.find(Tnode) == nodeSet.end()) {
				nodeSet.emplace(Tnode);
				nodeQue.push(Tnode);
			}
			if(nodeSet.find(Enode) == nodeSet.end()) {
				nodeSet.emplace(Enode);
				nodeQue.push(Enode);
			}
		}
		levelWidth = nodeQue.size();
		//cerr << "varlevel = " << i << " ; level width = " << levelWidth << endl;
		width = max(width, levelWidth);
	}
	if(isOutput)
		cout << "[GET WIDTH OF CONS BDD]:" << width << endl;
	
	return width;
}

template <typename dagType>
vector<int>* BcsAlgorithm<dagType>::getEgPath() {
	return &egPath;
}

template <typename dagType>
void BcsAlgorithm<dagType>::updateExpInfo(int &var, int &nodes, int &width, int& states, int &steps, long &traArcs, double &ddtime, double &searchtime, int &ddArcs) {
    var = numOfvars;
    nodes = numOfDDnodes;
    width = DDwidth;
    if(tableLsize == tableBsize + 1) {
        states = tableLsize;
    } else {
        cerr << "[#ERROR#] L size != B size + 1" << endl;
        cerr << "\t L size = " <<  tableLsize << endl;
        cerr << "\t B size = " <<  tableBsize << endl;
        exit(1);
    }
    steps = numOfsteps;
    traArcs = slideCount;
    ddtime = ddTime;
    searchtime = searchTime;
	ddArcs = (numOfDDnodes - 2) * 2;
}

#endif