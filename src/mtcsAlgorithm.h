#ifndef MTCSALGORITHM_H
#define MTCSALGORITHM_H
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <stack>
#include <cmath>
#include "dag.h"
#include "meddly.h"

using namespace MEDDLY;

struct fMddNode{
    int ng;
    int level;
    fMddNode* other;
    vector<fMddNode*> children;
};

struct sortMddNode {  
    bool operator()(fMddNode* a, fMddNode* b) const {
        return a->level < b->level;
    }  
};  

template <typename dagType>
class MtcsAlgorithm {
private:
    bool isShort;
    // node vertex(int) means node group
    DAG<dagType>* dag;
    vector<int> var2ng;
    vector<vector<int>> varsMddEdge2eg;
    vector<unordered_map<int, int>> varsEg2mddEdge;
    vector<unpacked_node*> mddNodeList;
    vector<int> mddNode2level;
    //dagProblem problemType;
    string consFileName;
    domain* d;
    expert_forest* f;
    policies p;
    dd_edge* consRoot;
    fMddNode* fconsRoot;
    fMddNode* fmddTnode;
    fMddNode* fmddFnode;
    vector<unordered_map<node_handle, int>> T_slow;
    vector<int> egPath;
    bool notDomainReduction;
    unordered_map<node_handle, fMddNode*> fmddVisited;
    vector<unordered_map<int, int>> varsEgisIn;
    bool isOutput;
    bool isGrouped;
    // exp info
    int numOfvars;
    int numOfDDnodes;
    int numOfarcs;
    int DDwidth;
    int tableLsize;
    int tableBsize;
    int numOfsteps; 
    long slideCount;
    double ddTime;
    double searchTime;

public:
    MtcsAlgorithm(DAG<dagType>* dag, 
        string consFileName, 
        bool isShort = true,
        bool isoutput = true,
        bool isGro = false,
        bool noDoReduce = false);

	void setIsShort(bool newBool);

	/*
		mtcs algorithm
		return the short/long est path
	*/
    int MTCS(bool isDetail = false);

    int MTCS_F(bool isDetail = false);

    double setCons(string consNewName);

	void printCons(int printOption = 1);

    double updateFastMdd();

    inline fMddNode* createFMddNode(int ng, int level);

    void destoryFastMdd();

	void writePictureConsMDD(string mddName = DEFAULT_MDD_PICTURE_NAME);

    vector<unordered_map<node_handle, int>>* getTtable();

    vector<int>* getEgPath();

    int getConsTreeWidth();

    void updateExpInfo(int &var, int &nodes, int &width, int& states, int &steps, long &traArcs, double &ddtime, double &searchtime, int &ddarcs);

private:
    inline node_handle followMTDD(node_handle n, int vsng, int vtng, int eg);

    inline fMddNode* followMTDD_F(fMddNode* curNode, int vsng, int eg);

    inline fMddNode* followMTDD_F_G(fMddNode* curNode, int vsng, int vtng, int eg);

    inline fMddNode* followMTDD_noDReduce(fMddNode* curNode, int vsng, int vtng, int eg);
};

/*
    method
*/
template <typename dagType>
MtcsAlgorithm<dagType>::MtcsAlgorithm(DAG<dagType>* g, string consFN, bool isS, bool isout, bool isGro, bool nDR) {
    dag = g;
    consFileName = consFN;
    isShort = isS;
    slideCount = 0;
    consRoot = nullptr;
    fconsRoot = nullptr;
    notDomainReduction = nDR;
    isOutput = isout;
    isGrouped = isGro;
    setCons(consFileName);
}

template <typename dagType>
void MtcsAlgorithm<dagType>::setIsShort(bool nb) {
    isShort = nb;
}

template <typename dagType>
vector<unordered_map<node_handle, int>>* MtcsAlgorithm<dagType>::getTtable() {
    return &T_slow;
}

template <typename dagType>
void MtcsAlgorithm<dagType>::printCons(int op) {
    cout << "[PRINT CONS MDD IN:" << consFileName << "]\n\t";
    ostream_output meddlyout(cout);
    consRoot->show(meddlyout, op);
}

template <typename dagType>
void MtcsAlgorithm<dagType>::writePictureConsMDD(string mddName) {
    string fileName = DEFAULT_PICTURE_DIR + mddName;
    consRoot->writePicture(fileName.c_str(), "pdf");
    cout << "[WRITE PICTURE MDD IN FILE:" << fileName << ".pdf]" << endl;
}

template <typename dagType>
inline node_handle MtcsAlgorithm<dagType>::followMTDD(node_handle n, int vsng, int vtng, int eg) {
    unpacked_node* curNode;
    if(f->isTerminalNode(n))
        return n;
#if !GENERAL
    int mddNg = var2ng[f->getNodeLevel(n) - 1];
    if(vsng < mddNg)
        return n;
    if(vsng > mddNg)
        return 0;
    int varLevel = f->getNodeLevel(n);
    int egPosNum;
    curNode = mddNodeList[n];
    auto egPos = varsEg2mddEdge[varLevel].find(eg);
    if(egPos == varsEg2mddEdge[varLevel].end()) {
        egPosNum = varsEg2mddEdge[varLevel][-10];
    } else {
        egPosNum = egPos->second;
    }
    slideCount++;
    n = curNode->d(egPosNum);
#else
    if(vsng == var2ng[f->getNodeLevel(n) - 1]) {
        slideCount++;
        // nodel handel 0 means F
        int egPosNum;
        int varLevel = f->getNodeLevel(n);
        auto egPos = varsEg2mddEdge[varLevel].find(eg);
        if(egPos == varsEg2mddEdge[varLevel].end()) {
            egPosNum = varsEg2mddEdge[varLevel][-10];
        } else {
            egPosNum = egPos->second;
        }
        curNode = mddNodeList[n];
        n = curNode->d(egPosNum);
        if(f->isTerminalNode(n))
            return n;
    }
    while(vtng > var2ng[f->getNodeLevel(n) - 1]) {
        slideCount++;
        curNode = mddNodeList[n];
        n = curNode->d( curNode->getSize()-1 );
        if(f->isTerminalNode(n))
            return n;
    }
#endif
    /*
    if(eg == 13011) {
        cerr << "; n after = " << n << endl;
    }
    */
    return n;
}

template <typename dagType>
inline fMddNode* MtcsAlgorithm<dagType>::followMTDD_F(fMddNode* n, int vsng, int eg) {
    if(n->ng < 1 || vsng < n->ng) return n;
    if(vsng > n->ng) return fmddFnode;
    // vsng == n->ng
    auto egPos = varsEg2mddEdge[n->level].find(eg);
    if(egPos == varsEg2mddEdge[n->level].end()) {
        n = n->other;
    } else {
        n = n->children[egPos->second];
    }
    slideCount++;
    return n;
}

template <typename dagType>
inline fMddNode* MtcsAlgorithm<dagType>::followMTDD_F_G(fMddNode* n, int vsng, int vtng, int eg) {
    if(n->ng < 1)
        return n;
    if(vsng == n->ng) {
        slideCount++;
        // nodel handel 0 means F
        auto egPos = varsEg2mddEdge[n->level].find(eg);
        if(egPos == varsEg2mddEdge[n->level].end()) {
            //n = n->children[0];
            n = n->other;
        } else {
            n = n->children[egPos->second];
        }
        /*
        if(varsEgisIn[n->level][eg]) {
            n = n->children[varsEg2mddEdge[n->level][eg]];
        } else {
            n = n->other;
        }
        */
       if(n->ng < 1) return n;
    }
    while(vtng > n->ng) {
        slideCount++;
        //n = n->children[0];
        n = n->other;
        if(n->ng < 1) return n;
    }
    return n;
}

template <typename dagType>
inline fMddNode* MtcsAlgorithm<dagType>::followMTDD_noDReduce(fMddNode* n, int vsng, int vtng, int eg) {
    if(n->ng < 1)
        return n;
    if(vsng == n->ng) {
        slideCount++;
        // nodel handel 0 means F
        n = n->children[varsEg2mddEdge[n->level][eg]];
        if(n->ng < 1) return n;
    }
    while(vtng > n->ng) {
        slideCount++;
        n = n->other;
        if(n->ng < 1) return n;
    }
    return n;
}

template <typename dagType>
double MtcsAlgorithm<dagType>::setCons(string consNewName) {
    struct timeval t1,t2,t3,t4;
    consFileName = consNewName;
    var2ng.clear();
    varsMddEdge2eg.clear();
    varsEg2mddEdge.clear();
    mddNodeList.clear();
    mddNode2level.clear();
    varsEgisIn.clear();
    unordered_map<int, int> ng2var;
    vector<vector<pair<int, int>>> consInfo;
    vector<int> consFileInfo;
    unordered_map<string, dd_edge> nodeCahe;
    //cleanup();
    //initialize();
        gettimeofday(&t1,NULL);
    
    string line;
	ifstream in(consFileName, ios::in);
	bool isfirstn = true;
	bool isfirstr = true;
    bool isfirsti = true;
    bool isfirstin = true;
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
    /*
		there we need now:
		var number and map ng to var
		var's bound, i.e. how many edges is constrained from the target.
	*/
    vector<int> clause = commonTool::sepnum(line);
	int numofcla = clause[3];
	int countOfindex = 0;
	int varLevel;
    // read all node group index and sort them
	for(int j = 0; j < numofcla; j++){
        //cerr << "[]is reading clause:" << j << endl;
        vector<pair<int, int>> rowInfo;
		getline(in, line);
		clause = commonTool::sepnum(line);
		for(int i : clause){
            //cerr << "\t[]is reading var:" << i << endl;
			if(i==0) break;
            //isEgInDomain[ abs(i) ] = true;
            consFileInfo.push_back( abs(i) );
			Edge e = dag->getEdgeFromGro( abs(i) );
			int sourNodeGro = dag->getGroupFromNode( e.source() );
            rowInfo.push_back( make_pair(i > 0 ? sourNodeGro : -sourNodeGro, abs(i)) );
            if(commonTool::findExistEG(sourNodeGro, &var2ng) == -1)
                var2ng.push_back(sourNodeGro);
		}
        consInfo.push_back(rowInfo);
	}
    sort(var2ng.begin(), var2ng.end(), greater<int>());
    for(int i = 0; i < var2ng.size(); i++) {
        ng2var[var2ng[i]] = i;
    }
    int bounds[ var2ng.size() ];
    varsMddEdge2eg = vector<vector<int>>(var2ng.size()+1);
    varsEg2mddEdge = vector<unordered_map<int, int>>(var2ng.size()+1);
    //varsEgisIn = vector<unordered_map<int, int>>(var2ng.size()+1);
    // here programa finish the work to give any vertice a var
    if(notDomainReduction) {
        for(int i = 0; i < var2ng.size(); i++) {
            int var = i + 1;
            int sourNodeGro = var2ng[i];
            int sourceNode = dag->getNodeFromGroup(sourNodeGro);
            if(!isGrouped) {
                varsMddEdge2eg[var].push_back(-10);
                varsEg2mddEdge[var][-10] = varsMddEdge2eg[var].size() - 1;

                vector<Edge>* edgeList = dag->getEdgeList(sourceNode);
                for(int j = 0; j < edgeList->size(); j++) {
                    Edge* e = &(edgeList->at(j));
                    varsMddEdge2eg[var].push_back(e->group());
                    varsEg2mddEdge[var][e->group()] = varsMddEdge2eg[var].size() - 1;
                    //varsEgisIn[var][e->group()] = 1;
                }
            } else {
                for(auto ng : dag->ng2itsEg[sourNodeGro]) {
                    varsMddEdge2eg[var].push_back(ng);
                    varsEg2mddEdge[var][ng] = varsMddEdge2eg[var].size() - 1;
                    //varsEgisIn[var][ng] = 1;
                }
            }
            bounds[i] = varsMddEdge2eg[i+1].size();
        }
    } else {
        /*
        for(int i = 0; i < var2ng.size(); i++) {
            int var = i + 1;
            int sourNodeGro = var2ng[i];
            int sourceNode = dag->getNodeFromGroup(sourNodeGro);
            vector<Edge>* edgeList = dag->getEdgeList(sourceNode);
            for(int j = 0; j < edgeList->size(); j++) {
                Edge* e = &(edgeList->at(j));
                varsEgisIn[var][e->group()] = 0;
            }
        }
        */
        for(int i : consFileInfo) {
            Edge e = dag->getEdgeFromGro( abs(i) );
            int sourNodeGro = dag->getGroupFromNode( e.source() );
            int var = ng2var[sourNodeGro] + 1;
            if(varsEg2mddEdge[var].count(i) != 0)
                continue;
            varsMddEdge2eg[var].push_back(i);
            varsEg2mddEdge[var][i] = varsMddEdge2eg[var].size() - 1;
            //varsEgisIn[var][i] = 1;
        }

        for(int i = 1; i < varsMddEdge2eg.size(); i++) {
            if(isGrouped) {
                int node = dag->getNodeFromGroup( var2ng[i - 1] );
                vector<Edge>* nodeEdgeList = dag->getEdgeList(node);
                int layerEgNum = nodeEdgeList->size();
                if(varsMddEdge2eg[i].size() == layerEgNum)
                    continue;
            }
            varsMddEdge2eg[i].insert(varsMddEdge2eg[i].begin(), -10);
            for (auto& pair : varsEg2mddEdge[i]) {
                pair.second += 1;
            }
            varsEg2mddEdge[i][-10] = 0;
        }
/*
        for(int i = 1; i < varsMddEdge2eg.size(); i++) {
#if !GENERAL
            int node = dag->getNodeFromGroup( var2ng[i - 1] );
            vector<Edge>* nodeEdgeList = dag->getEdgeList(node);
            int layerEgNum = nodeEdgeList->size();
            if(varsMddEdge2eg[i].size() == layerEgNum)
                continue;
#endif
            varsMddEdge2eg[i].push_back(-10);
            varsEg2mddEdge[i][-10] = varsMddEdge2eg[i].size() - 1;
        }
*/
        for(int i = 0; i < var2ng.size(); i++) {
            bounds[i] = varsMddEdge2eg[i+1].size();
        }
    }
    // debug start
    /*
    cerr << "[DEBUG INFO]:bounds\n\t";
    for(int i = 0; i < var2ng.size(); i++) {
        cerr << bounds[i] << " ";
    }
    cerr << endl;
    */
    // debug finish
    d = createDomain();
    d->createVariablesBottomUp(bounds, var2ng.size());
    p.useDefaults(false);
    p.setPessimistic();
    p.setFullyReduced();
    //p.setQuasiReduced();
    //p.setIdentityReduced();
    expert_forest* mdd = (expert_forest*)d->createForest(
      false,                    // this is not a relation
      range_type::BOOLEAN,          // terminals are either true or false
      edge_labeling::MULTI_TERMINAL,    // disables edge-labeling
      p
      );
    f = (expert_forest*)d->createForest(
      false,                    // this is not a relation
      range_type::BOOLEAN,          // terminals are either true or false
      edge_labeling::MULTI_TERMINAL,    // disables edge-labeling
      p
      );
#if 0
    cout << "[DEBUG INFO]PRINTING var2ng:";
    for(auto i : var2ng) {
        cout << i << " ";
    }
    cout << endl;

    cout << "[DEBUG INFO]PRINTING varsMddEdge2eg:";
    for(auto a : varsMddEdge2eg) {
        cout << "\n\t";
        for(auto b : a) {
            cout << b << " ";
        }
    }
    cout << endl;

    cout << "[DEBUG INFO]PRINTING varsEg2mddEdge:";
    for(auto a : varsEg2mddEdge) {
        cout << "\n\t";
        for(auto b : a) {
            cout  << b.first << ":" << b.second << " ";
        }
    }
    cout << endl;
    
    cout << "[DEBUG INFO]PRINTING consInfo:";
    for(auto a : consInfo) {
        cout << "\n\t";
        for(auto b : a) {
            cout << "[" << b.first << " " << b.second << "] ";
        }
    }
    cout << endl;
#endif
    ostream_output meddlyout(cout);
    dd_edge root(mdd);
    dd_edge logicTrue(mdd);
    dd_edge tmp(mdd);
    dd_edge iteratRoot(mdd);
    dd_edge iteratTmp(mdd);
    mdd->createEdge(true, logicTrue);
    // assign mdd node for every var
    for(auto row : consInfo) {
        for(auto info : row) {
            dd_edge nnode(mdd);
            string cacheKey = to_string(info.first) + " " + to_string(info.second);
            auto cacheInfo = nodeCahe.find(cacheKey);
            if(cacheInfo != nodeCahe.end()) {
                continue;
            } else {
                bool isNot = info.first < 0;
                int sour = isNot ? -info.first : info.first;
                int eg = info.second;
                if((varLevel=commonTool::findExistEG(sour, &var2ng)) == -1) {
                    cerr << "[ERROR] Unmarked var while struct mdd\n";
                    return 0;
                } 
                // leve 0 is termianl
                varLevel++;
                // contruct mdd node
                int edgeNum = varsMddEdge2eg[varLevel].size();
                int trueEdge = varsEg2mddEdge[varLevel][eg];
                bool* assign = new bool[edgeNum];
                for(int i = 0; i < edgeNum; i++) {
                    assign[i] = false;
                }
                assign[trueEdge] = true;
                mdd->createEdgeForVar(varLevel, false, assign, nnode);
                // var or -var
                if(isNot) apply(DIFFERENCE, logicTrue, nnode, nnode);
                //nnode.show(meddlyout, 2);
                /*
                try {
                    nnode.show(meddlyout, 2);
                } catch (error e) {
                    cerr << "\n[MEDDLY ERROR WHILE MAKE NODE]\n";
                    cerr << "\tError name: " << e.getName() << "\n";
                    cerr << "\tIn file: " << e.getFile() << "\n";
                    cerr << "\tIn line: " << e.getLine() << "\n";
                    exit(EXIT_FAILURE);
                }
                */
                nodeCahe[cacheKey] = nnode;
            }
        }
    }
    // struct MDD
    int interaCount = 0;
    int interaRowCount = 0;
    int interaLimit = 100;
//[old time]
        //gettimeofday(&t1,NULL);
    int countRow = 1;
    for(auto row : consInfo) {
        //cerr << "[]handeling row:" << countRow++ << "/" << consInfo.size() << endl;
        // int coutCol = 1;
        for(auto info : row) {
            //cerr << "\t[]handeling col:" << coutCol++ << "/" << row.size() << " info:" << info.first << " " << info.second <<endl;
            dd_edge nnode(mdd);
            string cacheKey = to_string(info.first) + " " + to_string(info.second);
            auto cacheInfo = nodeCahe.find(cacheKey);
            if(cacheInfo == nodeCahe.end()) {
                cerr << "[ERROR] Uassigned var while struct mdd\n";
                return 0; 
            }
            nnode = cacheInfo->second;
            if(isfirstin){
                //cerr << "\t\tis COPY" << endl;
				apply(COPY, nnode, iteratTmp);
                //tmp = nnode;
                //cerr << "\t\tCOPY finish" << endl;
				isfirstin = false;
			}else{
                //cerr << "\t\tis caluse handel" << endl;
                apply(isCnf ? UNION : INTERSECTION, nnode, iteratTmp, iteratTmp);
                //cerr << "\t\tcaluse handel finish" << endl;
			}
            interaRowCount++;

            if(interaRowCount < interaLimit) continue;
            isfirstin = true;
            interaRowCount = 0;
            if(isfirstn) {
                apply(COPY, iteratTmp, tmp);
                isfirstn = false;
            } else {
                apply(isCnf ? UNION : INTERSECTION, iteratTmp, tmp, tmp);
            }
        }
        if(row.size() < interaLimit) {
            tmp = iteratTmp;
        } else {
            apply(isCnf ? UNION : INTERSECTION, iteratTmp, tmp, tmp);
        }
        isfirstin = true;
        isfirstn = true;
        interaRowCount = 0;
        if(isfirsti){
            iteratRoot = tmp;
            isfirsti = false;
        }else{
            //cerr << "\t\tis itera root handel" << endl;
            apply(isCnf ? INTERSECTION : UNION, tmp, iteratRoot, iteratRoot);
            //cerr << "\t\troot itera handel finish" << endl;
        }
        interaCount++;
        if(interaCount < interaLimit) continue;

        isfirsti = true;
        interaCount = 0;
        if(isfirstr){
            root = iteratRoot;
            isfirstr = false;
        }else{
            //cerr << "\t\tis root handel" << endl;
            apply(isCnf ? INTERSECTION : UNION, iteratRoot, root, root);
            //cerr << "\t\troot handel finish" << endl;
        }
    }
    if(consInfo.size() < interaLimit) {
        root = iteratRoot;
    } else {
        apply(isCnf ? INTERSECTION : UNION, iteratRoot, root, root);
    }
        gettimeofday(&t2,NULL);
    double timeUse = (t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0)*1000;
    if(consRoot != nullptr)
        delete consRoot;
    consRoot = new dd_edge(f);
    try {
        apply(COPY, root, *consRoot);
    } catch (error e) {
        cerr << "[MEDDLY ERROR]\n";
        cerr << "\tError name: " << e.getName() << "\n";
        cerr << "\tIn file: " << e.getFile() << "\n";
        cerr << "\tIn line: " << e.getLine() << "\n";
        exit(EXIT_FAILURE);
    }
    // store mdd node information for efficiency
    unordered_set<node_handle> visited;
    int headNodeHandle = consRoot->getNode();
    stack<node_handle> nodeSt;
    mddNodeList = vector<unpacked_node*>(headNodeHandle + 1);
    mddNode2level = vector<int>(headNodeHandle + 1);
    visited.emplace(0);
    visited.emplace(-1);
    nodeSt.push(headNodeHandle);
    while(!nodeSt.empty()) {
        int curNodeHandel = nodeSt.top();
        nodeSt.pop();
        unpacked_node* curNode = f->newUnpacked(curNodeHandel, FULL_ONLY);
        mddNodeList[curNodeHandel] = curNode;
        mddNode2level[curNodeHandel] = f->getNodeLevel(curNodeHandel);
        visited.emplace(curNodeHandel);
        for(int i = 0; i < curNode->getNNZs(); i++) {
            node_handle nextNode = curNode->d(i);
            if(visited.count(nextNode) == 0)
                nodeSt.push(nextNode);
        }
    }
    // fast mdd
    gettimeofday(&t3,NULL);
    double fMddSize = updateFastMdd();
    gettimeofday(&t4,NULL);
    double fMddTimeUse = (t4.tv_sec - t3.tv_sec + (t4.tv_usec - t3.tv_usec)/1000000.0)*1000;
    if(isOutput) {
        cout << "[MTCS SET NEW CONS FROM FILE:" << consFileName << "]\n";
        cout << "\tmdd type:"  << (notDomainReduction ? "orignal" : "domain reduction") << "\n";
        cout << "\tmdd time:"  << timeUse << "ms" << "\n";
        cout << "\tmdd node number:" << visited.size() << "\n";
        cout << "\tfmdd time:"  << fMddTimeUse << "ms" << "\n";
        cout << "\tfmdd node number:"  << fmddVisited.size() << endl;
        cout << "\tfmdd size (number of arcs):"  << (int)fMddSize << endl;
    }
    for(int i = 1; i < varsMddEdge2eg.size(); i++) {
        varsEg2mddEdge[i].erase(-10);
    }
    numOfvars = var2ng.size();
    ddTime = timeUse;
    numOfDDnodes = fmddVisited.size();
    numOfarcs = (int)fMddSize;
    DDwidth = getConsTreeWidth();
    return timeUse;
}

template <typename dagType>
int MtcsAlgorithm<dagType>::MTCS(bool isDetail) {
    struct timeval t1,t2;
    int numOfvex = dag->getVertexNum();
	int step  = 0;
	int vt, wt, vsng, vtng, eg;
    egPath.clear();
    T_slow.clear();
    T_slow = vector<unordered_map<node_handle, int>>(numOfvex);
#if DEBUG
	vector<unordered_map<node_handle, pair<Edge*, node_handle>>> B(numOfvex);
#endif
    node_handle n = consRoot->getNode();
    node_handle np;
    //int debug_count = 0;
    slideCount = 0;
    T_slow[0][n] = 0;
        gettimeofday(&t1,NULL);
    for (int s = 0; s < numOfvex; s++) {
		// use link table to travel edge
        //cerr << "Is handeling node s : " << s << endl;
        //if(s == 550) cerr << "node " << s << "'s step number:" << step << endl;
		vector<Edge>* vexEdge = dag->getEdgeList(s);
        Edge* curEdge;
		for (int i = 0; i < vexEdge->size(); i++) {
			curEdge = &((*vexEdge)[i]);
			vt = curEdge->target();
			wt = curEdge->weight();
            eg = curEdge->group();
            vsng = dag->getGroupFromNode(s);
            vtng = dag->getGroupFromNode(vt);
            /*
            if(vsng == 901 && curEdge->group() == 21355) {
                cerr << "\tedge_wt = " << wt << " edge_gro = " << curEdge->group()  << " edge number = " << vexEdge->size() << endl;
            }
            */
			for (auto nl : T_slow[s]) {
                /*
                if(vsng == 901 && debug_count < 4) {
                    cerr << "vsng = " << vsng << " vt = " << vt << " mdd node = " << nl.first << endl;
                    debug_count++;
                }
                */
				step++;
				node_handle n = nl.first;
				int l = nl.second;
                //cerr << "\tfollowMTDD start" << endl;
				np = followMTDD(n, vsng, vtng, eg);
                //cerr << "\tfollowMTDD end" << endl;
				// if n = nonT continue
				if(np == 0) {
                    continue;
                }
				// if T[vt][n'] < l + wt
				if (T_slow[vt].count(np) != 0) {
					// find other path
					if ((isShort && T_slow[vt][np] > l + wt)||
							(!isShort && T_slow[vt][np] < l + wt)) {
						T_slow[vt][np] = l + wt;
#if DEBUG
						B[vt][np] = make_pair(curEdge, n);
#endif
					}
				} else {
					// haven't be traveled
					T_slow[vt][np] = l + wt;
#if DEBUG
					B[vt][np] = make_pair(curEdge, n);
#endif
				}
			}
		}
	}
		gettimeofday(&t2,NULL);
    double timeUse = (t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0)*1000;
	int res =  -1;
    if (T_slow[numOfvex-1].count(-1) != 0) {
        res =  T_slow[numOfvex-1][-1];
    }
	string pathType = isShort ? "Short" : "Long";
	cout << "[MTCS " << pathType << "est" << " Path's Length is:";
    if(res != -1) {
        cout << res;
    } else {
        cout << "UNREACHED";
    }
    cout << "]\n";
	cout << "\tstep number is:" << step << "\n";
	cout << fixed << setprecision(3) <<"\tmtcs time is:" << timeUse << "ms" << endl;
    if(slideCount != 0)
        cout << "\tmdd slide number is:" << slideCount << "\n";
	// print path of node
#if DEBUG
    if(res != -1) {
        stack<int> pathStack;
        stack<int> egStack;
        stack<int> ngStack;
        auto p = B[numOfvex - 1][-1];
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
#endif
    // debug
    /*
    int showTnode = 21;
    for(auto tp : T_slow[showTnode]){
		cout << "\tT["<< showTnode << "][" << tp.first << "] = " << tp.second << "\n";
	}
    */
    if(!isDetail) return res;
	// print the table of T[v][n]
	cout << "\tprint the T table:\n";
	for (int s = 0; s < numOfvex; s++){
        if(T_slow[s].size() == 0) continue;
		cout << "\t#\n";
		for(auto tp : T_slow[s]){
			cout << "\tT["<< s << "][" << tp.first << "] = " << tp.second << "\n";
		}
	}
    
	cout << endl;
	return res;
}

template <typename dagType>
int MtcsAlgorithm<dagType>::MTCS_F(bool isDetail) {
    struct timeval t1,t2;
        gettimeofday(&t1,NULL);
    int numOfvex = dag->getVertexNum();
	int step  = 0;
	int vt, wt, vsng, vtng, eg;
	vector<unordered_map<fMddNode*, int>> T(numOfvex);
	vector<unordered_map<fMddNode*, pair<Edge*, fMddNode*>>> B(numOfvex);
    fMddNode* n = fconsRoot;
    fMddNode* np;
    //int debug_count = 0;
    slideCount = 0;
    T[0][n] = 0;
    for (int s = 0; s < numOfvex; s++) {
		// use link table to travel edge
        //cerr << "Is handeling node s : " << s << endl;
        //if(s == 18287) cerr << "node " << s << "'s step number:" << step << endl;
		vector<Edge>* vexEdge = dag->getEdgeList(s);
        Edge* curEdge;
		for (int i = 0; i < vexEdge->size(); i++) {
			curEdge = &((*vexEdge)[i]);
			vt = curEdge->target();
			wt = curEdge->weight();
            eg = curEdge->group();
            vsng = dag->getGroupFromNode(s);
            vtng = dag->getGroupFromNode(vt);
            /*
            if(vsng == 901 && curEdge->group() == 21355) {
                cerr << "\tedge_wt = " << wt << " edge_gro = " << curEdge->group()  << " edge number = " << vexEdge->size() << endl;
            }
            */
			for (auto nl : T[s]) {
                /*
                if(vsng == 901 && debug_count < 4) {
                    cerr << "vsng = " << vsng << " vt = " << vt << " mdd node = " << nl.first << endl;
                    debug_count++;
                }
                */
				step++;
				n = nl.first;
				int l = nl.second;
                //cerr << "\tfollowMTDD start" << endl;
                //np = followMTDD_F_G(n, vsng, vtng, eg);
                
                if(notDomainReduction) {
                    np = followMTDD_noDReduce(n, vsng, vtng, eg);
                } else {
                    np = followMTDD_F_G(n, vsng, vtng, eg);
                }
                
                //cerr << "\tfollowMTDD end" << endl;
				// if n = nonT continue
				if(np->level == 0) {
                    continue;
                }
				// if T[vt][n'] < l + wt
				if (T[vt].count(np) != 0) {
					// find other path
					if ((isShort && T[vt][np] > l + wt)||
							(!isShort && T[vt][np] < l + wt)) {
						T[vt][np] = l + wt;
						B[vt][np] = make_pair(curEdge, n);
					}
				} else {
					// haven't be traveled
					T[vt][np] = l + wt;
					B[vt][np] = make_pair(curEdge, n);
				}
			}
		}
	}
	int res =  -1;
    fMddNode* termianlNode = nullptr;
    for(auto nl : T[numOfvex-1]) {
        termianlNode = nl.first;
        if(termianlNode->ng == -1) {
            res = nl.second;
            break;
        }
    }
    if(res != -1) {
        stack<int> pathStack;
        stack<int> egStack;
        auto p = B[numOfvex - 1][termianlNode];
        while (p.first->source() != 0) {
            egStack.push(p.first->group());
            p = B[p.first->source()][p.second];
        }
        egStack.push(p.first->group());
        /*
        cout << "\tprint path:"; 
        while (!egStack.empty()) {
            cout << egStack.top() << " ";
            egStack.pop();
        }
        */
    } 
    gettimeofday(&t2,NULL);
    double timeUse = (t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0)*1000;

	string pathType = isShort ? "Short" : "Long";
    if(isOutput) {
        cout << "[MTCS FAST " << pathType << "est" << " Path's Length is:";
        if(res != -1) {
            cout << res;
        } else {
            cout << "UNREACHED";
        }
        cout << "]\n";
        cout << "\tmdd type:"  << (notDomainReduction ? "orignal" : "domain reduction") << "\n"; 
        cout << "\tstep number is:" << step << "\n";
        cout << fixed << setprecision(3) <<"\tmtcs fast time is:" << timeUse << "ms" << endl;
        if(slideCount != 0)
            cout << "\tmdd slide number is:" << slideCount << "\n";
    }
    // update info
    numOfsteps = step;
    searchTime = timeUse;
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
    if(res != -1) {
        stack<int> pathStack;
        stack<int> egStack;
        auto p = B[numOfvex - 1][termianlNode];
        cout << "\tprint path:";
        pathStack.push(p.first->target());
        //cout << " " << p.first->target() << " ";
        while (p.first->source() != 0) {
            pathStack.push(p.first->source());
            egStack.push(p.first->group());
            //cout << p.first->source() << " ";
            p = B[p.first->source()][p.second];
        }
        pathStack.push(0);
        egStack.push(p.first->group());
        while (!pathStack.empty()) {
            cout << pathStack.top() << " ";
            pathStack.pop();
        }
        cout << "\n\tprint edge group:";
        while (!egStack.empty()) {
            cout << egStack.top() << " ";
            egStack.pop();
        }
        cout << "\n";
    } 
#endif
    if(!isDetail) return res;
	// print the table of T[v][n]
	cout << "\tprint the T table:\n";
	for (int s = 0; s < numOfvex; s++){
        if(T[s].size() == 0) continue;
		cout << "\t#\n";
		for(auto tp : T[s]){
			cout << "\tT["<< s << "][" << tp.first << "] = " << tp.second << "\n";
		}
	}
    
	cout << endl;
	return res;
}

template <typename dagType>
inline fMddNode* MtcsAlgorithm<dagType>::createFMddNode(int ng, int level) {
    fMddNode* newNode = new fMddNode();
    newNode->ng = ng;
    newNode->level = level;
    newNode->other = nullptr;
    return newNode;
}

template <typename dagType>
double MtcsAlgorithm<dagType>::updateFastMdd() {
    double size = 0;
    destoryFastMdd();
    //unordered_map<node_handle, fMddNode*> visited;
    fmddTnode = createFMddNode(-1, -1);
    fmddFnode = createFMddNode(0, 0);
    fmddVisited.emplace(-1, fmddTnode);
    fmddVisited.emplace(0, fmddFnode);
    stack<node_handle> nodeSt;
    
    int consNodeHandel = consRoot->getNode();
    int consNodeNg = var2ng[f->getNodeLevel(consNodeHandel) - 1];
    nodeSt.push(consNodeHandel);
    fconsRoot = createFMddNode(consNodeNg, f->getNodeLevel(consNodeHandel));
    fmddVisited.emplace(consNodeHandel, fconsRoot);
    while(!nodeSt.empty()) {
        fMddNode* curFMddNode;
        int curNodeHandel = nodeSt.top();
        nodeSt.pop();
        auto fMddInfo = fmddVisited.find(curNodeHandel);
        if(fMddInfo == fmddVisited.end()) {
            curFMddNode = createFMddNode(var2ng[f->getNodeLevel(curNodeHandel) - 1],  f->getNodeLevel(curNodeHandel));
            fmddVisited.emplace(curNodeHandel, curFMddNode);
        } else {
            curFMddNode = fMddInfo->second;
        }
        unpacked_node* curNode = f->newUnpacked(curNodeHandel, FULL_ONLY);
        fMddNode* nextFMddNode;
        for(int i = 0; i < curNode->getNNZs(); i++) {
            size++;
            node_handle nextNodeHandel = curNode->d(i);
            auto fMddInfo2 = fmddVisited.find(nextNodeHandel);
            if(fMddInfo2 == fmddVisited.end()) {
                nextFMddNode = createFMddNode(var2ng[f->getNodeLevel(nextNodeHandel) - 1], f->getNodeLevel(nextNodeHandel) );
                fmddVisited.emplace(nextNodeHandel, nextFMddNode);
                nodeSt.push(nextNodeHandel);
            } else {
                nextFMddNode = fMddInfo2->second;
            }
            curFMddNode->children.push_back(nextFMddNode);
        }
        curFMddNode->other = curFMddNode->children[0];
        //curFMddNode->other = curFMddNode->children.back();
#if GENERAL
        //curFMddNode->children.erase(curFMddNode->children.begin());
        //curFMddNode->children.pop_back();
#endif
    }
    return size;
}

template <typename dagType>
void MtcsAlgorithm<dagType>::destoryFastMdd() {
    for (auto& pair : fmddVisited) {
    // pair 是 std::pair<MEDDLY::node_handle, fMddNode*>
    // pair.second 是 fMddNode* 类型的指针
        delete pair.second;  // 释放指针指向的内存
    }
    fmddVisited.clear(); // 清空映射
}

template <typename dagType>
int MtcsAlgorithm<dagType>::getConsTreeWidth() {
    int levelNum = var2ng.size();
    int width = 0;
    priority_queue<fMddNode*, vector<fMddNode*>, sortMddNode> nodeQue;

    unordered_set<fMddNode*> nodeSet;
	nodeQue.push(fconsRoot);
	nodeSet.emplace(fconsRoot);  

    for(int i = levelNum; i > 0; i--) {
		int levelWidth;
		/*
		priority_queue<fMddNode*, vector<fMddNode*>, sortMddNode> nodeQue_copy = nodeQue;
		while(!nodeQue_copy.empty()) {
			cout << nodeQue_copy.top()->level << " ";
			nodeQue_copy.pop();
		}
		cout << endl;
        */
		while(nodeQue.top()->level == i) {
			//cerr << "\ttop level = " << nodeQue.top()->level << endl;
			fMddNode* Pnode = nodeQue.top();
			nodeQue.pop();
            for(auto child : Pnode->children) {
                if(nodeSet.find(child) == nodeSet.end()) {
                    nodeSet.emplace(child);
                    nodeQue.push(child);
			    }
            }
            if(nodeSet.find(Pnode->other) == nodeSet.end()) {
                    nodeSet.emplace(Pnode->other);
                    nodeQue.push(Pnode->other);
			}
		}
        // here
		levelWidth = nodeQue.size();
		//cerr << "varlevel = " << i << " ; level width = " << levelWidth << endl;
		width = max(width, levelWidth);
	}
    if(isOutput)
	    cout << "[GET WIDTH OF CONS MDD]:" << width << endl;
	
	return width;
}

template <typename dagType>
vector<int>* MtcsAlgorithm<dagType>::getEgPath() {
    return &egPath;
}

//int &var, int &nodes, int &width, int& states, int &steps, long &traArcs, double &ddtime, double &searchtime
template <typename dagType>
void MtcsAlgorithm<dagType>::updateExpInfo(int &var, int &nodes, int &width, int& states, int &steps, long &traArcs, double &ddtime, double &searchtime, int &ddarcs) {
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
    ddarcs = numOfarcs;
}

#endif