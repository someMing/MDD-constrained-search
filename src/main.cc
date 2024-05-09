#include "testData.h"
#include "bcsAlgorithm.h"
#include "mtcsAlgorithm.h"

using namespace std;

int main(int argc, char *argv[]) {

	string dirName;
	string dagName;
	string consName;

    if (argc != 6) {
        cerr << "Usage: " << argv[0] << " <DAGDir> <DAGFile> <constraint> <isGrouped> <isShortest>" << endl;
        return 1;
    }

	dirName = argv[1];
    dagName = argv[2];
    consName = argv[3];
    int isGrouped = stoi(argv[4]);
    int isShortestPath = stoi(argv[5]);

	initialize();

	TestData td;
	td.setDirName(dirName);
	td.setDagName(dagName);

	DAG<string>* dag;
	if (dagName.find("Pos") != string::npos) {
        dag = td.readPosDag<string>();
    } else {
        dag = td.readDag<string>();
    }

	if(isGrouped) {
		dag->updateGroupInfo();
	}

	bool isOutput = false;

	int bcs_vars;
    int bcs_nodes;
	int bcs_arcs;
    int bcs_width;
	int bcs_states;
	int bcs_steps;
	long bcs_traArcs;
	double bcs_ddTime;
	double bcs_searchTime;

	int mcs_vars;
    int mcs_nodes;
	int mcs_arcs;
    int mcs_width;
	int mcs_states;
	int mcs_steps;
	long mcs_traArcs;
	double mcs_ddTime;
	double mcs_searchTime;

	int Remcs_vars;
    int Remcs_nodes;
	int Remcs_arcs;
    int Remcs_width;
	int Remcs_states;
	int Remcs_steps;
	long Remcs_traArcs;
	double Remcs_ddTime;
	double Remcs_searchTime;


	BcsAlgorithm<string> bcs(dag, consName, isShortestPath, isOutput);
	MtcsAlgorithm<string> mcs(dag, consName, isShortestPath, isOutput, isGrouped, true);
	MtcsAlgorithm<string> re_mcs(dag, consName, isShortestPath, isOutput, isGrouped);

	//int &var, int &nodes, int &width, int& states, int &steps, long &traArcs, double &ddtime, double &searchtime
	bcs.setCons(consName);
	bcs.BCS();
	bcs.updateExpInfo(bcs_vars, bcs_nodes, bcs_width, bcs_states, bcs_steps, bcs_traArcs, bcs_ddTime, bcs_searchTime, bcs_arcs);

	mcs.setCons(consName);
	mcs.MTCS_F();
	mcs.updateExpInfo(mcs_vars, mcs_nodes, mcs_width, mcs_states, mcs_steps, mcs_traArcs, mcs_ddTime, mcs_searchTime, mcs_arcs);

	re_mcs.setCons(consName);
	re_mcs.MTCS_F();
	re_mcs.updateExpInfo(Remcs_vars, Remcs_nodes, Remcs_width, Remcs_states, Remcs_steps, Remcs_traArcs, Remcs_ddTime, Remcs_searchTime, Remcs_arcs);

	cout << "\n[DAG INFO] " << (isShortestPath ? "Shortest Path " : "Longest Path ") << "and " << (isGrouped ? "Grouped" : "Not Grouped") << "\n";

	// output
	ofstream file_bcs("result_bcs.csv", ofstream::trunc);
	if (file_bcs.is_open()) {
		file_bcs << "DAG," << "Constraint," << "Vars," << "Nodes," << "Arcs," << "Width," << "Entries," << "Steps," << "Traverse arcs,"
				<< "Compilation time (ms)," << "Search time (ms)," << "Total time (ms)," << endl;

        file_bcs << dagName << "," << consName << "," << bcs_vars << ',' << bcs_nodes << ',' << bcs_arcs << ',' << bcs_width << ',' << bcs_states << ',' << bcs_steps << ',' << bcs_traArcs << ',' 
			<< bcs_ddTime << ',' << bcs_searchTime << ',' << bcs_ddTime + bcs_searchTime << ',';
        file_bcs.close();
        cout << "[Data] bcs result written successfully." << endl;
    } else {
        cout << "[ERROR]Unable to open the file res_bcs.csv." << endl;
    }

	ofstream file_mcs("result_mcs.csv", ofstream::trunc);
	if (file_mcs.is_open()) {
		file_mcs << "DAG," << "Constraint," << "Vars," << "Nodes," << "Arcs," << "Width," << "Entries," << "Steps," << "Traverse arcs,"
				<< "Compilation time (ms)," << "Search time (ms)," << "Total time (ms)," << endl;

        file_mcs << dagName << "," << consName << "," << mcs_vars << ',' << mcs_nodes << ',' << mcs_arcs << ',' << mcs_width << ',' << mcs_states << ',' << mcs_steps << ',' 
			<< mcs_traArcs << ',' << mcs_ddTime << ',' << mcs_searchTime << ',' << mcs_ddTime + mcs_searchTime << ',';
        file_mcs.close();
        cout << "[Data] mcs without domain reduction result written successfully." << endl;
    } else {
        cout << "[ERROR]Unable to open the file res_mcs.csv." << endl;
    }

	ofstream file_remcs("result_remcs.csv", ofstream::trunc);
	if (file_remcs.is_open()) {
		file_remcs << "DAG," << "Constraint," << "Vars," << "Nodes," << "Arcs," << "Width," << "Entries," << "Steps," << "Traverse arcs,"
				<< "Compilation time (ms)," << "Search time (ms)," << "Total time (ms)," << endl;

        file_remcs << dagName << "," << consName << "," << Remcs_vars << ',' << Remcs_nodes << ',' << Remcs_arcs << ',' << Remcs_width << ',' << Remcs_states << ',' << Remcs_steps << ',' 
			<< Remcs_traArcs << ',' << Remcs_ddTime << ',' << Remcs_searchTime << ',' << Remcs_ddTime + Remcs_searchTime << ',';
        file_remcs.close();
        cout << "[Data] mcs with domain reduction result written successfully." << endl;
    } else {
        cout << "[ERROR]Unable to open the file res_remcs.csv." << endl;
    }
	
    return 0;
}