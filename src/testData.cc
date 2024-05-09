#include "testData.h"

TestData::TestData(string dirN, string dagN, string chkN, string disN) {
    dirName = dirN;
    dagFileName = dagN;
    chkFileName = chkN;
	disFileName = disN;
    problemType = unKnow;
    safeModel = true;
}

void TestData::setDirName(string newName) {
	dirName = newName;
}

void TestData::setDagName(string newName) {
	dagFileName = newName;
}


void TestData::setChkName(string newName) {
    chkFileName = newName;
}

void TestData::setDisName(string newName) {
    disFileName = newName;
}

string TestData::getDagFileName() {
	return dirName + dagFileName;
}

string TestData::getChkFileName() {
	return dirName + chkFileName + ".cnf";
}

string TestData::getDisFileName() {
	return dirName + disFileName + ".cnf";
}

void TestData::allowGenerate() {
    safeModel = false;
}