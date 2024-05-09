#ifndef DAG_H
#define DAG_H

#include <vector>
#include <iostream>
#include "edge.h"
#include "define.h"
using namespace std;

template <typename vertexType>
class DAG {

private: 
    int numVertex;
    int numEdge;
    int numEdgeGro;
    vector<vertexType> vertex;

public:
    vector<vector<Edge>> edgeList;
    vector<int> node2ng;
    vector<Edge> eg2edge;
    vector<int> ng2node;
    vector<unordered_set<int>> ng2itsEg;

public:
    DAG(int n, int eg = 400);

    // s-source t-target wt-weight gro-group
    void setEdge(int s, int t, int wt = 0, int gro = -1);

    // for debug
    bool isEdge(int s, int t);

    /*
        display the DAG graph
        vexPrint : whether the vex can be print directly
    */
    void show(bool vexPrint = false);

    void updateGroupInfo();

    inline void setVertex(vertexType vex, int v);

    inline vector<Edge>* getEdgeList(int v);

    inline Edge getEdgeFromGro(int edgeGro);

    inline void setNodeGro(int nodeIndex, int nodeGro);

    inline int getVertexNum() const;

    inline int getEdgeNum() const;

    inline int getEdgeGroupNum() const;

    inline int getGroupFromNode(int nodeIndex) const;

    inline int getNodeFromGroup(int ng) const;

    inline vertexType getVertex(int v) const;
};

/*
	DAG method
    for the use of template, all methods is implented here
*/
template <typename vertexType>
DAG<vertexType>::DAG(int n, int numEg) {
    numVertex = n;
    numEdge = 0;
    vertex = vector<vertexType>(n);
    edgeList = vector<vector<Edge>>(n, vector<Edge>());
    node2ng = vector<int>(n);
    eg2edge = vector<Edge>(numEg+1, Edge());
    ng2node = vector<int>(n+1, -1);
}

template <typename vertexType>
void DAG<vertexType>::setEdge(int s, int t, int wt, int gro) {
#if DEBUG
    if(isEdge(s, t)) {
        cerr << "[DAG CLASS WARNING]: Set a edge repeatedly!" << endl;
        return;
    }
#endif
    Edge e(s, t, wt, gro);
    edgeList[s].push_back(e);
    numEdge++;
    if(eg2edge[gro].source() == -1) 
        eg2edge[gro] = e;
#if 0
    cout << "is seting edge gro " << gro << "\n";
    for(int i = 0; i < 11; i++) {
		Edge e = this->getEdgeFromGro(i);
		cout << e.group() << " ";
	}
	cout << endl;
#endif
}

template <typename vertexType>
bool DAG<vertexType>::isEdge(int s, int t) {
    for(Edge e : edgeList[s]) {
        if(e.target() == t)
            return true;
    }
    return false;
}

template <typename vertexType>
void DAG<vertexType>::show(bool vexPrint) {
    cout << "\n[DAG SHOW vex num:" << numVertex << "; edge num:"
        << numEdge << "]\n";
    for(int i = 0; i < numVertex; i++) {
        cout << "\tvex index :" << i << ";";
        if(vexPrint) {
            cout << " vex:" << vertex[i] << ";";
        }
        cout << " out-edge num: " << edgeList[i].size() << "\n";
        for(Edge e : edgeList[i]) {
            cout << "\t  edge:" << e.source() << " " << e.target() << " "
                << e.weight() << " " << e.group() << "\n";
        }
    }
    cout << endl;
}

template<typename vertexType>
inline void DAG<vertexType>::setVertex(vertexType vex, int v) {
    vertex[v] = vex;
}

template <typename vertexType>
inline vector<Edge>* DAG<vertexType>::getEdgeList(int v) {
    return &edgeList[v];
}

template <typename vertexType>
inline int DAG<vertexType>::getVertexNum() const {
    return numVertex;
}

template <typename vertexType>
inline int DAG<vertexType>::getEdgeNum() const {
    return numEdge;
}

template <typename vertexType>
inline int DAG<vertexType>::getEdgeGroupNum() const {
    return eg2edge.size()-1;
}

template <typename vertexType>
inline vertexType DAG<vertexType>::getVertex(int v) const {
    return vertex[v];
}

template <typename vertexType>
inline void DAG<vertexType>::setNodeGro(int index, int gro) {
    node2ng[index] = gro;
    if(ng2node[gro] == -1 || this->getEdgeList(index)->size() > this->getEdgeList(ng2node[gro])->size())
        ng2node[gro] = index;
}

template <typename vertexType>
inline Edge DAG<vertexType>::getEdgeFromGro(int edgeGro) {
    return eg2edge[edgeGro];
}

template <typename vertexType>
inline int DAG<vertexType>::getGroupFromNode(int n) const {
    return node2ng[n];
}

template <typename vertexType>
inline int DAG<vertexType>::getNodeFromGroup(int n) const {
    return ng2node[n];
}

template <typename vertexType>
void DAG<vertexType>::updateGroupInfo() {
    unordered_set<int> ngSet;
    for(int i = 0; i < node2ng.size(); i++) {
        ngSet.emplace(node2ng[i]);
    }
    ng2itsEg = vector<unordered_set<int>>(ngSet.size());
    for(int i = 0; i < numVertex; i++) {
        int ng = node2ng[i];
        for(int j = 0; j < edgeList[i].size(); j++) {
            Edge* e = &edgeList[i][j];
            ng2itsEg[ng].emplace(e->group());
        }
    }
}

#endif