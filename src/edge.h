/*
	the class implement edge
	In BCS:
		give every edge a ser to express node (every BDD node express one edge)
		give every edge a group to express node (every BDD node express one group edge)
		gro is a int can express grop or ser
*/

#ifndef EDGE_H
#define EDGE_H
using namespace std;

class Edge {
		 
private:
	int sour, tar, wt, gro;	 

public:
	Edge();

	Edge(int s, int t, int w, int g = -1);

	int source() const;

	int target() const;

	int weight() const;

	int group() const;
};

/*
	inline Edge method
*/
inline int Edge::source() const {
	return sour;
}
inline int Edge::target() const {
	return tar;
}
inline int Edge::weight() const {
	return wt;
}
inline int Edge::group() const {
	return gro;
}

#endif
