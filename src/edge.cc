#include "edge.h"

Edge::Edge() {
    sour = tar = wt = gro = -1;
}

Edge::Edge(int s, int t, int w, int g) {
    sour = s;
    tar = t;
    wt = w;
    gro = g;
}