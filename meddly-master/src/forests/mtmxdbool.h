/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2009, Iowa State University Research Foundation, Inc.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MEDDLY_MTMXDBOOL_H
#define MEDDLY_MTMXDBOOL_H

#include "mtmxd.h"

namespace MEDDLY {
  class mt_mxd_bool;
};

// ******************************************************************

/**
    Forest for multi-terminal, mxd, boolean range.
*/
class MEDDLY::mt_mxd_bool : public mtmxd_forest {
  public:

    mt_mxd_bool(unsigned dsl, domain *d, const policies &p, int* level_reduction_rule=NULL, bool tv=false);
    ~mt_mxd_bool();

    virtual void createEdge(bool val, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge& e);
    virtual void createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a);
    virtual void evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;

    virtual node_handle unionOneMinterm(node_handle a,  int* from,  int* to, int level);

    virtual void showTerminal(output &s, node_handle tnode) const;
    virtual void writeTerminal(output &s, node_handle tnode) const;
    virtual node_handle readTerminal(input &s);

  protected:
    virtual const char* codeChars() const;
    node_handle unionOneMinterm_r(int index, int vh, node_handle nh,  int* from,  int* to);
    bool checkTerminalMinterm(node_handle a,  int* from,  int* to,int levelRead, node_handle& c);

};

#endif

