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

#ifndef MEDDLY_EVMDD_H
#define MEDDLY_EVMDD_H

#include "ev.h"
#include "../oper_binary.h"

namespace MEDDLY {
  class evmdd_forest;
};

class MEDDLY::evmdd_forest : public ev_forest {
  public:
    evmdd_forest(unsigned dsl, domain* d, range_type t, edge_labeling ev,
      const policies &p, int* level_reduction_rule=NULL);

    virtual void swapAdjacentVariables(int level);
    virtual void moveDownVariable(int high, int low);
    virtual void moveUpVariable(int low, int high);

  protected:
    template <class OPERATION, typename TYPE>
    inline void evaluateT(const dd_edge &f, const int* vlist, TYPE &val) const
    {
      if (f.getForest() != this) throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
      if (vlist == 0) throw error(error::INVALID_VARIABLE, __FILE__, __LINE__);

      // assumption: vlist does not contain any special values (-1, -2, etc).
      // vlist contains a single element.
      node_handle node = f.getNode();

      f.getEdgeValue(val);
      while (!isTerminalNode(node)) {
        TYPE ev;
        getDownPtr(node, vlist[getVarByLevel(getNodeLevel(node))], ev, node);
        val = (node) ? OPERATION::apply(val, ev) : ev;
      }
    }

  public:
    /// Special case for createEdge(), with only one minterm.
    template <class OPERATION, typename TYPE>
    inline void
    createEdgePath(int k, const int* vlist, TYPE &ev, node_handle &ed)
    {
        if (0==ed) return;
        for (int i=1; i<=k; i++) {
          if (DONT_CARE == vlist[i]) {
            // make a redundant node
            if (isFullyReduced()) continue;
            int sz = getLevelSize(i);
            unpacked_node* nb = unpacked_node::newFull(this, i, sz);
            nb->d_ref(0) = ed;
            nb->setEdge(0, ev);
            for (int v=1; v<sz; v++) {
              nb->d_ref(v) = linkNode(ed);
              nb->setEdge(v, ev);
            }
            createReducedNode(-1, nb, ev, ed);
          } else {
            // make a singleton node
            unpacked_node* nb = unpacked_node::newSparse(this, i, 1);
            nb->i_ref(0) = vlist[i];
            nb->d_ref(0) = ed;
            nb->setEdge(0, ev);
            createReducedNode(-1, nb, ev, ed);
          }
        } // for i
    }

};


//
// Helper class for createEdge
//

namespace MEDDLY {

  template <class OPERATION, typename T>
  class evmdd_edgemaker {
      evmdd_forest* F;
      const int* const* vlist;
      const T* values;
      int* order;
      int N;
      int K;
      binary_operation* unionOp;
    public:
      evmdd_edgemaker(evmdd_forest* f, const int* const* mt, const T* v,
        int* o, int n, int k, binary_operation* unOp)
      {
        F = f;
        vlist = mt;
        values = v;
        order = o;
        N = n;
        K = k;
        unionOp = unOp;
      }

      inline const int* unprimed(int i) const {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, N);
        return vlist[order[i]];
      }
      inline int unprimed(int i, int k) const {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, N);
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1, k, K+1);
        return vlist[order[i]][k];
      }
      inline T term(int i) const {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, N);
        return values ? values[order[i]]: 1;
      }
      inline void swap(int i, int j) {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, N);
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, j, N);
        MEDDLY::SWAP(order[i], order[j]);
      }

      inline void createEdge(T &ev, node_handle &ed) {
        createEdge(K, 0, N, ev, ed);
      }

      /**
          Recursive implementation of createEdge(),
          for use by evmdd_forest descendants.
      */
      void createEdge(int k, int start, int stop, T &ev, node_handle &ed) {
        MEDDLY_DCASSERT(k>=0);
        MEDDLY_DCASSERT(stop > start);
        //
        // Fast special case
        //
        if (1==stop-start) {
          ev = term(start);
          ed = bool_Tencoder::value2handle(true);
          F->createEdgePath<OPERATION, T>(k, unprimed(start), ev, ed);
          return;
        }
        //
        // Check terminal case
        //
        if (0==k) {
          ev = term(start);
          for (int i=start+1; i<stop; i++) {
            OPERATION::unionEq(ev, term(i));
          }
          ed = bool_Tencoder::value2handle(true);
          return;
        }

        // size of variables at level k
        unsigned lastV = unsigned(F->getLevelSize(k));
        // index of end of current batch
        int batchP = start;

        //
        // Move any "don't cares" to the front, and process them
        //
        unsigned nextV = lastV;
        for (int i=start; i<stop; i++) {
          if (DONT_CARE == unprimed(i, k)) {
            if (batchP != i) {
              swap(batchP, i);
            }
            batchP++;
          } else {
            MEDDLY_DCASSERT(unprimed(i, k) >= 0);
            nextV = MIN(nextV, unsigned(unprimed(i, k)));
          }
        }
        dd_edge dontcare(F);
        if (batchP > start) {
          T dc_val;
          node_handle dc_ptr;
          createEdge(k-1, start, batchP, dc_val, dc_ptr);
          dontcare.set(dc_ptr, dc_val);
        } else {
          OPERATION::makeEmptyEdge(dontcare);
        }

        //
        // Start new node at level k
        //
        unpacked_node* nb = unpacked_node::newSparse(F, k, lastV);
        unsigned z = 0; // number of nonzero edges in our sparse node

        //
        // For each value v,
        //  (1) move those values to the front
        //  (2) process them, if any
        //  (3) union with don't cares
        //
        dd_edge these(F);

        for (unsigned v = (dontcare.getNode()) ? 0 : nextV;
             v<lastV;
             v = (dontcare.getNode()) ? v+1 : nextV)
        {
          nextV = lastV;
          //
          // neat trick!
          // shift the array over, because we're done with the previous batch
          //
          start = batchP;

          //
          // (1) move anything with value v, to the "new" front
          //
          for (int i=start; i<stop; i++) {
            if (v == unprimed(i, k)) {
              if (batchP != i) {
                swap(batchP, i);
              }
              batchP++;
            } else {
              nextV = MIN(nextV, unsigned(unprimed(i, k)));
            }
          } // for i

          //
          // (2) recurse if necessary
          //
          if (batchP > start) {
            T these_val;
            node_handle these_ptr;
            createEdge(k-1, start, batchP, these_val, these_ptr);
            these.set(these_ptr, these_val);
          } else {
            OPERATION::makeEmptyEdge(these);
          }

          //
          // (3) union with don't cares
          //
          if (dontcare.getNode()) {
            if (these.getNode()) {
              // Need to do a union
              if (0==unionOp) {
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
              }
              MEDDLY_DCASSERT(unionOp);
              unionOp->computeTemp(dontcare, these, these);
            } else {
              these = dontcare;
            }
          }

          //
          // add to sparse node, unless empty
          //
          if (0==these.getNode()) continue;
          nb->i_ref(z) = v;
          nb->d_ref(z) = F->linkNode(these);
          T temp;
          these.getEdgeValue(temp);
          nb->setEdge(z, temp);
          z++;
        } // for v

        //
        // Cleanup
        //
        // F->unlinkNode(dc_ptr);
        nb->shrinkSparse(z);

        F->createReducedNode(-1, nb, ev, ed);

      }; // method createEdge

  }; // class evmdd_edgemaker

}; // namespace MEDDLY

#endif
