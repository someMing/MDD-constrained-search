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

#ifndef MEDDLY_UNIQUE_TABLE_H
#define MEDDLY_UNIQUE_TABLE_H

#include "defines.h"
#include "forest.h"

namespace MEDDLY {
    class unique_table;
};

/** Unique table for discovering duplicate nodes.

    This is now a stand-alone, non-template class
    designed specifically for expert_forests.
 */
class MEDDLY::unique_table {
    public:
        unique_table(expert_forest *ef);
        ~unique_table();

        /// Get the hash table size for a given variable
        inline unsigned getSize(int var) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, min_var, var, max_var+1);
            return tables[var].getSize();
        }

        /// Get the number of unique nodes for a given variable
        inline unsigned getNumEntries(int var) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, min_var, var, max_var+1);
            return tables[var].getNumEntries();
        }

        /// Get the memory used for a given variable's table
        inline unsigned getMemUsed(int var) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, min_var, var, max_var+1);
            return tables[var].getMemUsed();
        }

        /** If the table of the given variable contains key, move it to
            the front of the list.  Otherwise, do nothing.
            Returns the item if found, 0 otherwise.

            Class T must have the following methods:
                unsigned hash():    return the hash value for this item.
                bool equals(int p): return true iff this item equals node p.
        */
        template <typename T>
        inline node_handle find(const T &key, int var) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, min_var, var, max_var+1);
            return tables[var].find(key);
        }


        /** Add the item to the front of the list of the corresponding
            variable.  Used when we KNOW that the item is not in the
            unique table already.
        */
        inline void add(unsigned hash, node_handle item) {
            int level = parent->getNodeLevel(item);
            int var = parent->getVarByLevel(level);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, min_var, var, max_var+1);
            tables[var].add(hash, item);
        }


        /** If the table of the corresponding variable contains key,
            remove it and return it.  I.e., the exact key.
            Otherwise, return 0.
        */
        inline node_handle remove(unsigned hash, node_handle item) {
            int var = parent->getVarByLevel(parent->getNodeLevel(item));
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, min_var, var, max_var+1);
            return tables[var].remove(hash, item);
        }

        /**
            Clear the items in the table of the given variable and
            reset the state.
        */
        inline void clear(int var) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, min_var, var, max_var+1);
            return tables[var].clear();
        }


        /**
            Retrieve the items in the table of the given variable.
                @param  var     Variable of interest
                @param  NH      Array of node handles
                @param  sz      Dimension of NH array, and max items
                                to retrieve.

                @return         Number of items retrieved; will be
                                at most sz.
        */
        inline unsigned getItems(int var, node_handle* NH, unsigned sz) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, min_var, var, max_var+1);
            return tables[var].getItems(NH, sz);
        }



        /// Get the total size (sum over all variables)
        unsigned getSize() const;

        /// Get the total number of nodes (sum over all variables)
        unsigned getNumEntries() const;

        /// Get the total memory used (sum over all variables)
        unsigned getMemUsed() const;

        void reportStats(output &s, const char* pad, unsigned flags) const;

        /// For debugging
        void show(output &s) const;


    private:
        // ========================================================
        class subtable {
            public:
                subtable();
                ~subtable();

                inline unsigned getSize() const {
                    return size;
                }
                inline unsigned getNumEntries() const {
                    return num_entries;
                }
                inline unsigned getMemUsed() const {
                    return size * sizeof(node_handle);
                }

                void reportStats(output &s, const char* pad, unsigned flags)
                    const;

                /// For debugging
                void show(output& s) const;

                /**
                    Initialize the sub table. Must be called before use.
                */
                void init(expert_forest *ef);

                /** If table contains key, move it to the front of the list.
                    Otherwise, do nothing.
                    Returns the item if found, 0 otherwise.

                    Class T must have the following methods:
                        unsigned hash():    return the hash value for this item.
                        bool equals(int p): return true iff this item equals node p.
                */
                template <typename T> int find(const T &key) const
                {
                    unsigned h = key.hash() % size;
                    MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, h, size);
                    node_handle prev = 0;
                    for (node_handle ptr = table[h];
                            ptr != 0;
                            ptr = parent->getNext(ptr))
                    {
                        if (parent->areDuplicates(ptr, key)) {
                            // MATCH
                            if (ptr != table[h]) {
                                // Move to front
                                MEDDLY_DCASSERT(prev);
                                parent->setNext(prev, parent->getNext(ptr));
                                parent->setNext(ptr, table[h]);
                                table[h] = ptr;
                            }
                            MEDDLY_DCASSERT(table[h] == ptr);
                            return ptr;
                        }
                        prev = ptr;
                    }

                    return 0;
                }

                /** Add the item to the front of the list.
                    Used when we KNOW that the item is not
                    in the unique table already.
                */
                void add(unsigned hash, node_handle item);

                /** If table contains key, remove it and return it.
                    I.e., the exact key.
                    Otherwise, return 0.
                */
                int remove(unsigned hash, node_handle item);

                /**
                    Remove all the items in the table and reset the state.
                */
                void clear();

                /**
                    Retrieve the items in the table of the given variable.
                        @param  var     Variable of interest
                        @param  NH      Array of node handles
                        @param  sz      Dimension of NH array, and max items
                                        to retrieve.

                        @return         Number of items retrieved; will be
                                        at most sz.
                */
                unsigned getItems(node_handle* items, unsigned sz) const;

            private:  // helper methods
                /// Empty the hash table into a list; returns the list.
                node_handle convertToList();

                /// A series of inserts; doesn't check for duplicates or expand.
                void buildFromList(node_handle front);

                /// Expand the hash table (if possible)
                void expand();

                /// Shrink the hash table
                void shrink();

            private:
                expert_forest* parent;
                unsigned size;
                unsigned num_entries;
                unsigned next_expand;
                unsigned next_shrink;
                node_handle* table;

                static const unsigned MAX_SIZE = 1073741824;
                static const unsigned MIN_SIZE = 8;
        }; // class subtable
        // ========================================================


    private:
        expert_forest* parent;
        subtable* tables;
        int min_var, max_var;
};

#endif // #include guard
