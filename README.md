## Introduction
This is the source code and datasets of the paper "A Multi-Valued Decision Diagram-Based Approach to Constrained Optimal Path Problems over Directed Acyclic Graphs"

and Mingwei's graduation thesis.

## Code author
Mingwei Zhang - 13337649640@163.com

## Build
System requirements: `Ubuntu`

Lib requirements: have installed cudd-3.0.0 and meddly-master

Warning: You have to install meddly-master using the directory 'meddly-master', Otherwise the Makefile will meet a mistake.
	 If you cannot install meddly-master, contact the author, and I will give you a installed version of meddly-master.

open terminal in this directory

`make` - build main program

`make clean` - clean all * .o and the executable file



## Run code

`./main <DAGDir> <DAGFile> <constraint> <isGrouped> <isShortest>`

`DAGDir` - directed acyclic graph file directory

`DAGFile` - directed acyclic graph file name

`constraint` - constraint file (.cnf file or .dnf file) 

`isGrouped` - It Indicates whether the vertices of the DAG are grouped (0 or 1).  It should be set to 1 when you test the datasets Knapsack and Viterbi path.

`isShortest` - It Indicates It tells you whether we want to find the shortest path or the longest path (0 or 1). 



For example,  if we want to test the DAG `./data/experimentData/expPosDag` with `./data/experimentData/expPosCHK2.cnf` constraint:

```
./main ./data/experimentData/ expPosDag ./data/experimentData/expPosCHK2.cnf 1 0
```


The results of BCS, MCS without domain reduction, and MCS with domain reduction are saved in the file `result_bcs.csv`, `result_mcs.csv`, and `result_remcs`,  respectively.



## Datasets

All data sets in the ijcai-24 paper are saved in ```./data/experimentData/```, including  4 DAG files and 24 constraint files.

All data sets in the graduation thesis are saved in ```./data/experimentDataGraduate/```, including  4 DAG files and 64 constraint files.

The characters in the file name represent the following:

`Cit` - the shortest path problem over DAGs  (DAG)

`Knap` - the 0-1 knapsack problem (Knapsack)

`ED` -  the edit distance problem （Edit distance）

`Pos` - the Viberti path problem （Viterbi path）

`CHK` - Check constraint 

`DIS` - Dis constraint
