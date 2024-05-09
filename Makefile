CC = g++
CUDD_DIR = ./cudd-3.0.0
MDD_DIR = ./meddly-master
INCLUDE = -I./src -I$(CUDD_DIR)/cudd -I$(MDD_DIR)/include/meddly
MDD_LIB = $(MDD_DIR)/lib
CUDD_LIB = $(CUDD_DIR)/lib
LIB = ./lib
LIBS = $(MDD_LIB)/libmeddly.a $(CUDD_LIB)/libcudd.a $(LIB)/libgmp.a
TARGET = main
SRCS = ./$(wildcard *.cc)
OBJS = ./build/$(patsubst %cc, %o, %(SRCS))
SRC_WHERE = ./src/
BUILD_WHERE = ./build/
OType = -O3

run:clean $(TARGET)

$(TARGET): $(SRC_WHERE)main.cc $(BUILD_WHERE)edge.o $(BUILD_WHERE)testData.o
	$(CC) $(OType) $(INCLUDE) $^ -o $@ $(LIBS)

$(BUILD_WHERE)edge.o: $(SRC_WHERE)edge.cc
	$(CC) $(OType) $(INCLUDE) -c $^ -o $@

$(BUILD_WHERE)testData.o: $(SRC_WHERE)testData.cc
	$(CC) $(OType) $(INCLUDE) -c $^ -o $@

clean:
	rm -f *.o $(BUILD_WHERE)*.o main