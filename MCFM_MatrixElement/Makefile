CXX=g++
CC=gcc
CFLAGS=-O2 -Wall
LDFLAGS=-Llib
PRFFLAGS=-lProof
THRDFLAGS=-lThread
INS=-I$(ROOTSYS)/include/root
INS2=-I$(ROOFITSYS)/include
INSS=-I./include

LD1=-L$(ROOFITSYS)/lib

CFLAGS += `root-config --cflags`
LIBS += `root-config --glibs`

LDa=-lRooFitCore
LDb=-lRooFit

OBJ1=nlome_root.o
OBJ2=makeTxtFile.o
OBJ3=test_hzz_gridZZ.o
OBJ4=test_hzz_gridH.o
OBJ5=splitTrees.o

.PHONY: clean all main test

all: nlome splitTrees

nlome: nlome_root.o
	$(CXX) -o nlome_root.exe $(OBJ1) -lgfortran -lNLOME $(LIBS) -LNLOME 

test2: test_hzz_cpp.o
	$(CXX) -o test_hzz.exe test_hzz_cpp.o -lgfortran -LNLOME -lNLOME

txt: makeTxtFile.o
	$(CXX) -o makeTxtFile.exe $(OBJ2) $(LIBS)

gridZZ: test_hzz_gridZZ.o
	$(CXX) -o test_hzz_gridZZ.exe $(OBJ3) -lgfortran -LNLOME -lNLOME

gridH: test_hzz_gridH.o
	$(CXX) -o test_hzz_gridH.exe $(OBJ4) -lgfortran -LNLOME -lNLOME

splitTrees: splitTrees.o
	$(CXX) -o splitTrees.exe $(OBJ5) $(LIBS)


clean:
	@rm *.o test_hzz *.exe

##############RULES##############                                                                                                                           
.cc.o:
	$(CXX) $(CFLAGS) $(INS) -c $<
.cpp.o:
	$(CXX) $(CFLAGS) $(INS) -c $<
