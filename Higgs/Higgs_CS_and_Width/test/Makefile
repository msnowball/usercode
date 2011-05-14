CXX=g++
CC=gcc
CFLAGS=-O2 -Wall
LDFLAGS=-Llib
PRFFLAGS=-lProof
THRDFLAGS=-lThread
INS=-I$(ROOTSYS)/include/root
INSS=-I../../Higgs_CS_and_Width/include
INSSS=-I../../Higgs_CS_and_Width/src

CFLAGS += `root-config --cflags`
LIBS += `root-config --glibs`

OBJS=test.o
OBJSS=validationTest.o
OBJSSS=validationTest_pdf.o
OBJSSSS=validationTest_scale.o
OBJSSSSS=validationTest_BR.o

.PHONY: clean all main test

all: $(OBJS)

main: test.o
	$(CXX) -o test.exe $(OBJS) $(LIBS)

validate: validationTest.o
	$(CXX) -o validationTest.exe $(OBJSS) $(LIBS)

validatePdf: validationTest_pdf.o
	$(CXX) -o validationTest_pdf.exe $(OBJSSS) $(LIBS)

validateScale: validationTest_scale.o
	$(CXX) -o validationTest_scale.exe $(OBJSSSS) $(LIBS)

validateBR: validationTest_BR.o
	$(CXX) -o validationTest_BR.exe $(OBJSSSSS) $(LIBS)

test: all main
	./main

clean:
	@rm *.o *.exe *~ 

new:
	@rm *.root *.eps *.gif *~

##############RULES##############
.cc.o:
	$(CXX) $(CFLAGS) $(INSS) $(INSSS) -c $<
.cpp.o:
	$(CXX) $(CFLAGS) $(INS) -c $<


