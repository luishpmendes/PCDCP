INC = /opt/gurobi652/linux64/include/
CPP = g++
CARGS = -Wall -m64 -O3 -fno-inline -std=c++11
CPPLIB = -L/opt/gurobi652/linux64/lib/ -lgurobi_c++ -lgurobi65

clean: 
	rm -rf phase1 phase2

phase1: phase1.cpp
	$(CPP) $(CARGS) -o phase1 phase1.cpp -I$(INC) $(CPPLIB) -lpthread -lm

phase2: phase2.cpp
	$(CPP) $(CARGS) -o phase2 phase2.cpp -I$(INC) $(CPPLIB) -lpthread -lm
