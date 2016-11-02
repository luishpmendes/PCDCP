INC = /opt/gurobi652/linux64/include/
CPP = g++
CARGS = -Wall -m64 -O3 -fno-inline
CPPLIB = -L/opt/gurobi652/linux64/lib/ -lgurobi_c++ -lgurobi65

clean: 
	rm -rf grasp

grasp: grasp.cpp
	$(CPP) $(CARGS) -o grasp grasp.cpp -I$(INC) $(CPPLIB) -lpthread -lm
