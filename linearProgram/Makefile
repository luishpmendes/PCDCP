CPPC = g++
CFLAGS = -Wall
DEBUGFLAGS = -g -O0 -fno-inline

COMPILER = $(CPPC)
CFLAGS += $(DEBUGFLAGS)

clean: 
	rm -rf linearProgram

linearProgram: linearProgram.cpp
	$(COMPILER) $(CFLAGS) linearProgram.cpp -o linearProgram
