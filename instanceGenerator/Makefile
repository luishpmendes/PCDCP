CPPC = g++
CFLAGS = -Wall
DEBUGFLAGS = -g -O0 -fno-inline -std=c++11

COMPILER = $(CPPC)
CFLAGS += $(DEBUGFLAGS)

clean: 
	rm -rf instanceGenerator

instanceGenerator: instanceGenerator.cpp
	$(COMPILER) $(CFLAGS) instanceGenerator.cpp -o instanceGenerator