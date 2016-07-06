#!/bin/bash
make plotGenerator;
path="linearProgram";
for n in 50 100 250 500
do
    for d in 0.3 0.5 0.7
    do
        for k in 0 2 4 6 8 10
        do
            for t in 0 1
            do
                ./plotGenerator $path $n ${d//.} $k $t;
            done
        done
    done
done
make clean;