#!/bin/bash
rm -rf ../input;
mkdir -p ../input;
make instanceGenerator;
for n in 50 100 250 500
do
    for d in 0.3 0.5 0.7
    do
        for k in 0 2 4 6 8 10
        do
            for t in 0 1
            do
                ./instanceGenerator $n $d $k $t > "../input/instanceN"$n"D"${d//.}"K"$k"T"$t".in";
            done
        done
    done
done
make clean;
