#!/bin/bash
rm -rf ../input;
mkdir -p ../input;
make instanceGenerator;
for n in 10 20 50 100 200
do
    for d in 0.3 0.5 0.7
    do
        for k in 0 2 4 6 8 10
        do
            for t in 0 1
            do
                for r in 0 1 2
                do
                    for p in 0.25 0.5
                    do
                        ./instanceGenerator $n $d $k $t $r $p > "../input/instanceN"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}".in";
                    done
                done
            done
        done
    done
done
make clean;
