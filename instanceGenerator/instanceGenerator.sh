#!/bin/bash
rm -rf ../input;
mkdir -p ../input;
make instanceGenerator;
for n in 10 20 50 100
do
    for d in 0.3 0.5 0.7
    do
        for k in 0 10 20
        do
            for t in 0 1
            do
                for p in 0.2
                do
                    for i in 0
                    do
                        echo "N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i;
                        ./instanceGenerator $n $d $k $t $p > "../input/instanceN"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i".in";
                    done
                done
            done
        done
    done
done
make clean;
