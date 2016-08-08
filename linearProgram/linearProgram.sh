#!/bin/bash
rm -rf output;
mkdir -p output;
make linearProgram;
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
                    mkdir -p "output/N"$n"D"${d//.}"K"$k"T"$t"R"$r;
                    ./linearProgram < "../input/instanceN"$n"D"${d//.}"K"$k"T"$t"R"$r".in" > "output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"/result.out";
                done
            done
        done
    done
done
make clean;
