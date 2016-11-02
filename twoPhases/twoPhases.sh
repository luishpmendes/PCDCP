#!/bin/bash
rm -rf output;
mkdir -p output;
make phase1;
make phase2;
timeLimit = 10.0;
for n in 10 20 50 100
do
    for d in 0.3 0.5 0.7
    do
        for k in 0 5 10 20
        do
            for t in 0 1
            do
                for p in 0.1 0.5 1.0
                do
                    echo "phase1 - N"$n"D"${d//.}"K"$k"T"$t"P"${p//.};
                    mkdir -p "output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.};
                    ./phase1 $timeLimit < "../input/instanceN"$n"D"${d//.}"K"$k"T"$t"P"${p//.}".in" > "output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"/result.tmp";
                done
            done
        done
    done
done
for n in 10 20 50 100
do
    for d in 0.3 0.5 0.7
    do
        for k in 0 5 10 20
        do
            for t in 0 1
            do
                for p in 0.1 0.5 1.0
                do
                    echo "phase2 - N"$n"D"${d//.}"K"$k"T"$t"P"${p//.};
                    mkdir -p "output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.};
                    ./phase2 $timeLimit < "output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"/result.tmp" > "output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"/result.out";
                done
            done
        done
    done
done
make clean;
