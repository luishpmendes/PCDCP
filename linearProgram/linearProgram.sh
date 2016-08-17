#!/bin/bash
rm -rf output;
mkdir -p output;
make linearProgram;
for n in 10 20 50 100
do
    for d in 0.3 0.5 0.7
    do
        for k in 0 5 10 20
        do
            for t in 0 1
            do
                for r in 0 1 2
                do
                    for p in 0.1 0.5 1.0
                    do
                        echo "N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.};
                        mkdir -p "output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.};
                        ./linearProgram < "../input/instanceN"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}".in" > "output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/result.out";
                    done
                done
            done
        done
    done
done
make clean;
