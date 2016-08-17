#!/bin/bash
make solutionValidator;
path="linearProgram";
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
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/validation.txt";
                        ./solutionValidator $path $n ${d//.} $k $t $r ${p//.} > "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/validation.txt";
                    done
                done
            done
        done
    done
done
gnuplot generator.plt;
make clean;
