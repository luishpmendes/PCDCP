#!/bin/bash
make plotGenerator;
path="linearProgram";
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
                    rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"/vertices.txt";
                    rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"/edges.txt";
                    rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"/solutionVertices.txt";
                    rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"/solutionEdges.txt";
                    rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"/plot.eps";
                    ./plotGenerator $path $n ${d//.} $k $t $r;
                done
            done
        done
    done
done
gnuplot generator.plt;
make clean;
