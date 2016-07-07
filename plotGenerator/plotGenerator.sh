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
                rm "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"/vertices.txt";
                rm "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"/edges.txt";
                rm "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"/solutionVertices.txt";
                rm "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"/solutionEdges.txt";
                rm "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"/plot.eps";
                ./plotGenerator $path $n ${d//.} $k $t;
            done
        done
    done
done
gnuplot generator.plt;
make clean;
