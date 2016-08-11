#!/bin/bash
make plotGenerator;
path="linearProgram";
for n in 10 20 50 100
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
                        echo "N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.};
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/vertices.txt";
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/edges.txt";
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/solutionVertices.txt";
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/solutionEdges.txt";
                        rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"R"$r"P"${p//.}"/plot.eps";
                        ./plotGenerator $path $n ${d//.} $k $t $r;
                    done
                done
            done
        done
    done
done
gnuplot generator.plt;
make clean;
