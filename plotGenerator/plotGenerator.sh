#!/bin/bash
make plotGenerator;
for path in "linearProgram" "twoPhases"
do
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
                            echo $path" - N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i;
                            rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/vertices.txt";
                            rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/edges.txt";
                            rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/solutionVertices.txt";
                            rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/solutionEdges.txt";
                            rm -f "../"$path"/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"/result.pdf";
                            ./plotGenerator $path $n ${d//.} $k $t ${p//.} $i;
                            gnuplot -e "path = '${path}'; n='${n}'; d='${d}'; k='${k}'; t='${t}'; p='${p}'; D='${d//.}'; P='${p//.}'; i='$i';" generator.plt;
                        done
                    done
                done
            done
        done
    done
done
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
                        for a in 0.3 0.5 0.7
                        do
                            echo "grasp - N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${a//.};
                            rm -f "../grasp/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${a//.}"/vertices.txt";
                            rm -f "../grasp/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${a//.}"/edges.txt";
                            rm -f "../grasp/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${a//.}"/solutionVertices.txt";
                            rm -f "../grasp/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${a//.}"/solutionEdges.txt";
                            rm -f "../grasp/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"A"${a//.}"/result.pdf";
                            ./plotGenerator "grasp" $n ${d//.} $k $t ${p//.} $i ${a//.};
                            gnuplot -e "path = 'grasp'; n='${n}'; d='${d}'; k='${k}'; t='${t}'; p='${p}'; D='${d//.}'; P='${p//.}'; i='$i'; a='$a'; A='${a//.}';" generator.plt;
                        done
                    done
                done
            done
        done
    done
done
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
                        for ps in 10 50 100
                        do
                            for mr in 0.1 0.2 0.3
                            do
                                echo "geneticAlgorithm - N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"PS"$ps"MR"${mr//.};
                                rm -f "../geneticAlgorithm/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"PS"$ps"MR"${mr//.}"/vertices.txt";
                                rm -f "../geneticAlgorithm/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"PS"$ps"MR"${mr//.}"/edges.txt";
                                rm -f "../geneticAlgorithm/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"PS"$ps"MR"${mr//.}"/solutionVertices.txt";
                                rm -f "../geneticAlgorithm/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"PS"$ps"MR"${mr//.}"/solutionEdges.txt";
                                rm -f "../geneticAlgorithm/output/N"$n"D"${d//.}"K"$k"T"$t"P"${p//.}"I"$i"PS"$ps"MR"${mr//.}"/result.pdf";
                                ./plotGenerator "geneticAlgorithm" $n ${d//.} $k $t ${p//.} $i $ps ${mr//.};
                                gnuplot -e "path = 'geneticAlgorithm'; n='${n}'; d='${d}'; k='${k}'; t='${t}'; p='${p}'; D='${d//.}'; P='${p//.}'; i='$i'; ps='$ps'; MR='${mr//.}';" generator.plt;
                            done
                        done
                    done
                done
            done
        done
    done
done
make clean;
