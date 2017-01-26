#! /usr/bin/gnuplot
set terminal pdf color font "times"
unset key
set size ratio -1 1,1
set xrange [0:1]
set yrange [0:1]
if (path eq "grasp") {
    set output "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."A".A."/graph.pdf"
} else {
    set output "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."/graph.pdf"
}
set title "Instancia do PCDCP (N = ".n."; D = ".d."; K = ".k."; T = ".t."; P = ".p.")"
if (path eq "grasp") {
    vertices = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."A".A."/vertices.txt"
} else {
    vertices = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."/vertices.txt"
}
if (path eq "grasp") {
    edges = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."A".A."/edges.txt"
} else {
    edges = "../".path."/output/N".n."D".D."K".k."T".t."P".P."I".i."/edges.txt"
}
plot vertices using ($1):($2):($3) with points pointtype 7 pointsize variable linecolor rgb "black", \
	 edges with lines linetype 1 linewidth 1 linecolor rgb "forest-green"
