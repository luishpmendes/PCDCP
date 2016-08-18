#! /usr/bin/gnuplot
set terminal pdf color font "times"
unset key
set size ratio -1 1,1
set xrange [0:1]
set yrange [0:1]
set output "../".path."/output/N".n."D".D."K".k."T".t."P".P."/graph.pdf"
set title "Instancia do PCDCP (N = ".n."; D = ".d."; K = ".k."; T = ".t."; P = ".p.")"
vertices = "../".path."/output/N".n."D".D."K".k."T".t."P".P."/vertices.txt"
edges = "../".path."/output/N".n."D".D."K".k."T".t."P".P."/edges.txt"
plot vertices using ($1):($2):($3) with points pointtype 7 pointsize variable linecolor rgb "black", \
	 edges with lines linetype 1 linewidth 1 linecolor rgb "forest-green"
