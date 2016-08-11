set terminal postscript eps color font "times"
unset key
set size ratio -1 1,1
set xrange [0:1]
set yrange [0:1]
do for [n in "10 20 50 100 200"] {
	do for [d in "03 05 07"] {
		do for [k in "0 2 4 6 8 10"] {
			do for [t in "0 1"] {
				do for [r in "0 1 2"] {
					do for [p in "025 05"] {
						set output "../linearProgram/output/N".n."D".d."K".k."T".t."R".r."P".p."/plot.eps"
						set title "Instancia do PCDCP (N = ".n.";D = ".d.";K = ".k.";T = ".t."; R = ".r."; P = ".p.")"
						vertices = "../linearProgram/output/N".n."D".d."K".k."T".t."R".r."P".p."/vertices.txt"
						solutionEdges = "../linearProgram/output/N".n."D".d."K".k."T".t."R".r."P".p."/solutionEdges.txt"
						solutionVertices = "../linearProgram/output/N".n."D".d."K".k."T".t."R".r."P".p."/solutionVertices.txt"
						plot vertices using ($1):($2):($3) with points pointtype 7 pointsize variable linecolor rgb "black", \
							 solutionEdges with lines linetype 1 linewidth 1 linecolor rgb "forest-green", \
							 solutionVertices with circles linetype 1 linewidth 0.5 linecolor rgb "red"
					}
				}
			}
		}
	}
}
