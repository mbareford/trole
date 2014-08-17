set autoscale
set palette rgbformulae 33,13,10
set pm3d map
unset pm3d
unset ztics
unset key
unset title
unset label
set xrange [*:*]
set yrange [*:*]
set zrange [*:*]

path = '"./data/output'
path2 = '"./data/recompute'
set macros

set terminal postscript eps enhanced color solid 18

unset label
set xlabel "Flare Energy (dimensionless)"
set ylabel "Flare Count (normalised)"

set xrange [0.0:0.08]
set yrange [0.0:0.08]

set key at 0.07,0.07
set key font ",20"
set key spacing 1.25

set output @path/figs/wrpf_rx10e5_lle3.eps"
set title "10^2 loops undergoing 1000 relaxation events" 
plot @path2/o_wrbc_lle3.txt" using 1:($2/100000.0) with lines lt 3 lw 2 title "published", @path/o_wrbc_lle3.txt" using 1:($2/100000.0) with lines lt 1 lw 2 title "recomputated"
set output @path/figs/wrpf_rx10e5_lle2.eps"
set title "10^3 loops undergoing 100 relaxation events"
plot @path2/o_wrbc_lle2.txt" using 1:($2/100000.0) with lines lt 3 lw 2 title "published", @path/o_wrbc_lle2.txt" using 1:($2/100000.0) with lines lt 1 lw 2 title "recomputated"
set output @path/figs/wrpf_rx10e5_lle1.eps"
set title "10^4 loops undergoing 10 relaxation events"
plot @path2/o_wrbc_lle1.txt" using 1:($2/100000.0) with lines lt 3 lw 2 title "published", @path/o_wrbc_lle1.txt" using 1:($2/100000.0) with lines lt 1 lw 2 title "recomputated"
set output @path/figs/wrpf_rx10e5_lle0.eps"
set title "10^5 loops undergoing 1 relaxation event"
plot @path2/o_wrbc_lle0.txt" using 1:($2/100000.0) with lines lt 3 lw 2 title "published", @path/o_wrbc_lle0.txt" using 1:($2/100000.0) with lines lt 1 lw 2 title "recomputated"

unset label
