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

path = '"../../data/output'
set macros

set terminal postscript eps enhanced solid 18

set xrange [*:*]
set yrange [*:*]

unset label
set xlabel "Flare Energy (dimensionless)"
set ylabel "Flare Count (normalised)"

set xrange [*:*]
set yrange [*:*]

set output @path/figs/wrpf.eps"
plot @path/o_wrbc.txt" using 1:($2/100000.0) with lines lt -1
unset label
