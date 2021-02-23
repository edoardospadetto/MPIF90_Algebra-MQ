set terminal pdf size 7, 5 font "Latin Modern Math, 25"
set output "./gnuplot/times.pdf"

set xlabel "log(N)"
set ylabel "log(time) [s]"
set grid
set key top left
set format y "10^{%T}"
set logscale

plot "./results/times.txt" u 1:2 w lp notitle 
