set terminal pdf size 15, 5 font "Latin Modern Math, 25"
set encoding utf8
set output "h_i_111.pdf"
set multiplot layout 1,2
load 'plasma.pal'

NT = 13
jx = 1.0 
jy = 1.0
jz = 1.0
j = (jx > jy ? jx : jy)
f(x) = abs(x) < 2*(j-jz) ? -j - 0.25*x**2/(j-jz) : -abs(x)-jz   
g(x) = abs(x) < 2 ? - 0.25*x**2 -1 : -abs(x) 

set xlabel "λ"
set ylabel "e"
set grid

set rmargin at screen 0.43
set lmargin at screen 0.05
set bmargin at screen 0.25
set label 11 center at graph 0.5 ,char 1 "(a) Ising Model" font ",30"
unset key 
plot for [i=2:NT] 'results/datain.txt' u (($1==i)?$2:1/0):(($1==i)?$3:1/0) w l ls i/2 lw 1.5 title 'N = '.i ,\
     g(x) lc rgb 'black' dt 2 title 'MF' 

set rmargin at screen 0.87
set lmargin at screen 0.49
set bmargin at screen 0.25
set label 11 "(b) Heisenberg Model" font ",30"
set key at screen 1,screen 0.9
plot for [i=2:NT] 'results/datahn.txt' u (($1==i)?$2:1/0):(($1==i)?$3:1/0) w l ls i/2 lw 1.5 title 'N = '.i ,\
     f(x) lc rgb 'black' dt 2 title 'MF' 
