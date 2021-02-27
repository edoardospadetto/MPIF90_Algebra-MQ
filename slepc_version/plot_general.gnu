set terminal pdf size 9, 5 font "Latin Modern Math, 25"
set encoding utf8
set output "h_221.pdf"
load 'plasma.pal'

NT = 13
jx = 2.0 
jy = 2.0
jz = 1.0
j = (jx > jy ? jx : jy)
f(x) = abs(x) < 2*(j-jz) ? -j - 0.25*x**2/(j-jz) : -abs(x)-jz   
g(x) = abs(x) < 2 ? - 0.25*x**2 -1 : -abs(x) 

set xlabel "λ"
set ylabel "e"
set grid
set key outside 
plot for [i=2:NT] './results/h_221.txt' u (($1==i)?$2:1/0):(($1==i)?$3:1/0) w l ls i/2 lw 2 title 'N = '.i ,\
     f(x) lc rgb 'black' dt 2 title 'MF'

