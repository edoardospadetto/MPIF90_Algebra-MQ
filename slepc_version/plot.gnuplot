set terminal jpeg
set output 'data.jpg'
set xlabel 'lambda'
set ylabel 'e'

j = 1.5
jz = 1

NT = 12

f(x) = abs(x) < 2*(j-jz) ? -x**2/(4*(j-jz)) - j : -abs(x)-jz
plot for [i=2:NT] 'data.txt' u (($1==i)?$2:1/0):(($1==i)?$3:1/0) w l ls i/1.5 lw 1 title ''.i ,\
f(x) lc rgb 'black' title 'MF'

