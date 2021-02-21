set terminal pdf size 9, 5 font "Latin Modern Math, 25"
set encoding utf8

#set title "Ground state density for different N"
set xlabel "λ" 
set ylabel "e"
set grid
set key outside 

f(x) = 0 <= x && x <= 2 ? -1 - x**2 / 4 : -x    
g(x) = -1-x            

set output "gs_H.pdf"
plot "eig_H_2.txt" u 1:2 w l title "N=2", \
     "eig_H_3.txt" u 1:2 w l title "N=3", \
     "eig_H_4.txt" u 1:2 w l title "N=4", \
     "eig_H_5.txt" u 1:2 w l title "N=5", \
     "eig_H_6.txt" u 1:2 w l title "N=6", \
     "eig_H_7.txt" u 1:2 w l title "N=7", \
     "eig_H_8.txt" u 1:2 w l title "N=8", \
     g(x) title "Mean field" dt 2

set output "gs_I.pdf"
plot "eig_I_2.txt" u 1:2 w l title "N=2", \
     "eig_I_3.txt" u 1:2 w l title "N=3", \
     "eig_I_4.txt" u 1:2 w l title "N=4", \
     "eig_I_5.txt" u 1:2 w l title "N=5", \
     "eig_I_6.txt" u 1:2 w l title "N=6", \
     "eig_I_7.txt" u 1:2 w l title "N=7", \
     "eig_I_8.txt" u 1:2 w l title "N=8", \
     f(x) title "Mean field" dt 2
     