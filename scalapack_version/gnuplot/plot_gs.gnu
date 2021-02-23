set terminal pdf size 9, 5 font "Latin Modern Math, 25"
set encoding utf8

#set title "Ground state density for different N"
set xlabel "λ" 
set ylabel "e"
set grid
set key outside 

jx = 1.0 
jy = 1.0
jz = 1.0

j = (jx > jy ? jx : jy)


f(x) = abs(x) < 2*(j-jz) ? -j - 0.25*x**2/(j-jz) : -abs(x) - jz   
g(x) = abs(x) < 2*(j-jz) ? -0.25*x**2 - 1 : -abs(x) 

set output "./results/gs_H.pdf"
plot "./results/eig_H_2.txt" u 1:2 w l title "N=2", \
     "./results/eig_H_3.txt" u 1:2 w l title "N=3", \
     "./results/eig_H_4.txt" u 1:2 w l title "N=4", \
     "./results/eig_H_5.txt" u 1:2 w l title "N=5", \
     "./results/eig_H_6.txt" u 1:2 w l title "N=6", \
     "./results/eig_H_7.txt" u 1:2 w l title "N=7", \
     "./results/eig_H_8.txt" u 1:2 w l title "N=8", \
     "./results/eig_H_9.txt" u 1:2 w l title "N=9", \
     "./results/eig_H_10.txt" u 1:2 w l title "N=10", \
     f(x) title "Mean field" dt 2

set output "./results/gs_I.pdf"
plot "./results/eig_I_2.txt" u 1:2 w l title "N=2", \
     "./results/eig_I_3.txt" u 1:2 w l title "N=3", \
     "./results/eig_I_4.txt" u 1:2 w l title "N=4", \
     "./results/eig_I_5.txt" u 1:2 w l title "N=5", \
     "./results/eig_I_6.txt" u 1:2 w l title "N=6", \
     "./results/eig_I_7.txt" u 1:2 w l title "N=7", \
     "./results/eig_I_8.txt" u 1:2 w l title "N=8", \
     g(x) title "Mean field" dt 2
     
