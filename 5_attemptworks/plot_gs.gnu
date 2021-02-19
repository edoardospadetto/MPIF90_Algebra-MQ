set terminal pdf size 9, 5 font "Latin Modern Math, 25"
set output "gs.pdf"

set encoding utf8

#set title "Ground state density for different N"
set xlabel "λ" 
set ylabel "e"
set grid
set key outside 

f(x) = 0 <= x && x <= 2 ? -1 - x**2 / 4 : -x               

plot "data3.txt" u 1:2 w l title "N=2", \
     f(x) title "Mean field" dt 2 
    #  "eig_3.txt" u 1:2 w l title "N=3", \
    #  "eig_4.txt" u 1:2 w l title "N=4", \
    #  "eig_5.txt" u 1:2 w l title "N=5", \
    #  "eig_6.txt" u 1:2 w l title "N=6", \
    #  "eig_7.txt" u 1:2 w l title "N=7", \
    #  "eig_8.txt" u 1:2 w l title "N=8", \
    #  "eig_9.txt" u 1:2 w l title "N=9", \
    #  "eig_10.txt" u 1:2 w l title "N=10"
     