set terminal jpeg
set output 'data.jpg'
set xlabel 'N'
set ylabel '\l'
plot 'data.txt' u 1:2:3 w image notitle
