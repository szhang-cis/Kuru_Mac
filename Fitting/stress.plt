#GNUPLOT
set encoding utf8
set terminal postscript eps enhanced 'Helvetica' 7 color
set output 'stress.eps'
set multiplot
set size 0.33,0.33
#
set origin 0.0,0.5
set key top right
set title 'S11'
set xlabel 'Time[Months]'
set ylabel 'Cauchy Stress[kPa]'
plot [0:*][*:*] 'stress.dat' u 1:($2/1000) t '' w l
#
set origin 0.33,0.5
set key bottom right
set title 'S22'
set xlabel 'Time[Months]'
set ylabel 'Cauchy Stress[kPa]'
plot [0:*][*:*] 'stress.dat' u 1:($3/1000) t '' w l
#
set origin 0.66,0.5
set key bottom right
set title 'S33'
set xlabel 'Time[Months]'
set ylabel 'Cauchy Stress[kPa]'
plot [0:*][*:*] 'stress.dat' u 1:($4/1000) t '' w l
#
set origin 0.0,0.0
set key bottom left
set title 'S12'
set xlabel 'Time[Months]'
set ylabel 'Cauchy Stress[kPa]'
plot [0:*][*:*] 'stress.dat' u 1:($5/1000) t '' w l
#
set origin 0.33,0.0
set key bottom left
set title 'S13'
set xlabel 'Time[Months]'
set ylabel 'Cauchy Stress[kPa]'
plot [0:*][*:*] 'stress.dat' u 1:($6/1000) t '' w l
#
set origin 0.66,0.0
set key bottom left
set title 'S23'
set xlabel 'Time[Months]'
set ylabel 'Cauchy Stress[kPa]'
plot [0:*][*:*] 'stress.dat' u 1:($7/1000) t '' w l
