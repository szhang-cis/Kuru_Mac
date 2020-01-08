#GNUPLOT
set encoding utf8
set terminal postscript eps enhanced 'Helvetica' 10 color
set output 'analstre.eps'
set multiplot
set size 0.5,0.5
#
set origin 0.0,0.0
set key top left
set title 'Circumferential Stress'
set xlabel 'Steps'
set ylabel 'Cauchy Stress[kPa]'
plot 'analstrethin.dat' u 1:($6/1000) t 'Thin' w l,'analstrethick.dat' u 1:($6/1000) t 'Thick' w l,'analstrethinm.dat' u 1:($6/1000) t 'Thin-mod' w l

set origin 0.5,0.5
set key top left
set title 'Longitudinal Stress'
set xlabel 'Steps'
set ylabel 'Cauchy Stress[kPa]'
plot 'analstrethin.dat' u 1:($7/1000) t 'Thin' w l,'analstrethick.dat' u 1:($7/1000) t 'Thick' w l,'analstrethinm.dat' u 1:($7/1000) t 'Thin-mod' w l

set origin 0.5,0.0
set key top right
set title 'Radial Stress'
set xlabel 'Steps'
set ylabel 'Cauchy Stress[kPa]'
plot 'analstrethin.dat' u 1:($5/1000) t 'Thin' w l,'analstrethick.dat' u 1:($5/1000) t 'Thick' w l,'analstrethinm.dat' u 1:($5/1000) t 'Thin-mod' w l

set origin 0.0,0.5
set key top right
set title 'Radius'
set xlabel 'Steps'
set ylabel 'Radius[mm]'
plot 'analstrethin.dat' u 1:($4*10) t 'Thin' w l,'analstrethick.dat' u 1:($4*10) t 'Thick' w l,'analstrethinm.dat' u 1:($4*10) t 'Thin-mod' w l
