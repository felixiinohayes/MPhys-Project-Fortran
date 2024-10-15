set xrange [ -0.09 : 0.09]
set terminal pngcairo enhanced font "DejaVu,30" fontscale 1 size 2000, 2000
set output "band.png"
set border
unset xtics
unset ytics
set encoding iso_8859_1
set size ratio 0 1.0,1.0
set yrange [5.5: 6.5]
unset key
set mytics 2
set parametric
set palette defined (0 "#deebf7", 1 "#c6dbef", 2 "#9ecae1", 3 "#6baed6", 4 "#4292c6", 5 "#2171b5", 6 "#084594")
set trange [-10:10]
set multiplot
plot "super_H_top.dat" u 1:2:3 with l lt 1 lw 1 linecolor palette notitle
unset multiplot
