ef=  4.18903780
set xtics ( \
"L"   -0.837034, \
"A"    0.000000, \
"H"    0.966523 )
set xrange [   -0.837034:    0.966523]
set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 7.50in
set output "band.pdf"
set encoding iso_8859_1
set size ratio 0 1.0,1.0
set ylabel "E-E_{CBM} (eV)"
set yrange [ -2 : 2.0 ]
unset key
set ytics 1.0 scale 1 nomirror out
set mytics 2
set parametric
set trange [-10:10]
plot "band.dat" u 1:($2-ef) with l lt 1 lw 3,\
    0.000000,t with l lt 2  lc -1,\
t,0 with l lt 2  lc -1
