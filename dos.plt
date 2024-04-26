# Set the output file format and name (optional)
set terminal png
set output 'plot.png'

# Set up plot settings
set xlabel "X-axis"
set ylabel "Y-axis"
set title "Multiple lines plot"

# Define the plot style and colors
set style line 1 lc rgb 'blue' lw 2
set style line 2 lc rgb 'red' lw 2
set style line 3 lc rgb 'green' lw 2
# You can define more styles if you have more line colors

# Plot the data
plot 'ene_total.dat' using 2:3:1 with lines lc variable title "Lines"