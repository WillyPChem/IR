#!/usr/bin/gnuplot

set terminal png
set output "IR_HCl2.png"

set xlabel 'Wavenumber'
set ylabel 'Absorption Intensity'
plot 'Absorption_Spectrum.txt'
