#!/usr/bin/gnuplot
set term gif enhanced animate

set output 'Orb.gif'
set title "{/Symbol Y}(x,t)"

set xrange [-2:2]
set yrange [-2.8:2.8]
#set xrange [0:14000e-10]
#set yrange [0:14000e-10]
#set zrange [0:14000e-10]
set xlabel "Position (Bohr Radii)"
set ylabel "Amplitude"

i=0
n=999

load "orb.gnuplot"
