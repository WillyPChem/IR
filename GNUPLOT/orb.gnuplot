plot 'orb.txt' using 1:2 index i w l lw 3 title sprintf("Re({/Symbol Y}(x,t), t(fs): %f",(i+1)*0.02418), \
'orb.txt' using 1:3 index i w l lw 3 title sprintf("P(x,t), t(fs): %f",(i+1)*0.02418), \
'orb.txt' using 1:4 index i w l lw 3 title sprintf("V(x, t=0), t(fs): %f",(i+1)*0.02418), \
'orb.txt' using 1:($5+$4) index i w l lw 3 title sprintf("V(x,t), t(fs): %f",(i+1)*0.02418)
#'Orbital.txt' using 1:4 index i w l lw 3 title sprintf("P(x,t), t(fs): %f",(i+1)*0.144)

i=i+1

if(i<n) reread
