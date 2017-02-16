INCLUDE    = -I.
CFLAGS     = -O2 -g
EXEC       = Herm.exe
CXX        = g++

${EXEC}: hermite_polynomial_prb.cpp hermite_polynomial.o
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} hermite_polynomial_prb.cpp hermite_polynomial.o -o ${EXEC}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<

