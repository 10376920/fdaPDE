original:
	cd makevars/original/ &&\
    R CMD SHLIB ../../src/*.cpp ../../src/*.c -o ../../fdaPDE_original.so
schur:
	cd makevars/schur/ &&\
	R CMD SHLIB ../../src/*.cpp ../../src/*.c -o ../../fdaPDE_schur.so
mumps:
	cd makevars/mumps/ &&\
	R CMD SHLIB ../../src/*.cpp ../../src/*.c -o ../../fdaPDE_MUMPS.so
clean:
	rm src/*.o
