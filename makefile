original:
	cd makevars/original/ &&\
    R CMD SHLIB ../../src/*.cpp ../../src/*.c -o ../../fdaPDE_original.so
schur:
	cd makevars/schur/ &&\
	R CMD SHLIB ../../src/*.cpp ../../src/*.c -o ../../fdaPDE_schur.so
mumps_pord:
	cd makevars/mumps_pord/ &&\
	R CMD SHLIB ../../src/*.cpp ../../src/*.c -o ../../fdaPDE_MUMPS_PORD.so
mumps_scotch:
	cd makevars/mumps_scotch/ &&\
	R CMD SHLIB ../../src/*.cpp ../../src/*.c -o ../../fdaPDE_MUMPS_SCOTCH.so
mumps_whole:
	cd makevars/mumps_whole/ &&\
	R CMD SHLIB ../../src/*.cpp ../../src/*.c -o ../../fdaPDE_MUMPS_whole.so
stochastic:
	cd makevars/stochastic/ &&\
	R CMD SHLIB ../../src/*.cpp ../../src/*.c -o ../../fdaPDE_stochastic.so
temp:
	cd makevars/temp/ &&\
	R CMD SHLIB ../../src/*.cpp ../../src/*.c -o ../../fdaPDE_temp.so
woodbury_whole:
	MAKEFLAGS="CXX=mpic++  SHLIB_CXXLD=mpic++" &&\
	cd makevars/woodbury_whole/ &&\
	R CMD SHLIB ../../src/*.cpp ../../src/*.c ../../src/LinearSolvers/*.cpp -o ../../fdaPDE_woodbury_whole.so
woodbury_decomposeQ:
	cd makevars/woodbury_decomposeQ/ &&\
	R CMD SHLIB ../../src/*.cpp ../../src/*.c -o ../../fdaPDE_woodbury_decomposeQ.so
clean:
	rm src/*.o src/LinearSolvers/*.o 
