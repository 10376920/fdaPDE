# Completare la variabile MUMPS PAHT con il percorso di mumps (CON I)
# #da completare per utente
MUMPS_I = -I/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/mumps/5.0.1/include
MUMPS_L = -L/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/mumps/5.0.1/lib -L/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/scotch/6.0.4/lib
MUMPS_l =  -ldmumps -lptscotcherr

PKG_LIBS =$(MUMPS_L) $(MUMPS_l)
