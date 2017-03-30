# To be completed with the corresponding paths. Make sure to keep the -I/-L/-l at the beginning.

MUMPS_I = -I/u/geo2/sw/Toolchains/gcc+system/4.8/Packages/libs/mumps-mpi/4.10/include
MUMPS_L = -L/u/geo2/sw/Toolchains/gcc+system/4.8/Packages/libs/mumps-mpi/4.10/lib -L/u/geo2/sw/Toolchains/gcc+system/4.8/Packages/libs/scotch-mpi/5.1/lib
MUMPS_l =  -ldmumps -lptscotcherr
PKG_LIBS =$(MUMPS_L) $(MUMPS_l)
