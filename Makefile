FC	= gfortran-mp-11
FFLAGS	= -O3 -g -fbacktrace -fbounds-check
NC_INC	= $(shell nf-config --cflags)
NC_LIB	= $(shell nf-config --flibs)
AE_INC  = -I/Users/joakimkjellsson/Downloads/aerobulk/include -I/Users/joakimkjellsson/Downloads/aerobulk/mod/
AE_LIB  = -L/Users/joakimkjellsson/Downloads/aerobulk/lib/ -L/Users/joakimkjellsson/Downloads/aerobulk/mod/ -laerobulk 

all: main.exe

main.exe: main.F90
	$(FC) -o main.exe $(NC_INC) $(NC_LIB) $(AE_INC) $(AE_LIB) $(FFLAGS) main.F90
	chmod u+rx main.exe
clean:
	rm *.o *.mod *.exe
