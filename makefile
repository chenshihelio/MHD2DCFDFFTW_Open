CC=mpicc
# CC = gcc


sysname=$(uname -n)
UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
	OPTIONS = -O3 -I/usr/local/include -L/usr/local/lib -lfftw3
endif
ifeq ($(UNAME), Linux)
	OPTIONS = -O3 -I/opt/fftw/3.3.4/intel/mvapich2_ib/include -L/opt/fftw/3.3.4/intel/mvapich2_ib/lib -lfftw3
endif

# # for cygwin
# OPTIONS = -O3 -I/usr/include -L/usr/lib -lfftw3

OPTIONS = -O3 -lm -std=c99 -I/usr/local/include -L/usr/local/lib -lfftw3

EXE = mhd.exe 
MAIN = mhd.c

FILE01 = initialize.c mathlib.c usrInitialize.c writeOutput.c \
	timeAdvance.c parallel.c calcDerivs.c calcRHS.c filtering.c \
	divBcleaning.c

OBJMAIN = ${MAIN:.f90=.o}
OBJ01 = ${FILE01:.c=.o}
OBJ   = $(OBJMAIN) $(OBJ01)


$(EXE): $(OBJ)
	$(CC) -o $(EXE) $(OBJ) $(OPTIONS) 
$(OBJMAIN):
	$(CC) -c $(MAIN) $(OPTIONS)
$(OBJ01):
	$(CC) -c $(FILE01) $(OPTIONS)

clean:
	rm *.o

cleanData:
	rm *.dat rec log
