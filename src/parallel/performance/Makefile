# Link to BLAS and LAPACK
LIB_LIST = -L${BLASDIR} ${BLASLIB} -lpthread

F95 = gfortran  # This is gfortran again, since we run it on 1 processor only
OPT = -O3 -Wall -fimplicit-none

all: testblas

OBJS = testblas.o lublas2.o lublas3.o rsymm.o check_lu.o bsblas2.o


%.o: %.f90                           
	$(F95) -c $(OPT) $<

testblas: $(OBJS)
	$(F95) -o $@ $(OBJS) $(LIB_LIST)

# Note that the blas are only used in the linking step

clean: 
	rm -f *.o  core.* testblas  # creates a command to remove old .o and core files

