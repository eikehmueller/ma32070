# Link to BLAS and LAPACK
LIB_LIST = -L${BLASDIR} ${BLASLIB} -lpthread

F95 = mpif90  # This is gfortran again, since we run it on 1 processor only
OPT = -O3 -Wall -fimplicit-none

all: testblas testcomm

OBJS_BLAS = testblas.o lublas2.o lublas3.o rsymm.o check_lu.o bsblas2.o
OBJS_COMM = testcomm.o pingpong.o

%.o: %.f90                           
	$(F95) -c $(OPT) $<

testblas: $(OBJS_BLAS)
	$(F95) -o $@ $(OBJS_BLAS) $(LIB_LIST)

testcomm: $(OBJS_COMM)
	$(F95) -o $@ $(OBJS_COMM) $(LIB_LIST)

# Note that the blas are only used in the linking step

clean: 
	rm -f *.o  core.* testblas testcomm  # remove old .o and core files

