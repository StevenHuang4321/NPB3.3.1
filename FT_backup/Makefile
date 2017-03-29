SHELL=/bin/sh
BENCHMARK=ft
BENCHMARKU=FT

include ../config/make.def

ERROR_FLAGS =  -Wall -Werror
CFLAGS = $(ERROR_FLAGS) -D_REENTRANT -DLinux -D_GNU_SOURCE -O0 -g -mcmodel=large

include ../sys/make.common

OBJS = ft.o unimem.o ${COMMON}/${RAND}.o ${COMMON}/print_results.o ${COMMON}/timers.o

${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -I/home/kwu/mpi_project/2016winter/NPB3.3-MPI/Unimem -L/home/kwu/mpi_project/20\
16winter/NPB3.3-MPI/Unimem  -o ${PROGRAM} ${OBJS} ${FMPI_LIB} -lnuma -lnumap -lpthread -LUnimem 

.f.o:
	${FCOMPILE} $<

ft.o:             ft.f  mpinpb.h npbparams.h

unimem.o: ../Unimem/unimem.c
        gcc -c ../Unimem/unimem.c -mcmodel=large -I/home/kwu/mpi_project/2016winter/NPB3.3-MPI/Unimem -L/home/kwu/mpi_project/2016winter/NPB3.3-MPI/Unimem -lnuma -lnumap -lpthread

clean:
	- rm -f *.o *~ mputil*
	- rm -f ft npbparams.h core
