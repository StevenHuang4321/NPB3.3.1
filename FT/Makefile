SHELL=/bin/sh
BENCHMARK=ft
BENCHMARKU=FT

include ../config/make.def

include ../sys/make.common

OBJS = ft.o ${COMMON}/${RAND}.o ${COMMON}/print_results.o ${COMMON}/timers.o mylib.o


${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${FMPI_LIB} -L/home/kwu/mpi_project/2016winter/NPB3.3-MPI/Unimem -L/opt/pasa_lib2/quartz-master/build/src/lib -lnvmemul -lnuma -lpthread 

mylib.o : ../mylib/mylib.c
	gcc -Wall ${CMPI_INC} -I/opt/pasa_lib2/quartz-master/src/lib -L/opt/pasa_lib2/quartz-master/build/src/lib -c ../mylib/mylib.c $(CFLAGS) -lnvmemul  -lnuma -lpthread 

.f.o:
	${FCOMPILE} $<

ft.o:             ft.f  global.h mpinpb.h npbparams.h

clean:
	- rm -f *.o *~ mputil*
	- rm -f ft npbparams.h core
