SHELL=/bin/sh
BENCHMARK=mg
BENCHMARKU=MG

include ../config/make.def

OBJS = mg.o ${COMMON}/print_results.o  \
       ${COMMON}/${RAND}.o ${COMMON}/timers.o mylib.o

include ../sys/make.common


${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${FMPI_LIB} -I/opt/pasa_lib2/quartz-master/src/lib -L/opt/pasa_lib2/quartz-master/build -lnvmemul -lnuma -lpthread 

mylib.o :	../mylib/mylib.c
	gcc -Wall ${CMPI_INC} -I/opt/pasa_lib2/quartz-master/src/lib -L/opt/pasa_lib2/quartz-master/build/src/lib -c ../mylib/mylib.c $(CFLAGS) -lnvmemul  -lnuma -lpthread

mg.o:		mg.f  globals.h mpinpb.h npbparams.h
	${FCOMPILE} mg.f

clean:
	- rm -f *.o *~ 
	- rm -f npbparams.h core



