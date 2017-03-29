SHELL=/bin/sh
BENCHMARK=cg
BENCHMARKU=CG

include ../config/make.def

OBJS = cg.o mylib.o ${COMMON}/print_results.o  \
       ${COMMON}/${RAND}.o ${COMMON}/timers.o

include ../sys/make.common


${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${FMPI_LIB} -I/opt/pasa_lib2/quartz-master/src/lib -L/opt/pasa_lib2/quartz-master/build/src/lib -lnvmemul -lnuma -lpthread  

mylib.o: ../mylib/mylib.c  
	gcc -Wall  ${CMPI_INC} -I/opt/pasa_lib2/quartz-master/src/lib -L/opt/pasa_lib2/quartz-master/build/src/lib -c ../mylib/mylib.c $(CFLAGS) -lnvmemul  -lnuma -lpthread 

cg.o:		cg.f  mpinpb.h npbparams.h timing.h
	${FCOMPILE} cg.f

clean:
	- rm -f *.o *~ 
	- rm -f npbparams.h core



