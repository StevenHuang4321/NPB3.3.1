ERROR_FLAGS = -std=gnu99 -Wall -Werror
CFLAGS = $(ERROR_FLAGS) -D_REENTRANT -DLinux -D_GNU_SOURCE -O0 -g

HEADER_FILES=numap.h

all: dynamic_lib static_lib

static_lib: numap.c 
	gcc $(CFLAGS) -c numap.c numap_analyse.c 
	ar rcs ../libnumap.a numap.o numap_analyse.o

dynamic_lib: numap.c $(HEADER_FILES)
	gcc $(CFLAGS) -fpic -c numap.c numap_analyse.c 
	gcc -shared -o ../libnumap.so numap.o numap_analyse.o -lnuma -lpfm

clean:
	rm -f *.o *.a *.so* 
