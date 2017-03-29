objects : test_f.o	test_c.o

test : $(objects)
	gfortran -o test $(objects)

test_c.o :  test_c.c
	gcc -c -g test_c.c

test_f.f : test_f.f
	gfortran -c -g test_f.f -fcray-pointer

clean : rm test $(objects)
