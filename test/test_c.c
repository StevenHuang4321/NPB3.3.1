#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <numaif.h>
/*
#pragma pack(2)
extern struct really_type{
  //float x, y, z[6];
  //double ydbl;
  float *x;
  float *y;
  float *z;
  //double *ydbl;
};
#pragma pack()*/
//extern struct Really test_;
/*
void find_memory_node_for_addr(void* ptr) {
  int numa_node = -1;
  if(get_mempolicy(&numa_node, NULL, 0, ptr, MPOL_F_NODE | MPOL_F_ADDR) < 0) {
    printf("WARNING: get_mempolicy failed\n");
  }
  printf("numa_node = %d\n", numa_node);

  }*/


void cmalloc_(void **test) {
  //printf("%p\n", *test);
    *test = malloc(sizeof(float));
 
}

void cmove_(void **test, int *size) {
  void *new = malloc(*size);
  memcpy(new, *test, *size);
  *test = new;
}
