#include "stdlib.h"
#include "stdio.h"
#include <string.h>  
#include "pmalloc.h" 
/*
extern struct partiti_size{
  int naa; 
  int nzz;
  int npcols; 
  int nprows;
  int proc_col; 
  int proc_row;
  int firstrow;
  int lastrow;
  int firstcol;
  int lastcol;
  int exch_proc;
  int exch_recv_length;
  int send_start;
  int send_len;

  };*/

//extern struct partiti_size* partit_size_;

/*
void partitsize_(struct partiti_size * p) {
  printf("%ld\n", p);
  printf ("%d\n", *p);
  p = (struct partiti_size *)malloc(sizeof (struct partiti_size));
  printf("%ld\n", p);
  // printf("%d\n", sizeof(partiti_size));
  //extern struct partiti_size* partit_size_ = malloc(sizeof(str))
  //printf("sssss  %ld\n", partit_size_);
  //printf("%d\n", partit_size_. naa);
  //  printf("11111111 %ld\n", *naa);
  // free(naa);
  //naa = 4;
  //naa = (int *)malloc(sizeof(int));
  //printf("22222222 %ld\n", *naa);
  //partit_size_ = (struct partit_size *)malloc(sizeof(partit_size));
    //  partit_size_->naa = 4;
    //printf("%d\n", partit_size_->naa);
  //printf("%ld\n", &(partit_size.naa));
  //free(partit_size_.naa);
  //partit_size_.naa = malloc(sizeof(int));
}
*/
void unimem_malloc_(void **data, int *size, int *pos){
  //printf("size: %d\n", *size);
  
  //printf("1: %d\n",** data);
  //printf("1: %ld\n", *data);
  if (*pos == 1) {
    //    free(*data);
    *data =(void *) malloc(*size);
    printf("malloc\n");
  } else {
    //    free(*data);
      *data =(void *)pmalloc(*size);
      printf("pmalloc\n");
  }
  //printf("2: %ld\n", *data);
  //  **data = 99;

  //int *data = (int*)malloc(*size);
  //printf("%d\t %d \t %d\n", *size, *pos, sizeof(data));
  //for(i=0;i<n1*2;i++){data[i] = (float)i*i+1.0;}
  /* int i;   
  for (i = 0; i < *size/sizeof(int); i++) {
    (int)data[i]) = i;
    printf("%d ", data[i]);
    }*/
  //*pos = 20;
  //printf("%d\n", sizeof(data));
 // return(data);
}

void dpfree_(void **ptr, int *size, int *pos){
  if (*pos == 0) {
    free(*ptr);
  } else {
    pfree(*ptr, *size);
  }
  // free(*ptr);
}

void datamovement_(void ** old_dp_ptr, int *size, int *pos){

  //  if (*old_dp_ptr == NULL) {
  //break;
    //return NULL;
  // } else {
    
    printf("old add: %ld\n", *old_dp_ptr);
    void * new_dp_ptr;

  if (*pos == 1) {
 
    new_dp_ptr = malloc(*size);
    
  } else {
    
    new_dp_ptr = pmalloc(*size);
    
  }
  printf("allocate add: %ld\n", new_dp_ptr);                                                                                                   
  //mbind call or a memcpy                                                                                                                          
  memcpy(new_dp_ptr, *old_dp_ptr, *size);
  *old_dp_ptr = new_dp_ptr;
  printf("new add: %ld\n", *old_dp_ptr);
  //dp_free(old_dp_ptr);
  // }
  //  return new_dp_ptr;
}
