#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <asm/unistd.h>
#include <stddef.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <stdint.h>
#include <emmintrin.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>
#include <omp.h>
#include <numa.h>
#include <numaif.h>
#include "mpi.h"

//******************************************************************************
/*#define _mm_clflushopt(addr)						\
  asm volatile(".byte 0x66; clflush %0" : "+m" (*(volatile char *)addr)); */

#define FLUSH_ALIGN ((uintptr_t)64)

#define ALIGN_MASK (FLUSH_ALIGN - 1)

#define CHUNK_SIZE 128 /* 16*8 */
#define CHUNK_SHIFT 7
#define CHUNK_MASK (CHUNK_SIZE - 1)

#define DWORD_SIZE 4
#define DWORD_SHIFT 2
#define DWORD_MASK (DWORD_SIZE - 1)

#define MOVNT_SIZE 16
#define MOVNT_MASK (MOVNT_SIZE - 1)
#define MOVNT_SHIFT 4

#define MOVNT_THRESHOLD 256

static size_t Movnt_threshold = MOVNT_THRESHOLD;


// timers
double start_clflush = 0, elapsed_clflush = 0;
double start_memcpy = 0, elapsed_memcpy = 0;
double start_dram_cache = 0, elapsed_dram_cache = 0;
double start_memmove = 0, elapsed_memmove = 0;
double start_msync = 0, elapsed_msync = 0;
double start_write = 0, elapsed_write = 0;



void c_openmp_init_(int *n){
  omp_set_num_threads(*n);  
}

/*
static void predrain_fence_sfence(void) {
  _mm_sfence(); // ensure CLWB or CLFLUSHOPT completes 
}
*/
static void flush_clflush(const void *addr, size_t len){
  // printf("addr %p len %zu\n", addr, len);
   start_clflush = MPI_Wtime();

  //  printf("%p %d\n", addr, len);
  uintptr_t uptr;                                                                                                                          
   //Loop through cache-line-size (typically 64B) aligned chunks                                                                                  
   // covering the given range.                                                                                                                   
   
  // #pragma omp parallel for
  for (uptr = (uintptr_t)addr & ~(FLUSH_ALIGN - 1);
       uptr < (uintptr_t)addr + len; uptr += FLUSH_ALIGN){
    _mm_clflush((char *)uptr);
  }
  elapsed_clflush += MPI_Wtime() - start_clflush;

}

void pmem_flush(const void *addr, size_t len){
  //printf("addr %p len %zu\n", addr, len);

  //  VALGRIND_DO_CHECK_MEM_IS_ADDRESSABLE(addr, len);                                                                                            
  //  Func_flush(addr, len);                                                                                                                      
    flush_clflush(addr, len);
}

//void c_memmove_movnt_(void *pmemdest, const void *src, size_t len) {  
void *c_memmove_movnt_(void **pmemdest_, const void **src_, size_t* len_) {
  //printf("pmemdest %p src %p len %zu\n",*pmemdest_, *src_, *len_);
  start_memmove = MPI_Wtime();

  __m128i xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;
  size_t i;
  __m128i *d;
  __m128i *s;
  void *pmemdest = *pmemdest_;
  void *src = *src_;
  void *dest1 = pmemdest;
  size_t cnt;
  size_t len = *len_;

  if (len == 0 || src == pmemdest){
     return pmemdest;
    //    printf("111111111111111111111\n");
    //return;
  }
  if (len < Movnt_threshold) {
    memmove(pmemdest, src, len);
      pmem_flush(pmemdest, len);
      //      printf("11111\n");        
    return pmemdest;
    //    printf("222222222222222222222\n");
    //return;
  }

  if ((uintptr_t)dest1 - (uintptr_t)src >= len) {
    //printf("!!!!!2222222222222222222\n");
                                                                                                           
    // Copy the range in the forward direction.                                                                                                   
                                                                                                                                                
    // This is the most common, most optimized case, used unless                                                                                  
    // the overlap specifically prevents it.                                                                                                      
     

    // copy up to FLUSH_ALIGN boundary 
    cnt = (uint64_t)dest1 & ALIGN_MASK;
    // printf("cnt = %zu\n", cnt);
    if (cnt > 0) {
      cnt = FLUSH_ALIGN - cnt;

      // never try to copy more the len bytes 
      if (cnt > len)
        cnt = len;

      uint8_t *d8 = (uint8_t *)dest1;
      const uint8_t *s8 = (uint8_t *)src;
      for (i = 0; i < cnt; i++) {
        *d8 = *s8;
        d8++;
        s8++;
      }
      
	pmem_flush(dest1, cnt);
	//	printf("222222\n");

      dest1 = (char *)dest1 + cnt;
      src = (char *)src + cnt;
      len -= cnt;
    }

    d = (__m128i *)dest1;
    s = (__m128i *)src;

    cnt = len >> CHUNK_SHIFT;
    for (i = 0; i < cnt; i++) {
      xmm0 = _mm_loadu_si128(s);
      xmm1 = _mm_loadu_si128(s + 1);
      xmm2 = _mm_loadu_si128(s + 2);
      xmm3 = _mm_loadu_si128(s + 3);
      xmm4 = _mm_loadu_si128(s + 4);
      xmm5 = _mm_loadu_si128(s + 5);
      xmm6 = _mm_loadu_si128(s + 6);
      xmm7 = _mm_loadu_si128(s + 7);
      s += 8;
      _mm_stream_si128(d,xmm0);
      _mm_stream_si128(d + 1,xmm1);
      _mm_stream_si128(d + 2,xmm2);
      _mm_stream_si128(d + 3,xmm3);
      _mm_stream_si128(d + 4,xmm4);
      _mm_stream_si128(d + 5, xmm5);
      _mm_stream_si128(d + 6,xmm6);
      _mm_stream_si128(d + 7,xmm7);
      //VALGRIND_DO_FLUSH(d, 8 * sizeof(*d));                                                                                                     
      d += 8;
    }
    
    // copy the tail (<128 bytes) in 16 bytes chunks 
    len &= CHUNK_MASK;
    //    printf("len = %zu\n", len);

    if (len != 0) {
      cnt = len >> MOVNT_SHIFT;
      for (i = 0; i < cnt; i++) {
        xmm0 = _mm_loadu_si128(s);
        _mm_stream_si128(d, xmm0);
        //VALGRIND_DO_FLUSH(d, sizeof(*d));                                                                                                       
        s++;
        d++;
      }
    }

    // copy the last bytes (<16), first dwords then bytes 
    len &= MOVNT_MASK;
    //printf("len = %zu\n", len);

    if (len != 0) {
      cnt = len >> DWORD_SHIFT;
      int32_t *d32 = (int32_t *)d;
      int32_t *s32 = (int32_t *)s;
      for (i = 0; i < cnt; i++) {
        _mm_stream_si32(d32, *s32);
        //VALGRIND_DO_FLUSH(d32, sizeof(*d32));                                                                                                   
        d32++;
        s32++;
      }
      cnt = len & DWORD_MASK;
      uint8_t *d8 = (uint8_t *)d32;
      const uint8_t *s8 = (uint8_t *)s32;

      for (i = 0; i < cnt; i++) {
        *d8 = *s8;
        d8++;
        s8++;
      }
     
	 pmem_flush(d32, cnt);
	 // printf("3333333\n");

      }
  } else {
    
    // printf("!!!!!333333333333333\n");
    // Copy the range in the backward direction.                                                                                                  
                                                                                                                                                
    // This prevents overwriting source data due to an                                                                                            
     // overlapped destination range.                
     

    dest1 = (char *)dest1 + len;
    src = (char *)src + len;

    cnt = (uint64_t)dest1 & ALIGN_MASK;
    if (cnt > 0) {
      // never try to copy more the len bytes 
      if (cnt > len)
        cnt = len;

      uint8_t *d8 = (uint8_t *)dest1;
      const uint8_t *s8 = (uint8_t *)src;
      for (i = 0; i < cnt; i++) {
        d8--;
        s8--;
        *d8 = *s8;
      }
     
	pmem_flush(d8, cnt);
	//	printf("4444444\n");

      dest1 = (char *)dest1 - cnt;
      src = (char *)src - cnt;
      len -= cnt;
    }

    d = (__m128i *)dest1;
    s = (__m128i *)src;
    cnt = len >> CHUNK_SHIFT;
    for (i = 0; i < cnt; i++) {
      xmm0 = _mm_loadu_si128(s - 1);
      xmm1 = _mm_loadu_si128(s - 2);
      xmm2 = _mm_loadu_si128(s - 3);
      xmm3 = _mm_loadu_si128(s - 4);
      xmm4 = _mm_loadu_si128(s - 5);
      xmm5 = _mm_loadu_si128(s - 6);
      xmm6 = _mm_loadu_si128(s - 7);
      xmm7 = _mm_loadu_si128(s - 8);
      s -= 8;
      _mm_stream_si128(d - 1, xmm0);
      _mm_stream_si128(d - 2, xmm1);
      _mm_stream_si128(d - 3, xmm2);
      _mm_stream_si128(d - 4, xmm3);
      _mm_stream_si128(d - 5, xmm4);
      _mm_stream_si128(d - 6, xmm5);
      _mm_stream_si128(d - 7, xmm6);
      _mm_stream_si128(d - 8, xmm7);
      d -= 8;
      // VALGRIND_DO_FLUSH(d, 8 * sizeof(*d));                                                                                                    
    }

    // copy the tail (<128 bytes) in 16 bytes chunks 
    len &= CHUNK_MASK;
    if (len != 0) {
      cnt = len >> MOVNT_SHIFT;
      for (i = 0; i < cnt; i++) {
        d--;
        s--;
        xmm0 = _mm_loadu_si128(s);
        _mm_stream_si128(d, xmm0);
        //      VALGRIND_DO_FLUSH(d, sizeof(*d));                                                                                                 
      }
    }
    // copy the last bytes (<16), first dwords then bytes 
    len &= MOVNT_MASK;
    if (len != 0) {
      cnt = len >> DWORD_SHIFT;
      int32_t *d32 = (int32_t *)d;
      int32_t *s32 = (int32_t *)s;
      for (i = 0; i < cnt; i++) {
        d32--;
        s32--;
	_mm_stream_si32(d32, *s32);
        //      VALGRIND_DO_FLUSH(d32, sizeof(*d32));                                                                                             
      }

      cnt = len & DWORD_MASK;
      uint8_t *d8 = (uint8_t *)d32;
      const uint8_t *s8 = (uint8_t *)s32;

      for (i = 0; i < cnt; i++) {
        d8--;
        s8--;
        *d8 = *s8;
      }
	pmem_flush(d8, cnt);
	//	printf("5555555\n");

    }
  }

  // serialize non-temporal store instructions 
  //  predrain_fence_sfence();
  //  printf("3333333333333333333333333333\n");
  
  // return;
  elapsed_memmove += MPI_Wtime() - start_memmove;


  return pmemdest;
}


void* c_memmove_movnt_dram_cache(void *pmemdest, const void *src, size_t len) {                                                              
//void *c_memmove_movnt_dramcache(void *pmemdest, const void **src_, size_t* len_) {
  //printf("pmemdest %p src %p len %zu\n",*pmemdest, *src_, *len_);
  //  start_memmove = MPI_Wtime();                                                                                                                

  __m128i xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;
  size_t i;
  __m128i *d;
  __m128i *s;
  // void *pmemdest = *pmemdest_;
  // void *src = *src_;
  void *dest1 = pmemdest;
  size_t cnt;
  //  size_t len = *len_;


  if (len == 0 || src == pmemdest){
    return pmemdest;
    //    printf("111111111111111111111\n");                                                                                                      
    //return;                                                                                                                                     
  }
  if (len < Movnt_threshold) {
    memmove(pmemdest, src, len);
    // if (*clflush_check == 1) {                                                                                                                 
    //    pmem_flush(pmemdest, len);
    // }                                                                                                                                         
    return pmemdest;
    //    printf("222222222222222222222\n");                                                                                                      
    //return;                                                                                                                                     
  }
  if ((uintptr_t)dest1 - (uintptr_t)src >= len) {

    // Copy the range in the forward direction.                                                                                                   

    // This is the most common, most optimized case, used unless                                                                                  
    // the overlap specifically prevents it.                                                                                                      


    // copy up to FLUSH_ALIGN boundary                                                                                                            
    cnt = (uint64_t)dest1 & ALIGN_MASK;
    if (cnt > 0) {
      cnt = FLUSH_ALIGN - cnt;

      // never try to copy more the len bytes                                                                                                     
      if (cnt > len)
        cnt = len;

      uint8_t *d8 = (uint8_t *)dest1;
      const uint8_t *s8 = (uint8_t *)src;
      for (i = 0; i < cnt; i++) {
        *d8 = *s8;
        d8++;
        s8++;
      }
      //      if (*clflush_check == 1){                                                                                                           
      //      pmem_flush(dest1, cnt);
      //}                                                                                                                                        
      dest1 = (char *)dest1 + cnt;
      src = (char *)src + cnt;
      len -= cnt;
    }

    d = (__m128i *)dest1;
    s = (__m128i *)src;
    cnt = len >> CHUNK_SHIFT;
    for (i = 0; i < cnt; i++) {
      xmm0 = _mm_loadu_si128(s);
      xmm1 = _mm_loadu_si128(s + 1);
      xmm2 = _mm_loadu_si128(s + 2);
      xmm3 = _mm_loadu_si128(s + 3);
      xmm4 = _mm_loadu_si128(s + 4);
      xmm5 = _mm_loadu_si128(s + 5);
      xmm6 = _mm_loadu_si128(s + 6);
      xmm7 = _mm_loadu_si128(s + 7);
      s += 8;
      _mm_stream_si128(d,xmm0);
      _mm_stream_si128(d + 1,xmm1);
      _mm_stream_si128(d + 2,xmm2);
      _mm_stream_si128(d + 3,xmm3);
      _mm_stream_si128(d + 4,xmm4);
      _mm_stream_si128(d + 5, xmm5);
      _mm_stream_si128(d + 6,xmm6);
      _mm_stream_si128(d + 7,xmm7);
      //VALGRIND_DO_FLUSH(d, 8 * sizeof(*d));                                                                                                     
      d += 8;
    }
    // copy the tail (<128 bytes) in 16 bytes chunks                                                                                              
    len &= CHUNK_MASK;
    if (len != 0) {
      cnt = len >> MOVNT_SHIFT;
      for (i = 0; i < cnt; i++) {
        xmm0 = _mm_loadu_si128(s);
        _mm_stream_si128(d, xmm0);
        //VALGRIND_DO_FLUSH(d, sizeof(*d));                                                                                                       
        s++;
        d++;
      }
    }

    // copy the last bytes (<16), first dwords then bytes                                                                                         
    len &= MOVNT_MASK;
    if (len != 0) {
      cnt = len >> DWORD_SHIFT;
      int32_t *d32 = (int32_t *)d;
      int32_t *s32 = (int32_t *)s;
      for (i = 0; i < cnt; i++) {
        _mm_stream_si32(d32, *s32);
        //VALGRIND_DO_FLUSH(d32, sizeof(*d32));                                                                                                   
        d32++;
        s32++;
      }
      cnt = len & DWORD_MASK;
      uint8_t *d8 = (uint8_t *)d32;
      const uint8_t *s8 = (uint8_t *)s32;

      for (i = 0; i < cnt; i++) {
        *d8 = *s8;
        d8++;
        s8++;
      }
      //        if (*clflush_check == 1) {                                                                                                       
      //      pmem_flush(d32, cnt);
      // }                                                                                                                                      
    }
  } else {
    // Copy the range in the backward direction.                                                                                                  

    // This prevents overwriting source data due to an                                                                                            
    // overlapped destination range.                                                                                                              


    dest1 = (char *)dest1 + len;
    src = (char *)src + len;

    cnt = (uint64_t)dest1 & ALIGN_MASK;
    if (cnt > 0) {
      // never try to copy more the len bytes                                                                                                     
      if (cnt > len)
        cnt = len;

      uint8_t *d8 = (uint8_t *)dest1;
      const uint8_t *s8 = (uint8_t *)src;
      for (i = 0; i < cnt; i++) {
        d8--;
        s8--;
        *d8 = *s8;
      }
      //      if (*clflush_check == 1){                                                                                                           
      //      pmem_flush(d8, cnt);
      //}                                                                                                                                        
      dest1 = (char *)dest1 - cnt;
      src = (char *)src - cnt;
      len -= cnt;
    }
    d = (__m128i *)dest1;
    s = (__m128i *)src;
    cnt = len >> CHUNK_SHIFT;
    for (i = 0; i < cnt; i++) {
      xmm0 = _mm_loadu_si128(s - 1);
      xmm1 = _mm_loadu_si128(s - 2);
      xmm2 = _mm_loadu_si128(s - 3);
      xmm3 = _mm_loadu_si128(s - 4);
      xmm4 = _mm_loadu_si128(s - 5);
      xmm5 = _mm_loadu_si128(s - 6);
      xmm6 = _mm_loadu_si128(s - 7);
      xmm7 = _mm_loadu_si128(s - 8);
      s -= 8;
      _mm_stream_si128(d - 1, xmm0);
      _mm_stream_si128(d - 2, xmm1);
      _mm_stream_si128(d - 3, xmm2);
      _mm_stream_si128(d - 4, xmm3);
      _mm_stream_si128(d - 5, xmm4);
      _mm_stream_si128(d - 6, xmm5);
      _mm_stream_si128(d - 7, xmm6);
      _mm_stream_si128(d - 8, xmm7);
      d -= 8;
      // VALGRIND_DO_FLUSH(d, 8 * sizeof(*d));                                                                                                    
    }
    // copy the tail (<128 bytes) in 16 bytes chunks                                                                                              
    len &= CHUNK_MASK;
    if (len != 0) {
      cnt = len >> MOVNT_SHIFT;
      for (i = 0; i < cnt; i++) {
        d--;
        s--;
        xmm0 = _mm_loadu_si128(s);
        _mm_stream_si128(d, xmm0);
        //      VALGRIND_DO_FLUSH(d, sizeof(*d));                                                                                                 
      }
    }
    // copy the last bytes (<16), first dwords then bytes                                                                                         
    len &= MOVNT_MASK;
    if (len != 0) {
      cnt = len >> DWORD_SHIFT;
      int32_t *d32 = (int32_t *)d;
      int32_t *s32 = (int32_t *)s;
      for (i = 0; i < cnt; i++) {
        d32--;
        s32--;
	_mm_stream_si32(d32, *s32);
        //      VALGRIND_DO_FLUSH(d32, sizeof(*d32));                                                                                             
      }

      cnt = len & DWORD_MASK;
      uint8_t *d8 = (uint8_t *)d32;
      const uint8_t *s8 = (uint8_t *)s32;

      for (i = 0; i < cnt; i++) {
        d8--;
        s8--;
        *d8 = *s8;
      }
      //      if (*clflush_check == 1){                                                                                                           
      // pmem_flush(d8, cnt);
      // }                                                                                                                                       
    }
  }

  // serialize non-temporal store instructions                                                                                                    
  //  predrain_fence_sfence();                                                                                                                    
  //  printf("3333333333333333333333333333\n");                                                                                                   
  // return;                                                                                                                                      
  //elapsed_memmove += MPI_Wtime() - start_memmove;                                                                                               
  return pmemdest;
}

//******************************************************************************
double tmp = 0;
FILE *f;

void *data_cache;
void *dest_data_cache;
/*
inline void do_flush_cache(volatile double *p) {
  asm volatile ("clflush (%0)" :: "r"(p));
}
*/
void c_flush_(void *arr, int *size_type, int *size_total) {                                                                              
  /*void *p = arr;                                                                                                                              
  int len = *size_total;
  //#pragma omp parallel for
  for(; p < arr + len; p += 64){                                                                                                      
    asm volatile ("clflush (%0)" :: "r"((volatile void *)p));
    }*/
  // flush_clflush(arr, *size_total);
} 

void c_malloc_(void **ptr, int *size) {
  *ptr=(void *)malloc(*size);
  memset(*ptr, 'c', *size);
}

void c_memcpy_(void **new, void **old, int *size) {
  start_memcpy = MPI_Wtime();
   memcpy(*new, *old, *size);
  elapsed_memcpy += MPI_Wtime() - start_memcpy;
  //  printf("%p %p\n", *new, *old);
   
  flush_clflush(*new, *size);
  /* void *p = old;                                                                                                                               
  int len = *size_total;   
  //  #pragma omp parallel for                                                                                                                    
  for(; p < arr + len; p += 64){                                                                                                                  
    asm volatile ("clflush (%0)" :: "r"((volatile void *)p));                                                                                
    }*/                                                                                                             
  //  Dong comment: call clflush after memcpy? 
}


char* concat(const char *s1, const char *s2) {
  const size_t len1 = strlen(s1);
  const size_t len2 = strlen(s2);
  char *result = malloc(len1+len2+1);//+1 for the zero-terminator
  //in real code you would check for errors in malloc here
  memcpy(result, s1, len1);
  memcpy(result+len1, s2, len2+1);//+1 to copy the null-terminator
  return result;
}
/*
void c_memfileopen_(int *curr_rank){
  char * path = "/dev/shm/kw_copy.data";
  char * rank_c = (char *)malloc(5);
  sprintf(rank_c, "%d", *curr_rank);
  char * fname = concat(path, rank_c);
 
  //  printf("%s\n", fname);
  if (access(fname, F_OK != -1)) {
    remove(fname);
    f = fopen(fname, "wb");
  } else {
    f = fopen(fname, "wb");
  }

  // Dong comment: after mem-based file write is done, call clflush?
}

//Dong comment: this is individual I/O, which is fine
void c_memfileopendisk_(int *curr_rank) {
  char * path = "kw_copy.data";
  char * rank_c = (char *)malloc(5);
  sprintf(rank_c, "%d", *curr_rank);
  char * fname = concat(path, rank_c);

  //  printf("%s\n", fname);                                                                                                                      
  if (access(fname, F_OK != -1)) {
    remove(fname);
    f = fopen(fname, "wb");
  } else {
    f = fopen(fname, "wb");
  }
                                                                                                           
}

void c_memfileopendisk_local_(int *curr_rank) {
  char * path = "/tmp/kw_copy.data";
  char * rank_c = (char *)malloc(5);
  sprintf(rank_c, "%d", *curr_rank);
  char * fname = concat(path, rank_c);

  //  printf("%s\n", fname);                                                                                                                      
  if (access(fname, F_OK != -1)) {
    remove(fname);
    f = fopen(fname, "wb");
  } else {
    f = fopen(fname, "wb");
  }

  }*/
/*
void c_memfileclose_(){
  fclose(f);
}
*/
void c_memwrite_(void **old_, int *unit_size, int *total_size, int *curr_rank, int *type){
  
  start_write = MPI_Wtime();

  void *old = *old_;
  char *fname;
  if (*type == 1) {
    
    // memory based file write
    
    char * path = "/dev/shm/kw_copy.data";
    char * rank_c = (char *)malloc(5);
    sprintf(rank_c, "%d", *curr_rank);
    fname = concat(path, rank_c);
    //  printf("%s\n", fname);                                                                                                                
    
    if (access(fname, F_OK != -1)) {
      remove(fname);
      f = fopen(fname, "wb");
    } else {
      f = fopen(fname, "wb");
    }

  } else if (*type == 2){
    // disk based file write (remote home)

    char * path = "kw_copy.data";
    char * rank_c = (char *)malloc(5);
    sprintf(rank_c, "%d", *curr_rank);
    fname = concat(path, rank_c);

    //  printf("%s\n", fname);                                                                                                                   
    if (access(fname, F_OK != -1)) {
      remove(fname);
      f = fopen(fname, "wb");
    } else {
      f = fopen(fname, "wb");
    }

  } else if (*type == 3){
    
    // disk based file write (local /tmp)
    
    char * path = "/tmp/kw_copy.data";
    char * rank_c = (char *)malloc(5);
    sprintf(rank_c, "%d", *curr_rank);
    fname = concat(path, rank_c);

    //  printf("%s\n", fname);                                                                                                                  
    if (access(fname, F_OK != -1)) {
      remove(fname);
      f = fopen(fname, "wb");
    } else {
      f = fopen(fname, "wb");
    }

  } else {
    printf("open file error !\n");
  }
  fwrite(old, *unit_size, *total_size, f);
  elapsed_write += MPI_Wtime() - start_write;

  start_msync = MPI_Wtime();
  msync(old, *total_size, MS_SYNC);
  elapsed_msync += MPI_Wtime() - start_msync;
  //Dong comment: after checkpoint at the end of each iteration, call fclose for disk-based checkpoint
  fclose(f);
}


void c_memcmp_(void* old, void *new, int *len) {
  if (memcmp(old, new, *len) != 0) {
    printf("yes!!!!!!!!!\n");
  } else {
    printf("no!!!!!!!!!!!!!\n");
  }
}
/*
void c_msync_(void *addr, size_t *len) {
  msync(addr, *len, MS_SYNC);
}
*/

void c_lib_init_(int *curr_rank) {
  // if (*curr_rank == 0 || * curr_rank == 4 || *curr_rank == 8 || *curr_rank == 12) {
    data_cache = malloc(128000000);
    memset(data_cache, 128000000, '0');
    dest_data_cache = malloc(128000000);
    //}
}
void c_dram_cache_cp_(void **dest_data_cache, void **data_cache, int *len){
   
    start_dram_cache = MPI_Wtime();

    memcpy(*dest_data_cache, *data_cache, *len);

    elapsed_dram_cache += MPI_Wtime() - start_dram_cache;
}

void c_dram_cache_move_(void **dest_data_cache, void **data_cache, int *len){
                                                                                                   
    start_dram_cache = MPI_Wtime();
    c_memmove_movnt_dram_cache(*dest_data_cache, *data_cache, *len);
    elapsed_dram_cache += MPI_Wtime() - start_dram_cache;

}


void find_memory_node_for_addr(void* ptr) {
  int numa_node = -1;
  if(get_mempolicy(&numa_node, NULL, 0, ptr, MPOL_F_NODE | MPOL_F_ADDR) < 0) {
    printf("WARNING: get_mempolicy failed\n");
  }
  printf("numa_node = %d\n", numa_node);

}


void  unimem_malloc_(void **data, int *size, int *pos) {

  if (*pos == 1) {
    // bandwidth mode                                                                                                                           
    //allocate on DRAM                                                                                                                            
    numa_set_preferred(0);
    numa_run_on_node(0);
    //        *data =(void *) malloc(*size);                                                                                                      
    *data = (void *)numa_alloc_onnode(*size, 0);
    memset(*data, 'c', *size);
    find_memory_node_for_addr(*data);
    //printf("DRAM\n");                                                                                                                           

  } else {
    // bandwidth mode                                                                                                                          
    // allocate on NVM                                                                                                                            
    numa_set_preferred(1);
    numa_run_on_node(0);                                                                                            
    *data = (void *)numa_alloc_onnode(*size, 1);
    memset(*data, 'c', *size);
    find_memory_node_for_addr(*data);                                                                                                         
    //printf("NVM\n");                                                                                                                         
  }
}

void c_print_timer_(int *curr_rank) {
 
  printf("memcpy: %d %f\n", *curr_rank, elapsed_memcpy);
  printf("clflush: %d %f\n", *curr_rank, elapsed_clflush);
  printf("data_cache: %d %f\n", *curr_rank, elapsed_dram_cache);
  printf("memmove: %d %f\n\n\n", *curr_rank, elapsed_memmove - elapsed_clflush);
  printf("write: %d %f\n", *curr_rank, elapsed_write);
  printf("msync: %d %f\n", *curr_rank, elapsed_msync);

}

void c_memcpy_init(void **dest_data, void **data_src, int *len) {
  memcpy(*dest_data, *data_src, *len);
}

void c_memcpy_nek_(void *src, int *len){
  start_dram_cache = MPI_Wtime();
  // allocate on NVM                                                                                                                              
  numa_set_preferred(1);
  numa_run_on_node(0);
  void *dest1 = (void *)numa_alloc_onnode(*len, 1);

  memcpy(dest1, src, *len);
  find_memory_node_for_addr(dest1);

  elapsed_dram_cache += MPI_Wtime() - start_dram_cache;

  
  start_memcpy = MPI_Wtime();

  // allocate on NVM                                                                                                                          
  numa_set_preferred(1);
  numa_run_on_node(0);
  void *dest2 = (void *)numa_alloc_onnode(*len, 1);
  //  memset(dest, 'c', *size);
  memcpy(dest2, dest1, *len); 
  find_memory_node_for_addr(dest2);

  elapsed_memcpy += MPI_Wtime() - start_memcpy;

  flush_clflush(dest2, *len);

  
}
void c_switch_(void **a, void **b){
  void *tmp = *a;
  *a = *b;
  *b = tmp;
}

/*
void c_memcpy_nek_(int *len) {
  start_memcpy = MPI_Wtime();                                                                                                                     
     
  double size =*len/16;
  // allocate on NVM                                                                                                                              
  numa_set_preferred(1);                                                                                                                          
  numa_run_on_node(0);                                                                                                                            
  void *dest = (void *)numa_alloc_onnode(*len, 1);                                                                                                
  //  memset(dest, 'c', *size);                                                                                                                   
  memcpy(dest, src, *len);                                                                                                                        
  find_memory_node_for_addr(dest);                                                                                                                
                                                                                                                                                  
  elapsed_memcpy += MPI_Wtime() - start_memcpy;                                                                                                   

  flush_clflush(dest, *len);
  }*/
