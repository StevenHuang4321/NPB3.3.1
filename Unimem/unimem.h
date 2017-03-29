#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "pmalloc.h"
#include <sys/types.h>
#include <sys/syscall.h>
#include <numa.h>
#include <numaif.h>
#include <pthread.h>
//#include "numap_link.h"
#include "numap/numap.h"

#define DRAM_SIZE 1000  //MB   cg 250  ft 150 mg 100                                                                                              
#define MAX_PHASE_NUM 20
#define MAX_OBJ_NUM 10
//#define MAX_HW_EVENT 2                                                                                                        
#define BW_NVM 8224  //MB/S                                                                                                     
#define BW_DRAM 13203 //MB/S                                                                                                    
#define cache_line_size 64*1e-6    // 64B  * 10e-6                                                                             
#define num_process 1

#define LOOKAHEAD_WIN -1   //-1 means the whole loop                                                                            
#define MODEL_FREQ  20  //2                                                                                                    
#define  sampling_rate 2000 // 65000                                                                                            
#define num_ele 6  // cg 7 ft 6  mg 4 BT 16
#define CF 400 //BW/bw

//sampling variables                                                                                                            
struct numap_sampling_measure sm;
int res;  //error variables                                                                                                     

//phase and iteration count
int re_sampling = 0;                                                                                                    
unsigned long long exe_time_before_move = 0;
unsigned long long exe_time_after_move = 0; 
int phase_counter = 0;
int iter_num = 0;
int sampling_check = 0;
int phase_total = 0;

int curr_dram_used = 0;
//int c[num_ele + 1][DRAM_SIZE +1];   //Max benefit value                                                                    
   
long double dp[num_ele + 1][DRAM_SIZE + 1];                                                                                                       

long double cache = cache_line_size * sampling_rate;
unsigned long long overlap_time;
unsigned long long overlap_time_max;
long double overlap_;
long double benefit_local;
long double benefit_global;

int placement_method;   //data placement method   -1 global_view   1 local_view

//helper thread
pthread_t tid;
pthread_attr_t attr;


typedef struct PHASE_INFO{
  unsigned long long exe_time;

  //  hw_events events_counter[num_ele];
  obj_info objs_book_local[num_ele];
  //  data_movement_task move_tasks;                                                                                            
  //int phase_type;  //COMP=0; COMM=1                                                                                           
} phase_info;

phase_info phase_record[MAX_PHASE_NUM];
//obj_addr addr_book[num_ele];
obj_info objs_book_global[num_ele];

typedef struct MOVE_INFO{
  void *ptr;
  int size;
} move_info;

//function list                                                                                                                 
static __inline__ unsigned long long rdtsc(void) {
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

// compartor of qsort                                                                                                          
                                                                                                                                
int comparator (const void * a, const void * b) {

  unsigned long long tmp_a = ((obj_info*)a) -> benefit;
  unsigned long long tmp_b = ((obj_info*)b) -> benefit;

  return (tmp_b - tmp_a);

}

long double max(long double a, long double b) { return (a > b)? a : b; }

/*
int unimem_init_(void **ptr_colidx,  unsigned long long *colidx_size,
                 void **ptr_a, unsigned long long *a_size,
                 void **ptr_w, unsigned long long *w_size,
                 void **ptr_z, unsigned long long *z_size,
                 void **ptr_p, unsigned long long *p_size,
                 void **ptr_q, unsigned long long *q_size,
                 void **ptr_r, unsigned long long *r_size);*/
int unimem_init_(void ** ptr_colidx, unsigned long long* colidx_addr, unsigned long long *colidx_size,
                 void ** ptr_a, unsigned long long* a_addr, unsigned long long *a_size,
                 void ** ptr_w, unsigned long long* w_addr, unsigned long long *w_size,
                 void ** ptr_z, unsigned long long* z_addr, unsigned long long *z_size,
                 void ** ptr_p, unsigned long long* p_addr, unsigned long long *p_size,
                 void ** ptr_q, unsigned long long* q_addr, unsigned long long *q_size,
                 void ** ptr_r, unsigned long long* r_addr, unsigned long long *r_size,
		 int *init_pos);

void unimem_init_ft_(void ** ptr_u, unsigned long long* u_addr, unsigned long long *u_size,
		     void ** ptr_u0, unsigned long long* u0_addr, unsigned long long *u0_size,
		     void ** ptr_u1, unsigned long long* u1_addr, unsigned long long *u1_size,
		     void ** ptr_u2, unsigned long long* u2_addr, unsigned long long *u2_size,
		     void ** ptr_sums, unsigned long long* sums_addr, unsigned long long *sums_size,
		     void ** ptr_twiddle, unsigned long long* twiddle_addr, unsigned long long *twiddle_size,
		     int *init_pos);

void unimem_init_mg_(void ** ptr_buff, unsigned long long* buff_addr, unsigned long long * buff_size,
		     void ** ptr_u, unsigned long long*u_addr,unsigned long long * u_size,
                     void ** ptr_v, unsigned long long*v_addr,unsigned long long * v_size,
                     void ** ptr_r, unsigned long long*r_addr,unsigned long long * r_size,
                     int *init_pos);
void unimem_init_bt_(void ** ptr_u, unsigned long long* u_addr, unsigned long long * u_size,
		     void ** ptr_us, unsigned long long* us_addr, unsigned long long * us_size,
		     void ** ptr_vs, unsigned long long* vs_addr, unsigned long long * vs_size,
		     void ** ptr_ws, unsigned long long* ws_addr, unsigned long long * ws_size,
		     void ** ptr_qs, unsigned long long* qs_addr, unsigned long long * qs_size,
		     void ** ptr_rho_i, unsigned long long* rho_i_addr, unsigned long long * rho_i_size,
		     void ** ptr_square, unsigned long long* square_addr, unsigned long long * square_size,
		     void ** ptr_rhs, unsigned long long *rhs_addr, unsigned long long * rhs_size,
		     void ** ptr_forcing, unsigned long long *forcing_addr, unsigned long long *forcing_size,
		     void ** ptr_out_buffer, unsigned long long *out_buffer_addr, unsigned long long *out_buffer_size,
		     void ** ptr_in_buffer, unsigned long long *in_buffer_addr, unsigned long long *in_buffer_size,
		     void ** ptr_fjac, unsigned long long *fjac_addr, unsigned long long *fjac_size,
		     void ** ptr_njac, unsigned long long *njac_addr, unsigned long long *njac_size,
		     void ** ptr_lhsa, unsigned long long *lhsa_addr, unsigned long long *lhsa_size,
		     void ** ptr_lhsb, unsigned long long *lhsb_addr, unsigned long long *lhsb_size,
		     void ** ptr_lhsc, unsigned long long *lhsc_addr, unsigned long long *lhsc_size,
                     int * init_pos);



int numap_sampling_init();
int numap_read_begin_();
int numap_read_end_(int phase_counter);
void beginonephase_();
void endonephase_();
void begin_one_iteration(int* main_loop);
void end_one_iteration(int *main_loop);
long double objs_benefit_local(int phase_index, int phase_index_max);
long double objs_benefit_global();
long double knapSack(unsigned long max_W, obj_info objs_book[], int n);
unsigned long long overlap_gap(int phase_index, int item, int phase_index_max);
void  unimem_malloc_(void **data, int *size, int *pos);
void find_memory_node_for_addr(void* ptr);
void* unimem_move(int item_index);
void update_obj_info(int local_global, int item_index, int phase);
//void move_test_a(void **data_old);
void unimem_malloc_ft_(unsigned long long *ptr);
