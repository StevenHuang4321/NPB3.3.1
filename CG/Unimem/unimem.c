#include "unimem.h"

// unimem init api (store targeted data obejct address and numap counter init)
int unimem_init_(void ** ptr_colidx, unsigned long long* colidx_addr, unsigned long long *colidx_size,
		 void ** ptr_a, unsigned long long* a_addr, unsigned long long *a_size,
		 void ** ptr_w, unsigned long long* w_addr, unsigned long long *w_size,
		 void ** ptr_z, unsigned long long* z_addr, unsigned long long *z_size,
		 void ** ptr_p, unsigned long long* p_addr, unsigned long long *p_size,
		 void ** ptr_q, unsigned long long* q_addr, unsigned long long *q_size,
		 void ** ptr_r, unsigned long long* r_addr, unsigned long long *r_size,
                 int *init_pos){
  /*
  addr_book[0].addr_begin = *colidx_addr;
  addr_book[0].addr_end = *colidx_addr + *colidx_size;
  addr_book[1].addr_begin = *a_addr;
  addr_book[1].addr_end = *a_addr + *a_size;
  addr_book[2].addr_begin = *x_addr;  //????
  addr_book[2].addr_end = *x_addr + *x_size;
  addr_book[3].addr_begin = *z_addr;
  addr_book[3].addr_end = *z_addr + *z_size;
  addr_book[4].addr_begin = *p_addr;
  addr_book[4].addr_end = *p_addr + *p_size;
  addr_book[5].addr_begin = *q_addr;
  addr_book[5].addr_end = *q_addr + *q_size;
  addr_book[6].addr_begin = *r_addr;
  addr_book[6].addr_end = *r_addr + *r_size;
  */
  //printf("ptr = %llu\n", (unsigned long long)*ptr_colidx);
  //printf("loc(colidx): %llu\n",*colidx_addr);
   
  //init local_view data object info
  int i;
  //printf("current pos = %d\n",*init_pos);
  for (i = 0; i < 8; i++) {
      phase_record[i].objs_book_local[0].ptr = &(*ptr_colidx);
      phase_record[i].objs_book_local[0].addr_begin = *colidx_addr;
      phase_record[i].objs_book_local[0].addr_end = *colidx_addr + *colidx_size;
      phase_record[i].objs_book_local[0].real_size = *colidx_size;
      phase_record[i].objs_book_local[0].size = *colidx_size * 1.25e-7 ;
      phase_record[i].objs_book_local[0].curr_pos = *init_pos;       // -1 NVM DRAM 1
      phase_record[i].objs_book_local[0].flag =0;

      phase_record[i].objs_book_local[1].ptr = &(*ptr_a);
      phase_record[i].objs_book_local[1].addr_begin = *a_addr;
      phase_record[i].objs_book_local[1].addr_end = *a_addr + *a_size;
      phase_record[i].objs_book_local[1].real_size = *a_size;
      phase_record[i].objs_book_local[1].size = *a_size * 1.25e-7 ;
      phase_record[i].objs_book_local[1].curr_pos = *init_pos;
      phase_record[i].objs_book_local[1].flag =0;

      phase_record[i].objs_book_local[2].ptr = &(*ptr_w);
      phase_record[i].objs_book_local[2].addr_begin =  *w_addr;                                                                
      phase_record[i].objs_book_local[2].addr_end = *w_addr + *w_size;
      phase_record[i].objs_book_local[2].real_size = *w_size;
      phase_record[i].objs_book_local[2].size = *w_size * 1.25e-7 ;
      phase_record[i].objs_book_local[2].curr_pos = *init_pos;
      phase_record[i].objs_book_local[2].flag =0;

      phase_record[i].objs_book_local[3].ptr = &(*ptr_z);
      phase_record[i].objs_book_local[3].addr_begin = *z_addr;
      phase_record[i].objs_book_local[3].addr_end = *z_addr + *z_size;
      phase_record[i].objs_book_local[3].real_size = *z_size;
      phase_record[i].objs_book_local[3].size = *z_size * 1.25e-7 ;
      phase_record[i].objs_book_local[3].curr_pos = *init_pos;
      phase_record[i].objs_book_local[3].flag =0;

      phase_record[i].objs_book_local[4].ptr = &(*ptr_p);
      phase_record[i].objs_book_local[4].addr_begin = *p_addr;
      phase_record[i].objs_book_local[4].addr_end = *p_addr + *p_size;
      phase_record[i].objs_book_local[4].real_size = *p_size;
      phase_record[i].objs_book_local[4].size = *p_size * 1.25e-7 ;
      phase_record[i].objs_book_local[4].curr_pos = *init_pos;
      phase_record[i].objs_book_local[4].flag =0;

      phase_record[i].objs_book_local[5].ptr = &(*ptr_q);
      phase_record[i].objs_book_local[5].addr_begin = *q_addr;
      phase_record[i].objs_book_local[5].addr_end = *q_addr + *q_size;
      phase_record[i].objs_book_local[5].real_size = *q_size;
      phase_record[i].objs_book_local[5].size = *q_size * 1.25e-7 ;
      phase_record[i].objs_book_local[5].curr_pos = *init_pos;
      phase_record[i].objs_book_local[5].flag =0;

      phase_record[i].objs_book_local[6].ptr = &(*ptr_r);
      phase_record[i].objs_book_local[6].addr_begin = *r_addr;
      phase_record[i].objs_book_local[6].addr_end = *r_addr + *r_size;
      phase_record[i].objs_book_local[6].real_size = *r_size;
      phase_record[i].objs_book_local[6].size = *r_size * 1.25e-7;
      phase_record[i].objs_book_local[6].curr_pos = *init_pos;
      phase_record[i].objs_book_local[6].flag = 0;
  }
  // init obj usage for each phase                                                                       
                  
  // phase 0 
  phase_record[0].objs_book_local[0].used = 1;
  phase_record[0].objs_book_local[1].used = 0;
  phase_record[0].objs_book_local[2].used = 1;
  phase_record[0].objs_book_local[3].used = 0;
  phase_record[0].objs_book_local[4].used = 1;
  phase_record[0].objs_book_local[5].used = 0;
  phase_record[0].objs_book_local[6].used = 0;
  
   // phase 1                                       
                                                      
  phase_record[1].objs_book_local[0].used = 0;
  phase_record[1].objs_book_local[1].used = 0;
  phase_record[1].objs_book_local[2].used = 1;
  phase_record[1].objs_book_local[3].used = 0;
  phase_record[1].objs_book_local[4].used = 0;
  phase_record[1].objs_book_local[5].used = 1;
  phase_record[1].objs_book_local[6].used = 0;
    
 // phase 2                                       
                
  phase_record[2].objs_book_local[0].used = 0;
  phase_record[2].objs_book_local[1].used = 0;
  phase_record[2].objs_book_local[2].used = 1;
  phase_record[2].objs_book_local[3].used = 0;
  phase_record[2].objs_book_local[4].used = 0;
  phase_record[2].objs_book_local[5].used = 1;
  phase_record[2].objs_book_local[6].used = 0;
    
 // phase 3
  phase_record[3].objs_book_local[0].used = 0;
  phase_record[3].objs_book_local[1].used = 0;
  phase_record[3].objs_book_local[2].used = 1;
  phase_record[3].objs_book_local[3].used = 0;
  phase_record[3].objs_book_local[4].used = 1;
  phase_record[3].objs_book_local[5].used = 1;
  phase_record[3].objs_book_local[6].used = 0;
   
  // phase 4                                                                                                  
  phase_record[4].objs_book_local[0].used = 0;
  phase_record[4].objs_book_local[1].used = 0;
  phase_record[4].objs_book_local[2].used = 0;
  phase_record[4].objs_book_local[3].used = 0;
  phase_record[4].objs_book_local[4].used = 0;
  phase_record[4].objs_book_local[5].used = 0;
  phase_record[4].objs_book_local[6].used = 0;
  
   // phase 5                                                                                    
  phase_record[5].objs_book_local[0].used = 0;
  phase_record[5].objs_book_local[1].used = 0; 
  phase_record[5].objs_book_local[2].used = 0;
  phase_record[5].objs_book_local[3].used = 1;
  phase_record[5].objs_book_local[4].used = 1;
  phase_record[5].objs_book_local[5].used = 1;
  phase_record[5].objs_book_local[6].used = 1;
 
    // phase 6                                                                                   
  phase_record[6].objs_book_local[0].used = 0;
  phase_record[6].objs_book_local[1].used = 0;
  phase_record[6].objs_book_local[2].used = 0;
  phase_record[6].objs_book_local[3].used = 0;
  phase_record[6].objs_book_local[4].used = 0;
  phase_record[6].objs_book_local[5].used = 0;
  phase_record[6].objs_book_local[6].used = 0;
   
  // phase 7                                                                                         

  phase_record[7].objs_book_local[0].used = 0;
  phase_record[7].objs_book_local[1].used = 0; 
  phase_record[7].objs_book_local[2].used = 0;
  phase_record[7].objs_book_local[3].used = 0;
  phase_record[7].objs_book_local[4].used = 1;
  phase_record[7].objs_book_local[5].used = 0;
  phase_record[7].objs_book_local[6].used = 1;

  
  // init objects info for global_view
  for (i = 0; i < num_ele; i++) {
    objs_book_global[i].ptr = phase_record[0].objs_book_local[i].ptr;
    objs_book_global[i].addr_begin = phase_record[0].objs_book_local[i].addr_begin;
    objs_book_global[i].addr_end = phase_record[0].objs_book_local[i].addr_end;
    objs_book_global[i].real_size = phase_record[0].objs_book_local[i].real_size;
    objs_book_global[i].size = (phase_record[0].objs_book_local[i].addr_end - phase_record[0].objs_book_local[i].addr_begin) * 1.25e-7;
    objs_book_global[i].data_copy_time = objs_book_global[i].size * 2 / (double)BW_NVM;
    objs_book_global[i].curr_pos = *init_pos; // -1 NVM 1 DRAM
    objs_book_global[i].flag = 0;
    //printf("copy time = %Lf\n", objs_book_global[i].data_copy_time);
    //objs_book_global[i].data_copy_time = (addr_book[i].addr_end - addr_book[i].addr_begin) * 2 / (BW_NVM * 10e7);
  }

  // init helper thread 
  pthread_attr_init(&attr);
  
  // init numap library
  numap_sampling_init();
  //printf("numap init\n");
  return 0;
}

int numap_sampling_init() {
  // Init numap                                                                                                                
  res = numap_init();
  if(res < 0) {
    fprintf(stderr, "numap_init : %s\n", numap_error_message(res));
    return -1;
  }

  return 0;
}

int numap_read_begin_(){

  res = numap_sampling_init_measure(&sm, 1, sampling_rate, 64);
  if(res < 0) {
    fprintf(stderr, "numap_sampling_init error : %s\n", numap_error_message(res));
    return -1;
  }
  sm.tids[0] = syscall(SYS_gettid);

  // Start memory read access sampling                                                                                                            
  //  printf("\nStarting memory sampling\n");
  res = numap_sampling_wr_start(&sm);
  if(res < 0) {
    fprintf(stderr, "numap_sampling_start error : %s\n", numap_error_message(res));
    return -1;
  }
   
  return 0;
}

int numap_read_end_(int phase_counter){

  res = numap_sampling_wr_stop(&sm);
  if(res < 0) {
    printf("numap_sampling_stop error : %s\n", numap_error_message(res));
    return -1;
  }
  
  //printf("\nMemory read sampling results\n");

  //         numap_sampling_write_print(&sm, 0);
	      numap_sampling_wr_summary(&sm, num_ele,  phase_record[phase_counter].objs_book_local);
  
  // for ()
  //printf("", phase_record[phase_counter].objs_book_local[]);
  res = numap_sampling_end(&sm);
  if (res < 0) {
    printf("numap_sampling_close file error: %s\n", numap_error_message(res));
  }
  
  return 0;
}

void begin_one_phase_() {
  if (sampling_check == 1) {
    // printf("phase_counter %d\n", phase_counter);

   phase_record[phase_counter].exe_time = rdtsc();  //starting time                                             
                                                               
   //TODO: start collecting events
 
   numap_read_begin_();
  } else {
    //TODO: Read “phase_info_record[global_phase_counter+LOOKAHEAD_WIN]”;
    //if(data copy cost is smaller than the time between the two phases)
    if(1) {
	//Ask a helper thread to do data copy by putting a data copy request into the queue;
    }

  }
}

void end_one_phase_() {
  if (sampling_check == 1) {
    phase_record[phase_counter].exe_time = rdtsc() - phase_record[phase_counter].exe_time;  //exe time                        
  
    // printf("exe_time: %lld\n", phase_record[phase_counter].exe_time);                                           
    // phase_record[phase_counter].hw_events[x] = //end events                                                                  
    numap_read_end_(phase_counter);
    
  }
  phase_counter++;
}

void begin_one_iteration_(int* main_loop) {
  // printf("!!!!!!!!!!!start \n");

  phase_counter = 0;
  //  if(iter_num == 0  && *main_loop == 1) {
  if((iter_num % MODEL_FREQ == 0 || iter_num == 0 ) && *main_loop == 1) {
    sampling_check = 1;
  } else {
    sampling_check = 0;
  }
  //printf("sampling_check  =%d \n", sampling_check);
}


void end_one_iteration_(int *main_loop) {
  //printf("!!!!!!!!!!!end \n");
  if(sampling_check == 1 && *main_loop == 1) {
   
    //TODO: data placement algorithm
    
    //TODO: calculate cost of global view and local view
    benefit_global = benefit_local = 0;
    placement_method = -1;

    //local_view
    int phase_index, item_index;
    for (phase_index = 0 ; phase_index < phase_counter; phase_index++) {
      //printf("main phase = %d\n", i);
      benefit_local += objs_benefit_local(phase_index, phase_counter);
    }
    
    //global_view
    for (item_index = 0; item_index  < num_ele; item_index++) {
      objs_book_global[item_index].counts = 0;
      for (phase_index = 0; phase_index < phase_counter; phase_index++) {
	objs_book_global[item_index].counts += phase_record[phase_index].objs_book_local[item_index].counts;
      }
    }
    benefit_global = objs_benefit_global();     
  
    for (item_index = 0; item_index  < num_ele; item_index++) {
      //      printf("gloabl_view loc(%lld) size= %lu count(%d)  benefit: %Lf flag = %d\n", objs_book_global[item_index].addr_begin, objs_book_global[item_index].size, objs_book_global[item_index].counts, objs_book_global[item_index].benefit, objs_book_global[item_index].flag);
    }
    
    //printf("benefit_local = %Lf \n", benefit_local);
    //printf("benefit_gloabl = %Lf \n", benefit_global);
    
    if (benefit_global >= benefit_local) {
      placement_method = -1;
    } else {
      placement_method = 1;
    }
    
    //void *data_colidx;
    //void *data_a;
    // void *data_p;
    //numa_set_preferred(0);
    //  numa_run_on_node(0);                                                                                                            
    
    //data_colidx = (void *)numa_alloc_onnode(objs_book_global[0].real_size, 0);
    //memset(data_colidx, 'c', objs_book_global[0].real_size);
    //memcpy(data_colidx, *ptr_colidx, objs_book_global[0].real_size);
    //printf("11111111loc(a)%p\n",objs_book_global[1].ptr);
    //  move_test_a( objs_book_global[1].ptr);
    //    move_test_a(&(*ptr_a)); 
    //data_a = (void *)numa_alloc_onnode(objs_book_global[1].real_size, 0);
    // memset(data_a, 'c', objs_book_global[1].real_size);
    //memcpy(data_a, *ptr_a, objs_book_global[1].real_size);
    //data_p = (void *)numa_alloc_onnode(objs_book_global[4].real_size, 0);
    //memset(data_p, 'c', objs_book_global[4].real_size);
    //memcpy(data_p, *ptr_p, objs_book_global[4].real_size);
    //*ptr_colidx = data_colidx;
    // *ptr_a = data_a;
    //*ptr_p = data_p;
    
    
    //    move_test_(*ptr_colidx, (int)objs_book_global[0].size);
    // global_view to move data objects
     if (placement_method == -1) {
      // create a helper thread to excute data object move
      int i;
      for (i = 0; i < num_ele; i++) {
	//	if (objs_book_global[i].flag == 1) {
	  pthread_create(&tid, &attr, unimem_move, i);
	  pthread_join(tid, NULL);
	  //}
      }
      
      }
    // local_view to move data objects
    if (placement_method == 1) {

    }
    sampling_check = 0;
  }
  // if(iter_num == 0  && *main_loop == 1) {
  if ((iter_num % MODEL_FREQ == 0 || iter_num == 0 ) && *main_loop == 1) {
    phase_total = phase_counter;
  }
  if (*main_loop == 1) {
    iter_num++;
  }
  //  printf("current iter: %d\n", iter_num);
}



long double objs_benefit_local(int phase_index, int phase_index_max) {
  
  int item_index;
  for (item_index = 0; item_index < num_ele; item_index++) {
    
     // check if this data object is highly memory intensive
     if (phase_record[phase_index].objs_book_local[item_index]. counts > 0) {
       
       overlap_time_max = overlap_gap(phase_index, item_index, phase_index_max);
       overlap_ = (long double)overlap_time_max * 10e-9;
      
         // check if overlap time can cover data copy time * 2
	if (overlap_ >= objs_book_global[item_index].data_copy_time) {
	  //  printf("bigger!!!");
	   phase_record[phase_index].objs_book_local[item_index].benefit = (phase_record[phase_index].objs_book_local[item_index].counts 
                                                       * cache/ BW_NVM) 
                                                       - (phase_record[phase_index].objs_book_local[item_index].counts
                                                       * cache / BW_DRAM);
	 } else {
	  //printf("smaller!!!");
	   phase_record[phase_index].objs_book_local[item_index].benefit = (phase_record[phase_index].objs_book_local[item_index].counts
     							       * cache/ BW_NVM)
    	                                                      - (phase_record[phase_index].objs_book_local[item_index].counts
     	                                                      * cache / BW_DRAM)
     	                                                      - (objs_book_global[item_index].data_copy_time - overlap_);
	 }
     } else {
              phase_record[phase_index].objs_book_local[item_index].benefit = -999;
     }
     //printf("loc(%lld) count(%d) benefit: %Lf\n", phase_record[phase_index].objs_book_local[item_index].addr_begin, phase_record[phase_index].objs_book_local[item_index].counts, phase_record[phase_index].objs_book_local[item_index].benefit);
     phase_record[phase_index].objs_book_local[item_index].flag = 0;

  }

  //call knapSack to choose maximum benefit decision 
 
  /* printf("max_local = %Lf\n", knapSack((unsigned long)DRAM_SIZE, phase_record[phase_index].objs_book_local, num_ele));  
  int i;
   for (i = 0; i < num_ele; i++) {
     printf("loc(%lld) size= %lu count(%d)  benefit: %Lf flag = %d\n", phase_record[phase_index].objs_book_local[i].addr_begin, phase_record[phase_index].objs_book_local[i].size,  phase_record[phase_index].objs_book_local[i].counts, phase_record[phase_index].objs_book_local[i].benefit, phase_record[phase_index].objs_book_local[i].flag);
     }
  */
  return knapSack((unsigned long)DRAM_SIZE, phase_record[phase_index].objs_book_local, num_ele);
  
}

long double objs_benefit_global(){
    int i;
    for (i = 0; i < num_ele; i++) {
      
      if ( objs_book_global[i]. counts > 0) {
      
	objs_book_global[i].benefit = (objs_book_global[i].counts * cache / BW_NVM)
                                    -(objs_book_global[i].counts * cache / BW_DRAM);
      //printf("!!!gloabl loc(%lld) count (%d) benefit: %Lf\n", objs_book_global[i].addr_begin, objs_book_global[i].counts, objs_book_global[i].benefit);
      } else {
	objs_book_global[i].benefit = -999;
      }
        objs_book_global[i].flag = 0;
    
    }

    //call knapSack to choose maximum benefit decision
                                                                         
    //    printf("max = %Lf\n", knapSack((unsigned long)DRAM_SIZE, objs_book_global, num_ele));
   
    return knapSack((unsigned long)DRAM_SIZE, objs_book_global, num_ele);
}

long double knapSack(unsigned long max_W, obj_info objs_book[], int n) {
  //int knapSack(int W, int wt[], int val[], int n) {

  //int i, w;
  //int dp[n+1][max_W+1];
  unsigned long w;

  // printf("max_W = %lu\n", max_W);
  int i;
  for (i = 0; i <= n; i++) {
    for (w = 0; w <= max_W; w++) {
  
      if (i==0 || w==0){
	dp[i][w] = 0;
      } else if (objs_book[i-1].size <= w) {
	dp[i][w] = max(objs_book[i-1].benefit + dp[i-1][w - objs_book[i-1].size],  dp[i-1][w]);
      } else {
	dp[i][w] = dp[i-1][w];
      }
    }
  }
  // mark data objects in bag
    for (i = n, w = max_W; i > 0; i--) {

      if (dp[i][w] > dp[i - 1][w]) {
	objs_book[i - 1].flag = 1;
	w = w - objs_book[i].size;
      }
    }
    /*
  i = n, w = max_W;  
  while (dp[i][w] != 0) {  
    if (dp[i-1][w-objs_book[i].size]+objs_book[i].benefit == dp[i][w]) {  
      
      objs_book[i].flag = 1;  
      w -= objs_book[i].size;  
      --i;  
    }  
    } */

  return dp[n][max_W];

}


// find largest overlap time gap for each object after current phase
unsigned long long overlap_gap(int phase_index, int item, int phase_index_max){
    
  overlap_time = 0;
  overlap_time_max = 0;
  int i, move_phase = 0;
  for (i = phase_index + 1; i < phase_index_max; i++) {
    if (phase_record[i].objs_book_local[item].used == 0) {
      overlap_time += phase_record[i].exe_time;
      move_phase ++;
    } else {
      if (overlap_time > overlap_time_max) {
        overlap_time_max = overlap_time;
	overlap_time = 0;
	//find begin phase to move
	
	if (i - move_phase >= 0) {
	  phase_record[i - move_phase].objs_book_local[item].move = 1;
	} else {
	  phase_record[move_phase + 1 - i ].objs_book_local[item].move = 1;
	}
      }
    }
    
    if (i == phase_index_max - 1) {
      i = 0;
    } 
    if (i == phase_index) {
      break;
    }
  }
 
  //printf("phase =%d item = %d phase_Max = %d\n", phase_index, item, phase_index_max);
  return overlap_time_max;
}

//void dpmalloc_(void **data, int *size, int *pos) {
void  unimem_malloc_(void **data, int *size, int *pos) {
  
  if (*pos == 1) {
    printf("DRAM\N");
    //allocate on DRAM
        numa_set_preferred(0);
	//  numa_run_on_node(0);
	//    *data =(void *) malloc(*size);
	*data = (void *)numa_alloc_onnode(*size, 0);
	memset(*data, 'c', *size);
	//	find_memory_node_for_addr(*data);
       
  } else if (*pos == -1){
    printf("NVM\n");
    // allocate on NVM
    numa_set_preferred(1);
    
    //    numa_run_on_node(1);
    //    *data =(void *)pmalloc(*size);
	*data = (void *)numa_alloc_onnode(*size, 1);
	memset(*data, 'c', *size);
	//	find_memory_node_for_addr(*data);
  }
}


void find_memory_node_for_addr(void* ptr) {
  int numa_node = -1;
  if(get_mempolicy(&numa_node, NULL, 0, ptr, MPOL_F_NODE | MPOL_F_ADDR) < 0) {
     printf("WARNING: get_mempolicy failed\n");
  }
  printf("numa_node = %d\n", numa_node);
  
}

void unimem_move(int item_index) {

  //  printf("item_index = %d\n", item_index);
  //global_view to move data
  if (placement_method == -1){
    void **data_old = objs_book_global[item_index].ptr;
    //    printf("22222222222loc(a)%p\n", *data_old);
    //find_memory_node_for_addr(*data_old);
    void *data_new;
    //if (objs_book_global[item_index].flag == 1){
    if (objs_book_global[item_index].flag == 1 && objs_book_global[item_index].curr_pos == -1){
      
      // put in the bag and current on NVM 
      numa_set_preferred(0);
      data_new = (void *)numa_alloc_onnode(objs_book_global[item_index].real_size, 0);
      //memset(data_new, 'c', objs_book_global[item_index].real_size);      
      memcpy(data_new, *data_old, objs_book_global[item_index].real_size);
      *data_old = data_new;
      //  find_memory_node_for_addr(*data_old);
      objs_book_global[item_index].curr_pos = 1;
      // printf("333333333333loc(a)%p\n", *data_old);
      update_obj_info(-1, item_index, -1);

    } else if (objs_book_global[item_index].flag == 0 && objs_book_global[item_index].curr_pos == 1) {
      
      // not in the bag but current on DRAM
      numa_set_preferred(1);
      data_new = (void *)numa_alloc_onnode(objs_book_global[item_index].real_size, 1);
      //memset(data_new, 'c', objs_book_global[item_index].real_size);      
      memcpy(data_new, *data_old, objs_book_global[item_index].real_size);
      *data_old = data_new;
      //find_memory_node_for_addr(*data_old);
      //printf("333333333333loc(a)%p\n", *data_old);
      objs_book_global[item_index].curr_pos = -1;
      update_obj_info(-1, item_index, -1);

    }                                                    
   
  } else if (placement_method == 1) {
    //local_view to move data
    


    }
  
  //pthread_exit(0);
}

void update_obj_info(int local_global, int item_index, int phase) {
  // -1 updated global  1 updated local
  // global_view phase -1
  if (local_global == -1) {
    // update global info
    void ** ptr = objs_book_global[item_index].ptr;
    objs_book_global[item_index].addr_begin = (unsigned long long)*ptr;
    objs_book_global[item_index].addr_end = objs_book_global[item_index].addr_begin + objs_book_global[item_index].real_size;
    
    //update local info
    int i;
    for (i = 0; i < phase_counter; i++) {
      phase_record[i].objs_book_local[item_index].addr_begin = (unsigned long long)*ptr;
      phase_record[i].objs_book_local[item_index].addr_end = phase_record[i].objs_book_local[item_index].addr_begin + phase_record[i].objs_book_local[item_index].real_size;
    }
  } else if (local_global == 1) {
    void ** ptr =  phase_record[phase].objs_book_local[item_index].ptr;
    //update global info
    objs_book_global[item_index].ptr = phase_record[phase].objs_book_local[item_index].ptr;
    objs_book_global[item_index].addr_begin = (unsigned long long)*ptr;
    objs_book_global[item_index].addr_end = objs_book_global[item_index].addr_begin + objs_book_global[item_index].real_size;
    
    //update local info
    int i;
    for (i = 0; i < phase_counter; i++) {
      phase_record[i].objs_book_local[item_index].addr_begin = (unsigned long long)*ptr;
      phase_record[i].objs_book_local[item_index].addr_end = phase_record[i].objs_book_local[item_index].addr_begin + phase_record[i].objs_book_local[item_index].real_size;
    }
  }
}

/*
void move_test_a(void **data_old) {
  // printf("22222222222loc(a)%p\n", *data_old);
  find_memory_node_for_addr(*data_old);
  void *data_a;
  numa_set_preferred(0);
  data_a = (void *)numa_alloc_onnode(objs_book_global[1].real_size, 0);                                                                          
  memset(data_a, 'c', objs_book_global[1].real_size);                                                                                            
  memcpy(data_a, *data_old, objs_book_global[1].real_size);
  *data_old = data_a;
  find_memory_node_for_addr(*data_old);
  printf("333333333333loc(a)%p\n", *data_old);

}

*/
