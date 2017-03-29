/*typedef struct HW_EVENTS {                                                                                                                        
  unsigned long long addr_begin;                                                                                                                  
  unsigned long long addr_end;                                                                                                                    
  int counts;     // #of main memory access(last level cache miss)                                                                                
                                                                                                                                                 
} hw_events;                                                                                                                                      
                                                                                                                                                  
                                                                                                                                                  
typedef struct OBJ_ADDR {                                                                                                                     
   unsigned long long addr_begin;                                                                                                             
   unsigned long long addr_end;                                                                                                                   
} obj_addr;                                                                                                                                       
*/
typedef struct OBJ_INFO {
  void **ptr;
  unsigned long long addr_begin;
  unsigned long long addr_end;                                                                                                                
  unsigned long size;    //mb
  unsigned long real_size;   //bit 
  int counts;
  long double  benefit;
  long double data_copy_time;   //data copy to and back                                                                                           
  int flag;     // if in the bag 0 no 1 yes
  int used;    // 0 no 1 yes                                                                                                                      
  int move;    // 0 no 1 yes                          
  int curr_pos;   //-1 NVM DRAM 1
                                                                         
} obj_info;
