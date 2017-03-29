c---------------------------------------------------------------------
c  Parameter lm (declared and set in "npbparams.h") is the log-base2 of 
c  the edge size max for the partition on a given node, so must be changed 
c  either to save space (if running a small case) or made bigger for larger 
c  cases, for example, 512^3. Thus lm=7 means that the largest dimension 
c  of a partition that can be solved on a node is 2^7 = 128. lm is set 
c  automatically in npbparams.h
c  Parameters ndim1, ndim2, ndim3 are the local problem dimensions. 
c---------------------------------------------------------------------

      include 'npbparams.h'
    
      integer nm      ! actual dimension including ghost cells for communications
     >      , nv      ! size of rhs array
     >      , nr      ! size of residual array
     >      , nm2     ! size of communication buffer
     >      , maxlevel! maximum number of levels
      pointer(ptr_nm, nm)
      pointer(ptr_nv, nv)
      pointer(ptr_nr, nr)
      pointer(ptr_nm2, nm2)
      pointer(ptr_maxlevel, maxlevel)
  
      parameter( nm=2+2**lm, nv=(2+2**ndim1)*(2+2**ndim2)*(2+2**ndim3) )
      parameter( nm2=2*nm*nm, maxlevel=(lt_default+1) )
      parameter( nr = (8*(nv+nm**2+5*nm+14*lt_default-7*lm))/7 )
      integer maxprocs
      pointer(ptr_maxprocs, maxprocs)     

      parameter( maxprocs = 131072 )  ! this is the upper proc limit that 
                                      ! the current "nr" parameter can handle
c---------------------------------------------------------------------
      integer nbr(3,-1:1,maxlevel), msg_type(3,-1:1)
      pointer(ptr_nbr, nbr)
      pointer(ptr_msg_type, msg_type)

      integer  msg_id(3,-1:1,2),nx(maxlevel),ny(maxlevel),nz(maxlevel)
      pointer(ptr_msg_id, msg_id)
      pointer(ptr_nx, nx)
      pointer(ptr_ny, ny)
      pointer(ptr_nz, nz)
      common /mg3/ ptr_nbr,ptr_msg_type,ptr_msg_id,ptr_nx,ptr_ny,ptr_nz

      character class
      pointer(ptr_class, class)

      common /ClassType/ptr_class

      integer debug_vec(0:7)
      pointer(ptr_debug_vec, debug_vec)
      common /my_debug/ ptr_debug_vec

      integer ir(maxlevel), m1(maxlevel), m2(maxlevel), m3(maxlevel)
      pointer(ptr_ir, ir)
      pointer(ptr_m1, m1)
      pointer(ptr_m2, m2)
      pointer(ptr_m3, m3)

      integer lt, lb
      pointer(ptr_lt, lt)
      pointer(ptr_lb, lb)
      common /fap/ ptr_ir,ptr_m1,ptr_m2,ptr_m3,ptr_lt,ptr_lb

      logical dead(maxlevel), give_ex(3,maxlevel), take_ex(3,maxlevel)
      pointer(ptr_dead, dead)
      pointer(ptr_give_ex, give_ex)
      pointer(ptr_take_ex, take_ex)
      common /comm_ex/ ptr_dead, ptr_give_ex, ptr_take_ex

c---------------------------------------------------------------------
c  Set at m=1024, can handle cases up to 1024^3 case
c---------------------------------------------------------------------
      integer m
       pointer(ptr_m ,m)
c      parameter( m=1037 )
      parameter( m=nm+1 )

      double precision buff(nm2,4)
      pointer(ptr_buff, buff)
      common /buffer/ ptr_buff

c---------------------------------------------------------------------
      integer t_bench, t_init, t_psinv, t_resid, t_rprj3, t_interp, 
     >        t_norm2u3, t_comm3, t_rcomm, t_last
      pointer(ptr_t_bench, t_bench)
      pointer(ptr_t_init, t_init)
      pointer(ptr_t_psinv, t_psinv)
      pointer(ptr_t_resid, t_resid)
      pointer(ptr_t_rprj3, t_rprj3)
      pointer(ptr_t_interp, t_interp)
      pointer(ptr_t_norm2u3, t_norm2u3)
      pointer(ptr_t_comm3, t_comm3)
      pointer(ptr_t_rcomm, t_rcomm)
      pointer(ptr_t_last, t_last)

      parameter (t_bench=1, t_init=2, t_psinv=3, t_resid=4, t_rprj3=5,  
     >        t_interp=6, t_norm2u3=7, t_comm3=8, 
     >        t_rcomm=9, t_last=9)

      logical timeron
      pointer(ptr_timeron, timeron)
      common /timers/ ptr_timeron



