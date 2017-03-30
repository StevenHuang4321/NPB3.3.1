!-------------------------------------------------------------------------!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3         !
!                                                                         !
!                                   C G                                   !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is part of the NAS Parallel Benchmark 3.3 suite.      !
!    It is described in NAS Technical Reports 95-020 and 02-007           !
!                                                                         !
!    Permission to use, copy, distribute and modify this software         !
!    for any purpose with or without fee is hereby granted.  We           !
!    request, however, that all derived work reference the NAS            !
!    Parallel Benchmarks 3.3. This software is provided "as is"           !
!    without express or implied warranty.                                 !
!                                                                         !
!    Information on NPB 3.3, including the technical report, the          !
!    original specifications, source code, results and information        !
!    on how to submit new results, is available at:                       !
!                                                                         !
!           http://www.nas.nasa.gov/Software/NPB/                         !
!                                                                         !
!    Send comments or suggestions to  npb@nas.nasa.gov                    !
!                                                                         !
!          NAS Parallel Benchmarks Group                                  !
!          NASA Ames Research Center                                      !
!          Mail Stop: T27A-1                                              !
!          Moffett Field, CA   94035-1000                                 !
!                                                                         !
!          E-mail:  npb@nas.nasa.gov                                      !
!          Fax:     (650) 604-3957                                        !
!                                                                         !
!-------------------------------------------------------------------------!


c---------------------------------------------------------------------
c
c Authors: M. Yarrow
c          C. Kuszmaul
c          R. F. Van der Wijngaart
c          H. Jin
c
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c---------------------------------------------------------------------
      program cg
c---------------------------------------------------------------------
c---------------------------------------------------------------------


      implicit none

      include 'mpinpb.h'
      include 'timing.h'
      integer status(MPI_STATUS_SIZE), request, ierr

      include 'npbparams.h'


c---------------------------------------------------------------------
c  num_procs must be a power of 2, and num_procs=num_proc_cols*num_proc_rows.
c  num_proc_cols and num_proc_cols are to be found in npbparams.h.
c  When num_procs is not square, then num_proc_cols must be = 2*num_proc_rows.
c---------------------------------------------------------------------
      integer    num_procs 
      pointer (ptr_num_procs, num_procs)
c      call DPMALLOC(ptr_num_procs, sizeof(num_procs), 0)
      parameter( num_procs = num_proc_cols * num_proc_rows )
      pointer (ptr_num_pro_cols, num_proc_cols)
      pointer (ptr_num_proc_rows, num_proc_rows)


c---------------------------------------------------------------------
c  Class specific parameters: 
c  It appears here for reference only.
c  These are their values, however, this info is imported in the npbparams.h
c  include file, which is written by the sys/setparams.c program.
c---------------------------------------------------------------------

C----------
C  Class S:
C----------
CC       parameter( na=1400, 
CC      >           nonzer=7, 
CC      >           shift=10., 
CC      >           niter=15,
CC      >           rcond=1.0d-1 )
C----------
C  Class W:
C----------
CC       parameter( na=7000,
CC      >           nonzer=8, 
CC      >           shift=12., 
CC      >           niter=15,
CC      >           rcond=1.0d-1 )
C----------
C  Class A:
C----------
CC       parameter( na=14000,
CC      >           nonzer=11, 
CC      >           shift=20., 
CC      >           niter=15,
CC      >           rcond=1.0d-1 )
C----------
C  Class B:
C----------
CC       parameter( na=75000, 
CC      >           nonzer=13, 
CC      >           shift=60., 
CC      >           niter=75,
CC      >           rcond=1.0d-1 )
C----------
C  Class C:
C----------
CC       parameter( na=150000, 
CC      >           nonzer=15, 
CC      >           shift=110., 
CC      >           niter=75,
CC      >           rcond=1.0d-1 )
C----------
C  Class D:
C----------
CC       parameter( na=1500000, 
CC      >           nonzer=21, 
CC      >           shift=500., 
CC      >           niter=100,
CC      >           rcond=1.0d-1 )
C----------
C  Class E:
C----------
CC       parameter( na=9000000, 
CC      >           nonzer=26, 
CC      >           shift=1500., 
CC      >           niter=100,
CC      >           rcond=1.0d-1 )
      


      integer    nz
      pointer (ptr_nz, nz)
      parameter( nz = na*(nonzer+1)/num_procs*(nonzer+1)+nonzer
     >              + na*(nonzer+2+num_procs/256)/num_proc_cols )

c------------------------------------------------------------------
c  convert common block variable to regular variable
c------------------------------------------------------------------      
c      common / partit_size  /  naa, nzz, 
c     >                         npcols, nprows,
c     >                         proc_col, proc_row,
c     >                         firstrow, 
c     >                         lastrow, 
c     >                         firstcol, 
c     >                         lastcol,
c     >                         exch_proc,
c     >                         exch_recv_length,
c     >                         send_start,
c     >                         send_len
      integer                  naa, nzz, 
     >                         npcols, nprows,
     >                         proc_col, proc_row,
     >                         firstrow, 
     >                         lastrow, 
     >                         firstcol, 
     >                         lastcol,
     >                         exch_proc,
     >                         exch_recv_length,
     >                         send_start,
     >                         send_len
c      pointer (ptr_naa, naa)
      pointer (ptr_nzz, nzz)
      pointer (ptr_npcols, npcols)
      pointer (ptr_nprows, nprows)
      pointer (prt_proc_col, proc_col)
      pointer (ptr_proc_row, proc_row)
      pointer (ptr_firstrow, firstrow)
      pointer (ptr_lastrow, lastrow)
      pointer (ptr_firstcol, firstcol)
      pointer (ptr_lastcol, lastcol)
      pointer (ptr_exch_proc, exch_proc)
      pointer (ptr_exch_recv_length, exch_recv_length)
      pointer (ptr_send_start, send_start)
      pointer (ptr_send_len, send_len)

c      common / main_int_mem /  colidx,     rowstr,
c     >                         iv,         arow,     acol
      integer                  colidx(nz), rowstr(na+1),
     >                         iv(2*na+1), arow(nz), acol(nz)

      pointer (ptr_colidx, colidx)
      pointer (ptr_rowstr, rowstr)
      pointer (ptr_iv, iv)
      pointer (ptr_arow, arow)
      pointer (ptr_acol, acol)

c      common / main_flt_mem /  v,       aelt,     a,
c     >                         x,
c     >                         z,
c     >                         p,
c     >                         q,
c     >                         r,
c     >                         w
      double precision         v(na+1), aelt(nz), a(nz),
     >                         x(na/num_proc_rows+2),
     >                         z(na/num_proc_rows+2),
     >                         p(na/num_proc_rows+2),
     >                         q(na/num_proc_rows+2),
     >                         r(na/num_proc_rows+2),
     >                         w(na/num_proc_rows+2),
     >                         p_o(na/num_proc_rows+2),
     >                         p_e(na/num_proc_rows+2),
     >                         r_o(na/num_proc_rows+2),
     >                         r_e(na/num_proc_rows+2),
     >                         z_o(na/num_proc_rows+2),
     >                         z_e(na/num_proc_rows+2)

      pointer (ptr_v, v)
      pointer (ptr_aelt, aelt)
      pointer (ptr_a, a)
      pointer (ptr_x, x)
      pointer (ptr_z, z)
      pointer (ptr_p, p)
      pointer (ptr_q, q)
      pointer (ptr_r, r)
      pointer (ptr_w, w)

      pointer (ptr_p_n, p_e)
      pointer (ptr_p_o, p_o)
      pointer (ptr_r_n, r_e)
      pointer (ptr_r_o, r_o)
      pointer (ptr_z_n, z_e)
      pointer (ptr_z_o, z_o)

c      common /urando/          amult, tran
      double precision         amult, tran
      pointer (ptr_amult, amult)
      pointer (ptr_tran, tran)
c---------------------------------------------------------------
      integer            l2npcols
      pointer (ptr_l2npcols, l2npcols)
      integer            reduce_exch_proc(num_proc_cols)
      pointer (ptr_reduce_exch_proc,reduce_exch_proc)
      integer            reduce_send_starts(num_proc_cols)
      pointer (ptr_reduce_send_starts, reduce_send_starts)
      integer            reduce_send_lengths(num_proc_cols)
      pointer (ptr_reduce_send_lengt, reduce_send_lengths)
      integer            reduce_recv_starts(num_proc_cols)
      pointer (ptr_reduce_recv_starts, reduce_recv_starts)
      integer            reduce_recv_lengths(num_proc_cols)
      pointer (ptr_reduce_recv_lengt, reduce_recv_lengths)

      integer            i, j, k, it
      pointer (ptr_i, i)
      pointer (ptr_j, j)
      pointer (ptr_k, k)
      pointer (ptr_it, it)
      double precision   zeta, randlc
      pointer (ptr_zeta, zeta)
c      pointer (ptr_randlc, randlc)
      external           randlc
      double precision   rnorm
      pointer (ptr_rnorm, rnorm)
      double precision   norm_temp1(2), norm_temp2(2)
      pointer (ptr_norm_temp1, norm_temp1)
      pointer (ptr_norm_temp2, norm_temp2)
      
      double precision   t, tmax, mflops
      pointer (ptr_t, t)
      pointer (ptr_tmax, tmax)
      pointer (ptr_mflops, mflops)
      external           timer_read
      double precision   timer_read
      character          class
      pointer (ptr_class, class)
      logical            verified
      pointer(ptr_verified, verified)
      double precision   zeta_verify_value, epsilon, err
      pointer (ptr_zeta_verify_value, zeta_verify_value)
      pointer (ptr_epsilon, epsilon)
      pointer (ptr_err, err)
      double precision tsum(t_last+2), t1(t_last+2),
     >                 tming(t_last+2), tmaxg(t_last+2)
      pointer (ptr_tsum, tsum)
      pointer (ptr_t1, t1)
      pointer (ptr_tming, tming)
      pointer (ptr_tmaxg, tmaxg)
      
      character        t_recs(t_last+2)*8
      pointer (ptr_t_recs, t_recs)

      data t_recs/'total', 'conjg', 'rcomm', 'ncomm',
     >            ' totcomp', ' totcomm'/



c      integer testdata
c      pointer (ptr_testdata, testdata)
c--------------------------------------------------------------------
c     kai
c--------------------------------------------------------------------
      integer main_loop
      integer  data_pos
      integer curr_rank
      common /myrank/ curr_rank
      
      integer copy_size
      common /mytest/ copy_size

      double precision   p_copy(na/num_proc_rows+2),
     >                   x_copy(na/num_proc_rows+2),
     >                   z_copy(na/num_proc_rows+2),
     >                   x_copy2(na/num_proc_rows+2),
     >                   z_copy2(na/num_proc_rows+2),
     >                   p_copy2(na/num_proc_rows+2)

      pointer (ptr_p_copy, p_copy)
      pointer (ptr_x_copy, x_copy)
      pointer (ptr_z_copy, z_copy)
      pointer (ptr_p_copy2, p_copy2)
      pointer (ptr_x_copy2, x_copy2)
      pointer (ptr_z_copy2, z_copy2)

c---------------------------------------------------------------------
c        set data init posistion
c---------------------------------------------------------------------
      copy_size = sizeof(p)
      data_pos = -1
      print *, data_pos
c---------------------------------------------------------------------
c Re-malloc via C function
c---------------------------------------------------------------------     
c      call c_openmp_init(4)

      call unimem_malloc(ptr_p_copy, sizeof(p_copy), data_pos)
      call unimem_malloc(ptr_x_copy, sizeof(x_copy), data_pos)
      call unimem_malloc(ptr_z_copy, sizeof(z_copy), data_pos)
      call unimem_malloc(ptr_p_copy2, sizeof(p_copy2), data_pos)
      call unimem_malloc(ptr_x_copy2, sizeof(x_copy2), data_pos)
      call unimem_malloc(ptr_z_copy2, sizeof(z_copy2), data_pos)

c     call DPMALLOC(ptr_num_procs, sizeof(num_procs), 0)
c      call DPMALLOC(ptr_num_pro_cols, sizeof(num_pro_cols), 0)
c      call DPMALLOC(ptr_num_pro_rows, sizeof(num_pro_rows), 0)
c      print *, sizeof(naa)
c      pointer (naa, ptr_naa)
c      PRINT *, "naa old:", loc(naa)
c      naa = 3
c      call PARTITSIZE(naa)
c      PRINT *, "naa new:", loc(naa)
c      PRINT *, naa;
c      print *, "old: ", nz
c      CALL DPMALLOC(ptr_testdata, sizeof(testdata), data_pos)
      CALL unimem_malloc(ptr_nz, sizeof(nz), data_pos)
c      print *, "new: ", nz
c      call DPMALLOC(ptr_naa, sizeof(naa),0)
c      call DPMALLOC(ptr_colidx, sizeof(colidx), 0)
c------------------------------------------------------------------------      
c      call DPMALLOC(ptr_naa, sizeof(naa), 0)
      call unimem_malloc(ptr_nzz, sizeof(nzz), data_pos)
      call unimem_malloc(ptr_npcols, sizeof(npcols), data_pos)
      call unimem_malloc(ptr_nprows, sizeof(nprows), data_pos)
      call unimem_malloc(prt_proc_col, sizeof(proc_col), data_pos)
      call unimem_malloc(ptr_proc_row, sizeof(proc_row), data_pos)
      call unimem_malloc(ptr_firstrow, sizeof(firstrow), data_pos)
      call unimem_malloc(ptr_lastrow, sizeof(lastrow), data_pos)
      call unimem_malloc(ptr_firstcol, sizeof(firstcol), data_pos)
      call unimem_malloc(ptr_lastcol, sizeof(lastcol), data_pos)
      call unimem_malloc(ptr_exch_proc, sizeof(exch_proc), data_pos)
      call unimem_malloc(ptr_exch_recv_length, sizeof(exch_recv_length), 
     >              data_pos)
      call unimem_malloc(ptr_send_start, sizeof(send_start), data_pos)
      call unimem_malloc(ptr_send_len, sizeof(send_len), data_pos)
      
  
      call unimem_malloc(ptr_colidx, sizeof(colidx), 1)
c      call unimem_malloc(ptr_colidx, sizeof(colidx), data_pos)
C      print *, "colidx:", sizeof(colidx)

c      call unimem_malloc(ptr_rowstr, sizeof(rowstr), 1)
      call unimem_malloc(ptr_rowstr, sizeof(rowstr), data_pos)

c      call unimem_malloc(ptr_iv, sizeof(iv), 1)
      call unimem_malloc(ptr_iv, sizeof(iv), data_pos)
      call unimem_malloc(ptr_arow, sizeof(arow), data_pos)
      call unimem_malloc(ptr_acol, sizeof(acol), data_pos)

c      call unimem_malloc(ptr_v, sizeof(v), 1)                                                                                             
      call unimem_malloc(ptr_v, sizeof(v), data_pos)
      call unimem_malloc(ptr_aelt, sizeof(aelt), data_pos)
      
c      print *,"1: ", loc(a)
c      print *,"1: ", ptr_a
      call unimem_malloc(ptr_a, sizeof(a), 1)
c      call unimem_malloc(ptr_a, sizeof(a), data_pos)
c      print *,"2: ", loc(a)
c      print *,"2: ", ptr_a

c      print *, "a:", sizeof(a)
      call unimem_malloc(ptr_x, sizeof(x), 1)
c      call unimem_malloc(ptr_x, sizeof(x), data_pos)      
      call unimem_malloc(ptr_z, sizeof(z), 1)
      call unimem_malloc(ptr_z_n, sizeof(z), 1)
       call unimem_malloc(ptr_z_o, sizeof(z), 1)
c      call unimem_malloc(ptr_z, sizeof(z), data_pos)
c      print *, "z:", sizeof(z)
      call unimem_malloc(ptr_p, sizeof(p), 1)
      call unimem_malloc(ptr_p_n, sizeof(p), 1)
      call unimem_malloc(ptr_p_o, sizeof(p), 1)
c      call unimem_malloc(ptr_p, sizeof(p), data_pos)
c      print *, "p:", sizeof(p)
c      call unimem_malloc(ptr_q, sizeof(q), 1)
      call unimem_malloc(ptr_q, sizeof(q), data_pos)
c      print *, "q:", sizeof(q)
c      call unimem_malloc(ptr_r, sizeof(r), 1)
      call unimem_malloc(ptr_r, sizeof(r), 1)
      call unimem_malloc(ptr_r_n, sizeof(r), 1)
      call unimem_malloc(ptr_r_o, sizeof(r), 1)

c      print *, "r:", sizeof(r)
c      call unimem_malloc(ptr_w, sizeof(w), 1)
      call unimem_malloc(ptr_w, sizeof(w), data_pos)
c      print *, "w:", sizeof(w)
      call unimem_malloc(ptr_amult, sizeof(amult), data_pos)
      call unimem_malloc(ptr_tran, sizeof(tran),data_pos)
c------------------------------------------------------------------------
      call unimem_malloc(ptr_l2npcols, sizeof(l2npcols), data_Pos)
      call unimem_malloc(ptr_reduce_exch_proc, sizeof(reduce_exch_proc),
     >              data_pos)
      call unimem_malloc(ptr_reduce_send_starts,
     >              sizeof(reduce_send_starts),
     >              data_pos)
      call unimem_malloc(ptr_reduce_send_lengt,
     >              sizeof(reduce_send_lengths),
     >              data_pos)
      call unimem_malloc(ptr_reduce_recv_starts,
     >              sizeof(reduce_recv_starts),
     >              data_pos)
      call unimem_malloc(ptr_reduce_recv_lengt,
     >              sizeof(reduce_recv_lengths),
     >              data_pos)
      call unimem_malloc(ptr_i, sizeof(i), data_pos)
      call unimem_malloc(ptr_j, sizeof(j), data_pos)

      call unimem_malloc(ptr_k, sizeof(k),data_pos)
      
      call unimem_malloc(ptr_it, sizeof(it), data_pos)
      call unimem_malloc(ptr_zeta, sizeof(zeta), data_pos)
c      call DPMALLOC(ptr_randlc, sizeof(randlc), 0)
      call unimem_malloc(ptr_rnorm, sizeof(rnorm), data_pos)
      call unimem_malloc(ptr_norm_temp1, sizeof(norm_temp1), data_pos)
      call unimem_malloc(ptr_norm_temp2, sizeof(norm_temp2), data_pos)
      call unimem_malloc(ptr_t, sizeof(t), data_pos)
      call unimem_malloc(ptr_tmax, sizeof(tmax), data_pos)
      call unimem_malloc(ptr_mflops, sizeof(mflops), data_pos)
      call unimem_malloc(ptr_class, sizeof(class), data_pos)
      call unimem_malloc(ptr_verified, sizeof(verified), data_pos)
      call unimem_malloc(ptr_zeta_verify_value, 
     >              sizeof(zeta_verify_value), 
     >              data_pos)
      call unimem_malloc(ptr_epsilon, sizeof(epsilon), data_pos)
      call unimem_malloc(ptr_err, sizeof(err), data_pos)
      call unimem_malloc(ptr_tsum, sizeof(tsum), data_pos)
      call unimem_malloc(ptr_t1, sizeof(t1), data_pos)
      call unimem_malloc(ptr_tming, sizeof(tming), data_pos)
      call unimem_malloc(ptr_tmaxg, sizeof(tmaxg), data_pos)
      call unimem_malloc(ptr_t_recs, sizeof(t_recs), data_pos)

c-------------------------------------------------------------------
      print*, "colidx: add:",loc(colidx),"size:",sizeof(colidx)
      print*, "a: add:",loc(a),"size:",sizeof(a)
      print*, "w: add:",loc(w),"size:",sizeof(w)
      print*, "z: add:",loc(z),"size:",sizeof(z)
      print*, "p: add:",loc(p),"size:",sizeof(p)
      print*, "r: add:",loc(r),"size:",sizeof(r)
      print*, "q: add:",loc(q),"size:",sizeof(q)
      print*, "rowst: size:", sizeof(rowstr)
      print*, "x size:", sizeof(x)
      print*, "v size:", sizeof(v)
      print*, "aelt size:", sizeof(aelt)
      print*, "acol size:", sizeof(acol)
      print*, "arow size:", sizeof(arow)
      print*, "iv size:", sizeof(iv)
c---------------------------------------------------------------------                                                                            
c  Set up mpi initialization and number of proc testing                                                                                           
c---------------------------------------------------------------------       
      call initialize_mpi

c---------------kai------------------------------------
c      call c_lib_init(curr_rank)
c------------------------------------------------------
c      CALL unimem_init_cg(ptr_colidx, loc(colidx), sizeof(colidx),
c     >                 ptr_a, loc(a), sizeof(a),
c     >                 ptr_w, loc(w), sizeof(w),
c     >                 ptr_z, loc(z), sizeof(z),
c     >                 ptr_p, loc(p), sizeof(p),
c     >                 ptr_q, loc(q), sizeof(q),
c     >                 ptr_r, loc(r), sizeof(r),
c     >                 data_pos)

c      CALL unimem_init(ptr_colidx, sizeof(colidx),
c     >                 ptr_a, sizeof(a),
c     >                 ptr_w, sizeof(w),
c     >                 ptr_z, sizeof(z),
c     >                 ptr_p, sizeof(p),
c     >                 ptr_q, sizeof(q),
c     >                 ptr_r, sizeof(r))

c---------------------------------------------------------------------
      if( na .eq. 1400 .and. 
     &    nonzer .eq. 7 .and. 
     &    niter .eq. 15 .and.
     &    shift .eq. 10.d0 ) then
         class = 'S'
         zeta_verify_value = 8.5971775078648d0
      else if( na .eq. 7000 .and. 
     &         nonzer .eq. 8 .and. 
     &         niter .eq. 15 .and.
     &         shift .eq. 12.d0 ) then
         class = 'W'
         zeta_verify_value = 10.362595087124d0
      else if( na .eq. 14000 .and. 
     &         nonzer .eq. 11 .and. 
     &         niter .eq. 15 .and.
     &         shift .eq. 20.d0 ) then
         class = 'A'
         zeta_verify_value = 17.130235054029d0
      else if( na .eq. 75000 .and. 
     &         nonzer .eq. 13 .and. 
     &         niter .eq. 75 .and.
     &         shift .eq. 60.d0 ) then
         class = 'B'
         zeta_verify_value = 22.712745482631d0
      else if( na .eq. 150000 .and. 
     &         nonzer .eq. 15 .and. 
     &         niter .eq. 75 .and.
     &         shift .eq. 110.d0 ) then
         class = 'C'
         zeta_verify_value = 28.973605592845d0
      else if( na .eq. 1500000 .and. 
     &         nonzer .eq. 21 .and. 
     &         niter .eq. 100 .and.
     &         shift .eq. 500.d0 ) then
         class = 'D'
         zeta_verify_value = 52.514532105794d0
      else if( na .eq. 9000000 .and. 
     &         nonzer .eq. 26 .and. 
     &         niter .eq. 100 .and.
     &         shift .eq. 1.5d3 ) then
         class = 'E'
         zeta_verify_value = 77.522164599383d0
      else
         class = 'U'
      endif

      if( me .eq. root )then
         write( *,1000 ) 
         write( *,1001 ) na
         write( *,1002 ) niter
         write( *,1003 ) nprocs
         write( *,1004 ) nonzer
         write( *,1005 ) shift
 1000 format(//,' NAS Parallel Benchmarks 3.3 -- CG Benchmark', /)
 1001 format(' Size: ', i10 )
 1002 format(' Iterations: ', i5 )
 1003 format(' Number of active processes: ', i5 )
 1004 format(' Number of nonzeroes per row: ', i8)
 1005 format(' Eigenvalue shift: ', e8.3)
      endif

      if (.not. convertdouble) then
         dp_type = MPI_DOUBLE_PRECISION
      else
         dp_type = MPI_REAL
      endif


      naa = na
      nzz = nz
c      print *, "nz:", loc(nz)

c---------------------------------------------------------------------
c  Set up processor info, such as whether sq num of procs, etc
c---------------------------------------------------------------------
      call setup_proc_info( num_procs, 
     >                      num_proc_rows, 
     >                      num_proc_cols,
     >                      naa,nzz, npcols, nprows,
     >                      proc_col, proc_row, firstrow, lastrow,
     >                      exch_proc, exch_recv_length,
     >                      send_start, send_len)


c---------------------------------------------------------------------
c  Set up partition's submatrix info: firstcol, lastcol, firstrow, lastrow
c---------------------------------------------------------------------
      call setup_submatrix_info( l2npcols,
     >                           reduce_exch_proc,
     >                           reduce_send_starts,
     >                           reduce_send_lengths,
     >                           reduce_recv_starts,
     >                           reduce_recv_lengths,
     >                           naa, nzz, npcols, nprows,
     >                           proc_col, proc_row,firstrow,
     >                           lastrow, firstcol, lastcol,
     >                           exch_proc, exch_recv_length,
     >                           send_start, send_len )


      do i = 1, t_last
         call timer_clear(i)
      end do

c---------------------------------------------------------------------
c  Inialize random number generator
c---------------------------------------------------------------------
      tran    = 314159265.0D0
      amult   = 1220703125.0D0
      zeta    = randlc( tran, amult )

c---------------------------------------------------------------------
c  Set up partition's sparse random matrix for given class size
c---------------------------------------------------------------------
      call makea(naa, nzz, a, colidx, rowstr, nonzer,
     >           firstrow, lastrow, firstcol, lastcol, 
     >           rcond, arow, acol, aelt, v, iv, shift,amult, tran)



c---------------------------------------------------------------------
c  Note: as a result of the above call to makea:
c        values of j used in indexing rowstr go from 1 --> lastrow-firstrow+1
c        values of colidx which are col indexes go from firstcol --> lastcol
c        So:
c        Shift the col index vals from actual (firstcol --> lastcol ) 
c        to local, i.e., (1 --> lastcol-firstcol+1)
c---------------------------------------------------------------------
      do j=1,lastrow-firstrow+1
         do k=rowstr(j),rowstr(j+1)-1
            colidx(k) = colidx(k) - firstcol + 1
         enddo
      enddo

c---------------------------------------------------------------------
c  set starting vector to (1, 1, .... 1)
c---------------------------------------------------------------------
      do i = 1, na/num_proc_rows+1
         x(i) = 1.0D0
      enddo

      zeta  = 0.0d0

c---------------------------------------------------------------------
c---->
c  Do one iteration untimed to init all code and data page tables
c---->                    (then reinit, start timing, to niter its)
c---------------------------------------------------------------------
      main_loop = 0
      do it = 1, 1

c---------------------------------------------------------------------
c  The call to the conjugate gradient routine:
c---------------------------------------------------------------------
         call conj_grad ( 
     >  ptr_p_n,ptr_p_o, ptr_r_n,ptr_r_o,ptr_z_n,ptr_z_o,   
     >                    ptr_colidx,
     >                    rowstr,
     >                    ptr_x,
     >                    ptr_z,
     >                    ptr_a,
     >                    ptr_p,
     >                    ptr_q,
     >                    ptr_r,
     >                    ptr_w,
     >                    rnorm, 
     >                    l2npcols,
     >                    reduce_exch_proc,
     >                    reduce_send_starts,
     >                    reduce_send_lengths,
     >                    reduce_recv_starts,
     >                    reduce_recv_lengths,
     >                    naa, nzz,
     >                    npcols, nprows,
     >                    proc_col, proc_row,
     >                    firstrow,
     >                    lastrow,
     >                    firstcol,
     >                    lastcol,
     >                    exch_proc,
     >                    exch_recv_length,
     >                    send_start,
     >                    send_len,
     >                    main_loop,
     >                    ptr_p_copy,
     >                    ptr_x_copy,
     >                    ptr_z_copy,
     >                       ptr_p_copy2,
     >                       ptr_x_copy2,
     >                       ptr_z_copy2)

c---------------------------------------------------------------------
c  zeta = shift + 1/(x.z)
c  So, first: (x.z)
c  Also, find norm of z
c  So, first: (z.z)
c---------------------------------------------------------------------
         norm_temp1(1) = 0.0d0
         norm_temp1(2) = 0.0d0
         do j=1, lastcol-firstcol+1
            norm_temp1(1) = norm_temp1(1) + x(j)*z(j)
            norm_temp1(2) = norm_temp1(2) + z(j)*z(j)
         enddo

         do i = 1, l2npcols
            if (timeron) call timer_start(t_ncomm)
            call mpi_irecv( norm_temp2,
     >                      2, 
     >                      dp_type,
     >                      reduce_exch_proc(i),
     >                      i,
     >                      mpi_comm_world,
     >                      request,
     >                      ierr )
            call mpi_send(  norm_temp1,
     >                      2, 
     >                      dp_type,
     >                      reduce_exch_proc(i),
     >                      i,
     >                      mpi_comm_world,
     >                      ierr )
            call mpi_wait( request, status, ierr )
            if (timeron) call timer_stop(t_ncomm)

            norm_temp1(1) = norm_temp1(1) + norm_temp2(1)
            norm_temp1(2) = norm_temp1(2) + norm_temp2(2)
         enddo

         norm_temp1(2) = 1.0d0 / sqrt( norm_temp1(2) )


c---------------------------------------------------------------------
c  Normalize z to obtain x
c---------------------------------------------------------------------
         do j=1, lastcol-firstcol+1      
            x(j) = norm_temp1(2)*z(j)    
         enddo                           


      enddo                              ! end of do one iteration untimed


c---------------------------------------------------------------------
c  set starting vector to (1, 1, .... 1)
c---------------------------------------------------------------------
c
c  NOTE: a questionable limit on size:  should this be na/num_proc_cols+1 ?
c
      do i = 1, na/num_proc_rows+1
         x(i) = 1.0D0
      enddo

      zeta  = 0.0d0

c---------------------------------------------------------------------
c  Synchronize and start timing
c---------------------------------------------------------------------
      do i = 1, t_last
         call timer_clear(i)
      end do
      call mpi_barrier( mpi_comm_world,
     >                  ierr )

      call timer_clear( 1 )
      call timer_start( 1 )

c---------------------------------------------------------------------
c---->
c  Main Iteration for inverse power method
c---->
c---------------------------------------------------------------------
      main_loop = 1
      do it = 1, niter
c      do it = 1, 10
c---------------------------------------------------------------------
c  The call to the conjugate gradient routine:
c---------------------------------------------------------------------
         call conj_grad ( 
     >     ptr_p_n,ptr_p_o, ptr_r_n,ptr_r_o,ptr_z_n,ptr_z_o,
     >                    ptr_colidx,
     >                    rowstr,
     >                    ptr_x,
     >                    ptr_z,
     >                    ptr_a,
     >                    ptr_p,
     >                    ptr_q,
     >                    ptr_r,
     >                    ptr_w,
     >                    rnorm, 
     >                    l2npcols,
     >                    reduce_exch_proc,
     >                    reduce_send_starts,
     >                    reduce_send_lengths,
     >                    reduce_recv_starts,
     >                    reduce_recv_lengths,
     >                    naa, nzz,
     >                    npcols, nprows,
     >                    proc_col, proc_row,
     >                    firstrow,
     >                    lastrow,
     >                    firstcol,
     >                    lastcol,
     >                    exch_proc,
     >                    exch_recv_length,
     >                    send_start,
     >                    send_len,
     >                    main_loop,
     >                    ptr_p_copy,
     >                    ptr_x_copy,
     >                    ptr_z_copy,
     >                       ptr_p_copy2,
     >                       ptr_x_copy2,
     >                       ptr_z_copy2)


c---------------------------------------------------------------------
c  zeta = shift + 1/(x.z)
c  So, first: (x.z)
c  Also, find norm of z
c  So, first: (z.z)
c---------------------------------------------------------------------
         norm_temp1(1) = 0.0d0
         norm_temp1(2) = 0.0d0
         do j=1, lastcol-firstcol+1
            norm_temp1(1) = norm_temp1(1) + x(j)*z(j)
            norm_temp1(2) = norm_temp1(2) + z(j)*z(j)
         enddo

         do i = 1, l2npcols
            if (timeron) call timer_start(t_ncomm)
            call mpi_irecv( norm_temp2,
     >                      2, 
     >                      dp_type,
     >                      reduce_exch_proc(i),
     >                      i,
     >                      mpi_comm_world,
     >                      request,
     >                      ierr )
            call mpi_send(  norm_temp1,
     >                      2, 
     >                      dp_type,
     >                      reduce_exch_proc(i),
     >                      i,
     >                      mpi_comm_world,
     >                      ierr )
            call mpi_wait( request, status, ierr )
            if (timeron) call timer_stop(t_ncomm)

            norm_temp1(1) = norm_temp1(1) + norm_temp2(1)
            norm_temp1(2) = norm_temp1(2) + norm_temp2(2)
         enddo

         norm_temp1(2) = 1.0d0 / sqrt( norm_temp1(2) )


         if( me .eq. root )then
            zeta = shift + 1.0d0 / norm_temp1(1)
            if( it .eq. 1 ) write( *,9000 )
            write( *,9001 ) it, rnorm, zeta
         endif
 9000 format( /,'   iteration           ||r||                 zeta' )
 9001 format( 4x, i5, 7x, e20.14, f20.13 )

c---------------------------------------------------------------------
c  Normalize z to obtain x
c---------------------------------------------------------------------
         do j=1, lastcol-firstcol+1      
            x(j) = norm_temp1(2)*z(j)    
         enddo                           


      enddo                              ! end of main iter inv pow meth

      call timer_stop( 1 )
c-----------------kai-----------------------------                                                                                                
      call c_print_timer(curr_rank)
c------------------------------------------------- 
c---------------------------------------------------------------------
c  End of timed section
c---------------------------------------------------------------------

      t = timer_read( 1 )

      call mpi_reduce( t,
     >                 tmax,
     >                 1, 
     >                 dp_type,
     >                 MPI_MAX,
     >                 root,
     >                 mpi_comm_world,
     >                 ierr )

      if( me .eq. root )then
         write(*,100)
 100     format(' Benchmark completed ')

         epsilon = 1.d-10
         if (class .ne. 'U') then

            err = abs( zeta - zeta_verify_value )/zeta_verify_value
            if( err .le. epsilon ) then
               verified = .TRUE.
               write(*, 200)
               write(*, 201) zeta
               write(*, 202) err
 200           format(' VERIFICATION SUCCESSFUL ')
 201           format(' Zeta is    ', E20.13)
 202           format(' Error is   ', E20.13)
            else
               verified = .FALSE.
               write(*, 300) 
               write(*, 301) zeta
               write(*, 302) zeta_verify_value
 300           format(' VERIFICATION FAILED')
 301           format(' Zeta                ', E20.13)
 302           format(' The correct zeta is ', E20.13)
            endif
         else
            verified = .FALSE.
            write (*, 400)
            write (*, 401)
            write (*, 201) zeta
 400        format(' Problem size unknown')
 401        format(' NO VERIFICATION PERFORMED')
         endif


         if( tmax .ne. 0. ) then
            mflops = float( 2*niter*na )
     &                  * ( 3.+float( nonzer*(nonzer+1) )
     &                    + 25.*(5.+float( nonzer*(nonzer+1) ))
     &                    + 3. ) / tmax / 1000000.0
         else
            mflops = 0.0
         endif

c-----------------kai-----------------------------                                                                                                
c      call c_print_timer(curr_rank)
c------------------------------------------------- 
         call print_results('CG', class, na, 0, 0,
     >                      niter, nnodes_compiled, nprocs, tmax,
     >                      mflops, '          floating point', 
     >                      verified, npbversion, compiletime,
     >                      cs1, cs2, cs3, cs4, cs5, cs6, cs7)


      endif


      if (.not.timeron) goto 999

      do i = 1, t_last
         t1(i) = timer_read(i)
      end do
      t1(t_conjg) = t1(t_conjg) - t1(t_rcomm)
      t1(t_last+2) = t1(t_rcomm) + t1(t_ncomm)
      t1(t_last+1) = t1(t_total) - t1(t_last+2)

      call MPI_Reduce(t1, tsum,  t_last+2, dp_type, MPI_SUM, 
     >                0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(t1, tming, t_last+2, dp_type, MPI_MIN, 
     >                0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(t1, tmaxg, t_last+2, dp_type, MPI_MAX, 
     >                0, MPI_COMM_WORLD, ierr)

      if (me .eq. 0) then
         write(*, 800) nprocs
         do i = 1, t_last+2
            tsum(i) = tsum(i) / nprocs
            write(*, 810) i, t_recs(i), tming(i), tmaxg(i), tsum(i)
         end do
      endif

 800  format(' nprocs =', i6, 11x, 'minimum', 5x, 'maximum', 
     >       5x, 'average')
 810  format(' timer ', i2, '(', A8, ') :', 3(2x,f10.4))

 999  continue
      call mpi_finalize(ierr)




      end                              ! end main





c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine initialize_mpi
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none

      include 'mpinpb.h'
      include 'timing.h'

      integer   ierr, fstatus
c      pointer (ptr_ierr, ierr)
c      pointer (ptr_fstatus, fstatus)
c      call DPMALLOC(ptr_ierr, sizeof(ierr), 0)
c      call DPMALLOC(ptr_fstatus, sizeof(fstatus),0)
c---------------------------kai--------------------
      integer curr_rank
      common /myrank/ curr_rank
c-------------------------------------------------
      call mpi_init( ierr )
      call mpi_comm_rank( mpi_comm_world, me, ierr )
      call mpi_comm_size( mpi_comm_world, nprocs, ierr )
      root = 0
      
      curr_rank = me

      if (me .eq. root) then
         open (unit=2,file='timer.flag',status='old',iostat=fstatus)
         timeron = .false.
         if (fstatus .eq. 0) then
            timeron = .true.
            close(2)
         endif
      endif

      call mpi_bcast(timeron, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)

      return
      end



c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine setup_proc_info( num_procs, 
     >                            num_proc_rows, 
     >                            num_proc_cols,
     >                            naa,nzz, npcols, nprows,
     >                            proc_col, proc_row, firstrow, lastrow,
     >                            exch_proc, exch_recv_length, 
     >                            send_start, send_len)
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none

      include 'mpinpb.h'

c      common / partit_size  /  naa, nzz, 
c     >                         npcols, nprows,
c     >                         proc_col, proc_row,
c     >                         firstrow, 
c     >                         lastrow, 
c     >                         firstcol, 
c     >                         lastcol,
c     >                         exch_proc,
c     >                         exch_recv_length,
c     >                         send_start,
c     >                         send_len
      integer                  naa, nzz, 
     >                         npcols, nprows,
     >                         proc_col, proc_row,
     >                         firstrow, 
     >                         lastrow, 
     >                         firstcol, 
     >                         lastcol,
     >                         exch_proc,
     >                         exch_recv_length,
     >                         send_start,
     >                         send_len

      integer   num_procs, num_proc_cols, num_proc_rows
      integer   i, ierr
      integer   log2nprocs

c---------------------------------------------------------------------
c  num_procs must be a power of 2, and num_procs=num_proc_cols*num_proc_rows
c  When num_procs is not square, then num_proc_cols = 2*num_proc_rows
c---------------------------------------------------------------------
c  First, number of procs must be power of two. 
c---------------------------------------------------------------------
      if( nprocs .ne. num_procs )then
          if( me .eq. root ) write( *,9000 ) nprocs, num_procs
 9000     format(      /,'Error: ',/,'num of procs allocated   (', 
     >                 i4, ' )',
     >                 /,'is not equal to',/,
     >                 'compiled number of procs (',
     >                 i4, ' )',/   )
          call mpi_finalize(ierr)
          stop
      endif


      i = num_proc_cols
 100  continue
          if( i .ne. 1 .and. i/2*2 .ne. i )then
              if ( me .eq. root ) then  
                 write( *,* ) 'Error: num_proc_cols is ',
     >                         num_proc_cols,
     >                        ' which is not a power of two'
              endif
              call mpi_finalize(ierr)
              stop
          endif
          i = i / 2
          if( i .ne. 0 )then
              goto 100
          endif
      
      i = num_proc_rows
 200  continue
          if( i .ne. 1 .and. i/2*2 .ne. i )then
              if ( me .eq. root ) then 
                 write( *,* ) 'Error: num_proc_rows is ',
     >                         num_proc_rows,
     >                        ' which is not a power of two'
              endif
              call mpi_finalize(ierr)
              stop
          endif
          i = i / 2
          if( i .ne. 0 )then
              goto 200
          endif
      
      log2nprocs = 0
      i = nprocs
 300  continue
          if( i .ne. 1 .and. i/2*2 .ne. i )then
              write( *,* ) 'Error: nprocs is ',
     >                      nprocs,
     >                      ' which is not a power of two'
              call mpi_finalize(ierr)
              stop
          endif
          i = i / 2
          if( i .ne. 0 )then
              log2nprocs = log2nprocs + 1
              goto 300
          endif

CC       write( *,* ) 'nprocs, log2nprocs: ',nprocs,log2nprocs

      
      npcols = num_proc_cols
      nprows = num_proc_rows


      return
      end




c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine setup_submatrix_info( l2npcols,
     >                                 reduce_exch_proc,
     >                                 reduce_send_starts,
     >                                 reduce_send_lengths,
     >                                 reduce_recv_starts,
     >                                 reduce_recv_lengths, 
     >                                 naa, nzz, npcols, nprows, 
     >                                 proc_col, proc_row,firstrow,
     >                                 lastrow, firstcol, lastcol,
     >                                 exch_proc, exch_recv_length,
     >                                 send_start, send_len)                                 
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none

      include 'mpinpb.h'

      integer      col_size, row_size

c      common / partit_size  /  naa, nzz, 
c     >                         npcols, nprows,
c     >                         proc_col, proc_row,
c     >                         firstrow, 
c     >                         lastrow, 
c     >                         firstcol, 
c     >                         lastcol,
c     >                         exch_proc,
c     >                         exch_recv_length,
c     >                         send_start,
c     >                         send_len
      integer                  naa, nzz, 
     >                         npcols, nprows,
     >                         proc_col, proc_row,
     >                         firstrow, 
     >                         lastrow, 
     >                         firstcol, 
     >                         lastcol,
     >                         exch_proc,
     >                         exch_recv_length,
     >                         send_start,
     >                         send_len

      integer   reduce_exch_proc(*)
      integer   reduce_send_starts(*)
      integer   reduce_send_lengths(*)
      integer   reduce_recv_starts(*)
      integer   reduce_recv_lengths(*)

      integer   i, j
      integer   div_factor
      integer   l2npcols


      proc_row = me / npcols
      proc_col = me - proc_row*npcols



c---------------------------------------------------------------------
c  If naa evenly divisible by npcols, then it is evenly divisible 
c  by nprows 
c---------------------------------------------------------------------

      if( naa/npcols*npcols .eq. naa )then
          col_size = naa/npcols
          firstcol = proc_col*col_size + 1
          lastcol  = firstcol - 1 + col_size
          row_size = naa/nprows
          firstrow = proc_row*row_size + 1
          lastrow  = firstrow - 1 + row_size
c---------------------------------------------------------------------
c  If naa not evenly divisible by npcols, then first subdivide for nprows
c  and then, if npcols not equal to nprows (i.e., not a sq number of procs), 
c  get col subdivisions by dividing by 2 each row subdivision.
c---------------------------------------------------------------------
      else
          if( proc_row .lt. naa - naa/nprows*nprows)then
              row_size = naa/nprows+ 1
              firstrow = proc_row*row_size + 1
              lastrow  = firstrow - 1 + row_size
          else
              row_size = naa/nprows
              firstrow = (naa - naa/nprows*nprows)*(row_size+1)
     >                 + (proc_row-(naa-naa/nprows*nprows))
     >                     *row_size + 1
              lastrow  = firstrow - 1 + row_size
          endif
          if( npcols .eq. nprows )then
              if( proc_col .lt. naa - naa/npcols*npcols )then
                  col_size = naa/npcols+ 1
                  firstcol = proc_col*col_size + 1
                  lastcol  = firstcol - 1 + col_size
              else
                  col_size = naa/npcols
                  firstcol = (naa - naa/npcols*npcols)*(col_size+1)
     >                     + (proc_col-(naa-naa/npcols*npcols))
     >                         *col_size + 1
                  lastcol  = firstcol - 1 + col_size
              endif
          else
              if( (proc_col/2) .lt. 
     >                           naa - naa/(npcols/2)*(npcols/2) )then
                  col_size = naa/(npcols/2) + 1
                  firstcol = (proc_col/2)*col_size + 1
                  lastcol  = firstcol - 1 + col_size
              else
                  col_size = naa/(npcols/2)
                  firstcol = (naa - naa/(npcols/2)*(npcols/2))
     >                                                 *(col_size+1)
     >               + ((proc_col/2)-(naa-naa/(npcols/2)*(npcols/2)))
     >                         *col_size + 1
                  lastcol  = firstcol - 1 + col_size
              endif
CC               write( *,* ) col_size,firstcol,lastcol
              if( mod( me,2 ) .eq. 0 )then
                  lastcol  = firstcol - 1 + (col_size-1)/2 + 1
              else
                  firstcol = firstcol + (col_size-1)/2 + 1
                  lastcol  = firstcol - 1 + col_size/2
CC                   write( *,* ) firstcol,lastcol
              endif
          endif
      endif



      if( npcols .eq. nprows )then
          send_start = 1
          send_len   = lastrow - firstrow + 1
      else
          if( mod( me,2 ) .eq. 0 )then
              send_start = 1
              send_len   = (1 + lastrow-firstrow+1)/2
          else
              send_start = (1 + lastrow-firstrow+1)/2 + 1
              send_len   = (lastrow-firstrow+1)/2
          endif
      endif
          



c---------------------------------------------------------------------
c  Transpose exchange processor
c---------------------------------------------------------------------

      if( npcols .eq. nprows )then
          exch_proc = mod( me,nprows )*nprows + me/nprows
      else
          exch_proc = 2*(mod( me/2,nprows )*nprows + me/2/nprows)
     >                 + mod( me,2 )
      endif



      i = npcols / 2
      l2npcols = 0
      do while( i .gt. 0 )
         l2npcols = l2npcols + 1
         i = i / 2
      enddo


c---------------------------------------------------------------------
c  Set up the reduce phase schedules...
c---------------------------------------------------------------------

      div_factor = npcols
      do i = 1, l2npcols

         j = mod( proc_col+div_factor/2, div_factor )
     >     + proc_col / div_factor * div_factor
         reduce_exch_proc(i) = proc_row*npcols + j

         div_factor = div_factor / 2

      enddo


      do i = l2npcols, 1, -1

            if( nprows .eq. npcols )then
               reduce_send_starts(i)  = send_start
               reduce_send_lengths(i) = send_len
               reduce_recv_lengths(i) = lastrow - firstrow + 1
            else
               reduce_recv_lengths(i) = send_len
               if( i .eq. l2npcols )then
                  reduce_send_lengths(i) = lastrow-firstrow+1 - send_len
                  if( me/2*2 .eq. me )then
                     reduce_send_starts(i) = send_start + send_len
                  else
                     reduce_send_starts(i) = 1
                  endif
               else
                  reduce_send_lengths(i) = send_len
                  reduce_send_starts(i)  = send_start
               endif
            endif
            reduce_recv_starts(i) = send_start

      enddo


      exch_recv_length = lastcol - firstcol + 1


      return
      end




c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine conj_grad (
     > ptr_p_n,ptr_p_o, ptr_r_n,ptr_r_o,ptr_z_n,ptr_z_o,      
     >                       ptr_colidx,
     >                       rowstr,
     >                       ptr_x,
     >                       ptr_z,
     >                       ptr_a,
     >                       ptr_p,
     >                       ptr_q,
     >                       ptr_r,
     >                       ptr_w,
     >                       rnorm, 
     >                       l2npcols,
     >                       reduce_exch_proc,
     >                       reduce_send_starts,
     >                       reduce_send_lengths,
     >                       reduce_recv_starts,
     >                       reduce_recv_lengths, 
     >                       naa, nzz, 
     >                       npcols, nprows,
     >                       proc_col, proc_row,
     >                       firstrow, 
     >                       lastrow, 
     >                       firstcol, 
     >                       lastcol,
     >                       exch_proc,
     >                       exch_recv_length,
     >                       send_start,
     >                       send_len,
     >                       main_loop,
     >                       ptr_p_copy,
     >                       ptr_x_copy,
     >                       ptr_z_copy,
     >                       ptr_p_copy2,
     >                       ptr_x_copy2,
     >                       ptr_z_copy2)
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c  Floaging point arrays here are named as in NPB1 spec discussion of 
c  CG algorithm
c---------------------------------------------------------------------
 
      implicit none

      include 'mpinpb.h'
      include 'timing.h'

      integer status(MPI_STATUS_SIZE ), request


c      common / partit_size  /  naa, nzz, 
c     >                         npcols, nprows,
c     >                         proc_col, proc_row,
c     >                         firstrow, 
c     >                         lastrow, 
c     >                         firstcol, 
c     >                         lastcol,
c     >                         exch_proc,
c     >                         exch_recv_length,
c     >                         send_start,
c     >                         send_len
      integer                  naa, nzz, 
     >                         npcols, nprows,
     >                         proc_col, proc_row,
     >                         firstrow, 
     >                         lastrow, 
     >                         firstcol, 
     >                         lastcol,
     >                         exch_proc,
     >                         exch_recv_length,
     >                         send_start,
     >                         send_len



      double precision   x(*),
     >                   p(*),
     >                   r(*),
     >                   z(*), 
     >                   a(nzz)
      integer            colidx(nzz), rowstr(naa+1)

      double precision   
     >                   p_o(*),
     >                   p_e(*),
     >                   r_o(*),
     >                   r_e(*),
     >                   z_o(*),
     >                   z_e(*),
     >                   q(*),
     >                   w(*)                ! used as work temporary
      
      double precision   p_copy(*),
     >                   x_copy(*),
     >                   z_copy(*),
     >                   p_copy2(*),
     >                   x_copy2(*),
     >                   z_copy2(*)
      
c----------------------------------------------------------------
      pointer (ptr_a, a)
      pointer (ptr_colidx, colidx)
      pointer (ptr_q, q)
c      pointer (ptr_k, k)
      pointer (ptr_x, x)
      pointer (ptr_w, w)
      pointer (ptr_p,p)
      pointer (ptr_r,r)
      pointer (ptr_z,z)
      
      pointer (ptr_p_n, p_e)
      pointer (ptr_p_o, p_o)
      pointer (ptr_r_n, r_e)      
      pointer (ptr_r_o, r_o)
      pointer (ptr_z_n, z_e)
      pointer (ptr_z_o, z_o)  

      pointer (ptr_p_copy, p_copy)
      pointer (ptr_x_copy, x_copy)
      pointer (ptr_z_copy, z_copy)
      pointer (ptr_p_copy2, p_copy2)
      pointer (ptr_x_copy2, x_copy2)
      pointer (ptr_z_copy2, z_copy2)

c----------------------------------------------------------------
      integer   l2npcols
      integer   reduce_exch_proc(l2npcols)
      integer   reduce_send_starts(l2npcols)
      integer   reduce_send_lengths(l2npcols)
      integer   reduce_recv_starts(l2npcols)
      integer   reduce_recv_lengths(l2npcols)

      integer   i, j, k, ierr
      integer   cgit, cgitmax

      double precision   d, sum, rho, rho0, alpha, beta, rnorm

      external         timer_read
      double precision timer_read

      data      cgitmax / 25 /
c--------------------------------kai----------------------------------
      integer main_loop
      
      integer curr_rank
      common /myrank/ curr_rank
      integer copy_size
      common /mytest/ copy_size
      
c---------------------------------------------------------------------
c      print *,"3: ", loc(a)
c      print *,"3: ", ptr_a
      

c---------------------------------------------------------------------

      if (timeron) call timer_start(t_conjg)
c---------------------------------------------------------------------
c  Initialize the CG algorithm:
c---------------------------------------------------------------------
      do j=1,naa/nprows+1
         q(j) = 0.0d0
         z_e(j) = 0.0d0
         z_o(j) = 0.0d0
         r_o(j) = x(j)
         r_e(j) = x(j)
         p_o(j) = x(j)
         p_e(j) = x(j)
         w(j) = 0.0d0                 
      enddo


c---------------------------------------------------------------------
c  rho = r.r
c  Now, obtain the norm of r: First, sum squares of r elements locally...
c---------------------------------------------------------------------
      sum = 0.0d0
      do j=1, lastcol-firstcol+1
         sum = sum + r_e(j)*r_e(j)
      enddo

c---------------------------------------------------------------------
c  Exchange and sum with procs identified in reduce_exch_proc
c  (This is equivalent to mpi_allreduce.)
c  Sum the partial sums of rho, leaving rho on all processors
c---------------------------------------------------------------------
      do i = 1, l2npcols
         if (timeron) call timer_start(t_rcomm)
         call mpi_irecv( rho,
     >                   1,
     >                   dp_type,
     >                   reduce_exch_proc(i),
     >                   i,
     >                   mpi_comm_world,
     >                   request,
     >                   ierr )
         call mpi_send(  sum,
     >                   1,
     >                   dp_type,
     >                   reduce_exch_proc(i),
     >                   i,
     >                   mpi_comm_world,
     >                   ierr )
         call mpi_wait( request, status, ierr )
         if (timeron) call timer_stop(t_rcomm)

         sum = sum + rho
      enddo
      rho = sum


C      print *, "00000000000ptr_a: ", ptr_a
c---------------------------------------------------------------------
c---->
c  The conj grad iteration loop
c---->
c---------------------------------------------------------------------
      do cgit = 1, cgitmax
c          CALL begin_one_iteration(main_loop)
c---------------------phase 0 ----------------------------------------
c          CALL begin_one_phase
c---------------------------------------------------------------------
c  q = A.p
c  The partition submatrix-vector multiply: use workspace w
c---------------------------------------------------------------------
         call c_switch(ptr_p_n,ptr_p_o)
         call c_switch(ptr_r_n,ptr_r_o)
         call c_switch(ptr_z_n,ptr_z_o)

         do j=1,lastrow-firstrow+1
            sum = 0.d0
            do k=rowstr(j),rowstr(j+1)-1
               sum = sum + a(k)*p_e(colidx(k))
            enddo
            w(j) = sum
         enddo

c        CALL end_one_phase 
c-------------------------phase 1 ----------------------------------- 
c        CALL begin_one_phase
c---------------------------------------------------------------------
c  Sum the partition submatrix-vec A.p's across rows
c  Exchange and sum piece of w with procs identified in reduce_exch_proc
c---------------------------------------------------------------------
         do i = l2npcols, 1, -1
            if (timeron) call timer_start(t_rcomm)
            
          
            
            call mpi_irecv( q(reduce_recv_starts(i)),
     >                      reduce_recv_lengths(i),
     >                      dp_type,
     >                      reduce_exch_proc(i),
     >                      i,
     >                      mpi_comm_world,
     >                      request,
     >                      ierr )
            call mpi_send(  w(reduce_send_starts(i)),
     >                      reduce_send_lengths(i),
     >                      dp_type,
     >                      reduce_exch_proc(i),
     >                      i,
     >                      mpi_comm_world,
     >                      ierr )
            call mpi_wait( request, status, ierr )
            
           
            
            if (timeron) call timer_stop(t_rcomm)
           

            do j=send_start,send_start + reduce_recv_lengths(i) - 1
               w(j) = w(j) + q(j)
            enddo

           
         enddo

c         CALL end_one_phase
c-------------------------phase 2 -----------------------------------      
c         CALL begin_one_phase
c---------------------------------------------------------------------
c  Exchange piece of q with transpose processor:
c---------------------------------------------------------------------
         if( l2npcols .ne. 0 )then
            if (timeron) call timer_start(t_rcomm)

         
            call mpi_irecv( q,               
     >                      exch_recv_length,
     >                      dp_type,
     >                      exch_proc,
     >                      1,
     >                      mpi_comm_world,
     >                      request,
     >                      ierr )

            call mpi_send(  w(send_start),   
     >                      send_len,
     >                      dp_type,
     >                      exch_proc,
     >                      1,
     >                      mpi_comm_world,
     >                      ierr )
            call mpi_wait( request, status, ierr )

        

            if (timeron) call timer_stop(t_rcomm)
         else

            do j=1,exch_recv_length
               q(j) = w(j)
            enddo

         endif

c         CALL end_one_phase
c-------------------------phase 3------- -----------------------------
c         CALL begin_one_phase
c---------------------------------------------------------------------
c  Clear w for reuse...
c---------------------------------------------------------------------
      

         do j=1, max( lastrow-firstrow+1, lastcol-firstcol+1 )
            w(j) = 0.0d0
         enddo
         

c---------------------------------------------------------------------
c  Obtain p.q
c---------------------------------------------------------------------
         sum = 0.0d0
         do j=1, lastcol-firstcol+1
            sum = sum + p_e(j)*q(j)
         enddo

      
c         CALL end_one_phase
c----------------------phase 4-----------------------------------------
c         CALL begin_one_phase
c---------------------------------------------------------------------
c  Obtain d with a sum-reduce
c---------------------------------------------------------------------
         do i = 1, l2npcols
            if (timeron) call timer_start(t_rcomm)


            call mpi_irecv( d,
     >                      1,
     >                      dp_type,
     >                      reduce_exch_proc(i),
     >                      i,
     >                      mpi_comm_world,
     >                      request,
     >                      ierr )
            call mpi_send(  sum,
     >                      1,
     >                      dp_type,
     >                      reduce_exch_proc(i),
     >                      i,
     >                      mpi_comm_world,
     >                      ierr )

            call mpi_wait( request, status, ierr )


            if (timeron) call timer_stop(t_rcomm)


            sum = sum + d
            
         enddo

c         CALL end_one_phase
c--------------------------------phase 5 ------------------------
c         CALL begin_one_phase

         d = sum


c---------------------------------------------------------------------
c  Obtain alpha = rho / (p.q)
c---------------------------------------------------------------------
         alpha = rho / d

c---------------------------------------------------------------------
c  Save a temporary of rho
c---------------------------------------------------------------------
         rho0 = rho

c---------------------------------------------------------------------
c  Obtain z = z + alpha*p
c  and    r = r - alpha*q
c---------------------------------------------------------------------
         do j=1, lastcol-firstcol+1
            z_o(j) = z_e(j) + alpha*p_e(j)
            r_o(j) = r_e(j) - alpha*q(j)
         enddo
            
c---------------------------------------------------------------------
c  rho = r.r
c  Now, obtain the norm of r: First, sum squares of r elements locally...
c---------------------------------------------------------------------
         sum = 0.0d0
         do j=1, lastcol-firstcol+1
            sum = sum + r_o(j)*r_o(j)
         enddo
         
c         CALL end_one_phase
c-------------phase 6 -------------------------------------------------
c         CALL begin_one_phase
c---------------------------------------------------------------------
c  Obtain rho with a sum-reduce
c---------------------------------------------------------------------
         do i = 1, l2npcols
            if (timeron) call timer_start(t_rcomm)


            call mpi_irecv( rho,
     >                      1,
     >                      dp_type,
     >                      reduce_exch_proc(i),
     >                      i,
     >                      mpi_comm_world,
     >                      request,
     >                      ierr )
            call mpi_send(  sum,
     >                      1,
     >                      dp_type,
     >                      reduce_exch_proc(i),
     >                      i,
     >                      mpi_comm_world,
     >                      ierr )
            call mpi_wait( request, status, ierr )


            if (timeron) call timer_stop(t_rcomm)


            sum = sum + rho


         enddo
         
c         CALL end_one_phase
c------------------phase  7-------------------------------------------
c         CALL begin_one_phase

         rho = sum

c---------------------------------------------------------------------
c  Obtain beta:
c---------------------------------------------------------------------
         beta = rho / rho0

c---------------------------------------------------------------------
c  p = r + beta*p
c---------------------------------------------------------------------
         do j=1, lastcol-firstcol+1
            p_o(j) = r_o(j) + beta*p_e(j)
         enddo
         
c         CALL end_one_phase
c--------------kai------------------------
c         if (main_loop .eq. 1) then
c          call c_dram_cache_cp(ptr_p_copy, ptr_p, copy_size)
c          call c_dram_cache_cp(ptr_z_copy, ptr_z, copy_size)
c          call c_dram_cache_cp(ptr_x_copy, ptr_x, copy_size)
c          call c_dram_cache_move(ptr_p_copy, ptr_p, copy_size)
c          call c_dram_cache_move(ptr_z_copy, ptr_p, copy_size)
c          call c_dram_cache_move(ptr_x_copy, ptr_p, copy_size)


c          call c_memcpy(ptr_p_copy2, ptr_p_copy, copy_size)
c          call c_memcpy(ptr_z_copy2, ptr_z_copy, copy_size)
c          call c_memcpy(ptr_x_copy2, ptr_x_copy, copy_size)


c         call c_memmove_movnt(ptr_p_copy2, ptr_p_copy, copy_size)                                                                                
c         call c_memmove_movnt(ptr_z_copy2, ptr_z_copy, copy_size)                                                                                
c         call c_memmove_movnt(ptr_x_copy2, ptr_x_copy, copy_size)

c          call c_memwrite(ptr_p_copy, sizeof(sum), copy_size, 
c     >         curr_rank, 1)                                                                               
c          call c_memwrite(ptr_z_copy, sizeof(sum), copy_size, 
c     >         curr_rank, 1)                                                                               
c          call c_memwrite(ptr_x_copy, sizeof(sum), copy_size, 
c     >         curr_rank, 1)                                                                             
  
c          call c_memwrite(ptr_p_copy, sizeof(sum), copy_size, 
c     >         curr_rank, 2)                                                                               
c          call c_memwrite(ptr_z_copy, sizeof(sum), copy_size, 
c     >         curr_rank, 2)                                                                               
c          call c_memwrite(ptr_x_copy, sizeof(sum), copy_size, 
c     >         curr_rank, 2)                                                                             
  
c          call c_memwrite(ptr_p_copy, sizeof(sum), copy_size, 
c     >         curr_rank, 3)                                                                               
c          call c_memwrite(ptr_z_copy, sizeof(sum), copy_size, 
c     >         curr_rank, 3)                                                                               
c          call c_memwrite(ptr_x_copy, sizeof(sum), copy_size, 
c     >         curr_rank, 3)
c          end if
c----------------------------------------
      enddo                             ! end of do cgit=1,cgitmax
C      print *, "444444444ptr_a: ", ptr_a

c--ychuang
      do j=1,naa/nprows+1
         p(j) = p_o(j)
         r(j) = r_o(j)
         z(j) = z_o(j)                
      enddo

c---------------------------------------------------------------------
c  Compute residual norm explicitly:  ||r|| = ||x - A.z||
c  First, form A.z
c  The partition submatrix-vector multiply
c---------------------------------------------------------------------
      do j=1,lastrow-firstrow+1
         sum = 0.d0
         do k=rowstr(j),rowstr(j+1)-1
            sum = sum + a(k)*z(colidx(k))
         enddo
         w(j) = sum
      enddo


c---------------------------------------------------------------------
c  Sum the partition submatrix-vec A.z's across rows
c---------------------------------------------------------------------
      do i = l2npcols, 1, -1
         if (timeron) call timer_start(t_rcomm)
        
         call mpi_irecv( r(reduce_recv_starts(i)),
     >                   reduce_recv_lengths(i),
     >                   dp_type,
     >                   reduce_exch_proc(i),
     >                   i,
     >                   mpi_comm_world,
     >                   request,
     >                   ierr )
         call mpi_send(  w(reduce_send_starts(i)),
     >                   reduce_send_lengths(i),
     >                   dp_type,
     >                   reduce_exch_proc(i),
     >                   i,
     >                   mpi_comm_world,
     >                   ierr )
         call mpi_wait( request, status, ierr )
         if (timeron) call timer_stop(t_rcomm)

         do j=send_start,send_start + reduce_recv_lengths(i) - 1
            w(j) = w(j) + r(j)
         enddo
      enddo
      

c---------------------------------------------------------------------
c  Exchange piece of q with transpose processor:
c---------------------------------------------------------------------
      if( l2npcols .ne. 0 )then
         if (timeron) call timer_start(t_rcomm)
         call mpi_irecv( r,               
     >                   exch_recv_length,
     >                   dp_type,
     >                   exch_proc,
     >                   1,
     >                   mpi_comm_world,
     >                   request,
     >                   ierr )
   
         call mpi_send(  w(send_start),   
     >                   send_len,
     >                   dp_type,
     >                   exch_proc,
     >                   1,
     >                   mpi_comm_world,
     >                   ierr )
         call mpi_wait( request, status, ierr )
         if (timeron) call timer_stop(t_rcomm)
      else
         do j=1,exch_recv_length
            r(j) = w(j)
         enddo
      endif


c---------------------------------------------------------------------
c  At this point, r contains A.z
c---------------------------------------------------------------------
         sum = 0.0d0
         do j=1, lastcol-firstcol+1
            d   = x(j) - r(j)         
            sum = sum + d*d
         enddo
         
c---------------------------------------------------------------------
c  Obtain d with a sum-reduce
c---------------------------------------------------------------------
      do i = 1, l2npcols
         if (timeron) call timer_start(t_rcomm)
         call mpi_irecv( d,
     >                   1,
     >                   dp_type,
     >                   reduce_exch_proc(i),
     >                   i,
     >                   mpi_comm_world,
     >                   request,
     >                   ierr )
         call mpi_send(  sum,
     >                   1,
     >                   dp_type,
     >                   reduce_exch_proc(i),
     >                   i,
     >                   mpi_comm_world,
     >                   ierr )
         call mpi_wait( request, status, ierr )
         if (timeron) call timer_stop(t_rcomm)

         sum = sum + d
      enddo
      d = sum


      if( me .eq. root ) rnorm = sqrt( d )

      if (timeron) call timer_stop(t_conjg)


      return
      end                               ! end of routine conj_grad



c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine makea( n, nz, a, colidx, rowstr, nonzer,
     >                  firstrow, lastrow, firstcol, lastcol,
     >                  rcond, arow, acol, aelt, v, iv, shift,
     >                  amult, tran)
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit            none
      integer             n, nz, nonzer
      integer             firstrow, lastrow, firstcol, lastcol
      integer             colidx(nz), rowstr(n+1)
      integer             iv(2*n+1), arow(nz), acol(nz)
      double precision    v(n+1), aelt(nz)
      double precision    rcond, a(nz), shift
      double precision    amult, tran

c---------------------------------------------------------------------
c       generate the test problem for benchmark 6
c       makea generates a sparse matrix with a
c       prescribed sparsity distribution
c
c       parameter    type        usage
c
c       input
c
c       n            i           number of cols/rows of matrix
c       nz           i           nonzeros as declared array size
c       rcond        r*8         condition number
c       shift        r*8         main diagonal shift
c
c       output
c
c       a            r*8         array for nonzeros
c       colidx       i           col indices
c       rowstr       i           row pointers
c
c       workspace
c
c       iv, arow, acol i
c       v, aelt        r*8
c---------------------------------------------------------------------

      integer i, nnza, iouter, ivelt, ivelt1, irow, nzv, jcol

c---------------------------------------------------------------------
c      nonzer is approximately  (int(sqrt(nnza /n)));
c---------------------------------------------------------------------

      double precision  size, ratio, scale
      external          sparse, sprnvc, vecset

      size = 1.0D0
      ratio = rcond ** (1.0D0 / dfloat(n))
      nnza = 0

c---------------------------------------------------------------------
c  Initialize iv(n+1 .. 2n) to zero.
c  Used by sprnvc to mark nonzero positions
c---------------------------------------------------------------------

      do i = 1, n
           iv(n+i) = 0
      enddo
      do iouter = 1, n
         nzv = nonzer
         call sprnvc( n, nzv, v, colidx, iv(1), iv(n+1),amult,tran )
         call vecset( n, v, colidx, nzv, iouter, .5D0 )
         do ivelt = 1, nzv
              jcol = colidx(ivelt)
              if (jcol.ge.firstcol .and. jcol.le.lastcol) then
                 scale = size * v(ivelt)
                 do ivelt1 = 1, nzv
                    irow = colidx(ivelt1)
                    if (irow.ge.firstrow .and. irow.le.lastrow) then
                       nnza = nnza + 1
                       if (nnza .gt. nz) goto 9999
                       acol(nnza) = jcol
                       arow(nnza) = irow
                       aelt(nnza) = v(ivelt1) * scale
                    endif
                 enddo
              endif
         enddo
         size = size * ratio
      enddo


c---------------------------------------------------------------------
c       ... add the identity * rcond to the generated matrix to bound
c           the smallest eigenvalue from below by rcond
c---------------------------------------------------------------------
        do i = firstrow, lastrow
           if (i.ge.firstcol .and. i.le.lastcol) then
              iouter = n + i
              nnza = nnza + 1
              if (nnza .gt. nz) goto 9999
              acol(nnza) = i
              arow(nnza) = i
              aelt(nnza) = rcond - shift
           endif
        enddo


c---------------------------------------------------------------------
c       ... make the sparse matrix from list of elements with duplicates
c           (v and iv are used as  workspace)
c---------------------------------------------------------------------
      call sparse( a, colidx, rowstr, n, arow, acol, aelt,
     >             firstrow, lastrow,
     >             v, iv(1), iv(n+1), nnza )
      return

 9999 continue
      write(*,*) 'Space for matrix elements exceeded in makea'
      write(*,*) 'nnza, nzmax = ',nnza, nz
      write(*,*) ' iouter = ',iouter

      stop
      end
c-------end   of makea------------------------------

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine sparse( a, colidx, rowstr, n, arow, acol, aelt,
     >                   firstrow, lastrow,
     >                   x, mark, nzloc, nnza )
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit           none
      integer            colidx(*), rowstr(*)
      integer            firstrow, lastrow
      integer            n, arow(*), acol(*), nnza
      double precision   a(*), aelt(*)

c---------------------------------------------------------------------
c       rows range from firstrow to lastrow
c       the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
c---------------------------------------------------------------------
      integer            nzloc(n), nrows
      double precision   x(n)
      logical            mark(n)

c---------------------------------------------------
c       generate a sparse matrix from a list of
c       [col, row, element] tri
c---------------------------------------------------

      integer            i, j, jajp1, nza, k, nzrow
      double precision   xi

c---------------------------------------------------------------------
c    how many rows of result
c---------------------------------------------------------------------
      nrows = lastrow - firstrow + 1

c---------------------------------------------------------------------
c     ...count the number of triples in each row
c---------------------------------------------------------------------
      do j = 1, n
         rowstr(j) = 0
         mark(j) = .false.
      enddo
      rowstr(n+1) = 0

      do nza = 1, nnza
         j = (arow(nza) - firstrow + 1) + 1
         rowstr(j) = rowstr(j) + 1
      enddo

      rowstr(1) = 1
      do j = 2, nrows+1
         rowstr(j) = rowstr(j) + rowstr(j-1)
      enddo


c---------------------------------------------------------------------
c     ... rowstr(j) now is the location of the first nonzero
c           of row j of a
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c     ... do a bucket sort of the triples on the row index
c---------------------------------------------------------------------
      do nza = 1, nnza
         j = arow(nza) - firstrow + 1
         k = rowstr(j)
         a(k) = aelt(nza)
         colidx(k) = acol(nza)
         rowstr(j) = rowstr(j) + 1
      enddo


c---------------------------------------------------------------------
c       ... rowstr(j) now points to the first element of row j+1
c---------------------------------------------------------------------
      do j = nrows, 1, -1
          rowstr(j+1) = rowstr(j)
      enddo
      rowstr(1) = 1


c---------------------------------------------------------------------
c       ... generate the actual output rows by adding elements
c---------------------------------------------------------------------
      nza = 0
      do i = 1, n
          x(i)    = 0.0
          mark(i) = .false.
      enddo

      jajp1 = rowstr(1)
      do j = 1, nrows
         nzrow = 0

c---------------------------------------------------------------------
c          ...loop over the jth row of a
c---------------------------------------------------------------------
         do k = jajp1 , rowstr(j+1)-1
            i = colidx(k)
            x(i) = x(i) + a(k)
            if ( (.not. mark(i)) .and. (x(i) .ne. 0.D0)) then
             mark(i) = .true.
             nzrow = nzrow + 1
             nzloc(nzrow) = i
            endif
         enddo

c---------------------------------------------------------------------
c          ... extract the nonzeros of this row
c---------------------------------------------------------------------
         do k = 1, nzrow
            i = nzloc(k)
            mark(i) = .false.
            xi = x(i)
            x(i) = 0.D0
            if (xi .ne. 0.D0) then
             nza = nza + 1
             a(nza) = xi
             colidx(nza) = i
            endif
         enddo
         jajp1 = rowstr(j+1)
         rowstr(j+1) = nza + rowstr(1)
      enddo
CC       write (*, 11000) nza
      return
11000   format ( //,'final nonzero count in sparse ',
     1            /,'number of nonzeros       = ', i16 )
      end
c-------end   of sparse-----------------------------


c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine sprnvc( n, nz, v, iv, nzloc, mark, amult, tran)
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit           none
      double precision   v(*)
      integer            n, nz, iv(*), nzloc(n), nn1
      integer mark(n)
c      common /urando/    amult, tran
      double precision   amult, tran


c---------------------------------------------------------------------
c       generate a sparse n-vector (v, iv)
c       having nzv nonzeros
c
c       mark(i) is set to 1 if position i is nonzero.
c       mark is all zero on entry and is reset to all zero before exit
c       this corrects a performance bug found by John G. Lewis, caused by
c       reinitialization of mark on every one of the n calls to sprnvc
c---------------------------------------------------------------------

        integer            nzrow, nzv, ii, i, icnvrt

        external           randlc, icnvrt
        double precision   randlc, vecelt, vecloc


        nzv = 0
        nzrow = 0
        nn1 = 1
 50     continue
          nn1 = 2 * nn1
          if (nn1 .lt. n) goto 50

c---------------------------------------------------------------------
c    nn1 is the smallest power of two not less than n
c---------------------------------------------------------------------

100     continue
        if (nzv .ge. nz) goto 110
         vecelt = randlc( tran, amult )

c---------------------------------------------------------------------
c   generate an integer between 1 and n in a portable manner
c---------------------------------------------------------------------
         vecloc = randlc(tran, amult)
         i = icnvrt(vecloc, nn1) + 1
         if (i .gt. n) goto 100

c---------------------------------------------------------------------
c  was this integer generated already?
c---------------------------------------------------------------------
         if (mark(i) .eq. 0) then
            mark(i) = 1
            nzrow = nzrow + 1
            nzloc(nzrow) = i
            nzv = nzv + 1
            v(nzv) = vecelt
            iv(nzv) = i
         endif
         goto 100
110      continue
      do ii = 1, nzrow
         i = nzloc(ii)
         mark(i) = 0
      enddo
      return
      end
c-------end   of sprnvc-----------------------------


c---------------------------------------------------------------------
c---------------------------------------------------------------------
      function icnvrt(x, ipwr2)
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit           none
      double precision   x
      integer            ipwr2, icnvrt

c---------------------------------------------------------------------
c    scale a double precision number x in (0,1) by a power of 2 and chop it
c---------------------------------------------------------------------
      icnvrt = int(ipwr2 * x)

      return
      end
c-------end   of icnvrt-----------------------------


c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine vecset(n, v, iv, nzv, i, val)
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit           none
      integer            n, iv(*), nzv, i, k
      double precision   v(*), val

c---------------------------------------------------------------------
c       set ith element of sparse vector (v, iv) with
c       nzv nonzeros to val
c---------------------------------------------------------------------

      logical set

      set = .false.
      do k = 1, nzv
         if (iv(k) .eq. i) then
            v(k) = val
            set  = .true.
         endif
      enddo
      if (.not. set) then
         nzv     = nzv + 1
         v(nzv)  = val
         iv(nzv) = i
      endif
      return
      end
c-------end   of vecset-----------------------------

