!-------------------------------------------------------------------------!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3         !
!                                                                         !
!                                   S P                                   !
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
c Authors: R. F. Van der Wijngaart
c          W. Saphir
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       program MPSP
c---------------------------------------------------------------------

       include  'header.h'
       include  'mpinpb.h'
      
       integer          i, niter, step, c, error, fstatus
       pointer(ptr_i, i)
       pointer(ptr_niter, niter)
       pointer(ptr_step, step)
       pointer(ptr_c, c)
       pointer(ptr_error, error)
       pointer(ptr_fstatus, fstatus)

       external timer_read
       double precision mflops, t, tmax, timer_read
       pointer(ptr_mflops, mflops)
       pointer(ptr_t, t)
       pointer(ptr_tmax, tmax)
       
       logical          verified
       pointer(ptr_verified, verified)
       character        class
       pointer(ptr_class, class)
       double precision tsum(t_last+2), t1(t_last+2),
     >                  tming(t_last+2), tmaxg(t_last+2)
       pointer(ptr_tsum, tsum)
       pointer(ptr_t1, t1)
       pointer(ptr_tming, tming)
       pointer(ptr_tmaxg, tmaxg)

       character        t_recs(t_last+2)*8
       pointer(ptr_t_recs, t_recs)

       data t_recs/'total', 'rhs', 'xsolve', 'ysolve', 'zsolve', 
     >             'bpack', 'exch', 'xcomm', 'ycomm', 'zcomm',
     >             ' totcomp', ' totcomm'/

       integer data_pos, main_loop
       integer curr_rank
       common /myrank/curr_rank

       data_pos = -1
       main_loop = -1

c       call c_openmp_init(4)
       
c--------------------------------------kai----------------------------
       call unimem_malloc(ptr_u_copy, sizeof(u_copy), data_pos)
       call unimem_malloc(ptr_u_copy2, sizeof(u_copy2), data_pos)

       call unimem_malloc(ptr_ncells, sizeof(ncells), data_pos)
       call unimem_malloc(ptr_grid_points, sizeof(grid_points),data_pos)

       call unimem_malloc(ptr_tx1, sizeof(tx1), data_pos)
       call unimem_malloc(ptr_tx2, sizeof(tx2), data_pos)
       call unimem_malloc(ptr_tx3, sizeof(tx3), data_pos)
       call unimem_malloc(ptr_ty1, sizeof(ty1), data_pos)
       call unimem_malloc(ptr_ty2, sizeof(ty2), data_pos)
       call unimem_malloc(ptr_ty3, sizeof(ty3), data_pos)
       call unimem_malloc(ptr_tz1, sizeof(tz1), data_pos)
       call unimem_malloc(ptr_tz2, sizeof(tz2), data_pos)
       call unimem_malloc(ptr_tz3, sizeof(tz3), data_pos)

       call unimem_malloc(ptr_dx1, sizeof(dx1), data_pos)
       call unimem_malloc(ptr_dx2, sizeof(dx2), data_pos)
       call unimem_malloc(ptr_dx3, sizeof(dx3), data_pos)
       call unimem_malloc(ptr_dx4, sizeof(dx4), data_pos)
       call unimem_malloc(ptr_dx5, sizeof(dx5), data_pos)
       call unimem_malloc(ptr_dy1, sizeof(dy1), data_pos)
       call unimem_malloc(ptr_dy2, sizeof(dy2), data_pos)
       call unimem_malloc(ptr_dy3, sizeof(dy3), data_pos)
       call unimem_malloc(ptr_dy4, sizeof(dy4), data_pos)
       call unimem_malloc(ptr_dy5, sizeof(dy5), data_pos)
       call unimem_malloc(ptr_dz1, sizeof(dz1), data_pos)
       call unimem_malloc(ptr_dz2, sizeof(dz2), data_pos)
       call unimem_malloc(ptr_dz3, sizeof(dz3), data_pos)
       call unimem_malloc(ptr_dz4, sizeof(dz4), data_pos)
       call unimem_malloc(ptr_dz5, sizeof(dz5), data_pos)
       call unimem_malloc(ptr_dssp, sizeof(dssp), data_pos)
       call unimem_malloc(ptr_dt, sizeof(dt), data_pos)
       call unimem_malloc(ptr_ce, sizeof(ce), data_pos)

       call unimem_malloc(ptr_dxmax, sizeof(dxmax), data_pos)
       call unimem_malloc(ptr_dymax, sizeof(dymax), data_pos)
       call unimem_malloc(ptr_dzmax, sizeof(dzmax), data_pos)

       call unimem_malloc(ptr_xxcon1, sizeof(xxcon1), data_pos)
       call unimem_malloc(ptr_xxcon2, sizeof(xxcon2), data_pos)
       call unimem_malloc(ptr_xxcon3, sizeof(xxcon3), data_pos)
       call unimem_malloc(ptr_xxcon4, sizeof(xxcon4), data_pos)
       call unimem_malloc(ptr_xxcon5, sizeof(xxcon5), data_pos)

       call unimem_malloc(ptr_dx1tx1, sizeof(dx1tx1), data_pos)
       call unimem_malloc(ptr_dx2tx1, sizeof(dx2tx1), data_pos)
       call unimem_malloc(ptr_dx3tx1, sizeof(dx3tx1), data_pos)
       call unimem_malloc(ptr_dx4tx1, sizeof(dx4tx1), data_pos)
       call unimem_malloc(ptr_dx5tx1, sizeof(dx5tx1), data_pos)

       call unimem_malloc(ptr_yycon1, sizeof(yycon1), data_pos)
       call unimem_malloc(ptr_yycon2, sizeof(yycon2), data_pos)
       call unimem_malloc(ptr_yycon3, sizeof(yycon3), data_pos)
       call unimem_malloc(ptr_yycon4, sizeof(yycon4), data_pos)
       call unimem_malloc(ptr_yycon5, sizeof(yycon5), data_pos)

       call unimem_malloc(ptr_dy1ty1, sizeof(dy1ty1), data_pos)
       call unimem_malloc(ptr_dy2ty1, sizeof(dy2ty1), data_pos)
       call unimem_malloc(ptr_dy3ty1, sizeof(dy3ty1), data_pos)
       call unimem_malloc(ptr_dy4ty1, sizeof(dy4ty1), data_pos)
       call unimem_malloc(ptr_dy5ty1, sizeof(dy5ty1), data_pos)
       call unimem_malloc(ptr_zzcon1, sizeof(zzcon1), data_pos)
       call unimem_malloc(ptr_zzcon2, sizeof(zzcon2), data_pos)
       call unimem_malloc(ptr_zzcon3, sizeof(zzcon3), data_pos)
       call unimem_malloc(ptr_zzcon4, sizeof(zzcon4), data_pos)
       call unimem_malloc(ptr_zzcon5, sizeof(zzcon5), data_pos)

       call unimem_malloc(ptr_dz1tz1, sizeof(dz1tz1), data_pos)
       call unimem_malloc(ptr_dz2tz1, sizeof(dz2tz1), data_pos)
       call unimem_malloc(ptr_dz3tz1, sizeof(dz3tz1), data_pos)
       call unimem_malloc(ptr_dz4tz1, sizeof(dz4tz1), data_pos)
       call unimem_malloc(ptr_dz5tz1, sizeof(dz5tz1), data_pos)

       call unimem_malloc(ptr_dnxm1, sizeof(dnxm1), data_pos)
       call unimem_malloc(ptr_dnym1, sizeof(dnym1), data_pos)
       call unimem_malloc(ptr_dnzm1, sizeof(dnzm1), data_pos)
       call unimem_malloc(ptr_c1c2, sizeof(c1c2), data_pos)
       call unimem_malloc(ptr_c1c5, sizeof(c1c5), data_pos)
       call unimem_malloc(ptr_c3c4, sizeof(c3c4), data_pos)
       call unimem_malloc(ptr_c1345, sizeof(c1345), data_pos)
       call unimem_malloc(ptr_conz1, sizeof(conz1), data_pos)

       call unimem_malloc(ptr_c1, sizeof(c1), data_pos)
       call unimem_malloc(ptr_c1, sizeof(c1), data_pos)
       call unimem_malloc(ptr_c2, sizeof(c2), data_pos)
       call unimem_malloc(ptr_c3, sizeof(c3), data_pos)
       call unimem_malloc(ptr_c4, sizeof(c4), data_pos)
       call unimem_malloc(ptr_c5, sizeof(c5), data_pos)
       call unimem_malloc(ptr_c4dssp, sizeof(c4dssp), data_pos)
       call unimem_malloc(ptr_c5dssp, sizeof(c5dssp), data_pos)
       call unimem_malloc(ptr_dtdssp, sizeof(dtdssp), data_pos)
       call unimem_malloc(ptr_dttx1, sizeof(dttx1), data_pos)
       call unimem_malloc(ptr_bt, sizeof(bt), data_pos)
       call unimem_malloc(ptr_dttx2, sizeof(dttx2), data_pos)
       call unimem_malloc(ptr_dtty1, sizeof(dtty1), data_pos)
       call unimem_malloc(ptr_dtty2, sizeof(dtty2), data_pos)
       call unimem_malloc(ptr_dttz1, sizeof(dttz1), data_pos)
       call unimem_malloc(ptr_dttz2, sizeof(dttz2), data_pos)
       call unimem_malloc(ptr_c2dttx1, sizeof(c2dttx1), data_pos)
       call unimem_malloc(ptr_c2dtty1, sizeof(c2dtty1), data_pos)
       call unimem_malloc(ptr_c2dttz1, sizeof(c2dttz1), data_pos)
       call unimem_malloc(ptr_comz1, sizeof(comz1), data_pos)
       call unimem_malloc(ptr_comz4, sizeof(comz4), data_pos)
       call unimem_malloc(ptr_comz5, sizeof(comz5), data_pos)
       call unimem_malloc(ptr_comz6, sizeof(comz6), data_pos)
       call unimem_malloc(ptr_c3c4tx3, sizeof(c3c4tx3), data_pos)
       call unimem_malloc(ptr_c3c4ty3, sizeof(c3c4ty3), data_pos)
       call unimem_malloc(ptr_c3c4tz3, sizeof(c3c4tz3), data_pos)
       call unimem_malloc(ptr_c2iv, sizeof(c2iv), data_pos)
       call unimem_malloc(ptr_con43, sizeof(con43), data_pos)
       call unimem_malloc(ptr_con16, sizeof(con16), data_pos)

       call unimem_malloc(ptr_EAST, sizeof(EAST), data_pos)
       call unimem_malloc(ptr_SOUTH, sizeof(SOUTH), data_pos)
       call unimem_malloc(ptr_NORTH, sizeof(NORTH), data_pos)
       call unimem_malloc(ptr_WEST, sizeof(WEST), data_pos)
       call unimem_malloc(ptr_BOTTOM, sizeof(BOTTOM), data_pos)
       call unimem_malloc(ptr_TOP, sizeof(TOP), data_pos)

       call unimem_malloc(ptr_cell_coord, sizeof(cell_coord), data_pos)
       call unimem_malloc(ptr_cell_low, sizeof(cell_low), data_pos)
       call unimem_malloc(ptr_cell_high, sizeof(cell_high), data_pos)
       call unimem_malloc(ptr_cell_size, sizeof(cell_size), data_pos)
       call unimem_malloc(ptr_grid_size, sizeof(grid_size), data_pos)
       call unimem_malloc(ptr_successor, sizeof(successor), data_pos)
       call unimem_malloc(ptr_predecessor, sizeof(predecessor),data_pos)
       call unimem_malloc(ptr_slice, sizeof(slice), data_pos)
       call unimem_malloc(ptr_start, sizeof(start), data_pos)
       call unimem_malloc(ptr_end, sizeof(end), data_pos)

       call unimem_malloc(ptr_IMAX, sizeof(IMAX), data_pos)
       call unimem_malloc(ptr_JMAX, sizeof(JMAX), data_pos)
       call unimem_malloc(ptr_KMAX, sizeof(KMAX), data_pos)
       call unimem_malloc(ptr_MAX_CELL_DIM, sizeof(MAX_CELL_DIM),
     >                    data_pos)
       call unimem_malloc(ptr_BUF_SIZE, sizeof(BUF_SIZE), data_pos)

c       call unimem_malloc(ptr_u, sizeof(u), data_pos)
       call unimem_malloc(ptr_u, sizeof(u), 1)

       call unimem_malloc(ptr_us, sizeof(us), data_pos)
       call unimem_malloc(ptr_vs, sizeof(vs), data_pos)
       call unimem_malloc(ptr_ws, sizeof(ws), data_pos)
       call unimem_malloc(ptr_qs, sizeof(qs), data_pos)
       call unimem_malloc(ptr_ainv, sizeof(ainv), data_pos)
       call unimem_malloc(ptr_rho_i, sizeof(rho_i), data_pos)
       call unimem_malloc(ptr_speed, sizeof(speed), data_pos)
       call unimem_malloc(ptr_square, sizeof(square), data_pos)
c       call unimem_malloc(ptr_forcing, sizeof(forcing), data_pos)
       call unimem_malloc(ptr_forcing, sizeof(forcing), 1)

c       call unimem_malloc(ptr_rhs, sizeof(rhs), data_pos)
       call unimem_malloc(ptr_rhs, sizeof(rhs), 1)

c       call unimem_malloc(ptr_lhs, sizeof(lhs), data_pos)
       call unimem_malloc(ptr_lhs, sizeof(lhs), 1)

       call unimem_malloc(ptr_in_buffer, sizeof(in_buffer), data_pos)
       call unimem_malloc(ptr_out_buffer, sizeof(out_buffer), data_pos)

       call unimem_malloc(ptr_cv, sizeof(cv), data_pos)
       call unimem_malloc(ptr_rhon, sizeof(rhon), data_pos)
       call unimem_malloc(ptr_rhos, sizeof(rhos), data_pos)
       call unimem_malloc(ptr_rhoq, sizeof(rhoq), data_pos)
       call unimem_malloc(ptr_cuf, sizeof(cuf), data_pos)
       call unimem_malloc(ptr_q, sizeof(q), data_pos)
       call unimem_malloc(ptr_ue, sizeof(ue), data_pos)
       call unimem_malloc(ptr_buf, sizeof(buf), data_pos)

       call unimem_malloc(ptr_west_size, sizeof(west_size), data_pos)
       call unimem_malloc(ptr_east_size, sizeof(east_size), data_pos)
       call unimem_malloc(ptr_bottom_size, sizeof(bottom_size),data_pos)
       call unimem_malloc(ptr_top_size, sizeof(top_size), data_pos)
       call unimem_malloc(ptr_north_size, sizeof(north_size), data_pos)
       call unimem_malloc(ptr_south_size, sizeof(south_size), data_pos)
       call unimem_malloc(ptr_start_send_west, sizeof(start_send_west),
     >                    data_pos)
       call unimem_malloc(ptr_start_send_east, sizeof(start_send_east),
     >                    data_pos)
       call unimem_malloc(ptr_start_send_south,sizeof(start_send_south),
     >                    data_pos)
       call unimem_malloc(ptr_start_send_north,sizeof(start_send_north),
     >                    data_pos)
       call unimem_malloc(ptr_start_send_bottom,
     >                    sizeof(start_send_bottom), data_pos)
       call unimem_malloc(ptr_start_send_top, sizeof(start_send_top),
     >                    data_pos)
       call unimem_malloc(ptr_start_recv_west, sizeof(start_recv_west),
     >                    data_pos)
       call unimem_malloc(ptr_start_recv_east, sizeof(start_recv_east),
     >                    data_pos)
       call unimem_malloc(ptr_start_recv_south,sizeof(start_recv_south),
     >                    data_pos)
       call unimem_malloc(ptr_start_recv_north,sizeof(start_recv_north),
     >                    data_pos)
       call unimem_malloc(ptr_start_recv_bottom,
     >                    sizeof(start_recv_bottom), data_pos)
       call unimem_malloc(ptr_start_recv_top, sizeof(start_recv_top),
     >                    data_pos)

      call unimem_malloc(ptr_t_total, sizeof(t_total), data_pos)
      call unimem_malloc(ptr_t_rhs, sizeof(t_rhs), data_pos)
      call unimem_malloc(ptr_t_xsolve, sizeof(t_xsolve), data_pos)
      call unimem_malloc(ptr_t_ysolve, sizeof(t_ysolve), data_pos)
      call unimem_malloc(ptr_t_zsolve, sizeof(t_zsolve), data_pos)
      call unimem_malloc(ptr_t_bpack, sizeof(t_bpack), data_pos)
      call unimem_malloc(ptr_t_exch, sizeof(t_exch), data_pos)
      call unimem_malloc(ptr_t_xcomm, sizeof(t_xcomm), data_pos)
      call unimem_malloc(ptr_t_ycomm, sizeof(t_ycomm), data_pos)
      call unimem_malloc(ptr_t_zcomm, sizeof(t_zcomm), data_pos)
      call unimem_malloc(ptr_t_last, sizeof(t_last), data_pos)
      call unimem_malloc(ptr_timeron, sizeof(timeron), data_pos)

      call unimem_malloc(ptr_i, sizeof(i), data_pos)
      call unimem_malloc(ptr_niter, sizeof(niter), data_pos)
      call unimem_malloc(ptr_step, sizeof(step), data_pos)
      call unimem_malloc(ptr_c, sizeof(c), data_pos)
      call unimem_malloc(ptr_error, sizeof(error), data_pos)
      call unimem_malloc(ptr_fstatus, sizeof(fstatus), data_pos)

      call unimem_malloc(ptr_mflops, sizeof(mflops), data_pos)
      call unimem_malloc(ptr_t, sizeof(t), data_pos)
      call unimem_malloc(ptr_tmax, sizeof(tmax), data_pos)

      call unimem_malloc(ptr_verified, sizeof(verified), data_pos)
      call unimem_malloc(ptr_class, sizeof(class), data_pos)

      call unimem_malloc(ptr_tsum, sizeof(tsum), data_pos)
      call unimem_malloc(ptr_t1, sizeof(t1), data_pos)
      call unimem_malloc(ptr_tming, sizeof(tming), data_pos)
      call unimem_malloc(ptr_tmaxg, sizeof(tmaxg), data_pos)
      call unimem_malloc(ptr_t_recs, sizeof(t_recs), data_pos)

c---------------------------------------------------------------------
      print *, "sizeof(u)", sizeof(u)
      print *, "sizeof(us)", sizeof(us)
      print *, "sizeof(vs)", sizeof(vs)
      print *, "sizeof(ws)", sizeof(ws)
      print *, "sizeof(qs)", sizeof(qs)
      print *, "sizeof(rho_i)", sizeof(rho_i)
      print *, "sizeof(square)", sizeof(square)
      print *, "sizeof(rhs)", sizeof(rhs)
      print *, "sizeof(forcing)", sizeof(forcing)
      print *, "sizeof(out_buffer)", sizeof(out_buffer)
      print *, "sizeof(in_buffer)", sizeof(in_buffer)
      print *, "sizeof(lhs)", sizeof(lhs)




       call setup_mpi
       if (.not. active) goto 999

c---------------------------------------------------------------------
c      Root node reads input file (if it exists) else takes
c      defaults from parameters
c---------------------------------------------------------------------
       if (node .eq. root) then
          
          write(*, 1000)

          open (unit=2,file='timer.flag',status='old',iostat=fstatus)
          timeron = .false.
          if (fstatus .eq. 0) then
             timeron = .true.
             close(2)
          endif

          open (unit=2,file='inputsp.data',status='old', iostat=fstatus)
c
          if (fstatus .eq. 0) then
            write(*,233) 
 233        format(' Reading from input file inputsp.data')
            read (2,*) niter
            read (2,*) dt
            read (2,*) grid_points(1), grid_points(2), grid_points(3)
            close(2)
          else
            write(*,234) 
            niter = niter_default
            dt    = dt_default
            grid_points(1) = problem_size
            grid_points(2) = problem_size
            grid_points(3) = problem_size
          endif
 234      format(' No input file inputsp.data. Using compiled defaults')

          write(*, 1001) grid_points(1), grid_points(2), grid_points(3)
          write(*, 1002) niter, dt
          if (no_nodes .ne. total_nodes) write(*, 1004) total_nodes
          if (no_nodes .ne. maxcells*maxcells) 
     >        write(*, 1005) maxcells*maxcells
          write(*, 1003) no_nodes

 1000 format(//,' NAS Parallel Benchmarks 3.3 -- SP Benchmark',/)
 1001     format(' Size: ', i4, 'x', i4, 'x', i4)
 1002     format(' Iterations: ', i4, '    dt: ', F11.7)
 1004     format(' Total number of processes: ', i5)
 1005     format(' WARNING: compiled for ', i5, ' processes ')
 1003     format(' Number of active processes: ', i5, /)

       endif

       call mpi_bcast(niter, 1, MPI_INTEGER, 
     >                root, comm_setup, error)

       call mpi_bcast(dt, 1, dp_type, 
     >                root, comm_setup, error)

       call mpi_bcast(grid_points(1), 3, MPI_INTEGER, 
     >                root, comm_setup, error)

       call mpi_bcast(timeron, 1, MPI_LOGICAL, 
     >                root, comm_setup, error)


       call make_set

       do  c = 1, ncells
          if ( (cell_size(1,c) .gt. IMAX) .or.
     >         (cell_size(2,c) .gt. JMAX) .or.
     >         (cell_size(3,c) .gt. KMAX) ) then
             print *,node, c, (cell_size(i,c),i=1,3)
             print *,' Problem size too big for compiled array sizes'
             goto 999
          endif
       end do

       do  i = 1, t_last
          call timer_clear(i)
       end do

       call set_constants

       call initialize

       call lhsinit

       call exact_rhs

       call compute_buffer_size(5)

c---------------------------------------------------------------------
c      do one time step to touch all code, and reinitialize
c---------------------------------------------------------------------
       call adi
       call initialize

c---------------------------------------------------------------------
c      Synchronize before placing time stamp
c---------------------------------------------------------------------
       do  i = 1, t_last
          call timer_clear(i)
       end do
       call mpi_barrier(comm_setup, error)

       call timer_clear(1)
       call timer_start(1)

c       do  step = 1, niter
       do step = 1, 40
          if (node .eq. root) then
             if (mod(step, 20) .eq. 0 .or. 
     >           step .eq. 1) then
                write(*, 200) step
 200            format(' Time step ', i4)
              endif
          endif

          call adi

c----------------kai-----------------------------                                                                                                 
          call c_dram_cache_cp(ptr_u_copy, ptr_u, sizeof(u))
c          call c_dram_cache_move(ptr_u_copy, ptr_u, sizeof(u))

c         call c_memwrite(ptr_u_copy, sizeof(t), sizeof(u_copy), 
c     >         curr_rank, 1)                                                                                  
c         call c_memwrite(ptr_u_copy, sizeof(t), sizeof(u_copy),
c     >         curr_rank, 2)
         call c_memwrite(ptr_u_copy, sizeof(t), sizeof(u_copy),
     >         curr_rank, 3)
c          call c_memcpy(ptr_u_copy2, ptr_u_copy, sizeof(u_copy))
c          call c_memmove_movnt(ptr_u_copy2, ptr_u_copy, sizeof(u_copy))                                                                          
c------------------------------------------------- 

       end do

       call timer_stop(1)
       t = timer_read(1)
       
c------------kai---------
       call c_print_timer(curr_rank)
c-----------------------
       call verify(niter, class, verified)

       call mpi_reduce(t, tmax, 1, 
     >                 dp_type, MPI_MAX, 
     >                 root, comm_setup, error)

       if( node .eq. root ) then
          if( tmax .ne. 0. ) then
             mflops = (881.174*float( problem_size )**3
     >                -4683.91*float( problem_size )**2
     >                +11484.5*float( problem_size )
     >                -19272.4) * float( niter ) / (tmax*1000000.0d0)
          else
             mflops = 0.0
          endif

         call print_results('SP', class, grid_points(1), 
     >     grid_points(2), grid_points(3), niter, maxcells*maxcells, 
     >     total_nodes, tmax, mflops, '          floating point', 
     >     verified, npbversion,compiletime, cs1, cs2, cs3, cs4, cs5, 
     >     cs6, '(none)')
       endif

       if (.not.timeron) goto 999

       do i = 1, t_last
          t1(i) = timer_read(i)
       end do
       t1(t_xsolve) = t1(t_xsolve) - t1(t_xcomm)
       t1(t_ysolve) = t1(t_ysolve) - t1(t_ycomm)
       t1(t_zsolve) = t1(t_zsolve) - t1(t_zcomm)
       t1(t_last+2) = t1(t_xcomm)+t1(t_ycomm)+t1(t_zcomm)+t1(t_exch)
       t1(t_last+1) = t1(t_total)  - t1(t_last+2)

       call MPI_Reduce(t1, tsum,  t_last+2, dp_type, MPI_SUM, 
     >                 0, comm_setup, error)
       call MPI_Reduce(t1, tming, t_last+2, dp_type, MPI_MIN, 
     >                 0, comm_setup, error)
       call MPI_Reduce(t1, tmaxg, t_last+2, dp_type, MPI_MAX, 
     >                 0, comm_setup, error)

       if (node .eq. 0) then
          write(*, 800) total_nodes
          do i = 1, t_last+2
             tsum(i) = tsum(i) / total_nodes
             write(*, 810) i, t_recs(i), tming(i), tmaxg(i), tsum(i)
          end do
       endif
 800   format(' nprocs =', i6, 11x, 'minimum', 5x, 'maximum', 
     >        5x, 'average')
 810   format(' timer ', i2, '(', A8, ') :', 3(2x,f10.4))

 999   continue
       call mpi_barrier(MPI_COMM_WORLD, error)
       call mpi_finalize(error)

       end
