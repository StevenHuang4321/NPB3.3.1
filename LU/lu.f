!-------------------------------------------------------------------------!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3         !
!                                                                         !
!                                   L U                                   !
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
c Authors: S. Weeratunga
c          V. Venkatakrishnan
c          E. Barszcz
c          M. Yarrow
c
c---------------------------------------------------------------------

c---------------------------------------------------------------------
      program applu
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c   driver for the performance evaluation of the solver for
c   five coupled parabolic/elliptic partial differential equations.
c
c---------------------------------------------------------------------

      implicit none

      include 'mpinpb.h'
      include 'applu.incl'
      character class
      logical verified
      double precision mflops, timer_read
      integer i, ierr
      double precision tsum(t_last+2), t1(t_last+2),
     >                 tming(t_last+2), tmaxg(t_last+2)
      character        t_recs(t_last+2)*8

      data t_recs/'total', 'rhs', 'blts', 'buts', 'jacld', 'jacu', 
     >            'exch', 'lcomm', 'ucomm', 'rcomm',
     >            ' totcomp', ' totcomm'/

      integer data_pos, main_loop
      integer curr_rank
      common /myrank/ curr_rank

      data_pos = -1
      main_loop = -1
c---------------------------------kai--------------------------------

c      call c_openmp_init(4)                                                                                                                        

      call unimem_malloc(ptr_u_copy, sizeof(u_copy), data_pos)
      call unimem_malloc(ptr_u_copy2, sizeof(u_copy2), data_pos)



      call unimem_malloc(ptr_ipr_default, sizeof(ipr_default), data_pos)
      call unimem_malloc(ptr_omega_default, sizeof(omega_default), 
     >                   data_pos)

      call unimem_malloc(ptr_tolrsd1_def, sizeof(tolrsd1_def), data_pos)
      call unimem_malloc(ptr_tolrsd2_def, sizeof(tolrsd2_def), data_pos)
      call unimem_malloc(ptr_tolrsd3_def, sizeof(tolrsd3_def), data_pos)
      call unimem_malloc(ptr_tolrsd4_def, sizeof(tolrsd4_def), data_pos)
      call unimem_malloc(ptr_tolrsd5_def, sizeof(tolrsd5_def), data_pos)

      call unimem_malloc(ptr_c1, sizeof(c1), data_pos)
      call unimem_malloc(ptr_c2, sizeof(c2), data_pos)
      call unimem_malloc(ptr_c3, sizeof(c3), data_pos)
      call unimem_malloc(ptr_c4, sizeof(c4), data_pos)
      call unimem_malloc(ptr_c5, sizeof(c5), data_pos)
 
      call unimem_malloc(ptr_nx, sizeof(nx), data_pos)
      call unimem_malloc(ptr_ny, sizeof(ny), data_pos)
      call unimem_malloc(ptr_nz, sizeof(nz), data_pos)
      call unimem_malloc(ptr_nx0, sizeof(nx0), data_pos)
      call unimem_malloc(ptr_ny0, sizeof(ny0), data_pos)
      call unimem_malloc(ptr_nz0, sizeof(nz0), data_pos)
      call unimem_malloc(ptr_ipt, sizeof(ipt), data_pos)
      call unimem_malloc(ptr_ist, sizeof(ist), data_pos)
      call unimem_malloc(ptr_iend, sizeof(iend), data_pos)
      call unimem_malloc(ptr_jpt, sizeof(jpt), data_pos)
      call unimem_malloc(ptr_jst, sizeof(jst), data_pos)
      call unimem_malloc(ptr_jend, sizeof(jend), data_pos)
      call unimem_malloc(ptr_ii1, sizeof(ii1), data_pos)
      call unimem_malloc(ptr_ii2, sizeof(ii2), data_pos)
      call unimem_malloc(ptr_ji1, sizeof(ji1), data_pos)
      call unimem_malloc(ptr_ji2, sizeof(ji2), data_pos)     
      call unimem_malloc(ptr_ki1, sizeof(ki1), data_pos)
      call unimem_malloc(ptr_ki2, sizeof(ki2), data_pos)
      call unimem_malloc(ptr_dxi, sizeof(dxi), data_pos)
      call unimem_malloc(ptr_deta, sizeof(deta), data_pos)
      call unimem_malloc(ptr_dzeta, sizeof(dzeta), data_pos)
      call unimem_malloc(ptr_tx1, sizeof(tx1), data_pos)
      call unimem_malloc(ptr_tx2, sizeof(tx2), data_pos)
      call unimem_malloc(ptr_tx3, sizeof(tx3), data_pos)
      call unimem_malloc(ptr_ty1, sizeof(ty1), data_pos)
      call unimem_malloc(ptr_ty2, sizeof(ty2), data_pos)
      call unimem_malloc(ptr_ty3, sizeof(ty3), data_pos)
      call unimem_malloc(ptr_tz1, sizeof(ty1), data_pos)
      call unimem_malloc(ptr_tz2, sizeof(ty2), data_pos)
      call unimem_malloc(ptr_tz3, sizeof(ty3), data_pos)

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

c      call unimem_malloc(ptr_u, sizeof(u), data_pos)
c      call unimem_malloc(ptr_u, sizeof(u), 1)
      call unimem_malloc(ptr_u_n, sizeof(u), 1)
      call unimem_malloc(ptr_u_o, sizeof(u), 1)

c      call unimem_malloc(ptr_rsd, sizeof(rsd), data_pos)
      call unimem_malloc(ptr_rsd, sizeof(rsd), 1)

c      call unimem_malloc(ptr_frct, sizeof(frct), data_pos)
      call unimem_malloc(ptr_frct, sizeof(frct), 1)

c      call unimem_malloc(ptr_flux, sizeof(flux), data_pos)
      call unimem_malloc(ptr_flux, sizeof(flux), 1)

      call unimem_malloc(ptr_ipr, sizeof(ipr), data_pos)
      call unimem_malloc(ptr_inorm, sizeof(inorm), data_pos)

      call unimem_malloc(ptr_itmax, sizeof(itmax), data_pos)
      call unimem_malloc(ptr_invert, sizeof(invert), data_pos)
      call unimem_malloc(ptr_dt, sizeof(dt), data_pos)
      call unimem_malloc(ptr_omega, sizeof(omega), data_pos)
      call unimem_malloc(ptr_tolrsd, sizeof(tolrsd), data_pos)
      call unimem_malloc(ptr_rsdnm, sizeof(rsdnm), data_pos)
      call unimem_malloc(ptr_errnm, sizeof(errnm), data_pos)
      call unimem_malloc(ptr_frc, sizeof(frc), data_pos)
      call unimem_malloc(ptr_ttotal, sizeof(ttotal), data_pos)

c      call unimem_malloc(ptr_a, sizeof(a), data_pos)
      call unimem_malloc(ptr_a, sizeof(a), 1)

c      call unimem_malloc(ptr_b, sizeof(b), data_pos)
      call unimem_malloc(ptr_b, sizeof(b), 1)

c      call unimem_malloc(ptr_c, sizeof(c), data_pos)
      call unimem_malloc(ptr_c, sizeof(c), 1)

c      call unimem_malloc(ptr_d, sizeof(d), data_pos)
      call unimem_malloc(ptr_d, sizeof(d), 1)

      call unimem_malloc(ptr_ce,sizeof(ce), data_pos)
      call unimem_malloc(ptr_id, sizeof(id), data_pos)
      call unimem_malloc(ptr_ndim, sizeof(ndim), data_pos)
      call unimem_malloc(ptr_num, sizeof(num), data_pos)
      call unimem_malloc(ptr_xdim, sizeof(xdim), data_pos)
      call unimem_malloc(ptr_ydim, sizeof(ydim), data_pos)
      call unimem_malloc(ptr_row, sizeof(row), data_pos)
      call unimem_malloc(ptr_col, sizeof(col), data_pos)

      call unimem_malloc(ptr_north, sizeof(north), data_pos)
      call unimem_malloc(ptr_south, sizeof(south), data_pos)
      call unimem_malloc(ptr_east, sizeof(east), data_pos)
      call unimem_malloc(ptr_west, sizeof(west), data_pos)
      
      call unimem_malloc(ptr_from_s, sizeof(from_s), data_pos)
      call unimem_malloc(ptr_from_n, sizeof(from_n), data_pos)
      call unimem_malloc(ptr_from_e, sizeof(from_e), data_pos)
      call unimem_malloc(ptr_from_w, sizeof(from_w), data_pos)

c      call unimem_malloc(ptr_buf, sizeof(buf), data_pos)
      call unimem_malloc(ptr_buf, sizeof(buf), 1)

c      call unimem_malloc(ptr_buf1, sizeof(buf1), data_pos)
      call unimem_malloc(ptr_buf1, sizeof(buf1), 1)


      print *,"sizeof(u)", sizeof(u)
      print *,"sizeof(rsd)", sizeof(rsd)
      print *,"sizeof(frct)", sizeof(frct)
      print *,"sizeof(flux)", sizeof(flux)
      print *,"sizeof(a)", sizeof(a)
      print *,"sizeof(b)", sizeof(b)
      print *,"sizeof(c)", sizeof(c)
      print *,"sizeof(d)", sizeof(d)
      print *,"sizeof(buf)", sizeof(buf)
      print *,"sizeof(buf1)", sizeof(buf1)
      print *, "sizeof(ce)", sizeof(ce)

c---------------------------------------------------------------------
c   initialize communications
c---------------------------------------------------------------------
      call init_comm()

c---------------------------------------------------------------------
c   read input data
c---------------------------------------------------------------------
      call read_input()

      do i = 1, t_last
         call timer_clear(i)
      end do

c---------------------------------------------------------------------
c   set up processor grid
c---------------------------------------------------------------------
      call proc_grid()

c---------------------------------------------------------------------
c   determine the neighbors
c---------------------------------------------------------------------
      call neighbors()

c---------------------------------------------------------------------
c   set up sub-domain sizes
c---------------------------------------------------------------------
      call subdomain()

c---------------------------------------------------------------------
c   set up coefficients
c---------------------------------------------------------------------
      call setcoeff()

c---------------------------------------------------------------------
c   set the boundary values for dependent variables
c---------------------------------------------------------------------
      call setbv(u_n)

c---------------------------------------------------------------------
c   set the initial values for dependent variables
c---------------------------------------------------------------------
      call setiv(u_n)

c---------------------------------------------------------------------
c   compute the forcing term based on prescribed exact solution
c---------------------------------------------------------------------
      call erhs(u_n)

c---------------------------------------------------------------------
c   perform one SSOR iteration to touch all data and program pages 
c---------------------------------------------------------------------
      call ssor(1)

c---------------------------------------------------------------------
c   reset the boundary and initial values
c---------------------------------------------------------------------
      call setbv(u_n)
      call setiv(u_n)

c---------------------------------------------------------------------
c   perform the SSOR iterations
c---------------------------------------------------------------------
      call ssor(itmax)
c      call ssor(25)
c      call ssor(30)
c---------------------------------------------------------------------
c   compute the solution error
c---------------------------------------------------------------------
      call error(u_n)

c---------------------------------------------------------------------
c   compute the surface integral
c---------------------------------------------------------------------
      call pintgr(u_n)


c---------kai----------

      call c_print_timer(curr_rank)
c---------------------------------------------------------------------
c   verification test
c---------------------------------------------------------------------
      IF (id.eq.0) THEN
         call verify ( rsdnm, errnm, frc, class, verified )
         mflops = float(itmax)*(1984.77*float( nx0 )
     >        *float( ny0 )
     >        *float( nz0 )
     >        -10923.3*(float( nx0+ny0+nz0 )/3.)**2 
     >        +27770.9* float( nx0+ny0+nz0 )/3.
     >        -144010.)
     >        / (maxtime*1000000.)

         call print_results('LU', class, nx0,
     >     ny0, nz0, itmax, nnodes_compiled,
     >     num, maxtime, mflops, '          floating point', verified, 
     >     npbversion, compiletime, cs1, cs2, cs3, cs4, cs5, cs6, 
     >     '(none)')

      END IF

      if (.not.timeron) goto 999

      do i = 1, t_last
         t1(i) = timer_read(i)
      end do
      t1(t_rhs) = t1(t_rhs) - t1(t_exch)
      t1(t_last+2) = t1(t_lcomm)+t1(t_ucomm)+t1(t_rcomm)+t1(t_exch)
      t1(t_last+1) = t1(t_total) - t1(t_last+2)

      call MPI_Reduce(t1, tsum,  t_last+2, dp_type, MPI_SUM, 
     >                0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(t1, tming, t_last+2, dp_type, MPI_MIN, 
     >                0, MPI_COMM_WORLD, ierr)
      call MPI_Reduce(t1, tmaxg, t_last+2, dp_type, MPI_MAX, 
     >                0, MPI_COMM_WORLD, ierr)

      if (id .eq. 0) then
         write(*, 800) num
         do i = 1, t_last+2
            tsum(i) = tsum(i) / num
            write(*, 810) i, t_recs(i), tming(i), tmaxg(i), tsum(i)
         end do
      endif
 800  format(' nprocs =', i6, 11x, 'minimum', 5x, 'maximum', 
     >       5x, 'average')
 810  format(' timer ', i2, '(', A8, ') :', 3(2x,f10.4))

 999  continue
      call mpi_finalize(ierr)
      end


