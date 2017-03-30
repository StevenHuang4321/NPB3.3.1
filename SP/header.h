
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none

c---------------------------------------------------------------------
c The following include file is generated automatically by the
c "setparams" utility. It defines 
c      maxcells:      the square root of the maximum number of processors
c      problem_size:  12, 64, 102, 162 (for class T, A, B, C)
c      dt_default:    default time step for this problem size if no
c                     config file
c      niter_default: default number of iterations for this problem size
c---------------------------------------------------------------------

      include 'npbparams.h'

      integer           ncells, grid_points(3)
      pointer(ptr_ncells, ncells)
      pointer(ptr_grid_points, grid_points)
      common /global/   ptr_ncells, ptr_grid_points

      double precision  tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3, 
     >                  dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, 
     >                  dy5, dz1, dz2, dz3, dz4, dz5, dssp, dt, 
     >                  ce(5,13), dxmax, dymax, dzmax, xxcon1, xxcon2, 
     >                  xxcon3, xxcon4, xxcon5, dx1tx1, dx2tx1, dx3tx1,
     >                  dx4tx1, dx5tx1, yycon1, yycon2, yycon3, yycon4,
     >                  yycon5, dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
     >                  zzcon1, zzcon2, zzcon3, zzcon4, zzcon5, dz1tz1, 
     >                  dz2tz1, dz3tz1, dz4tz1, dz5tz1, dnxm1, dnym1, 
     >                  dnzm1, c1c2, c1c5, c3c4, c1345, conz1, c1, c2, 
     >                  c3, c4, c5, c4dssp, c5dssp, dtdssp, dttx1, bt,
     >                  dttx2, dtty1, dtty2, dttz1, dttz2, c2dttx1, 
     >                  c2dtty1, c2dttz1, comz1, comz4, comz5, comz6, 
     >                  c3c4tx3, c3c4ty3, c3c4tz3, c2iv, con43, con16
      pointer(ptr_tx1, tx1)
      pointer(ptr_tx2, tx2)
      pointer(ptr_tx3, tx3)
      pointer(ptr_ty1, ty1)
      pointer(ptr_ty2, ty2)
      pointer(ptr_ty3, ty3)
      pointer(ptr_tz1, tz1)
      pointer(ptr_tz2, tz2)
      pointer(ptr_tz3, tz3)
      pointer(ptr_dx1, dx1)
      pointer(ptr_dx2, dx2)
      pointer(ptr_dx3, dx3)
      pointer(ptr_dx4, dx4)
      pointer(ptr_dx5, dx5)
      pointer(ptr_dy1, dy1)
      pointer(ptr_dy2, dy2)
      pointer(ptr_dy3, dy3)
      pointer(ptr_dy4, dy4)
      pointer(ptr_dy5, dy5)
      pointer(ptr_dz1, dz1)
      pointer(ptr_dz2, dz2)
      pointer(ptr_dz3, dz3)
      pointer(ptr_dz4, dz4)
      pointer(ptr_dz5, dz5)
      pointer(ptr_dssp, dssp)
      pointer(ptr_dt, dt)
      pointer(ptr_ce, ce)
      pointer(ptr_dxmax, dxmax)
      pointer(ptr_dymax, dymax)
      pointer(ptr_dzmax, dzmax)
      pointer(ptr_xxcon1, xxcon1)
      pointer(ptr_xxcon2, xxcon2)
      pointer(ptr_xxcon3, xxcon3)
      pointer(ptr_xxcon4, xxcon4)
      pointer(ptr_xxcon5, xxcon5)
      pointer(ptr_dx1tx1, dx1tx1)
      pointer(ptr_dx2tx1, dx2tx1)
      pointer(ptr_dx3tx1, dx3tx1)
      pointer(ptr_dx4tx1, dx4tx1)
      pointer(ptr_dx5tx1, dx5tx1)
      pointer(ptr_yycon1, yycon1)
      pointer(ptr_yycon2, yycon2)
      pointer(ptr_yycon3, yycon3)
      pointer(ptr_yycon4, yycon4)
      pointer(ptr_yycon5, yycon5)
      pointer(ptr_dy1ty1, dy1ty1)
      pointer(ptr_dy2ty1, dy2ty1)
      pointer(ptr_dy3ty1, dy3ty1)
      pointer(ptr_dy4ty1, dy4ty1)
      pointer(ptr_dy5ty1, dy5ty1)
      pointer(ptr_zzcon1, zzcon1)
      pointer(ptr_zzcon2, zzcon2)
      pointer(ptr_zzcon3, zzcon3)
      pointer(ptr_zzcon4, zzcon4)
      pointer(ptr_zzcon5, zzcon5)
      pointer(ptr_dz1tz1, dz1tz1)
      pointer(ptr_dz2tz1, dz2tz1)
      pointer(ptr_dz3tz1, dz3tz1)
      pointer(ptr_dz4tz1, dz4tz1)
      pointer(ptr_dz5tz1, dz5tz1)
      pointer(ptr_dnxm1, dnxm1)
      pointer(ptr_dnym1, dnym1)
      pointer(ptr_dnzm1, dnzm1)
      pointer(ptr_c1c2, c1c2)
      pointer(ptr_c1c5, c1c5)
      pointer(ptr_c3c4, c3c4)
      pointer(ptr_c1345, c1345)
      pointer(ptr_conz1, conz1)
      pointer(ptr_c1, c1)
      pointer(ptr_c2, c2)
      pointer(ptr_c3, c3)
      pointer(ptr_c4, c4)
      pointer(ptr_c5, c5)
      pointer(ptr_c4dssp, c4dssp)
      pointer(ptr_c5dssp, c5dssp)
      pointer(ptr_dtdssp, dtdssp)
      pointer(ptr_dttx1, dttx1)
      pointer(ptr_bt, bt)
      pointer(ptr_dttx2, dttx2)
      pointer(ptr_dtty1, dtty1)
      pointer(ptr_dtty2, dtty2)
      pointer(ptr_dttz1, dttz1)
      pointer(ptr_dttz2, dttz2)
      pointer(ptr_c2dttx1, c2dttx1)
      pointer(ptr_c2dtty1, c2dtty1)
      pointer(ptr_c2dttz1, c2dttz1)
      pointer(ptr_comz1, comz1)
      pointer(ptr_comz4, comz4)
      pointer(ptr_comz5, comz5)
      pointer(ptr_comz6, comz6)
      pointer(ptr_c3c4tx3, c3c4tx3)
      pointer(ptr_c3c4ty3, c3c4ty3)
      pointer(ptr_c3c4tz3, c3c4tz3)
      pointer(ptr_c2iv, c2iv)
      pointer(ptr_con43, con43)
      pointer(ptr_con16, con16)

       common /constants/ ptr_tx1, ptr_tx2, ptr_tx3, ptr_ty1,
     >                   ptr_ty2, ptr_ty3, ptr_tz1, ptr_tz2, 
     >                  ptr_tz3,
     >                  ptr_dx1, ptr_dx2, ptr_dx3, ptr_dx4,
     >                  ptr_dx5, ptr_dy1, ptr_dy2, ptr_dy3, 
     >                  ptr_dy4,
     >                  ptr_dy5, ptr_dz1, ptr_dz2, ptr_dz3, 
     >                  ptr_dz4,
     >                  ptr_dz5, ptr_dssp, ptr_dt,
     >                  ptr_ce, ptr_dxmax, ptr_dymax,
     >                  ptr_dzmax, ptr_xxcon1, ptr_xxcon2,
     >                  ptr_xxcon3, ptr_xxcon4, ptr_xxcon5,
     >                  ptr_dx1tx1, ptr_dx2tx1, ptr_dx3tx1,
     >                  ptr_dx4tx1, ptr_dx5tx1, ptr_yycon1,
     >                  ptr_yycon2, ptr_yycon3, ptr_yycon4,
     >                  ptr_yycon5, ptr_dy1ty1, ptr_dy2ty1,
     >                  ptr_dy3ty1, ptr_dy4ty1, ptr_dy5ty1,
     >                  ptr_zzcon1, ptr_zzcon2, ptr_zzcon3,
     >                  ptr_zzcon4, ptr_zzcon5, ptr_dz1tz1,
     >                  ptr_dz2tz1, ptr_dz3tz1, ptr_dz4tz1,
     >                  ptr_dz5tz1, ptr_dnxm1, ptr_dnym1,
     >                  ptr_dnzm1, ptr_c1c2, ptr_c1c5, ptr_c3c4,
     >                  ptr_c1345, ptr_conz1, ptr_c1, ptr_c2,
     >                  ptr_c3, ptr_c4, ptr_c5, ptr_c4dssp,
     >                  ptr_c5dssp, ptr_dtdssp, ptr_dttx1, ptr_bt,
     >                  ptr_dttx2, ptr_dtty1, ptr_dtty2,
     >                  ptr_dttz1, ptr_dttz2, ptr_c2dttx1,
     >                  ptr_c2dtty1, ptr_c2dttz1, ptr_comz1,
     >                  ptr_comz4, ptr_comz5, ptr_comz6,
     >                  ptr_c3c4tx3, ptr_c3c4ty3, ptr_c3c4tz3,
     >                  ptr_c2iv, ptr_con43, ptr_con16

c      common /constants/ tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3,
c     >                  dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, 
c     >                  dy5, dz1, dz2, dz3, dz4, dz5, dssp, dt, 
c     >                  ce, dxmax, dymax, dzmax, xxcon1, xxcon2, 
c     >                  xxcon3, xxcon4, xxcon5, dx1tx1, dx2tx1, dx3tx1,
c     >                  dx4tx1, dx5tx1, yycon1, yycon2, yycon3, yycon4,
c     >                  yycon5, dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
c     >                  zzcon1, zzcon2, zzcon3, zzcon4, zzcon5, dz1tz1, 
c     >                  dz2tz1, dz3tz1, dz4tz1, dz5tz1, dnxm1, dnym1, 
c     >                  dnzm1, c1c2, c1c5, c3c4, c1345, conz1, c1, c2, 
c     >                  c3, c4, c5, c4dssp, c5dssp, dtdssp, dttx1, bt,
c     >                  dttx2, dtty1, dtty2, dttz1, dttz2, c2dttx1, 
c     >                  c2dtty1, c2dttz1, comz1, comz4, comz5, comz6, 
c     >                  c3c4tx3, c3c4ty3, c3c4tz3, c2iv, con43, con16

      integer           EAST, WEST, NORTH, SOUTH, 
     >                  BOTTOM, TOP

      pointer(ptr_EAST, EAST)
      pointer(ptr_WEST, WEST)
      pointer(ptr_NORTH, NORTH)
      pointer(ptr_SOUTH, SOUTH)
      pointer(ptr_BOTTOM, BOTTOM)
      pointer(ptr_TOP, TOP)

      parameter (EAST=2000, WEST=3000,      NORTH=4000, SOUTH=5000,
     >           BOTTOM=6000, TOP=7000)

      integer cell_coord (3,maxcells), cell_low (3,maxcells), 
     >        cell_high  (3,maxcells), cell_size(3,maxcells),
     >        predecessor(3),          slice    (3,maxcells),
     >        grid_size  (3),          successor(3),
     >        start      (3,maxcells), end      (3,maxcells)
      pointer(ptr_cell_coord, cell_coord)
      pointer(ptr_cell_low, cell_low)
      pointer(ptr_cell_high, cell_high)
      pointer(ptr_cell_size, cell_size)
      pointer(ptr_predecessor, predecessor)
      pointer(ptr_slice, slice)
      pointer(ptr_grid_size, grid_size)
      pointer(ptr_successor, successor)
      pointer(ptr_start, start)
      pointer(ptr_end, end)
 
      common /partition/ ptr_cell_coord, ptr_cell_low, 
     >                   ptr_cell_high, ptr_cell_size,
     >                   ptr_grid_size, ptr_successor,
     >                   ptr_predecessor, ptr_slice,
     >                   ptr_start, ptr_end

      integer IMAX, JMAX, KMAX, MAX_CELL_DIM, BUF_SIZE, IMAXP, JMAXP
      pointer(ptr_IMAX, IMAX)
      pointer(ptr_JMAX, JMAX)
      pointer(ptr_KMAX, KMAX)
      pointer(ptr_MAX_CELL_DIM, MAX_CELL_DIM)
      pointer(ptr_BUF_SIZE, BUF_SIZE)


      parameter (MAX_CELL_DIM = (problem_size/maxcells)+1)

      parameter (IMAX=MAX_CELL_DIM,JMAX=MAX_CELL_DIM,KMAX=MAX_CELL_DIM)
      parameter (IMAXP=IMAX/2*2+1,JMAXP=JMAX/2*2+1)

c---------------------------------------------------------------------
c +1 at end to avoid zero length arrays for 1 node
c---------------------------------------------------------------------
      parameter (BUF_SIZE=MAX_CELL_DIM*MAX_CELL_DIM*(maxcells-1)*60*2+1)

      double precision 
c     >   u       (-2:IMAXP+1,-2:JMAXP+1,-2:KMAX+1, 5,maxcells),
     >   u_n     (-2:IMAXP+1,-2:JMAXP+1,-2:KMAX+1, 5,maxcells),
     >   u_o     (-2:IMAXP+1,-2:JMAXP+1,-2:KMAX+1, 5,maxcells),
     >   us      (-1:IMAX,   -1:JMAX,   -1:KMAX,     maxcells),
     >   vs      (-1:IMAX,   -1:JMAX,   -1:KMAX,     maxcells),
     >   ws      (-1:IMAX,   -1:JMAX,   -1:KMAX,     maxcells),
     >   qs      (-1:IMAX,   -1:JMAX,   -1:KMAX,     maxcells),
     >   ainv    (-1:IMAX,   -1:JMAX,   -1:KMAX,     maxcells),
     >   rho_i   (-1:IMAX,   -1:JMAX,   -1:KMAX,     maxcells),
     >   speed   (-1:IMAX,   -1:JMAX,   -1:KMAX,     maxcells),
     >   square  (-1:IMAX,   -1:JMAX,   -1:KMAX,     maxcells),
     >   rhs     ( 0:IMAXP-1, 0:JMAXP-1, 0:KMAX-1, 5,maxcells),
     >   forcing ( 0:IMAXP-1, 0:JMAXP-1, 0:KMAX-1, 5,maxcells),
     >   lhs     ( 0:IMAXP-1, 0:JMAXP-1, 0:KMAX-1,15,maxcells),
     >   in_buffer(BUF_SIZE), out_buffer(BUF_SIZE),
     >   u_copy  (-2:IMAXP+1,-2:JMAXP+1,-2:KMAX+1, 5,maxcells),
     >   u_copy2 (-2:IMAXP+1,-2:JMAXP+1,-2:KMAX+1, 5,maxcells)

c      pointer(ptr_u, u)
      pointer(ptr_u_n, u_o)
      pointer(ptr_u_o, u_n)
      pointer(ptr_us, us)
      pointer(ptr_vs, vs)
      pointer(ptr_ws, ws)
      pointer(ptr_qs, qs)
      pointer(ptr_ainv, ainv)
      pointer(ptr_rho_i, rho_i)
      pointer(ptr_speed, speed)
      pointer(ptr_square, square)
      pointer(ptr_rhs, rhs)
      pointer(ptr_forcing, forcing)
      pointer(ptr_lhs, lhs)
      pointer(ptr_in_buffer, in_buffer)
      pointer(ptr_out_buffer, out_buffer)
      pointer(ptr_u_copy, u_copy)
      pointer(ptr_u_copy2, u_copy2)


      common /fields/  ptr_u_n, ptr_u_o,ptr_us, ptr_vs, ptr_ws, 
     >                 ptr_qs, ptr_ainv, ptr_rho_i, 
     >                 ptr_speed, ptr_square, 
     >                 ptr_rhs, ptr_forcing, ptr_lhs, 
     >                 ptr_in_buffer, ptr_out_buffer,
     >                 ptr_u_copy, ptr_u_copy2

      double precision cv(-2:MAX_CELL_DIM+1),   rhon(-2:MAX_CELL_DIM+1),
     >                 rhos(-2:MAX_CELL_DIM+1), rhoq(-2:MAX_CELL_DIM+1),
     >                 cuf(-2:MAX_CELL_DIM+1),  q(-2:MAX_CELL_DIM+1),
     >                 ue(-2:MAX_CELL_DIM+1,5), buf(-2:MAX_CELL_DIM+1,5)
      pointer(ptr_cv, cv)
      pointer(ptr_rhon, rhon)
      pointer(ptr_rhos, rhos)
      pointer(ptr_rhoq, rhoq)
      pointer(ptr_cuf, cuf)
      pointer(ptr_q, q)
      pointer(ptr_ue, ue)
      pointer(ptr_buf, buf)
      common /work_1d/ ptr_cv, ptr_rhon, ptr_rhos, 
     >                 ptr_rhoq, ptr_cuf, ptr_q, ptr_ue, ptr_buf

      integer  west_size, east_size, bottom_size, top_size,
     >         north_size, south_size, start_send_west, 
     >         start_send_east, start_send_south, start_send_north,
     >         start_send_bottom, start_send_top, start_recv_west,
     >         start_recv_east, start_recv_south, start_recv_north,
     >         start_recv_bottom, start_recv_top
      pointer(ptr_west_size, west_size)
      pointer(ptr_east_size, east_size)
      pointer(ptr_bottom_size, bottom_size)
      pointer(ptr_top_size, top_size)
      pointer(ptr_north_size, north_size)
      pointer(ptr_south_size, south_size)
      pointer(ptr_start_send_west, start_send_west)
      pointer(ptr_start_send_east, start_send_east)
      pointer(ptr_start_send_south, start_send_south)
      pointer(ptr_start_send_north, start_send_north)
      pointer(ptr_start_send_bottom, start_send_bottom)
      pointer(ptr_start_send_top, start_send_top)
      pointer(ptr_start_recv_west, start_recv_west)
      pointer(ptr_start_recv_east, start_recv_east)
      pointer(ptr_start_recv_south, start_recv_south)
      pointer(ptr_start_recv_north, start_recv_north)
      pointer(ptr_start_recv_bottom, start_recv_bottom)
      pointer(ptr_start_recv_top, start_recv_top)

       common /box/ ptr_west_size,ptr_east_size, ptr_bottom_size,
     >             ptr_top_size, ptr_north_size, ptr_south_size, 
     >             ptr_start_send_west, ptr_start_send_east, 
     >             ptr_start_send_south, ptr_start_send_north, 
     >             ptr_start_send_bottom, ptr_start_send_top,
     >             ptr_start_recv_west, ptr_start_recv_east, 
     >             ptr_start_recv_south, ptr_start_recv_north, 
     >             ptr_start_recv_bottom, ptr_start_recv_top

      integer t_total, t_rhs, t_xsolve, t_ysolve, t_zsolve, t_bpack, 
     >        t_exch, t_xcomm, t_ycomm, t_zcomm, t_last
      pointer(ptr_t_total, t_total)
      pointer(ptr_t_rhs, t_rhs)
      pointer(ptr_t_xsolve, t_xsolve)
      pointer(ptr_t_ysolve, t_ysolve)
      pointer(ptr_t_zsolve, t_zsolve)
      pointer(ptr_t_bpack, t_bpack)
      pointer(ptr_t_exch, t_exch)
      pointer(ptr_t_xcomm, t_xcomm)
      pointer(ptr_t_ycomm, t_ycomm)
      pointer(ptr_t_zcomm, t_zcomm)
      pointer(ptr_t_last, t_last)

      parameter (t_total=1, t_rhs=2, t_xsolve=3, t_ysolve=4, 
     >        t_zsolve=5, t_bpack=6, t_exch=7, t_xcomm=8, 
     >        t_ycomm=9, t_zcomm=10, t_last=10)
      logical timeron
      pointer(ptr_timeron, timeron)
      common /tflags/ ptr_timeron
