c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine  adi(u_new, u_old)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
        include 'header.h'
      double precision
     >u_new       (5,  -2:IMAX+1,-2:JMAX+1,-2:KMAX+1, maxcells),
     >u_old       (5,  -2:IMAX+1,-2:JMAX+1,-2:KMAX+1, maxcells)

      call copy_faces(u_old)

      call x_solve(u_old)

      call y_solve(u_old)

      call z_solve(u_old)

      call add(u_new, u_old)

      return
      end

