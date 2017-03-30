
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine  adi(un,uo)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      include 'header.h'
      double precision
     >un       (5,  -2:IMAX+1,-2:JMAX+1,-2:KMAX+1, maxcells),
     >uo       (5,  -2:IMAX+1,-2:JMAX+1,-2:KMAX+1, maxcells)

       call copy_faces(uo)

       call txinvr

       call x_solve(uo)

       call y_solve(uo)

       call z_solve(uo)

       call add(un,uo)

       return
       end

