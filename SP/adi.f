
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine  adi(un,uo)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      include 'header.h'
      double precision
     >un       (-2:IMAXP+1,-2:JMAXP+1,-2:KMAX+1, 5,maxcells),
     >uo       (-2:IMAXP+1,-2:JMAXP+1,-2:KMAX+1, 5,maxcells),

       call copy_faces(uo)

       call txinvr

       call x_solve(uo)

       call y_solve(uo)

       call z_solve(uo)

       call add(un,uo)

       return
       end

