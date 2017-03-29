c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine  add(u_new, u_old)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     addition of update to the vector u
c---------------------------------------------------------------------

      include 'header.h'
      double precision
     >u_new       (5,  -2:IMAX+1,-2:JMAX+1,-2:KMAX+1, maxcells),
     >u_old       (5,  -2:IMAX+1,-2:JMAX+1,-2:KMAX+1, maxcells)
      integer  c, i, j, k, m

      do     c = 1, ncells
         do     k = start(3,c), cell_size(3,c)-end(3,c)-1
            do     j = start(2,c), cell_size(2,c)-end(2,c)-1
               do     i = start(1,c), cell_size(1,c)-end(1,c)-1
                  do    m = 1, 5
                     u_new(m,i,j,k,c) = u_old(m,i,j,k,c) + 
     >                                  rhs(m,i,j,k,c)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end
