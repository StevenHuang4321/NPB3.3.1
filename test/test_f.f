      program test_f
      implicit none
      REAL(4) x,y,z
      pointer(ptr_x, x)
      pointer(ptr_y, y)
      pointer(ptr_z, z)
c      REAL(8) ydbl
c      pointer(ptr_ydbl, ydbl)
      COMMON / Really / ptr_x, ptr_y, ptr_z

c      x = 999
c      y = 100
c      z = 100
      print *, ptr_x
c      print *, loc(x),loc(y),loc(z),loc(ydbl)
      call cmalloc(ptr_x)
      call cmalloc(ptr_y)
      call cmalloc(ptr_z)
      print *, ptr_x, ptr_y, ptr_z

      x= 99
      y =100
      z =101
      print *, x,y,z
c      call cmove(ptr_x, sizeof(x))
c      call cmove(ptr_y, sizeof(y))
c      call cmove(ptr_z, sizeof(z))
      print *, x,y,z
      print *, ptr_x, ptr_y, ptr_z
      call foo
      end program test_f
      
      subroutine foo
      implicit none
      COMMON / Really / ptr_x, ptr_y, ptr_z
      REAL(4) x,y,z
      pointer(ptr_x, x)
      pointer(ptr_y, y)
      pointer(ptr_z, z)
      print *, ptr_x
      print *, x
      end subroutine
