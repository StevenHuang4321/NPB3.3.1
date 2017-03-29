c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  work_lhs.h
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      double precision fjac(5, 5, -2:MAX_CELL_DIM+1),
     >                 njac(5, 5, -2:MAX_CELL_DIM+1),
     >                 lhsa(5, 5, -1:MAX_CELL_DIM),
     >                 lhsb(5, 5, -1:MAX_CELL_DIM),
     >                 tmp1, tmp2, tmp3
      pointer(ptr_fjac, fjac)
      pointer(ptr_njac, njac)
      pointer(ptr_lhsa, lhsa)
      pointer(ptr_lhsb, lhsb)
      pointer(ptr_tmp1, tmp1)
      pointer(ptr_tmp2, tmp2)
      pointer(ptr_tmp3, tmp3)
      common /work_lhs/ ptr_fjac, ptr_njac, ptr_lhsa, 
     >                 ptr_lhsb, ptr_tmp1, ptr_tmp2, ptr_tmp3
