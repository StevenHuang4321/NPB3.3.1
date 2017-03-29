
#!/bin/bash

echo "yes!!!"
mpirun -np 4 -host n03,n04,n05,n07 ./bin/cg.C.4 2>&1 | tee output_bt_16lat.log
mpirun -np 4 -host n03,n04,n05,n07 ./bin/bt.C.4 2>&1 | tee output_bt_16lat.log
mpirun -np 4 -host n03,n04,n05,n07 ./bin/lu.C.4 2>&1 | tee output_lu_16lat.log
mpirun -np 4 -host n03,n04,n05,n07 ./bin/sp.C.4 2>&1 | tee output_sp_16lat.log
mpirun -np 4 -host n03,n04,n05,n07 ./bin/ft.C.4 2>&1 | tee output_ft_16lat.log
mpirun -np 4 -host n03,n04,n05,n07 ./bin/mg.C.4 2>&1 | tee output_mg_16lat.log
echo "Done!!!"
