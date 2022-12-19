echo 'NAME,RESULT' > results.csv
parallel ./implementations/tsp_run ::: ./problems/ch130.tsp ./problems/d198.tsp ./problems/eil76.tsp ./problems/fl1577.tsp ./problems/pr439.tsp ./problems/rat783.tsp ./problems/u1060.tsp ./problems/kroA100.tsp ./problems/lin318.tsp ./problems/pcb442.tsp
#./implementations/tsp_run ./problems/d198.tsp
#./implementations/tsp_run ./problems/eil76.tsp
#./implementations/tsp_run ./problems/fl1577.tsp
#./implementations/tsp_run ./problems/kroA100.tsp
#./implementations/tsp_run ./problems/lin318.tsp
#./implementations/tsp_run ./problems/pcb442.tsp
#./implementations/tsp_run ./problems/pr439.tsp
#./implementations/tsp_run ./problems/rat783.tsp
#./implementations/tsp_run ./problems/u1060.tsp