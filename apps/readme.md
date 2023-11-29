#ATPGraph

——————————依赖项———————————————
gcc >= 7.5.0 

openmpi || mpich

apt-get install openmpi-bin libopenmpi-dev libopenmpi-doc

——————————————数据转换——————————————
./bench2graph  <input bench file>  //output:  edgefile  gatetypefile

./bench2graph  b18_C.bench

./SNAPtoAdj  <input  edgefile>  <output csr filename> // edgefile convert to adj

./SNAPtoAdj  b18_C_graph.txt  b18_C_adj 

./preprocess.sh <bench name>

./preprocess.sh  b18_C

——————————————atalanta——————————————
./atalanta  <bench filename>

./atalanta  in002.bench

./atalanta  in003.bench

./atalanta  b17_C.bench

./atalanta  b18_C.bench

./atalanta  b19_C.bench

./atalanta  b20_C.bench

./atalanta  b21_C.bench

./atalanta  b22_C.bench

——————————————单机ATPGraph——————————————
./ATPG_FAN7 -t <gatetype filename> <adj filename>

./ATPG_FAN7 -t test17_gatetype.txt test17_adj

./ATPG_FAN7 -t in002_gatetype.txt in002_adj

./ATPG_FAN7 -t in003_gatetype.txt in003_adj

./ATPG_FAN7 -t in004_gatetype.txt in004_adj

./ATPG_FAN7 -t b17_C_gatetype.txt b17_C_adj 

./ATPG_FAN7 -t b18_C_gatetype.txt b18_C_adj 

./ATPG_FAN7 -t b19_C_gatetype.txt b19_C_adj 

./ATPG_FAN7 -t b20_C_gatetype.txt b20_C_adj 

./ATPG_FAN7 -t b21_C_gatetype.txt b21_C_adj 

./ATPG_FAN7 -t b22_C_gatetype.txt b22_C_adj 

./ATPG_FAN7 -t netcard_C_rmfg_gatetype.txt netcard_C_rmfg_adj

./ATPG_FAN7 -t b19c_10m_gatetype.txt b19c_10m_adj

——————————————分布式ATPGraph——————————————

mpirun -np <进程数> -hostfile <host filename> --allow-run-as-root(root用户下使用) ./ATPG_FAN_MPI12 -t <gatetype filename> <adj filename>

mpirun -np 16 -hostfile hostfile16 --allow-run-as-root ./ATPG_FAN_MPI12 -t b18_C_gatetype.txt b18_C_adj

mpirun -np 16 -hostfile hostfile16 --allow-run-as-root ./ATPG_FAN_MPI11 -t b18_C_gatetype.txt b18_C_adj

mpirun -np 16 -hostfile hostfile16 --allow-run-as-root ./ATPG_FAN_MPI8 -t b18_C_gatetype.txt b18_C_adj

——————————————单机多线程ATPGraph——————————————

mpirun -np <进程数>  --allow-run-as-root(root用户下使用) ./ATPG_FAN_MPI12 -t <gatetype filename> <adj filename>

mpirun -np 16  --allow-run-as-root ./ATPG_FAN_MPI8 -t b18_C_gatetype.txt b18_C_adj

——————————————性能分析——————————————
gprof ATPG_FAN gmon.out -q > time.txt

gprof atalanta gmon.out -q > time.txt

valgrind --tool=memcheck --leak-check=full ./ATPG_FAN7 -t in002_gatetype.txt in002_adj

valgrind --tool=memcheck --leak-check=full ./ATPG_FAN7 -t b18_C_gatetype.txt b18_C_adj 

valgrind --tool=massif ./ATPG_FAN -t in002_gatetype.txt in002_adj

——————————————其他算法——————————————

./Loop_detection -t test14_loop_gatetype.txt test14_loop_adj 

./ATPG_EST -t in002_gatetype.txt in002_adj
