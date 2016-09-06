  #!/bin/sh

 # Get current and core paths
current_path=$(pwd)
core_path="$(dirname "$(dirname "$current_path")")"

#Compilation
# Include the files to not cause a comnpilation error
g++-4.9 -Wall -I/usr/local/inlcude -I/usr/local/Cellar/open-mpi/1.10.2/include -I"$core_path" -c Main_Tester.cpp --std=c++11

g++-4.9 -Wall -I/usr/local/inlcude -I/usr/local/Cellar/open-mpi/1.10.2/include -c "$core_path/config.cpp" --std=c++11
g++-4.9 -Wall -I/usr/local/inlcude -I/usr/local/Cellar/open-mpi/1.10.2/include -c "$core_path/Clone.cpp" --std=c++11
g++-4.9 -Wall -I/usr/local/inlcude -I/usr/local/Cellar/open-mpi/1.10.2/include -c "$core_path/clonalFun.cpp" --std=c++11
g++-4.9 -Wall -I/usr/local/inlcude -I/usr/local/Cellar/open-mpi/1.10.2/include -c "$core_path/inout_funs.cpp" --std=c++11
g++-4.9 -Wall -I/usr/local/inlcude -I/usr/local/Cellar/open-mpi/1.10.2/include -c "$core_path/ClonalExpansion.cpp" --std=c++11 
g++-4.9 -Wall -I/usr/local/inlcude -I/usr/local/Cellar/open-mpi/1.10.2/include -c "$core_path/Random.cpp" --std=c++11 

g++-4.9 -L/usr/local/opt/libevent/lib -L/usr/local/Cellar/open-mpi/1.10.2/lib -lmpi_cxx -lmpi -lgsl -lgsl -lgslcblas  Main_Tester.o config.o Clone.o clonalFun.o inout_funs.o ClonalExpansion.o Random.o -o Main_Tester 

 rm config.o
 rm Clone.o
 rm clonalFun.o
 rm inout_funs.o
 rm ClonalExpansion.o
 rm Random.o
 rm Main_Tester.o
 

 mpirun -n 2 ./Main_Tester _Test

rm Main_Tester
