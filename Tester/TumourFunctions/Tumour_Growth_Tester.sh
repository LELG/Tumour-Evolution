  #!/bin/sh

 # Get current and core paths
current_path=$(pwd)
core_path="$(dirname "$(dirname "$current_path")")"

#Compilation
# Include the files to not cause a comnpilation error
g++ -Wall -shared -I/usr/local/inlcude -I"$core_path" -c Tumour_Growth_Tester.cpp --std=c++11

g++ -Wall -shared -I/usr/local/inlcude -c "$core_path/config.cpp" --std=c++11
g++ -Wall -shared -I/usr/local/inlcude -c "$core_path/Clone.cpp" --std=c++11
g++ -Wall -shared -I/usr/local/inlcude -c "$core_path/clonalFun.cpp" --std=c++11
g++ -Wall -shared -I/usr/local/inlcude -c "$core_path/inout_funs.cpp" --std=c++11
g++ -Wall -shared -I/usr/local/inlcude -c "$core_path/ClonalExpansion.cpp" --std=c++11 
g++ -Wall -shared -I/usr/local/inlcude -c "$core_path/Random.cpp" --std=c++11 


g++ -lgsl  Tumour_Growth_Tester.o config.o Clone.o clonalFun.o inout_funs.o ClonalExpansion.o Random.o -o Tumour_Growth_Tester 

 rm config.o
 rm Clone.o
 rm clonalFun.o
 rm inout_funs.o
 rm ClonalExpansion.o
 rm Random.o
 rm Tumour_Growth_Tester.o


 ./Tumour_Growth_Tester

rm Tumour_Growth_Tester