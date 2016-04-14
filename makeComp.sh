 g++ -Wall -I/usr/local/inlcude -c config.cpp --std=c++11
 g++ -Wall -I/usr/local/inlcude -c Clone.cpp --std=c++11
 g++ -Wall -I/usr/local/inlcude -c clonalFun.cpp --std=c++11
 g++ -Wall -I/usr/local/inlcude -c inout_funs.cpp --std=c++11
 g++ -Wall -I/usr/local/inlcude -c main.cpp --std=c++11
 g++ -Wall -I/usr/local/inlcude -c ClonalExpansion.cpp --std=c++11 
 g++ -Wall -I/usr/local/inlcude -c Random.cpp --std=c++11 

 g++ -lgsl config.o Clone.o clonalFun.o inout_funs.o ClonalExpansion.o Random.o main.o -o TumourEvolution 

 rm config.o
 rm Clone.o
 rm clonalFun.o
 rm inout_funs.o
 rm ClonalExpansion.o
 rm Random.o
 rm main.o

 ./TumourEvolution

 rm TumourEvolution
 