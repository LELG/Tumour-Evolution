 g++ -c config.cpp --std=c++11
 g++ -c Clone.cpp --std=c++11
 g++ -c clonalFun.cpp --std=c++11
  g++ -c inout_funs.cpp --std=c++11
 g++ -c main.cpp --std=c++11

 g++ config.o Clone.o clonalFun.o inout_funs.o main.o -o TumourEvolution

 rm config.o
 rm Clone.o
 rm clonalFun.o
 rm inout_funs.o
 rm main.o

 ./TumourEvolution

 rm TumourEvolution