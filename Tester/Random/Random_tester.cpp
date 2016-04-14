#include "Clone.h"
#include "ClonalExpansion.h"
#include "clonalFun.h"
#include "config.h"
#include "Random.h"
#include <iostream>         // std::cout,std::endl
#include <fstream>
#include <stdlib.h>         // EXIT_SUCCESS, EXIT_FAILURE
#include <stdexcept>        // std::exception, std::runtime_error
#include <memory>           // std::unique_ptr
#include <string>
#include <unistd.h>
#include <libgen.h>

using namespace std;

void do_Recovery_After_Replication();
void do_Poisson();
void do_Z();
void do_Binomial_Dy(int N, double dr);
void do_Binomial_NB(int N, double pr);
void do_Binomial_Mut(int N, double mr);
void do_uniform_mr(double mr);
void do_update_prolif(double pr);


char CWD[1024];
 
const string PATH = string(
                            dirname(
                                    dirname(
                                            getcwd( CWD, sizeof(CWD) )
                                            )
                                  )
                          ) + "/Test_Data/Random/";


//"./Test_Data/Random/";



int main(int argc, char**argv)
{
  cout << "SAVING FILES AT: " << PATH  << endl;
  do_Recovery_After_Replication(); 
  do_Poisson();
  do_Z();
  do_Binomial_Dy(200, 0.02);
  do_Binomial_NB(200, 0.03);
  do_Binomial_Mut(200,0.001);
  do_uniform_mr(0.001);
  do_update_prolif(0.01);


  return 0;

}

void do_Recovery_After_Replication()
{
  Random r;
  ofstream myfile;
  myfile.open (PATH+"Recovery_After_Replication.txt");

  for(int i = 0; i < 10000; i++)
    myfile << r.Recovery_After_Replication() << "\n";

  myfile.close();

}


void do_Poisson()
{
  Random r;
  ofstream myfile;
  myfile.open (PATH+"Poisson.txt");

  for(int i = 0; i < 10000; i++)
    myfile << r.Poisson() << "\n";

  myfile.close();

}

void do_Z()
{
  Random r;
  ofstream myfile;
  myfile.open (PATH+"Z.txt");

  for(int i = 0; i < 10000; i++)
    myfile << r.Z() << "\n";

  myfile.close();
}

void do_Binomial_Dy(int N, double dr)
{
  Random r;
  ofstream myfile;
  myfile.open (PATH+"Binomaial_Dying.txt");

  for(int i = 0; i < 10000; i++)
    myfile << r.Binomial_dying(N, dr) << "\n";

  myfile.close();
}

void do_Binomial_NB(int N, double pr)
{
  Random r;
  ofstream myfile;
  myfile.open (PATH+"Binomaial_Newborn.txt");

  for(int i = 0; i < 10000; i++)
    myfile << r.Binomial_dying(N, pr) << "\n";

  myfile.close();
}

void do_Binomial_Mut(int N, double mr)
{
  Random r;
  ofstream myfile;
  myfile.open (PATH+"Binomaial_Mut.txt");

  for(int i = 0; i < 10000; i++)
    myfile << r.Binomial_dying(N, mr) << "\n";

  myfile.close();
}

void do_uniform_mr( double mr)
{
  Random r;
  ofstream myfile;
  myfile.open (PATH+"Uniform_MR.txt");

  for(int i = 0; i < 10000; i++)
    myfile << r.Uniform_Mutation_Rate(mr) << "\n";

  myfile.close();
}

 void do_update_prolif(double pr)
 {
    Random r;
    ofstream myfile;
    myfile.open (PATH+"Update_Proliferation_Rate.txt");

  for(int i = 0; i < 10000; i++)
    myfile << r.Update_Proliferation_Rate(pr) << "\n";

  myfile.close();
 }









