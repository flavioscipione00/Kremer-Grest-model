#include "./sim.hpp"
#include <fstream>
#include <chrono>  //useful to check times 

int main(int argc, char **argv)
{

   mcsim<monomer> mc;

//setting the Random Number Generator to a casual seed
   mc.init_rng(-1);

//setting  the stating configuration:
// 1) assign per particle properties
// 2) seed and grow the polymers on a SC lattice

   mc.prepare_initial_conf();

//starting the run, 0 is the suffix of the ".txt" snapshot files
   mc.run(0);

  
  return 0;

};
