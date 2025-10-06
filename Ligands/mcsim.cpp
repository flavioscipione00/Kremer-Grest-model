#include "./sim.hpp"
#include <fstream>
#include <chrono>

int main(int argc, char **argv)
{

  mcsim<monomer> mc;

  mc.init_rng(-1);
  mc.prepare_initial_conf_stacked();
  mc.run(2);

  return 0;

};