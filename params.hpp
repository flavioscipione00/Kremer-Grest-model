#ifndef _PARAMS_
#define _PARAMS_
#include <fstream>
#include <iostream>
#include <string>
#include <map> // For storing parameters
#include <sstream>
#include <unordered_map>
#include <functional>
#include "./pvector.hpp"

class simpars
{
  using ntype=  double;
public:
  int nx, ny, nz; /* nx*ny*nz particelle */
  ntype T, P, mu; // temperature and pressure
  int Np; // numero di particelle
  long int maxadjstps, eqstps, adjstps, save_mgl_snapshot;
  long int savemeasure, outstps, totsteps; // Nsteps = simulations steps, outstps steps print something on current simulation status
  double rho, rc; // density
  int seed; // -1 means random
  pvector<ntype, 3> L; // box
  ntype sigma, epsilon, mass; // Lennard-Jones parameters
  ntype k, drmax; // FENE parameters
  ntype  deltra, deltra_cls,deltra_rot, deltra_2; // parameter of MC moves
  int Nx, Ny;
  ntype theta;

  simpars()
    {
      nx = 2;  // number of particles along each direction
      ny = 2;
      nz = 2;
      sigma=1.0;
      epsilon=1.0;
      rho = 0.5;
      rc = 4.0;
      seed=0;
      mass=1.0;
      adjstps = 1;
      maxadjstps = 70;
      eqstps=500;
      totsteps = 2000;
      save_mgl_snapshot = 100; 
      savemeasure=20;
      outstps=200;
      T = 2.0;
      P = 3.838; //se P*=\beta*P*v0, 1 < P* < 10 dove v0 è il volume di una particella 
      deltra = 0.2;
      deltra_2 = 0.2;
      deltra_cls = 0.02;
      deltra_rot = 0.314;
      k = 30.0;
      drmax = 1.5;
      theta = 1.0;

    }


//creating a loading function to be used in sim.hpp 
 void loadParameters(const std::string& filename) {
        // Map parameter names to their references
        std::unordered_map<std::string, std::function<void(const std::string&)>> paramMap = {
            {"nx", [&](const std::string& value) { nx = std::stoi(value); }},
            {"ny", [&](const std::string& value) { ny = std::stoi(value); }},
            {"nz", [&](const std::string& value) { nz = std::stoi(value); }},
            {"sigma", [&](const std::string& value) { sigma = std::stod(value); }},
            {"epsilon", [&](const std::string& value) { epsilon = std::stod(value); }},
            {"rho", [&](const std::string& value) { rho = std::stod(value); }},
            {"rc", [&](const std::string& value) { rc = std::stod(value); }},
            {"seed", [&](const std::string& value) { seed = std::stoi(value); }},
            {"mass", [&](const std::string& value) { mass = std::stod(value); }},
            {"adjstps", [&](const std::string& value) { adjstps = std::stoi(value); }},
            {"maxadjstps", [&](const std::string& value) { maxadjstps = std::stoi(value); }},
            {"eqstps", [&](const std::string& value) { eqstps = std::stoi(value); }},
            {"totsteps", [&](const std::string& value) { totsteps = std::stoi(value); }},
            {"save_mgl_snapshot", [&](const std::string& value) { save_mgl_snapshot = std::stoi(value); }},
            {"savemeasure", [&](const std::string& value) { savemeasure = std::stoi(value); }},
            {"outstps", [&](const std::string& value) { outstps = std::stoi(value); }},
            {"T", [&](const std::string& value) { T = std::stod(value); }},
            {"P", [&](const std::string& value) { P = std::stod(value); }},
            {"deltra", [&](const std::string& value) { deltra = std::stod(value); }},
            {"deltra_2", [&](const std::string& value) { deltra_2 = std::stod(value);}},
            {"mu", [&](const std::string& value) { mu = std::stod(value); }},
            {"Nx", [&](const std::string& value) { Nx = std::stod(value); }},
            {"Ny", [&](const std::string& value) { Ny = std::stod(value); }},
            {"k", [&](const std::string& value) { k = std::stod(value); }},
            {"drmax", [&](const std::string& value) { drmax = std::stod(value); }},
            {"theta", [&](const std::string& value) { theta = std::stod(value); }},
            {"deltra_cls", [&](const std::string& value) { deltra_cls = std::stod(value); }},
            {"deltra_rot", [&](const std::string& value) { deltra_rot = std::stod(value); }},
        

        };

        std::ifstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return;
        }

        std::string line;
        while (std::getline(file, line)) {
            std::istringstream ss(line);
            std::string key, value;
            if (std::getline(ss, key, '=') && std::getline(ss, value)) {
                key.erase(key.find_last_not_of(" \t") + 1); // Trim spaces
                value.erase(0, value.find_first_not_of(" \t")); // Trim spaces
                auto it = paramMap.find(key);
                if (it != paramMap.end()) {
                    it->second(value); // Call the corresponding lambda to set the value
                } else {
                    std::cerr << "Warning: Unknown parameter '" << key << "'" << std::endl;
                }
            }
        }
    }




};
#endif
