#ifndef _PARTICLE_
#define _PARTICLE_
/*
 * MC
 *
 * [ ] vij
 * [ ] store/restore methods
 * [ ] tra_move()
161.825 -139.632 182.49 0 1 2054
160.388 -140.698 183.384 0 1 2055 *
 */

constexpr double PI = 3.14159265358979323846;
#include "./pvector.hpp"
#include "./params.hpp"
#include <cmath>
#include <chrono>
#include <iomanip> // for setprecision()
using ntype =  double;
class particle
{
  // to store particle's position

public:
  ntype sigma, epsilon, rc, m;
  ntype k, drmax;
  ntype theta;
  ntype u, r_sq, delta;
  pvector<ntype,3> r,rold, v, f;
  int bound_to = -1,bound_to_old = -1; // particle's position and velocity

  // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
  // this method should be called when initializing particles
  void set_vcut(void);

  void tra_move(pvector<ntype,3> delr)
    {
      r+=delr;
    }
  void store()
    {
      rold = r;
    }
  void restore()
    {
      r = rold;
    }

  void bond_store()
  {

    bound_to_old = bound_to;


  }

  void bond_restore()
  {
    bound_to = bound_to_old;
  }


  void set_mass(ntype mass)
   {
     m = mass;
   }
  void set_sigma(ntype sig)
    {
      sigma = sig;
    }
  void set_epsilon(ntype eps)
    {
      epsilon = eps;
    }
  void set_rcut(ntype rcut)
    {
      rc = rcut;
    }


  void set_k(ntype kappa)
    {
      k = kappa;
  }

  void set_drmax(ntype Ro)
    {
      drmax = Ro;
  }

  void set_delta(ntype del)
    {
      delta = del;
  }

  void set_theta(ntype thet)
    {
      theta = thet;
  }

  void set_u(ntype u_)
    {
      u = u_;
    }

  void set_r_sq(ntype r_)
    {
      r_sq = r_;
    }

  void set_bound_to(int i)
    {

     bound_to = i;

    }


  particle()
    {
      sigma=1.0;
      epsilon=1.0;
      k = 30.0;
      drmax = 1.0;
      rc=2.5;
      m = 1.0;
      theta = 6.0;
      u = 1.0;
      r_sq = 1.0;
      bound_to = -1;
      delta = 0.5;
    }

  };

class monomer: public particle
{

  public:

  int polymer_id, polymer_type;
  ntype vijLJ(monomer  P, pvector<ntype,3> L)
  {
       ntype ene;
       ntype sig = (sigma + P.sigma)/2.0;
       pvector<ntype,3> Dr;
       ntype dist;
       ntype r2, r6;
       ntype rij = pow(2.0, 1.0 / 6.0) * sig;

       Dr = r - P.r;
       // MINIMUM IMAGE CONVENTION
       Dr -= L.mulcw(rint(Dr.divcw(L)));
       dist = Dr.norm();




     if(dist<rij)
     {
      // r2 =5.54 * 5.54 *  sigma/dist * sigma/dist;
       r2 = (sig/dist) * (sig/dist);
       r6 = r2 * r2 * r2;
       ene = 4.0 * epsilon * (r6 * r6 - r6) + epsilon;
     }
     else
       ene = 0.0;

   return ene;

}



 ntype vijFENE_sym(monomer  P, pvector<ntype,3> L)
  {
       ntype ene;
       pvector<ntype,3> Dr;
       ntype dist;
     //  ntype  eq = 6.55;
       ntype eq = 0.3;
       ntype dr = 0.25;
       ntype kappa = 50;

       Dr = r - P.r;
       // MINIMUM IMAGE CONVENTION
       Dr -= L.mulcw(rint(Dr.divcw(L)));
       dist = Dr.norm();

     if(fabs(dist-eq) < dr)
     {
       // std::cout << "dist = " << dist << " eq = " << eq << " dr = " << dr << "\n";
       ene = - 0.5 * kappa * dr * dr * log(1.0 - ((dist- eq)/dr) * ((dist-eq)/dr));
       return ene;
     }
     else
      return NAN;

 }



ntype vijFENE_zero(monomer P, pvector<ntype,3> L) {
    pvector<ntype,3> Dr = r - P.r;
    Dr -= L.mulcw(rint(Dr.divcw(L)));  // MIC
    ntype dist = Dr.norm();

    ntype eq = 0.0;       // distanza di equilibrio
    ntype drmax = 2.5;     // estensione massima (solo in allungamento)
    ntype k = 20.0;        // costante del FENE

    if (dist <= eq)
        return 0.0;

    ntype x = (dist - eq) / drmax;

    if (x >= 1.0)
        return NAN;  // oppure: return 1e10;

    return -0.5 * k * drmax * drmax * log(1.0 - x * x);
}



 ntype vijBend(monomer  P1, monomer P2,  pvector<ntype,3> L)
 {
  ntype ene;
  pvector<ntype,3> Dr1;
  pvector<ntype,3> Dr2;
  ntype dist1;
  ntype dist2;

  Dr1 = r - P1.r;
  Dr2 = P2.r - r;
  // MINIMUM IMAGE CONVENTION
  Dr1 -= L.mulcw(rint(Dr1.divcw(L)));
  Dr2 -= L.mulcw(rint(Dr2.divcw(L)));
  dist1 = Dr1.norm();
  dist2 = Dr2.norm();
  ntype dist;

  ntype cost = Dr1*Dr2/(dist1*dist2);
  ntype n = 1.0;
  //theta /= 2.0 * n;
  cost = std::clamp(cost, -1.0, 1.0);

  ene = theta * (1.0 - cost);
  return ene;
 }


ntype vij_square(monomer P, pvector<ntype,3> L, int i, int j)
{
    pvector<ntype,3> Dr = r - P.r;
    Dr -= L.mulcw(rint(Dr.divcw(L)));
    ntype dist = Dr.norm();

    // Hard core repulsion
    if (dist < r_sq)
        return NAN;

    // Determina chi è il ligando
    const bool is_i_lig = (polymer_type == 1);
    const bool is_j_lig = (P.polymer_type == 1);

    // Caso illegale: il ligando è già legato ad un altro G4 → shell vietata
    if (is_i_lig && bound_to != -1 && bound_to != j && dist < r_sq + delta)
        return NAN;

    if (is_j_lig && P.bound_to != -1 && P.bound_to != i && dist < r_sq + delta)
        return NAN;

    // Caso legame valido: ligando libero o già legato reciprocamente
    if (
        (is_i_lig && (bound_to == -1 || bound_to == j)) ||
        (is_j_lig && (P.bound_to == -1 || P.bound_to == i))
    )
    {
        if (dist < r_sq + delta)
            return -u;
    }

    // Tutto il resto → energia nulla
    return 0.0;
}



};
#endif