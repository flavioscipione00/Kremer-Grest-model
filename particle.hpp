#ifndef _PARTICLE_
#define _PARTICLE_
/*
 * MC
 *
 * [ ] vij
 * [ ] store/restore methods
 * [ ] tra_move()
 *
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
protected:
  ntype vcut;
public:
  ntype sigma, epsilon, rc, m;
  ntype k, drmax;
  ntype theta;
  pvector<ntype,3> r,rold, v, f; // particle's position and velocity

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


  void set_theta(ntype thet)
    {
      theta = thet;
  }


  particle()
    {
      sigma=1.0;
      epsilon=1.0;
      k = 30.0;
      drmax = 1.0;
      rc=2.5;
      vcut = 0.0;
      m = 1.0;
      theta = 6.0;
    }

  };

class monomer: public particle
{

  public:

  int polymer_id, polymer_type;
  ntype vijLJ(monomer  P, pvector<ntype,3> L, int i)
  {
       ntype ene;
       pvector<ntype,3> Dr;
       ntype dist;
       ntype r2, r6;
       ntype rcut = 0.0;

      rcut = rc; // Use the class member rc



       Dr = r - P.r;
       // MINIMUM IMAGE CONVENTION
       Dr -= L.mulcw(rint(Dr.divcw(L)));
       dist = Dr.norm();




     if(dist<rcut)
     {
       r2 = sigma/dist * sigma/dist;
       r6 = r2 * r2 * r2;
       ene = 4.0 * epsilon * (r6 * r6 - r6) + epsilon;
     }
     else
       ene = 0.0;

   return ene;

}

   ntype vijStack(monomer  P, pvector<ntype,3> L)
   {
        ntype ene;
        pvector<ntype,3> Dr;
        ntype dist;
        ntype r2, r6;

        Dr = r - P.r;
        // MINIMUM IMAGE CONVENTION
        Dr -= L.mulcw(rint(Dr.divcw(L)));
        dist = Dr.norm();




      if(dist<2.5)
      {
        r2 = sigma/dist * sigma/dist;
        r6 = r2 * r2 * r2;
        ene = 4.0 * epsilon * (r6 * r6 - r6) - 4.0* epsilon * (pow(sigma/2.5,12)-pow(sigma/2.5,6));
        ene *= 0.4;
      }
      else
        ene = 0.0;

    return ene;

   }


ntype vijFENE(monomer P, pvector<ntype,3> L) {
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
  cost = std::clamp(cost, -1.0, 1.0);
 // ntype t = acos(cost);

 // if(cost < -0.5 || (cost > 0.71 && cost < 0.9))  return NAN;

  ene = theta * (1.0 - cost);

  return ene;
 }


  void set_vcut(void)
    {
      vcut = 4.0*epsilon*(pow (5.54 *sigma/rc,12)-pow(5.54 * sigma/rc,6));
    }



};



#endif