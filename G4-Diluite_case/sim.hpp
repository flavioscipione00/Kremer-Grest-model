#ifndef _SIMCLASS_
#define _SIMCLASS_

#include <vector>
#include <string>
#include <fstream>
#include "./params.hpp"
#include "./particle.hpp"
#include "./pmatrix.hpp"
#include "./randnumgen.hpp"
#include <iomanip>

using ntype=double;

template<typename particle_type>
class sim
{
  using simp = simpars;
protected:
  simp pars;  // parameters
  std::vector<monomer> parts; //std::vector of particles
  std::vector<std::vector<int>> pol; //std::vector that contains the index of particles j that belong to polymer i
  int cell_nx , cell_ny, cell_nz , n_cell;
  pvector <ntype,3> cell_size;
  int long_polymers, short_polymers, total_polymers;


  pmatrixq<ntype,3> Identity = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  std::vector<int> head;
  std::vector<int> linked_list;
  std::vector<int> backup_linked_list;
  std::vector<int> head_backup;

 //if opt=1 calculate energies only for i < j
  ntype calcenergyi(int i, int opt=0)
    {
      int j;
      ntype enei=0.0;
      for (j=0; j < parts.size(); j++) // pars.Np it's the total number of particles
        {
          if (opt==1 && i >= j)
            continue;
          if (i==j)
            continue;
          // la classe particelle deve implementare un metodo vij per il calcolo dell'energia d'interazione
          enei += parts[i].vijLJ(parts[j], pars.L, 0);
          if (parts[i].polymer_id == parts[j].polymer_id && abs(i-j) == 1) enei += parts[i].vijFENE(parts[j], pars.L);

          // pars.L è un vettore con i lati del box
        }

      enei += bending_energy(i);

      return enei;
    }


    ntype bending_energy(int i)
    {
        if (i < 0 || i >= parts.size()) return 0.0;

        ntype ene = 0.0;

        // (i-2, i-1, i)
        if (i >= 2 &&
            parts[i].polymer_id == parts[i - 1].polymer_id &&
            parts[i - 1].polymer_id == parts[i - 2].polymer_id)
        {
            ene += parts[i - 1].vijBend(parts[i - 2], parts[i], pars.L);
        }

        // (i-1, i, i+1)
        if (i >= 1 && i + 1 < parts.size() &&
            parts[i].polymer_id == parts[i - 1].polymer_id &&
            parts[i].polymer_id == parts[i + 1].polymer_id)
        {
            ene += parts[i].vijBend(parts[i - 1], parts[i + 1], pars.L);
        }

        // (i, i+1, i+2)
        if (i + 2 < parts.size() &&
            parts[i].polymer_id == parts[i + 1].polymer_id &&
            parts[i + 1].polymer_id == parts[i + 2].polymer_id)
        {
            ene += parts[i + 1].vijBend(parts[i], parts[i + 2], pars.L);
        }

        return ene;
    }


  ntype totenergy()
    {
      ntype ene=0.0;
      for (auto i=0; i < parts.size(); i++)
        {
          ene+=calcenergyi(i, 1);
        }

      return ene;
    }

  void pbc(int i)
    {

      parts[i].r -= pars.L.mulcw(rint(parts[i].r.divcw(pars.L)));

    }

  void save_mgl_snapshot(long int t, int c)
    {
       std::fstream f;
       std::string s;

       s = "pos-" + std::to_string(t) + ".txt";
       f.open(s, std::ios::out|std::ios::trunc);
       for (int i=0; i < pars.Np; i++)
         {
           f << setprecision(10) <<parts[i].r(0) << " " << parts[i].r(1) << " " <<
             parts[i].r(2) << " " << parts[i].polymer_id << "\n";
         }
       f.close();
    }

public:


  void init_rng(int n)
    {
      if(n<0)
       {
        rng.rseed();
       }
      else rng.seed(n);
    }


  void init_cells(void) {
    cell_nx = 20;
    cell_ny = 20;
    cell_nz = 20;
    n_cell = cell_nx * cell_ny * cell_nz;


    cell_size(0) = pars.L(0) / cell_nx;
    cell_size(1) = pars.L(1) / cell_ny;
    cell_size(2) = pars.L(2) / cell_nz;

    if (n_cell <= 0) {
        throw std::runtime_error("Invalid number of cells: " + std::to_string(n_cell));
    }


    head.assign(n_cell, -1);
    linked_list.assign(pars.Np, -1);
    backup_linked_list.assign(pars.Np, -1);
    head_backup.assign(n_cell, -1);

  }





  void particle_to_cells(void)
   {

     int x = 0, y = 0, z = 0;
     std::cout << "cell_nx: " << cell_nx << ", cell_ny: " << cell_ny << ", cell_nz: " << cell_nz << "\n";
    for( int i = 0; i < parts.size(); i++)
    {
      x = floor((parts[i].r(0) + 0.5 * pars.L(0))/cell_size(0));
      y = floor((parts[i].r(1)+ 0.5 * pars.L(1))/cell_size(1));
      z = floor((parts[i].r(2)+ 0.5 * pars.L(2))/cell_size(2));

      x = std::max(0, std::min(cell_nx - 1, x));
      y = std::max(0, std::min(cell_ny - 1, y));
      z = std::max(0, std::min(cell_nz - 1, z));

       int index = x + y * cell_nx + z * cell_ny * cell_nx;

        if (index < 0 || index >= head.size()) {
            std::cerr << "Invalid cell index for particle " << i << ": " << index << "\n";
            continue;
        }


    linked_list[i] = head[index];
        head[index] = i;
    }

      for (int c = 0; c < n_cell; c++) {
    int loop_counter = 0;
    for (int p = head[c]; p >= 0; p = linked_list[p]) {
        loop_counter++;
        if (p < 0 || p >= pars.Np) {
            std::cerr << "Invalid particle index " << p << " in cell " << c << "\n";
            throw std::runtime_error("Corrupt linked list structure");
        }
        if (loop_counter > pars.Np) {
            std::cerr << "Infinite loop detected in cell " << c << "\n";
            throw std::runtime_error("Infinite loop in linked list");
        }
    }
}
}

ntype gauss(void)
{
  double x1, x2;
  do
    {
      x1 = rng.ranf();
      x2 = rng.ranf();
    }
  while (x1==0 || x2==0);
  return cos(2*M_PI*x2)*sqrt(-2.0*log(x1));
}



void update_cells(int i) {
    ntype shifted_old_x = parts[i].rold(0) + 0.5 * pars.L(0);
    ntype shifted_old_y = parts[i].rold(1) + 0.5 * pars.L(1);
    ntype shifted_old_z = parts[i].rold(2) + 0.5 * pars.L(2);

    ntype shifted_new_x = parts[i].r(0) + 0.5 * pars.L(0);
    ntype shifted_new_y = parts[i].r(1) + 0.5 * pars.L(1);
    ntype shifted_new_z = parts[i].r(2) + 0.5 * pars.L(2);


    int old_x = floor(shifted_old_x / cell_size(0));
    int old_y = floor(shifted_old_y / cell_size(1));
    int old_z = floor(shifted_old_z / cell_size(2));

    int new_x = floor(shifted_new_x / cell_size(0));
    int new_y = floor(shifted_new_y / cell_size(1));
    int new_z = floor(shifted_new_z / cell_size(2));

    new_x = std::max(0, std::min(cell_nx - 1, new_x));
    new_y = std::max(0, std::min(cell_ny - 1, new_y));
    new_z = std::max(0, std::min(cell_nz - 1, new_z));
    old_x = std::max(0, std::min(cell_nx - 1, old_x));
    old_y = std::max(0, std::min(cell_ny - 1, old_y));
    old_z = std::max(0, std::min(cell_nz - 1, old_z));


    int old_index = old_x + old_y * cell_nx + old_z * cell_nx * cell_ny;
    int new_index = new_x + new_y * cell_nx + new_z * cell_nx * cell_ny;


    if (old_index < 0 || old_index >= n_cell || new_index < 0 || new_index >= n_cell) {
        std::cerr << "Invalid index for particle " << i << ": old_index=" << old_index << ", new_index=" << new_index << "\n";
        return;
    }

    // Se la particella non cambia cella, non fare nulla
    if (old_index == new_index) return;

    // Rimuovi la particella dalla cella vecchia
    if (head[old_index] == i) {
        head[old_index] = linked_list[i];
    } else {
        int prev = head[old_index];
        while (linked_list[prev] != i) {
            prev = linked_list[prev];
        }
        linked_list[prev] = linked_list[i];
    }

    // Aggiungi la particella alla nuova cella
    linked_list[i] = head[new_index];
    head[new_index] = i;

}


ntype energy_linked_cells(int i, ntype rc, int opt = 0)
{
    int x = 0, y = 0, z = 0;
    int index = 0;
    ntype ene = 0.0;

    // Neighbor search grid index
    x = floor((parts[i].r(0) + 0.5 * pars.L(0)) / cell_size(0));
    y = floor((parts[i].r(1) + 0.5 * pars.L(1)) / cell_size(1));
    z = floor((parts[i].r(2) + 0.5 * pars.L(2)) / cell_size(2));

    // Loop over neighbor cells
    for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++)
            for (int dz = -1; dz <= 1; dz++)
            {
                int cx = (x + dx + cell_nx) % cell_nx;
                int cy = (y + dy + cell_ny) % cell_ny;
                int cz = (z + dz + cell_nz) % cell_nz;

                if (cx < 0 || cy < 0 || cz < 0 || cx >= cell_nx || cy >= cell_ny || cz >= cell_nz) continue;

                index = cx + cy * cell_nx + cz * cell_ny * cell_nx;
                if (index < 0 || index >= head.size()) continue;

                for (int j = head[index]; j >= 0; j = linked_list[j])
                {
                    if (i == j) continue;

                    // Evita doppio conteggio se opt == 1
                    if (opt == 1 && j >= i) continue;


                    // Determina se sono nella stessa catena
                    bool same_chain = (parts[i].polymer_id == parts[j].polymer_id);

                    // Per opt == 2 (cluster move): ignora interazioni intra-catena
                    if (opt == 2 && same_chain) continue;

                    // MIC
                    pvector<ntype, 3> Dr = parts[i].r - parts[j].r;
                    Dr -= pars.L.mulcw(rint(Dr.divcw(pars.L)));
                    ntype dist_sq = Dr.norm() * Dr.norm();

                    if( dist_sq > pars.rc * pars.rc *pars.sigma * pars.sigma) continue;



                    // === WCA ===
                    if (!same_chain || (abs(i - j) > 1 && opt != 2))
                        ene += parts[i].vijLJ(parts[j], pars.L, 0);

                    // === stacking intracatena (solo se non opt == 2)
                //     if (opt != 2 && same_chain && abs(i - j) > 1)
                  //      ene += parts[i].vijStack(parts[j], pars.L);
                }

            }

    // ✅ Add bonded FENE energy exactly once
    if (opt != 2)
    {
        if (i > 0 && parts[i].polymer_id == parts[i - 1].polymer_id) {
            ene += parts[i].vijFENE(parts[i - 1], pars.L);
            ene += parts[i].vijLJ(parts[i - 1], pars.L, 0);
            //e += parts[i].vijStack(parts[i - 1], pars.L);
        }
        if (i < parts.size() - 1 && parts[i].polymer_id == parts[i + 1].polymer_id) {
            ene += parts[i].vijFENE(parts[i + 1], pars.L);
            ene += parts[i].vijLJ(parts[i + 1], pars.L, 0);
          //  ene += parts[i].vijStack(parts[i + 1], pars.L);
        }

        ene += bending_energy(i);
    }
    return ene;

}

  ntype totenergy_linked_cells()
  {
    ntype ene = 0.0;
    for(int i = 0; i < parts.size(); i++)
      {
        ene += energy_linked_cells(i, pars.drmax, 1);

      }

    return ene;

  }

  void store_cell_lists(void)
  {
    for(int i = 0; i < pars.Np; i++)
     {
      backup_linked_list[i] = linked_list[i];
     }


   for( int i = 0; i< n_cell; i++)
   {

    head_backup[i] = head[i];

   }

  }

void restore_cell_lists(void)
  {
    for(int i = 0; i < pars.Np; i++)
     {
      linked_list[i] = backup_linked_list[i];

     }


   for( int i = 0; i< n_cell; i++)
   {

    head[i] = head_backup[i];

   }

}

  void set_poly(ntype alpha, ntype  L, int count, std::vector<monomer> &parts, int x)
{
    pvector<ntype,3> dr = {0.0, 0.0, 0.0};

    // First bond: Set poly[1] to a random unit vector added to poly[0]
    parts[count+1].r.random_orient();
    parts[count+1].r /= parts[count + 1].r.norm();
    parts[count+1].r = parts[count].r + L * parts[count +1].r;
    parts[count+1].polymer_id = parts[count].polymer_id;  // Ensure first bond length is L
    parts[count+1].set_epsilon(pars.epsilon);
    parts[count+1].set_sigma(pars.sigma);
    parts[count+1].set_mass(pars.mass);
    parts[count+1].set_rcut(pars.rc * pars.sigma);
    parts[count+1].set_k(pars.k);
    parts[count+1].set_drmax(pars.drmax);
    pbc(count+1);
    parts[count+1].polymer_type = parts[count].polymer_type;
    parts[count+1].set_theta(pars.theta);

    // Initial direction
    pvector<ntype,3> prev_dr = parts[count+1].r - parts[count].r;
    prev_dr /= prev_dr.norm();  // Normalize to unit length

    // Semi-flexible polymer growth
    for(int i = 2; i < x; i++)
    {
        // Generate a new random direction
        pvector<ntype,3> new_dr;
        /*new_dr(0) = gauss() * 0.05;
        new_dr(1) = gauss() * 0.05;
        new_dr(2) = gauss() * 0.05;
        new_dr /= new_dr.norm();  // Normalize to unit step */
        new_dr.random_orient();

        // Blend with previous direction to introduce stiffness
        dr = alpha * prev_dr + (1.0 - alpha) * new_dr;
        dr /= dr.norm();  // Normalize to unit step
        dr *= L;          // Scale to maintain constant bond length

        // Assign new position
        parts[count+i].r = parts[count+i-1].r + dr;
        parts[count+i].polymer_id = parts[count+i-1].polymer_id;
        parts[count+i].set_epsilon(pars.epsilon);
        parts[count+i].set_sigma(pars.sigma);
        parts[count+i].set_mass(pars.mass);
        parts[count+i].set_rcut(pars.rc * pars.sigma);
        parts[count+i].set_k(pars.k);
        parts[count+i].set_drmax(pars.drmax);
        pbc(count+i);
        parts[count+i].set_theta(pars.theta);
        parts[count+i].polymer_type = parts[count].polymer_type;
        // Update previous direction
        prev_dr = dr / L;  // Normalize to avoid drift over multiple steps
    }

}

void rotate_initial_chain(int idx)
{
    /* Rotazione casuale attorno al centro di massa della catena idx */
    ntype theta, sinw, cosw, ox, oy, oz;
    pmatrixq<ntype,3> Omega, OmegaSq, M;
    pvector<ntype,3> o, r0 = {0.0, 0.0, 0.0};

    // Genera angolo e asse casuali
    theta = PI * (rng.ranf() - 0.5);  // da -π/2 a +π/2
    sinw = sin(theta);
    cosw = (1.0 - cos(theta));
    o.random_orient();
    ox = o(0);
    oy = o(1);
    oz = o(2);

    // Calcola centro di massa della catena idx
    for (int j = 0; j < pol[idx].size(); j++)
        r0 += parts[pol[idx][j]].r;
    r0 /= ntype(pol[idx].size());

    // Costruisci matrice Omega (skew-symmetric)
    Omega(0,0) = 0;     Omega(0,1) = -oz; Omega(0,2) = oy;
    Omega(1,0) = oz;    Omega(1,1) = 0;   Omega(1,2) = -ox;
    Omega(2,0) = -oy;   Omega(2,1) = ox;  Omega(2,2) = 0;

    // Omega^2
    OmegaSq(0,0) = -oy*oy - oz*oz; OmegaSq(0,1) = ox*oy;          OmegaSq(0,2) = ox*oz;
    OmegaSq(1,0) = ox*oy;          OmegaSq(1,1) = -ox*ox - oz*oz; OmegaSq(1,2) = oy*oz;
    OmegaSq(2,0) = ox*oz;          OmegaSq(2,1) = oy*oz;          OmegaSq(2,2) = -ox*ox - oy*oy;

    // Matrice di rotazione
    M = -sinw * Omega + cosw * OmegaSq;
    M += Identity;

    // Applica la rotazione e reinserisci con PBC
    for (int j = 0; j < pol[idx].size(); j++)
    {
        parts[pol[idx][j]].r -= r0;
        parts[pol[idx][j]].r = M.transpose() * parts[pol[idx][j]].r;
        parts[pol[idx][j]].r += r0;
        pbc(pol[idx][j]);
    }
}




  void prepare_initial_conf(void)
    {
      int ix, iy, iz;
      int count=0, p_id = 0;
      pars.loadParameters("parametri.txt");
      pars.Np = pars.nx*pars.ny*pars.nz;
      pol.resize(pars.Np);
     // long_polymers =round(pars.Np/((pars.Nx/pars.Ny)+1));
     long_polymers = pars.Np;
      std::cout << "long_polymers: " << long_polymers << "\n";
      short_polymers = pars.Np - long_polymers;
      std::cout << "short_polymers: " << short_polymers << "\n";
      total_polymers = long_polymers + short_polymers;
      pars.Np = long_polymers *pars.Nx + short_polymers*pars.Ny;
      parts.resize(pars.Np);
      pars.L ={ntype(pars.nx), ntype(pars.ny), ntype(pars.nz)};
      pars.L *= 1.5  * cbrt(pars.Nx/10.0) * 10.0 * pars.sigma;
      std::vector<int> polymer_types(total_polymers, pars.Ny); // Tutti inizialmente polimeri corti (3 monomeri)
      std::fill(polymer_types.begin(), polymer_types.begin() + long_polymers, pars.Nx); // Assegniamo 100 lunghi (10 monomeri)


std::random_device rd;
std::mt19937 g(rd());
std::shuffle(polymer_types.begin(), polymer_types.end(), g);
int index = 0;

for(int i = 1; i < pars.nx + 1 ; i++)
{
 for(int j = 1; j < pars.nx +1; j++)
 {
  for(int k = 1; k < pars.nz +1 ; k++)
  {

    int polymer_length = polymer_types[index]; // Lunghezza attuale del polimero

    parts[count].r = {double(i),double(j),double(k)};
    pol[index].resize(polymer_length);
    parts[count].r *= 1.5 * 10.0 * cbrt(pars.Nx/10.0) * pars.sigma;
    for(int l = 0; l < polymer_length; l++)
     {
      pol[index][l] = count + l;
     }
    parts[count].r -= pars.L*0.5;
    p_id++;

    //std::cout << p_id << "\n";
    parts[count].polymer_id = p_id;
    parts[count].polymer_type = polymer_length >= 10 ? 0 : 1;
  //  std::cout << ligand[count].polymer_id << " " << "\n";

   parts[count].set_epsilon(pars.epsilon);
   parts[count].set_sigma(pars.sigma);
   parts[count].set_mass(pars.mass);
   parts[count].set_rcut(pars.rc * pars.sigma);
   parts[count].set_k(pars.k);
   parts[count].set_drmax(pars.drmax);
   parts[count].set_theta(pars.theta);


    set_poly(0.3, 1.0, count, parts, polymer_length);
    count += polymer_length;

    index++;


  }
 }
}

  }

  void run(void)
    {
      // intentionally void
    };
};

template<typename particle_type>
class mcsim: public sim<particle_type>
{
  // for calc_acceptance_and_adjust: total trial moves and accepted ones
  // for calculating acceptance rates.
  using bc=sim<particle_type>;
  using bc::parts, bc::pars, bc::pbc, bc::pol, bc::long_polymers, bc::short_polymers, bc::total_polymers,/*bc::head,bc::linked_lists,bc::head_backup,bc::backup_linked_lists,*/
        bc::save_mgl_snapshot, bc::Identity,
        bc::particle_to_cells, bc::energy_linked_cells, bc::totenergy_linked_cells,bc::init_cells,bc::update_cells,bc::store_cell_lists,bc::restore_cell_lists,
        bc::totenergy, bc::calcenergyi;

  // counters used to calculate acceptance rates
  long int tot_tra, tra_rej, tra_rej_cls, tot_tra_cls, tot_tra_rot_cls,tra_rej_rot_cls;



  void rot_cls_move(int i)
    {

       ntype ene = 0.0;
       store_cell_lists();


       for(int j = 0; j < pol[i].size(); j++)
       {
        ene += energy_linked_cells(pol[i][j], pars.rc, 2);
        parts[pol[i][j]].store();
       }

       alpha__rot_cls(i);
     //  std::cout << "cls_move" << "\n";
       acc_rot_cls(i,ene);

  }


  void alpha__rot_cls(int i)
  {
      /* random rotation around random orientation o */
      ntype theta, sinw, cosw, ox, oy, oz;
      pmatrixq<ntype,3> Omega, OmegaSq, Ro, M;
      pvector<ntype,3> o, r0;


      theta = pars.deltra_rot * 2.0 * (rng.ranf() - 0.5);
      sinw = sin(theta);
      cosw = (1.0 - cos(theta));

      pvector<ntype,3> Dr;
      ntype dist = 0.0;
      int count = 0;


          std::vector<pvector<ntype,3>> unwrapped;
          unwrapped.resize(pol[i].size());

          // Unwrap: primo monomero come origine
          unwrapped[0] = parts[pol[i][0]].r;

          for (int j = 1; j < pol[i].size(); j++)
          {
              pvector<ntype,3> d = parts[pol[i][j]].r - parts[pol[i][j - 1]].r;
              d -= pars.L.mulcw(rint(d.divcw(pars.L))); // MIC tra monomeri consecutivi
              unwrapped[j] = unwrapped[j - 1] + d;
          }

          // Calcola il centro di massa della catena unwrapped
          for (int j = 0; j < pol[i].size(); j++)
              r0 += unwrapped[j];
          r0 /= ntype(pol[i].size());


      o.random_orient();
      ox = o(0);
      oy = o(1);
      oz = o(2);

      // Skew-symmetric matrix (Omega)
      Omega(0,0) = 0;
      Omega(0,1) = -oz;
      Omega(0,2) = oy;
      Omega(1,0) = oz;
      Omega(1,1) = 0;
      Omega(1,2) = -ox;
      Omega(2,0) = -oy;
      Omega(2,1) = ox;
      Omega(2,2) = 0;

      // Omega squared
      OmegaSq(0,0) = -oy*oy - oz*oz;
      OmegaSq(0,1) = ox*oy;
      OmegaSq(0,2) = ox*oz;
      OmegaSq(1,0) = ox*oy;
      OmegaSq(1,1) = -ox*ox - oz*oz;
      OmegaSq(1,2) = oy*oz;
      OmegaSq(2,0) = ox*oz;
      OmegaSq(2,1) = oy*oz;
      OmegaSq(2,2) = -ox*ox - oy*oy;

      // Rotation matrix
      M = -sinw * Omega + cosw * OmegaSq;
      M += Identity;
      pmatrixq<ntype,3> Rt = M.transpose();

      // Rotate all particles in the cluster
      for(int j = 0; j < pol[i].size(); j++)
      {  // std::cout << parts[pol[i][j]].r -r0 << "\n";
        pvector<ntype, 3> r_rel = unwrapped[j] - r0;
        parts[pol[i][j]].r = M.transpose() * r_rel + r0;
       pbc(pol[i][j]);
       update_cells(pol[i][j]);
      }

  }



void acc_rot_cls(int i, ntype eno)
{
 ntype ene_new = 0.0;

 for(int j = 0; j < pol[i].size(); j++)
  {
   ene_new += energy_linked_cells(pol[i][j],pars.rc,2);
  }
 ntype deltaU = ene_new-eno;
 ntype xi=rng.ranf();

 if ((deltaU>0.0 && xi >= exp(-deltaU/pars.T)) || std::isnan(deltaU))
   {
    for(int j = 0; j < pol[i].size(); j++)
      {
       parts[pol[i][j]].restore();
      }
    tra_rej_rot_cls++;


    restore_cell_lists();
   }


}

  void alpha(int i)
    {


      pvector<ntype,3> delr;
      delr ={pars.deltra*2.0*(rng.ranf()-0.5),  pars.deltra*2.0*(rng.ranf()-0.5),
           pars.deltra*2.0*(rng.ranf()-0.5)};
      parts[i].tra_move(delr);
      pbc(i);

    }



  void acc(int i, ntype eno)
    {
   
     ntype ene_new = energy_linked_cells(i,pars.drmax,0);
    // std::cout << "ene_new: " << ene_new << "\n";
     ntype deltaU = ene_new-eno;
     ntype xi=rng.ranf();
     if ((deltaU>0.0 && xi >= exp(-deltaU/pars.T)) || std::isnan(deltaU))
       {

        parts[i].restore();

        tra_rej++;
       }

     else
       {
        update_cells(i);
       }

    }
    void move_NTV(int i)
    {
     ntype ene;
    // std::cout << "move" << "\n";
     ene = energy_linked_cells(i,pars.drmax,0);
  //   std::cout << "ene: " << ene << "\n";
     parts[i].store();
     alpha(i);
     acc(i,ene);

    }


   void cls_move(int i)
    {
     ntype ene = 0.0, ene_new = 0.0;
     pvector<ntype,3> delr;
     delr = {pars.deltra_cls*2.0*(rng.ranf()-0.5),  pars.deltra_cls*2.0*(rng.ranf()-0.5),
           pars.deltra_cls*2.0*(rng.ranf()-0.5)};

    store_cell_lists();


     for(int j = 0; j < pol[i].size(); j++)
     {
      ene += energy_linked_cells(pol[i][j], pars.rc, 2);
      parts[pol[i][j]].store();
      alpha_cls(pol[i][j], delr);
      ene_new += energy_linked_cells(pol[i][j], pars.rc, 2);
     }

     acc_cls(i,ene,ene_new);
    }


    void alpha_cls(int i, pvector<ntype,3> delr = {0.0, 0.0, 0.0})
    {

       parts[i].tra_move(delr);
       pbc(i);
      update_cells(i);

    }

  void acc_cls(int i, ntype eno, ntype ene_new = 0.0)
    {

     ntype deltaU = ene_new-eno;
     ntype xi=rng.ranf();

     if ((deltaU>0.0 && xi >= exp(-deltaU/pars.T)) || std::isnan(deltaU))
       {
        for(int j = 0; j < pol[i].size(); j++)
          {
           parts[pol[i][j]].restore();

          }

        restore_cell_lists();

        tra_rej_cls++;
       }

    }


  void calc_acceptance_and_adjust(void)
    {

      double acceptance, acceptance_cls;
      if(tot_tra> 0)
       {
       acceptance = (((double)(tot_tra-tra_rej))/tot_tra);
          if (acceptance > 0.5)
            {
              pars.deltra *= 1.1;
            }
          else if (acceptance < 0.35)
            {
              pars.deltra /= 1.1;
            }

      tot_tra=tra_rej=0;
       }

     if(tot_tra_cls > 0 && pars.deltra_cls <= pars.L(0)/2)
       {
       acceptance_cls = (((double)(tot_tra_cls-tra_rej_cls))/tot_tra_cls);
          if (acceptance_cls > 0.5)
            {
              pars.deltra_cls *= 1.5;
            }
          else if (acceptance_cls < 0.35)
            {
              pars.deltra_cls /= 1.5;
            }

      tot_tra_cls=tra_rej_cls=0;
       }

      if(tot_tra_rot_cls> 0 pars.deltra_rot <= PI)
        {
        double acceptance_rot_cls = (((double)(tot_tra_rot_cls-tra_rej_rot_cls))/tot_tra_rot_cls);
           if (acceptance_rot_cls > 0.5)
             {
               pars.deltra_rot *= 1.3;
             }
           else if (acceptance_rot_cls < 0.35)
             {
               pars.deltra_rot /= 1.3;
             }

       tot_tra_rot_cls=tra_rej_rot_cls=0;
        }

}



 public:
  void run(int c)
    {
      // ciclo sui passi MC
      int i, t, ip;
      tot_tra = tra_rej = tra_rej_cls =tot_tra_cls = tot_tra_rot_cls = tra_rej_rot_cls = 0;
      init_measures();
      init_cells();
      particle_to_cells();

      for (t = 0; t < pars.totsteps; t++)
        {

         for( int i = 0; i < pars.Np; i++)
            {
              if(rng.ranf() < 0.9)
              {
               ip = rng.ranf()*(pars.Np);
               move_NTV(ip);
               tot_tra++;
              }

             else
              {
               ip = rng.ranf()*(pars.nx * pars.ny * pars.nz);
              if(rng.ranf() < 0.5)
                {

                 cls_move(ip);
                 tot_tra_cls++;
                }
              else
               {
                rot_cls_move(ip);
                tot_tra_rot_cls++;
             }
           }
          }


         if ( t>= 0 && t%100 == 0)
            {
             save_mgl_snapshot(t, c);
            }



          if (t > 0 && pars.adjstps > 0 && t % pars.adjstps == 0 &&
              t < pars.maxadjstps)
            {
              calc_acceptance_and_adjust();
            }
        }

     //std::cout << "end to end distance: " << Ree() << "\n";
    }
};

#endif