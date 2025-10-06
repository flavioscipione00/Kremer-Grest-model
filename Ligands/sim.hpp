#ifndef _SIMCLASS_
#define _SIMCLASS_
#ifdef CLUSTER_MODE
  #define COUT(x) do {} while (0)
#else
  #define COUT(x) std::cout << x
#endif

#include <vector>
#include <string>
#include <fstream>
#include "./params.hpp"
#include "./particle.hpp"
#include "./pmatrix.hpp"
#include "./randnumgen.hpp"
#include <iomanip>
#include <unordered_set>


using ntype=double;

template<typename particle_type>
class sim
{
  using simp = simpars;
protected:
  simp pars;  // parametri per composizione
  std::vector<monomer> parts;
  std::vector<std::vector<int>> pol;
  int cell_nx , cell_ny, cell_nz , n_cell ;
  pvector <ntype,3> cell_size;
  int long_polymers, short_polymers, total_polymers;
  pmatrixq<ntype,3> Identity = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  std::vector<int> head;
  std::vector<int> linked_list;
  std::vector<int> backup_linked_list;
  std::vector<int> head_backup;

  std::unordered_map<int, std::unordered_set<int>>  g4_bond_map;
  std::unordered_map<int, std::unordered_set<int>>  g4_map_old;

 //if opt=1 calculate energies only for i < j
 /* ntype calcenergyi(int i, int opt=0)
    {
      int j;
      ntype enei=0.0;
      for (j=0; j < parts.size(); j++) // pars.Np è il numero totale di particelle
        {
          if (opt==1 && i >= j)
            continue;
          if (i==j)
            continue;
          // la classe particelle deve implementare un metodo vij per il calcolo dell'energia d'interazione
          enei += parts[i].vijLJ(parts[j], pars.L);
          if (parts[i].polymer_id == parts[j].polymer_id && abs(i-j) == 1) enei += parts[i].vijFENE(parts[j], pars.L);

          // pars.L è un vettore con i lati del box
        }

      enei += bending_energy(i);

      return enei;
    } */


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






 /* ntype totenergy()
    {
      ntype ene=0.0;
      for (auto i=0; i < parts.size(); i++)
        {
          ene+=calcenergyi(i, 1);
        }

      return ene;
    } */

  void pbc(int i)
    {

      parts[i].r -= pars.L.mulcw(rint(parts[i].r.divcw(pars.L)));

    }

      void save_mgl_snapshot(long int t, int c)
    {
       std::fstream f;
       std::string s;

       s = "cnf-" + std::to_string(t) + "_" + std::to_string(c) + ".txt";
       f.open(s, std::ios::out|std::ios::trunc);
       for (int i=0; i < pars.Np; i++)
         {
           f << parts[i].r(0) << " " << parts[i].r(1) << " " <<
             parts[i].r(2) << " " << parts[i].polymer_id << " " << parts[i].polymer_type <<" " << parts[i].bound_to << "\n";
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
     COUT("cell_nx: " << cell_nx << ", cell_ny: " << cell_ny << ", cell_nz: " << cell_nz << "\n");
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

void update_neighboring_ligands(int g4_index)
{
    if (g4_index < 0 || g4_index >= parts.size()) return;

    const double r_min = pars.r_sq;
    const double r_cutoff = r_min + parts[g4_index].delta;


    int x = floor((parts[g4_index].r(0) + 0.5 * pars.L(0)) / cell_size(0));
    int y = floor((parts[g4_index].r(1) + 0.5 * pars.L(1)) / cell_size(1));
    int z = floor((parts[g4_index].r(2) + 0.5 * pars.L(2)) / cell_size(2));

    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++)
    {
        int cx = (x + dx + cell_nx) % cell_nx;
        int cy = (y + dy + cell_ny) % cell_ny;
        int cz = (z + dz + cell_nz) % cell_nz;

        int index = cx + cy * cell_nx + cz * cell_ny * cell_nz;
        if (index < 0 || index >= head.size()) continue;

        for (int i = head[index]; i >= 0; i = linked_list[i])
        {
            if (i < 0 || i >= parts.size()) continue;
            if (parts[i].polymer_type != 1) continue; // solo ligandi


            pvector<ntype, 3> Dr = parts[i].r - parts[g4_index].r;
            Dr -= pars.L.mulcw(rint(Dr.divcw(pars.L)));
            double dist = Dr.norm();

            if (dist < r_cutoff && dist > r_min)
            {
                update_binding_state(i);
            }
        }
    }
}


void update_binding_state(int i)
{
    // ⛔ Ignora i non-ligandi
    if (parts[i].polymer_type != 1) return;

    const double r_min = pars.r_sq;
    const double r_cutoff = r_min + parts[i].delta;

    // === UNBINDING ===
    int j_old = parts[i].bound_to;
    if (j_old != -1)
    {
        pvector<ntype, 3> Dr = parts[i].r - parts[j_old].r;
        Dr -= pars.L.mulcw(rint(Dr.divcw(pars.L)));
        double dist = Dr.norm();

        if (dist > r_cutoff)
        {
            // Rimuovi il legame
            parts[i].bound_to = -1;

            // Rimuovi i dalla lista di j_old
           g4_bond_map[j_old].erase(i);
           if (g4_bond_map[j_old].empty()) g4_bond_map.erase(j_old);
        }
        else
        {
            return;  // rimane legato, non cerchiamo altri G4
        }
    }

    // === BINDING ===
    int j_best = -1;
    double r_best = r_cutoff + 1.0;

    int x = floor((parts[i].r(0) + 0.5 * pars.L(0)) / cell_size(0));
    int y = floor((parts[i].r(1) + 0.5 * pars.L(1)) / cell_size(1));
    int z = floor((parts[i].r(2) + 0.5 * pars.L(2)) / cell_size(2));

    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++)
    {
        int cx = (x + dx + cell_nx) % cell_nx;
        int cy = (y + dy + cell_ny) % cell_ny;
        int cz = (z + dz + cell_nz) % cell_nz;

        int index = cx + cy * cell_nx + cz * cell_ny * cell_nx;
        if (index < 0 || index >= head.size()) continue;

        for (int j = head[index]; j >= 0; j = linked_list[j])
        {
            if (j == i) continue;
            if (parts[j].polymer_type != 0) continue; // solo G4

            pvector<ntype, 3> Dr = parts[i].r - parts[j].r;
            Dr -= pars.L.mulcw(rint(Dr.divcw(pars.L)));
            double dist = Dr.norm();

            if (dist < r_best && dist < r_cutoff && dist > r_min)
            {
                j_best = j;
                r_best = dist;
            }
        }
    }

    if (j_best != -1)
    {
        //parts[j_best].bond_store();  // opzionale
        parts[i].bound_to = j_best;
        g4_bond_map[j_best].insert(i);

    }
}



ntype energy_linked_cells(int i, int opt = 0, std::unordered_set<int>* cluster_set = nullptr)
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
                    bool same_type = (parts[i].polymer_type == parts[j].polymer_type);


                    // Per opt == 2 (cluster move): ignora interazioni intra-catena
                    if (opt == 2 && same_chain) continue;


                    if(!same_type && !same_chain) ene += parts[i].vij_square(parts[j], pars.L,i,j);


                    // === WCA ===
                    if (!same_chain || (abs(i - j) > 1 && opt != 2))
                        ene += parts[i].vijLJ(parts[j], pars.L);

                }

            }

    // ✅ Add bonded FENE energy exactly once
    if (opt != 2)
    {
       if(parts[i].polymer_type == 1) // Solo per i monomeri di tipo 1
        {
            if (i > 0 && parts[i].polymer_id == parts[i - 1].polymer_id) {
                ene += parts[i].vijFENE_sym(parts[i - 1], pars.L);
                ene += parts[i].vijLJ(parts[i - 1], pars.L);
            }
            if (i < parts.size() - 1 && parts[i].polymer_id == parts[i + 1].polymer_id) {
                ene += parts[i].vijFENE_sym(parts[i + 1], pars.L);
                ene += parts[i].vijLJ(parts[i + 1], pars.L);
            }
        }


      else if(parts[i].polymer_type == 0) // Solo per i monomeri di tipo 0
        {
         if (i > 0 && parts[i].polymer_id == parts[i - 1].polymer_id) {
            ene += parts[i].vijFENE_zero(parts[i - 1], pars.L);
            ene += parts[i].vijLJ(parts[i - 1], pars.L);
        }
        if (i < parts.size() - 1 && parts[i].polymer_id == parts[i + 1].polymer_id) {
            ene += parts[i].vijFENE_zero(parts[i + 1], pars.L);
            ene += parts[i].vijLJ(parts[i + 1], pars.L);
        }

         ene += bending_energy(i);
        }




    }
    return ene;

}

  ntype totenergy_linked_cells()
  {
    ntype ene = 0.0;
    for(int i = 0; i < parts.size(); i++)
      {
        ene += energy_linked_cells(i, 1);

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



void set_poly(ntype alpha, ntype L, int count, std::vector<monomer> &parts, int x)
{
    pvector<ntype,3> dr = {0.0, 0.0, 0.0};
    pvector<ntype,3> prev_dr = {0.0, 0.0, 1.0};
   // ntype sg = rng.ranf();
   // sg >= 0.5? prev_dr *= -1.0 : prev_dr *= 1.0;

    // Primo bond verticale fisso
    parts[count+1].r = parts[count].r + prev_dr * L;
    parts[count+1].polymer_id = parts[count].polymer_id;
    pbc(count+1);
    if(parts[count].polymer_type == 1)
        parts[count+1].polymer_type = parts[count].polymer_type;
    else
        parts[count+1].polymer_type = parts[count].polymer_type;

    // Costruzione stacking rigido
    for(int i = 2; i < x; i++)
    {
        dr = L * prev_dr;
        parts[count+i].r = parts[count+i-1].r + dr;
        parts[count+i].polymer_id = parts[count].polymer_id;

        parts[count+i].polymer_type = parts[count].polymer_type;
        pbc(count+i);
        prev_dr = dr / L;
    }
}


void assign_properties_to_all(std::vector<monomer>& parts) {
  for (int i = 0; i < parts.size(); ++i) {
      int ptype = parts[i].polymer_type;
      parts[i].set_epsilon(pars.epsilon);
      if(ptype == 0)
          parts[i].set_sigma(pars.sigma_1);
      else if(ptype == 1)
          parts[i].set_sigma(pars.sigma_2);
      else
          parts[i].set_sigma(0.0); // Per altri tipi di polimero, se necessario

      parts[i].set_mass(pars.mass);
      parts[i].set_rcut((ptype == 0) ? pars.rc_1 : pars.rc_2);
      parts[i].set_k(pars.k);
      parts[i].set_drmax(pars.drmax);
      parts[i].set_theta((ptype == 0) ? pars.theta_1 : pars.theta_2);
      parts[i].set_r_sq(pars.r_sq);
      parts[i].set_u(pars.u);
      parts[i].set_delta(pars.delta);
      parts[i].set_bound_to(-1);
      parts[i].bound_to_old = -1;


  }
}




void prepare_initial_conf_stacked()
{
    int count = 0, p_id = 0;
    pars.loadParameters("parametri.txt");

    int mult = round((pars.Nx / pars.Ny));
    COUT("mult: " << mult << "\n");

    int total_sites = pars.nx * pars.ny * pars.nz;
    int long_sites = total_sites / 2;
    int short_blocks = total_sites - long_sites;

    long_polymers = long_sites;
    short_polymers = short_blocks * mult;
    total_polymers = long_polymers + short_polymers;

    pars.Np = long_polymers * pars.Nx + short_polymers * pars.Ny;
    parts.resize(pars.Np);
    pol.resize(total_polymers);


    // === RESCALING AUTOMATICO DELLO SPACING ===

    ntype spacing_0 = 1.15;            // spacing di riferimento
    ntype volume_scaling = 100.0;        // QUI puoi mettere il tuo fattore di volume
    int N_site = total_sites;

    ntype spacing = spacing_0 * pow(volume_scaling, 1.0/3.0);
    COUT("spacing: " << spacing << "\n");

    ntype lb = pars.lb; // lunghezza di bond usata per il set_poly
    ntype lb_2 = pars.lb_2;

    // === Calcolo dimensione box ===
    pars.L = {1.0, 1.0, 1.0};
    pars.L *= spacing * pars.Nx ;  // Questo resta come nel tuo codice

    std::vector<int> polymer_types(total_sites, pars.Nx);
    std::fill(polymer_types.begin() + long_sites, polymer_types.end(), -pars.Ny);

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(polymer_types.begin(), polymer_types.end(), g);

    int index = 0;

    for (int i = 0; i < pars.nx; i++) {
        for (int j = 0; j < pars.ny; j++) {
            for (int k = 0; k < pars.nz; k++) {
                pvector<ntype, 3> pos = {ntype(i), ntype(j), ntype(k)};
                pos(0) *= spacing;
                pos(1) *= spacing;
                pos(2) *= spacing;


                pos -= pars.L * 0.5;
                pos(2) += rng.ranf() * pars.L(0) * 0.5;

                int type = polymer_types[index];

                if (type > 0) {
                    // Polimero lungo
                    int polymer_length = type;
                    pol[p_id].resize(polymer_length);
                    parts[count].r = pos;
                    parts[count].polymer_id = p_id;
                    parts[count].polymer_type = 0;
                    pbc(count);
                    set_poly(1.0, lb, count, parts, polymer_length);
                    for (int l = 0; l < polymer_length; l++) {
                        pol[p_id][l] = count + l;
                    }
                    count += polymer_length;
                    p_id++;
                }
                else {
                    // m polimeri corti distinti
                    int polymer_length = -type;
                    int a = 0;
                    ntype sg = rng.ranf();
                    //sg >= 0.5 ? a = 1 : a = -1;

                    for (int m = 0; m < mult; m++) {
                        int current_pid = p_id + m;
                        pol[current_pid].resize(polymer_length);

                        pvector<ntype, 3> offset = {0.0, 0.0,  pars.Ny * +m * lb};
                        parts[count].r = pos + offset;
                        parts[count].polymer_id = current_pid;
                        parts[count].polymer_type = 1;
                        pbc(count);

                        set_poly(1.0, lb_2, count, parts, polymer_length);

                        for (int l = 0; l < polymer_length; l++) {
                            pol[current_pid][l] = count + l;
                        }
                        count += polymer_length;
                    }
                    p_id += mult;
                }
                index++;
            }
        }
    }

    assign_properties_to_all(parts);
}

void store(int i)
{
    parts[i].store();

    if (parts[i].polymer_type == 1)  // Ligando
    {
        parts[i].bond_store();

        int g4 = parts[i].bound_to;
        if (g4 != -1 && g4_bond_map.count(g4))
            g4_map_old[g4] = g4_bond_map[g4];  // Copia l'intero set
    }
    else  // G4
    {
        if (g4_bond_map.count(i))
        {
            g4_map_old[i] = g4_bond_map[i];  // Copia l'intero set

            for (int lig : g4_bond_map[i])
                parts[lig].bond_store();
        }
    }
}



void restore(int i)
{
    parts[i].restore();

    if (parts[i].polymer_type == 1)  // Ligando
    {
        int g4_new = parts[i].bound_to;
        int g4_old = parts[i].bound_to_old;

        parts[i].bond_restore();  // ripristina bound_to = g4_old

        // Rimuovi il ligando dalla mappa del G4 a cui si era legato (g4_new)
        if (g4_new != -1 && g4_new != g4_old && g4_bond_map.count(g4_new))
        {
            g4_bond_map[g4_new].erase(i);

            if (g4_bond_map[g4_new].empty())
                g4_bond_map.erase(g4_new);
        }

        // Ripristina la lista vecchia se esiste
        if (g4_old != -1 && g4_map_old.count(g4_old))
        {
            g4_bond_map[g4_old] = g4_map_old[g4_old];
        }
    }

    else  // G4
    {
        if (g4_map_old.count(i))
        {
            g4_bond_map[i] = g4_map_old[i];

            for (int lig : g4_map_old[i])
                parts[lig].bond_restore();
        }
    }
}



  void run(void)
    {
      // intentionally void
    }
};

template<typename particle_type>
class mcsim: public sim<particle_type>
{
  // for calc_acceptance_and_adjust: total trial moves and accepted ones
  // for calculating acceptance rates.
  using bc=sim<particle_type>;
  using bc::parts, bc::pars, bc::pbc, bc::pol, bc::long_polymers, bc::short_polymers, bc::total_polymers,
        bc::save_mgl_snapshot, bc::Identity,
        bc::particle_to_cells, bc::energy_linked_cells, bc::totenergy_linked_cells,bc::init_cells,bc::update_cells,bc::store_cell_lists,bc::restore_cell_lists, bc::cell_size,bc::linked_list,bc::head,bc::cell_ny,bc::cell_nx,bc::cell_nz,
         bc::update_binding_state, bc::update_neighboring_ligands, bc::store, bc::restore,bc::g4_bond_map,bc::g4_map_old;

  // counters used to calculate acceptance rates
  long int tot_tra, tra_rej,tot_tra_2, tra_rej_2, tra_rej_cls, tot_tra_cls, tot_tra_rot_cls,tra_rej_rot_cls;


 void rot_cls_move(int i)
{
    ntype ene = 0.0;
    ntype ene_new = 0.0;

    store_cell_lists();  // per sicurezza totale sulle celle
    bool is_ligand_cluster = (parts[pol[i][0]].polymer_type == 1);

    // Backup energia, posizioni, binding, celle
    for(int j = 0; j < pol[i].size(); j++)
    {
        ene += energy_linked_cells(pol[i][j], 2);
        store(pol[i][j]);
    }

    // Trial move (rotazione del cluster)
    alpha_rot_cls(i);

    // Aggiorna celle per TUTTI i monomeri
for (int j = 0; j < pol[i].size(); j++)
{
    pbc(pol[i][j]);
    update_cells(pol[i][j]);

    if (is_ligand_cluster)
    {
        update_binding_state(pol[i][j]);  // ligando → aggiorna sé stesso
    }
    else
    {
auto it = g4_bond_map.find(pol[i][j]);
if (it != g4_bond_map.end()) {
    std::vector<int> ligands(it->second.begin(), it->second.end());
    for (int k : ligands) {
        if (k < 0 || k >= parts.size()) {
            std::cerr << "❌ cls_move: indice k invalido = " << k
                      << " da G4 = " << pol[i][j]
                      << " (parts.size() = " << parts.size() << ")\n";
            continue;  // evita crash
        }
        update_binding_state(k);
    }
}
        update_neighboring_ligands(pol[i][j]);  // Aggiorna i ligandi vicini
    }

    ene_new += energy_linked_cells(pol[i][j], 2);
}


    acc_rot_cls(i, ene, ene_new);
}



  void alpha_rot_cls(int i)
  {
      /* random rotation around random orientation o */
      ntype theta, sinw, cosw, ox, oy, oz;
      pmatrixq<ntype,3> Omega, OmegaSq, Ro, M;
      pvector<ntype,3> o, r0;


      //theta = PI/2.0;
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
        parts[pol[i][j]].r = Rt * r_rel + r0;
      }

  }



void acc_rot_cls(int i, ntype eno, ntype ene_new)
{

 ntype deltaU = ene_new-eno;
 ntype xi=rng.ranf();

 if ((deltaU>0.0 && xi >= exp(-deltaU/pars.T)) || std::isnan(deltaU))
   {

    restore_cell_lists();
    for(int j = 0; j < pol[i].size(); j++)
      {
       restore(pol[i][j]);

      if(parts[pol[i][0]].polymer_type == 0)
        {
          update_neighboring_ligands(pol[i][j]);  // Aggiorna i ligandi vicini
        }
      }
    tra_rej_rot_cls++;



   }


}


  // //sto implementando solo la trial move (o-->n) per NTV
  void alpha(int i)
    {

   //   std::cout << "alpha" << "\n";
      pvector<ntype,3> delr;
      delr ={pars.deltra*2.0*(rng.ranf()-0.5),  pars.deltra*2.0*(rng.ranf()-0.5),
           pars.deltra*2.0*(rng.ranf()-0.5)};
      parts[i].tra_move(delr);
      pbc(i);

if (parts[i].polymer_type == 1)
{
    update_binding_state(i);  // ligando → aggiorna sé stesso
}
else
{
auto it = g4_bond_map.find(i);
if (it != g4_bond_map.end()) {
    std::vector<int> vec_copy(it->second.begin(), it->second.end());
    for (int k : vec_copy) {
        if (k < 0 || k >= parts.size()) {
            std::cerr << "❌ ERRORE: indice k fuori range: " << k << "\n";
            continue;
        }
        update_binding_state(k);
    }
}


    update_neighboring_ligands(i);  // Aggiorna i ligandi vicini
   }

    }



  void acc(int i, ntype eno)
    {
     //std::cout << "acc" << "\n";
     ntype ene_new = energy_linked_cells(i,0);
    // std::cout << "ene_new: " << ene_new << "\n";
     ntype deltaU = ene_new-eno;
     ntype xi=rng.ranf();
     if ((deltaU>0.0 && xi >= exp(-deltaU/pars.T)) || std::isnan(deltaU))
       {

        restore(i);
        if(parts[i].polymer_type == 0)
          {
            update_neighboring_ligands(i);  // Aggiorna i ligandi vicini
          }

       if(parts[i].polymer_type == 0)   
        {
          tra_rej++;
        }

       else
         {
          tra_rej_2++;
        }

       }

     else
       {
        update_cells(i);
       }
    }
    void move_NTV(int i)
    {
     ntype ene = 0.0;

     ene = energy_linked_cells(i,0);

     store_cell_lists();
     store(i);


     alpha(i);
     acc(i,ene);

    }


void cls_move(int i)
{
    ntype ene = 0.0;
    ntype ene_new = 0.0;
    bool is_ligand_cluster = (parts[pol[i][0]].polymer_type == 1);

    store_cell_lists();  // per sicurezza totale sulle celle complete

    // Backup energia, posizioni, binding, celle
    for(int j = 0; j < pol[i].size(); j++)
    {
        ene += energy_linked_cells(pol[i][j], 2);
        store(pol[i][j]);
    }



    // Trial move (traslazione)
    alpha_cls(i);

    // Aggiorna celle per TUTTI i monomeri
    for(int j = 0; j < pol[i].size(); j++)
    {
        pbc(pol[i][j]);
        update_cells(pol[i][j]);

   if (is_ligand_cluster)
    {
        update_binding_state(pol[i][j]);  // ligando → aggiorna sé stesso
    }
    else
    {
auto it = g4_bond_map.find(pol[i][j]);
if (it != g4_bond_map.end()) {
    std::vector<int> ligands(it->second.begin(), it->second.end());
    for (int k : ligands) {
        if (k < 0 || k >= parts.size()) {
            std::cerr << "❌ cls_move: indice k invalido = " << k
                      << " da G4 = " << pol[i][j]
                      << " (parts.size() = " << parts.size() << ")\n";
            continue;  // evita crash
        }
        update_binding_state(k);
    }
}
              update_neighboring_ligands(pol[i][j]);  // Aggiorna i ligandi vicini
    }




        ene_new += energy_linked_cells(pol[i][j], 2);
    }

    acc_cls(i,ene, ene_new);


}


    void alpha_cls(int i)
    {
     pvector<ntype,3> delr;
     delr = {pars.deltra_cls*2.0*(rng.ranf()-0.5),  pars.deltra_cls*2.0*(rng.ranf()-0.5),
           pars.deltra_cls*2.0*(rng.ranf()-0.5)};

     for(int j = 0; j < pol[i].size(); j++)
      {
       parts[pol[i][j]].tra_move(delr);

      }
    }

  void acc_cls(int i, ntype eno, ntype ene_new)
    {
     ntype deltaU = ene_new-eno;
     ntype xi=rng.ranf();

     if ((deltaU>0.0 && xi >= exp(-deltaU/pars.T)) || std::isnan(deltaU))
       {
        for(int j = 0; j < pol[i].size(); j++)
          {

           restore_cell_lists();
           restore(pol[i][j]);
           if(parts[pol[i][0]].polymer_type == 0)
            {
              update_neighboring_ligands(pol[i][j]);  // Aggiorna i ligandi vicini
            }
          }



        tra_rej_cls++;

          }


       }




  void calc_acceptance_and_adjust(void)
    {

      double acceptance, acceptance_2, acceptance_cls,acceptance_or;

      if(tot_tra> 0)
       {
       acceptance = (((double)(tot_tra-tra_rej))/tot_tra);
       COUT("rate tra: " << acceptance << " deltra=" << pars.deltra << "\n");
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

    if(tot_tra_2> 0)
       {
       acceptance_2 = (((double)(tot_tra_2-tra_rej_2))/tot_tra_2);
       COUT("rate tra: " << acceptance_2 << " deltra_2=" << pars.deltra_2 << "\n");
          if (acceptance_2 > 0.5)
            {
              pars.deltra_2 *= 1.1;
            }
          else if (acceptance < 0.35)
            {
              pars.deltra_2 /= 1.1;
            }

      tot_tra_2=tra_rej_2=0;
       }


     if(tot_tra_cls > 0 && (pars.deltra_cls < pars.L(0)/2.0))
       {
       acceptance_cls = (((double)(tot_tra_cls-tra_rej_cls))/tot_tra_cls);
       COUT("rate cls: " << acceptance_cls << " deltra_cls=" << pars.deltra_cls << "\n");
          if (acceptance_cls > 0.5)
            {
              pars.deltra_cls *=1.5;
            }
          else if (acceptance_cls < 0.30)
            {
              pars.deltra_cls /= 1.5;
            }

      tot_tra_cls=tra_rej_cls=0;
       }

    if( tot_tra_rot_cls > 0 && (pars.deltra_rot < PI))
       {
       acceptance_or = (((double)(tot_tra_rot_cls-tra_rej_rot_cls))/tot_tra_rot_cls);
       COUT("rate rot: " << acceptance_or << " deltra_rot_cls=" << pars.deltra_rot << "\n");
          if (acceptance_or > 0.5)
            {
              pars.deltra_rot *=1.1;
            }
          else if (acceptance_or < 0.30)
            {
              pars.deltra_rot /= 1.1;
            }

}
    }

  void init_measures(void)
    {
     std::ofstream file("output.txt", std::ios::out | std::ios::trunc);
     file.close();
    }



 public:
  void run(int c)
    {
      // ciclo sui passi MC
      int i, t, ip;
      tot_tra = tra_rej = tot_tra_2 = tra_rej_2 = tra_rej_cls =tot_tra_cls = tot_tra_rot_cls = tra_rej_rot_cls = 0;
      init_measures();
      init_cells();
      particle_to_cells();

      COUT(pars.L(0) << " " << pars.L(1) << " " << pars.L(2) << "\n");


for (t = 0; t < pars.totsteps; t++)
{

for (int i = 0; i < pars.Np; i++)
{
        double chi = rng.ranf();
        if (chi < 0.90)
        {
            ip = rng.ranf() * pars.Np;
            move_NTV(ip);
            if(parts[i].polymer_type ==0 )
             {
            tot_tra++;
             }
            else
             {
            tot_tra_2++;
             }
        }
        else if (chi < 0.95)
        {
            ip = rng.ranf() * total_polymers;
            cls_move(ip);
            tot_tra_cls++;
        }
        else
        {
            ip = rng.ranf() * total_polymers;
            rot_cls_move(ip);
            tot_tra_rot_cls++;
        }
    }

    if (t > 0 && pars.adjstps > 0 && t % pars.adjstps == 0 &&
              t < pars.maxadjstps)
            {
            // std::cout << "Adjusting...\n";
              calc_acceptance_and_adjust();
            }



         if ( t>=10000 && t%pars.savemeasure == 0)
         {
          //std::cout << "t: " << t << "\n";
           save_mgl_snapshot(t,c);

        }


        }
    }
};
#endif