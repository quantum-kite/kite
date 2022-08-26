/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

#include "Generic.hpp"
#include <Eigen/Dense>
#include <random>
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
template<typename T, unsigned D>
class Simulation;
#include "Global.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"
//#include "queue.hpp"
#include "Simulation.hpp"


template <typename T>
KPM_Vector<T,2u>::KPM_Vector(int mem, Simulation<T,2> & sim) :
  KPM_VectorBasis<T,2>(mem, sim),
  r(sim.r),
  tile{r.Ld[0], 1},
  tile_ghosts{1, r.Ld[0]},
  transf_max{r.ld[1], r.Ld[0]},
  h(sim.h),
  x(r.Ld),
  std(x.basis[1])
  {
    unsigned d;
    Coordinates <std::size_t, 3>     z(r.Ld);
    Coordinates <int, 3> x(r.nd), dist(r.nd);
    
    mult_t1_ghost_cor = new T**[r.Orb];
    for(unsigned io = 0; io < r.Orb; io++)
      {
	mult_t1_ghost_cor[io] = new T*[h.hr.NHoppings(io)];
	for(unsigned ib = 0; ib < h.hr.NHoppings(io); ib++)
	  mult_t1_ghost_cor[io][ib] = new T[TILE];
      }
    
    for(unsigned d = 0; d < 2; d++)
      for(unsigned b = 0; b < 2; b++)
	{
	  MemIndBeg[d][b] = new std::size_t[r.Orb];
	  MemIndEnd[d][b] = new std::size_t[r.Orb];
	}
    
    
    for(unsigned d = 0; d < 2u; d++)
      for(unsigned edge = 0; edge < 2; edge++)
	transf_bound[d][edge] = (r.boundary[d][edge] == true ? r.ld[D - 1 - d] : 0);
    
    for(unsigned edge = 0 ; edge < 2; edge++)
      if(r.boundary[1][edge] == 1) // Botton/Top edge Increase the size due left and right corners
	transf_bound[1][edge] += (r.boundary[0][0] == 1 ? NGHOSTS : 0 ) + (r.boundary[0][1] == 1 ? NGHOSTS : 0 ) ;
    
    for(std::size_t io = 0; io < r.Orb; io++)
      {
	d = 0;
	// Position of initial corner to copy the source
	MemIndBeg[d][0][io] = z.set({std::size_t(NGHOSTS),               std::size_t(NGHOSTS), io}).index;   
	MemIndBeg[d][1][io] = z.set({std::size_t(r.Ld[0] - 2 * NGHOSTS), std::size_t(NGHOSTS), io}).index;   
	// Position of initial corner to copy to the destiny 
	MemIndEnd[d][0][io] = z.set({std::size_t(0),                     std::size_t(NGHOSTS), io}).index;   
	MemIndEnd[d][1][io] = z.set({std::size_t(r.Ld[0] - NGHOSTS),     std::size_t(NGHOSTS), io}).index;   
	
	d = 1;
	// Bottom edge 
	MemIndBeg[d][0][io] = z.set({std::size_t(NGHOSTS), std::size_t(NGHOSTS),             io}).index;
	MemIndEnd[d][0][io] = z.set({std::size_t(NGHOSTS), std::size_t(0),                   io}).index;          
	if(r.boundary[1][0] == 1 && r.boundary[0][0] == 1) // Add the Left bottom Corner
	  { 
	    MemIndBeg[1][0][io] = z.set({std::size_t(0), std::size_t(NGHOSTS),             io}).index;
	    MemIndEnd[1][0][io] = z.set({std::size_t(0), std::size_t(0),                   io}).index;
	  } 
	// Top Edge 
	MemIndBeg[d][1][io] = z.set({std::size_t(NGHOSTS), std::size_t(r.Ld[1] - 2*NGHOSTS), io}).index;
	MemIndEnd[d][1][io] = z.set({std::size_t(NGHOSTS), std::size_t(r.Ld[1] - NGHOSTS),   io}).index;
	if(r.boundary[1][1] == 1 && r.boundary[0][0] == 1) // Add the Left bottom Corner
	  {
	    MemIndBeg[1][1][io] = z.set({std::size_t(0), std::size_t(r.Ld[1] - 2*NGHOSTS), io}).index;
	    MemIndEnd[1][1][io] = z.set({std::size_t(0), std::size_t(r.Ld[1] - NGHOSTS),   io}).index;
	  }
      }
    
    for(d = 0 ; d < D; d++)
      for(unsigned b  = 0 ; b < 2; b++)
	{
	  dist.set({0,0,0});
	  dist.coord[d] = int(b) * 2 - 1;
	  block[d][b] = x.set_coord( int(r.thread_id) ).add(dist).index;
	}
    
    // JPPP Calcular os vectores das fases devido às condicoes fronteira
    for(unsigned i = 0; i < D; i ++)
      for(unsigned j = 0; j < 3; j ++)
	Fact_Bnd[i][j] = new T[r.Ld[i]];
    
    initiate_vector();
  }

template <typename T>
KPM_Vector <T, 2>::~KPM_Vector(void){
  for(unsigned d = 0; d < 2; d++)
    for(unsigned b = 0; b < 2; b++)
      {
        delete MemIndBeg[d][b];
        delete MemIndEnd[d][b];
      }
  for(unsigned io = 0; io < r.Orb;io++)
    {	
      for(unsigned ib = 0; ib < h.hr.NHoppings(io); ib++)
        delete mult_t1_ghost_cor[io][ib];
      delete mult_t1_ghost_cor[io];
    }
  delete mult_t1_ghost_cor;
}

template <typename T>
void KPM_Vector <T, 2>::initiate_vector() {
  index = 0;
  
  // Check if the SEED variable is set to deterministic
  char *env;
  std::string seed("");
  env = getenv("SEED");
  if(env != NULL) seed = std::string(env);

  if(seed=="deterministic"){
      // For testing purposes only. Set only one site (0,0,0) to one, all else to zero
      // with PBC and translation invariance, this must yield the exact value of the trace
      // restricted to orbital 0

      v.setZero();
      Coordinates<std::size_t, 3> coord_Lt(r.Lt); // global sample, without ghosts
      Coordinates<std::size_t, 3> coord_Ld(r.Ld); // domain, with ghosts

      coord_Ld.set({NGHOSTS, NGHOSTS,0});
      r.convertCoordinates(coord_Lt, coord_Ld);

      unsigned x,y,o;
      x = coord_Lt.coord[0]; 
      y = coord_Lt.coord[1]; 
      o = coord_Lt.coord[2]; 
      if(x == 0 && y == 0 && o == 0) v(coord_Ld.index,index) = 1; 

  } else if(seed == "ones"){
      Coordinates<std::size_t, 3> x(r.Ld);
      for(std::size_t io = 0; io < r.Orb; io++)
        for(std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1++)
          for(std::size_t i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0++)
            v(x.set({i0,i1,io}).index, index) = 1.0/static_cast<value_type>(sqrt(value_type(r.Sizet - r.SizetVacancies)));


  } else {
      // Proceed as normal
      // The RNG is called inside the Random.cpp file, and there the code checks for a SEED again
      Coordinates<std::size_t, 3> x(r.Ld);
      for(std::size_t io = 0; io < r.Orb; io++)
        for(std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1++)
          for(std::size_t i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0++)
            v(x.set({i0,i1,io}).index, index) = simul.rnd.init()/static_cast<value_type>(sqrt(value_type(r.Sizet - r.SizetVacancies)));
  }




  for(unsigned i = 0; i < r.NStr; i++)
    {
      auto & vv = h.hV.position.at(i); 
      for(unsigned j = 0; j < vv.size(); j++)
        v(vv.at(j), index ) = 0. ;
    }
  initiate_phases();
}

template <typename T>
void KPM_Vector <T, 2>::initiate_phases() {
  // JPPP Calcular os vectores das fases devido às condicoes fronteira
  Coordinates<std::size_t, 3> x(r.Ld), z(r.Lt);
#pragma omp critical
  {
    for(unsigned i = NGHOSTS; i < r.Ld[0] - NGHOSTS; i++)
      {
	x.set({i, 0, 0});
	r.convertCoordinates(z,x);
	Fact_Bnd[0][2][x.coord[0]] = exp(ComplexTraits<T>::assign_value(0.0,-1.0*(h.BoundTwist[0] * int((int(z.coord[0]) + 1)/r.Lt[0]))));
	Fact_Bnd[0][1][x.coord[0]] = ComplexTraits<T>::assign_value(1.0,0.0);
	Fact_Bnd[0][0][x.coord[0]] = exp(ComplexTraits<T>::assign_value(0.0,h.BoundTwist[0] * int((int(r.Lt[0] - z.coord[0]))/r.Lt[0])));
      }
    for(unsigned i = NGHOSTS; i < r.Ld[1] - NGHOSTS; i++)
      {
	x.set({0, i, 0});
	r.convertCoordinates(z,x);
	Fact_Bnd[1][2][x.coord[1]] = exp(ComplexTraits<T>::assign_value(0.0,-1.0*(h.BoundTwist[1] * int((int(z.coord[1]) + 1)/r.Lt[1]))));
	Fact_Bnd[1][1][x.coord[1]] = ComplexTraits<T>::assign_value(1.0,0.0);
	Fact_Bnd[1][0][x.coord[1]] = exp(ComplexTraits<T>::assign_value(0.0,h.BoundTwist[1] * int((int(r.Lt[1] - z.coord[1]))/r.Lt[1])));
      }
  }

  // Testing initiate phases
  unsigned i=2;
	x.set({0, i, 0});
	r.convertCoordinates(z,x);

}

template <typename T>
T KPM_Vector <T, 2>::get_point()
{
  T point;
  if(r.thread_id == 0)
    {
      Coordinates<std::size_t, 3> x(r.Ld);
      std::size_t indice = x.set({std::size_t(r.Ld[0]/2),std::size_t(r.Ld[1]/2), std::size_t(0)}).index;
      point = v(indice,index);
    }
  else
    point = assign_value(0,0);
  return point;
}


template <typename T>
void KPM_Vector <T, 2>::build_site(unsigned long pos){
    // Builds an initial vector which is zero everywhere except
    // for a single site, where it is one
    Coordinates<unsigned long, 3> thread_coords(r.ld);
    Coordinates<unsigned long, 3> thread_coords_gh(r.Ld);
    Coordinates<unsigned long, 3> thread(r.nd);
    Coordinates<unsigned long, 3> total_coords(r.Lt);
    bool correct_thread;
    unsigned long T_thread[3]; // index of the thread
    unsigned long x_thread[3]; // position within the thread

    index = 0;
#pragma omp critical
{
    total_coords.set_coord(pos);
    for(unsigned d = 0; d < 2; d++){
        T_thread[d] = total_coords.coord[d]/r.ld[d];
        x_thread[d] = total_coords.coord[d]%r.ld[d];
    }
    T_thread[2] = 0;
    x_thread[2] = total_coords.coord[2];
    thread_coords.set_index(x_thread);
    thread.set_index(T_thread);

    //convert to coordinates with ghosts
    r.convertCoordinates(thread_coords_gh, thread_coords); 
 
    // check if the site is in the current thread
    correct_thread = thread.index == r.thread_id;
}
 
#pragma omp barrier
 v.setZero();
 v(thread_coords_gh.index,0) = T(correct_thread);
}

template <typename T>
void KPM_Vector <T, 2>::build_planewave(Eigen::Matrix<double,-1,1> & k, Eigen::Matrix<T,-1,1> & weight){
    // Builds an initial vector which is a plane wave with a specific value of k
    // weight is the weight of each orbital for this plane wave 
    // |k> = sum_{r,R} w(R) exp(i k.r + i k.R) |r,R>
    // r = lattice vector
    // R = orbital vector
    // w = weight (only depends on orbital)

    if(weight.rows() != r.Orb)
        std::cout << "Warning in build_planewave: the weight matrix must have the same number"
            " of elements as the number of orbitals.";


    index = 0;    // sets the KPM index to 0
    Coordinates<std::size_t, 3> local_coords(r.Ld), global_coords(r.Lt);

    Eigen::Map<Eigen::Matrix<std::size_t, 2, 1>> position(global_coords.coord); // spacial part of the coord vector in global_coords
    auto orb_a_coords = r.rLat.inverse() * r.rOrb;          // column i is the position of the i-th orbital in basis a
                                                            // r.rLat.inverse() each row is a bi / 2*M_PI 
    Eigen::Array<T, -1, 1> exp_R;
    exp_R = Eigen::Array<T, -1, 1>::Zero(r.Orb, 1);   // exponential related to the orbital
    T exp_r;                                                // exponential related to the lattice site

    // Calculate the exponential related to the orbital exp(i k R) w(R)
    // It is already divided by the norm, which is the number of total sites r.Nt
    for(std::size_t io = 0; io < r.Orb; io++)
        exp_R(io) = weight(io)*exp(assign_value(0, 2.0*M_PI*orb_a_coords.col(io).transpose()*k))/T(sqrt(r.Nt));

    // Calculate the exponential related to each unit cell
    for(std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1++)
      for(std::size_t i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0++)
        {
          local_coords.set({i0,i1,std::size_t(0)});
          r.convertCoordinates(global_coords, local_coords);          // Converts the coordinates within a thread to global coordinates
          
          exp_r = exp(assign_value(0,  2.0*M_PI * position.cast<double>().transpose()*k ));
          
          for(std::size_t io = 0; io < r.Orb; io++){
            local_coords.set({i0,i1,io});
            v(local_coords.index, 0) = exp_r*exp_R(io);
          }
        }

    KPM_VectorBasis<T,2u>::build_defect_planewave(k, weight);
}


template <typename T>
void KPM_Vector <T, 2>::build_wave_packet(Eigen::Matrix<double,-1,-1> & k, Eigen::Matrix<T,-1,-1> & psi0, double & sigma,
                                          Eigen::Matrix<double, 1, 2> & vb)
{
  index = 0;
  Coordinates<std::size_t, 3> x(r.Ld), z(r.Lt);
  Eigen::Matrix<T,-1,-1>  sum(r.Orb, 1);
  Eigen::Map<Eigen::Matrix<std::size_t,2, 1>> vv(z.coord);
  Eigen::Matrix<double, 1, 2> va;
  Eigen::Matrix<double, -1,-1> phase(1, psi0.cols());
  T soma = assign_value(0,0);
  auto bbs = r.rLat.inverse(); // each row is a bi / 2*M_PI 
  auto vOrb = bbs * r.rOrb;    // each columns is a vector in a1 basis

  double a00 = r.rLat.col(0).transpose()*r.rLat.col(0);
  double a11 = r.rLat.col(1).transpose()*r.rLat.col(1);
  double a01 = r.rLat.col(0).transpose()*r.rLat.col(1);
    
#pragma omp master 
  {
    simul.Global.mu = Eigen::Matrix<T,-1,-1>(psi0.cols(), 1);      
    for(int ik = 0; ik < psi0.cols(); ik++)
      simul.Global.mu(ik, 0) = assign_value( simul.rnd.get(), 0 );
  }
#pragma omp barrier

#pragma omp critical
  for(int ik = 0; ik < psi0.cols(); ik++)
    phase(0,ik)  = std::real(simul.Global.mu(ik, 0)) ;
#pragma omp barrier
  for(std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1++)
    for(std::size_t i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0++)
      {
        x.set({i0,i1,std::size_t(0)});
        r.convertCoordinates(z,x);
        double n1 = (double(z.coord[1]) - double(vb(0, 1)))/sigma;
        double n0 = (double(z.coord[0]) - double(vb(0, 0)))/sigma;
        double gauss = n0*n0 * a00 + n1*n1 * a11 + 2*n0*n1*a01;
        sum.setZero();
        va = vv.cast<double>().transpose();
	  
        for(int ik = 0; ik < psi0.cols(); ik++)
          {
            for(unsigned io = 0; io < r.Orb; io++)
              {
                double xx = (va  + vOrb.col(io).transpose() ) * k.col(ik);
                sum(io, 0) += psi0(io,ik) * exp(assign_value(-0.5*gauss,  2.*M_PI * (xx + 0.*phase(0,ik) )));
              }
          }
	  
        for(std::size_t io = 0; io < r.Orb; io++)
          {
            x.set({i0,i1,io});
            v(x.index, 0) = sum(io,0);
            soma += assign_value( std::real( sum(io,0) * std::conj(sum(io,0)) ), 0);
          }
      }
    
#pragma omp master
  simul.Global.soma = assign_value(0,0);
#pragma omp barrier
#pragma omp critical
  simul.Global.soma += soma;
#pragma omp barrier
  soma = simul.Global.soma;
  value_type soma2 = sqrt(std::real(soma));
    
  for(std::size_t io = 0; io < r.Orb; io++)
    for(std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1++)
      for(std::size_t i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0++)
        {
          x.set({i0,i1,io});
          v(x.index, 0) /= soma2;
        }
}

template <typename T>
template < unsigned MULT,bool VELOCITY> 
void KPM_Vector <T, 2>::build_regular_phases(int i1, unsigned axis)
{
  unsigned l[2 + 1], count;
  Coordinates<std::ptrdiff_t, 3>  global(r.Lt);
  Coordinates<std::ptrdiff_t, 3> local1(r.Ld);
  Coordinates<std::size_t,3> x(r.Ld);
  std::fill_n(l, 2, 3);
  l[2]  = r.Orb;
  Coordinates<std::ptrdiff_t, 2 + 1> b3(l);    
  Eigen::Map<Eigen::Matrix<std::ptrdiff_t,2, 1>> vee(b3.coord); // Column vector
    
    
  for(unsigned io = 0; io < r.Orb; io++)
    {
      const std::size_t ip = io * x.basis[2] ;
      const std::size_t j0 = ip + 0 + i1 * std;
      const std::size_t j1 = j0 + TILE * std;

      for(unsigned ib = 0; ib < h.hr.NHoppings(io); ib++)
        {
          T tt  = value_type(MULT + 1) * h.hr.hopping(ib, io);
          b3.set_coord(h.hr.dist(ib,io));        
          vee.array() -= 1;
          count = 0;
          if (VELOCITY)  tt  *=  h.hr.v.at(axis)(ib,io);
          for(std::size_t j = j0; j < j1; j += std )
            {
              r.convertCoordinates(global, local1.set_coord(j));
              value_type phase = vee(0)*global.coord[1]*r.ghost_pot(0,1);
              mult_t1_ghost_cor[io][ib][count] =  tt * multEiphase(phase);
              count++;
            }
        }
    }
}

template <typename T>
template < unsigned MULT> 
void KPM_Vector <T, 2>::initiate_stride(std::size_t & istr)
{
  std::size_t i0, i1;
  const std::size_t std = x.basis[1];
  // Periodic component of the Hamiltonian + Anderson disorder
  i0 = ((istr) % (r.lStr[0]) ) * TILE + NGHOSTS;
  i1 = ((istr) / r.lStr[0] ) * TILE + NGHOSTS;
  
  for(std::size_t io = 0; io < r.Orb; io++)
    {
      const std::size_t ip = io * x.basis[2];
      const std::size_t j0 = ip + i0 + i1 * std;
      const std::size_t j1 = j0 + TILE * std; //j0 and j1 define the limits of the for cycle


      for(std::size_t j = j0; j < j1; j += std )
        for(std::size_t i = j; i < j + TILE ; i++)
          phi0[i] = - value_type(MULT) * phiM2[i];
    }
}

template <typename T>
template < unsigned MULT> 
void inline KPM_Vector <T, 2>::mult_local_disorder(const  std::size_t & j0, const  std::size_t & io)
{
  const std::size_t j1 = j0 + TILE * std;
  const std::ptrdiff_t dd = (h.Anderson_orb_address[io] - std::ptrdiff_t(io))*r.Nd;
  // Anderson disorder
  if( h.Anderson_orb_address[io] >= 0)
    {
      for(std::size_t j = j0; j < j1; j += std )
        for(std::size_t i = j; i < j + TILE ; i++)
          phi0[i] += value_type(MULT + 1) * phiM1[i] * h.U_Anderson.at(i + dd);
    }
  else if (h.Anderson_orb_address[io] == - 1)
    {
      for(std::size_t j = j0; j < j1; j += std )
        for(std::size_t i = j; i < j + TILE ; i++)
          phi0[i] += value_type(MULT + 1) * phiM1[i] * h.U_Orbital.at(io);
    }

  // Custom local disorder
  if(h.is_custom_local_set){
    for(std::size_t j = j0; j < j1; j += std )
      for(std::size_t i = j; i < j + TILE ; i++)
        phi0[i] += value_type(MULT + 1) * phiM1[i] * h.custom_local(i);
  }
}

template <typename T>
void inline KPM_Vector <T, 2>::mult_regular_hoppings(const  std::size_t & j0, const  std::size_t & io){
  std::size_t count;
  const std::size_t j1 = j0 + TILE * std;
  std::size_t rr[2], hop[2], x, y; // Variables for TBC
  // Hoppings
  for(unsigned ib = 0; ib < h.hr.NHoppings(io); ib++){
      rr[0] = (j0 % r.Ld[0]);
      rr[1] = (j0 % (r.Ld[0] * r.Ld[1]))/(r.Ld[0]);
      
      const std::ptrdiff_t d1 = h.hr.distance(ib, io);
      // Determine the Supercell Jumps
      const std::size_t i_f = j0 + d1;
      hop[0] = (i_f % r.Ld[0] ) - rr[0] +1 ;
      hop[1] = (i_f % (r.Ld[0] * r.Ld[1]))/(r.Ld[0]) - rr[1] + 1;
      
      count = 0;
      y = 0;
      for(std::size_t j = j0; j < j1; j += std ) {
	const T t1 = mult_t1_ghost_cor[io][ib][count++] * Fact_Bnd[1][hop[1]][rr[1]+y];
	x = 0;
	for(std::size_t i = j; i < j + TILE ; i++) {
	  
	  phi0[i] += t1 * phiM1[i + d1] * Fact_Bnd[0][hop[0]][rr[0]+x];
	  x++;
	};
	y++;
      };
  }
}

template <typename T>
template <unsigned MULT, bool VELOCITY>
void KPM_Vector <T, 2>::KPM_MOTOR(KPM_Vector<T,2> *kpm_final, unsigned axis)
{
  std::size_t i0, i1;    
  phi0 = kpm_final->v.col(kpm_final->index).data();
  phiM1 = v.col( (memory - 1 + index) % memory ).data();
  phiM2 = v.col( (memory - 2 + index) % memory ).data();
    
  // Initialize tiles that have deffects connecting elements of a previous tile
  for(auto istr = h.cross_mozaic_indexes.begin(); istr != h.cross_mozaic_indexes.end() ; istr++)
    initiate_stride<MULT>(*istr);
    
  for( i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; i1 += TILE  ){
      build_regular_phases<MULT,VELOCITY>(i1, axis);
		
      for( i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; i0 += TILE ){
		    
          std::size_t istr = (i1 - NGHOSTS) / TILE * r.lStr[0] + (i0 - NGHOSTS) / TILE;
          if(h.cross_mozaic.at(istr))
            initiate_stride<MULT>(istr);
          // These four lines pertrain only to the magnetic field
          for(std::size_t io = 0; io < r.Orb; io++)
            {
              const std::size_t ip = io * x.basis[2];
              const std::size_t j0 = ip + i0 + i1 * std;
		
              // Local Energy
              if(!VELOCITY) mult_local_disorder<MULT>(j0, io);
		
              // Hoppings
              mult_regular_hoppings(j0, io);
            }

          KPM_VectorBasis<T,2u>::template multiply_defect<MULT, VELOCITY>(istr, phi0, phiM1, axis);
	  	    
          // Empty the vacancies in the tile
          auto & hV = h.hV.position.at(istr);
          for(auto k = hV.begin(); k != hV.end(); k++)
            phi0[*k] = 0.;

        }
    }

  for(auto vc =  h.hV.vacancies_with_defects.begin(); vc != h.hV.vacancies_with_defects.end(); vc++)
    phi0[*vc] = 0.;

    
  /* 
     Broken Imputirities:
     The bulk domain will receive contributions from the broken defects 
     located on the neighbour domains.
     We already subtract the vacancies from these contributions 
  */
    
  for(auto id = h.hd.begin(); id != h.hd.end(); id++)
    id->template multiply_broken_defect<MULT,VELOCITY>(phi0, phiM1,axis);
	  
  Exchange_Boundaries();


}
template <typename T>
void KPM_Vector <T, 2>::measure_wave_packet(T * bra, T * ket, T * results)  
{
  Coordinates<std::size_t,3> ad(r.Ld), at(r.Lt);
  ad.set({std::size_t(NGHOSTS), std::size_t(NGHOSTS), std::size_t(0)});
  r.convertCoordinates(at,ad);
  for(unsigned i = 0; i < 4; i++)
    results[i] *= 0.;

  for(unsigned io = 0; io < r.Orb; io++)
    {
      value_type deltax = r.rOrb(0,io);
      value_type deltay = r.rOrb(1,io);
	
      for(unsigned i1 = 0; i1 < r.ld[1]; i1++)
        {
          std::size_t ind = ad.set({std::size_t(NGHOSTS),std::size_t(i1 + NGHOSTS), std::size_t(io)}).index;
          value_type z1 = at.coord[1] + i1;
          value_type x0 = at.coord[0]*r.rLat(0,0) + z1 * r.rLat(0,1) + deltax;
          value_type y0 = at.coord[0]*r.rLat(1,0) + z1 * r.rLat(1,1) + deltay;	    

          T xl1 = assign_value(0., 0.);
          T xl2 = assign_value(0., 0.);
          T yl1 = assign_value(0., 0.);
          T yl2 = assign_value(0., 0.);
	    
          for(unsigned i0 = 0; i0 < r.ld[0]; i0++)
            {
              std::size_t j0 = ind + i0;
              value_type x = x0 + i0 * r.rLat(0,0);
              value_type y = y0 + i0 * r.rLat(1,0);
              T p = myconj(*(bra + j0)) * (*(ket + j0));
              xl1 += p * x;
              xl2 += p * x * x;
              yl1 += p * y;
              yl2 += p * y * y;
            }
    
          results[0] += xl1;
          results[1] += xl2;
          results[2] += yl1;
          results[3] += yl2;
        }
    }
  // Periodic component of the Hamiltonian + Anderson disorder
}
template <typename T>
void KPM_Vector <T, 2>::Exchange_Boundaries() {
  /*
    I have four boundaries to exchange with the other threads.
    First I will copy the lines along the a[1] direction to a consecutive shared vector
  */
#pragma omp barrier    
  Coordinates<std::size_t,3u> x(r.Ld), z(r.Lt);
  T  *phi = v.col(index).data();

    
  for(unsigned d = 0; d < 2; d++)
    {
      std::size_t BSize = r.Orb * transf_max[d] * NGHOSTS;
      T * ghosts_left = & simul.ghosts[0];
      T * ghosts_right = & simul.ghosts[BSize];

      for(std::size_t io = 0; io < r.Orb; io++)
        {
          std::size_t il = MemIndBeg[d][0][io];
          std::size_t ir = MemIndBeg[d][1][io];
	    
          for(std::size_t i = 0; i < transf_bound[d][0]; i++)
            {
              for(unsigned ig = 0; ig < NGHOSTS; ig++)
                ghosts_left [i + (ig + NGHOSTS*io) * transf_bound[d][0] ] = phi[il + ig * tile_ghosts[d]];
              il += tile[d];
            }
	    
          for(std::size_t i = 0; i < transf_bound[d][1]; i++)
            {
              for(unsigned ig = 0; ig < NGHOSTS; ig++)
                ghosts_right[i + (ig + NGHOSTS*io) * transf_bound[d][1] ] = phi[ir + ig * tile_ghosts[d]];
              ir += tile[d];
            }
        }
	
      // Copy the boundaries to the shared memory
      std::copy( ghosts_left, ghosts_left + 2*BSize, simul.Global.ghosts.begin() + 2*BSize * r.thread_id );	  
#pragma omp barrier
      auto neigh_left = simul.Global.ghosts.begin() + 2 * block[d][0] * BSize;
      auto neigh_right  = simul.Global.ghosts.begin() + 2 * block[d][1] * BSize;
      std::copy(neigh_right,         neigh_right + BSize , ghosts_right );     // From the left to the right
      std::copy(neigh_left + BSize,  neigh_left + 2*BSize, ghosts_left  )  ;   // From the right to the left
	
#pragma omp barrier	
      for(std::size_t io = 0; io < r.Orb; io++)
        {
          std::size_t il = MemIndEnd[d][0][io];
          std::size_t ir = MemIndEnd[d][1][io];
	    
          for(std::size_t i = 0; i < transf_bound[d][0]; i++)
            {
              for(int ig = 0; ig < NGHOSTS; ig++)
                phi[il + ig * tile_ghosts[d]] = ghosts_left [i + (ig + NGHOSTS * io) * transf_bound[d][0]];
              il += tile[d];
            }
	    
          for(std::size_t i = 0; i < transf_bound[d][1]; i++)
            {
              for(int ig = 0; ig < NGHOSTS; ig++)
                phi[ir + ig * tile_ghosts[d]] = ghosts_right[i + (ig + NGHOSTS * io) * transf_bound[d][1]];
              ir += tile[d];
            }
        }
    }
}

template <typename T>
void KPM_Vector <T, 2>::test_boundaries_system() {

  /*
    This  function tests if the boudaries exchange are well implemented 
  */

  Coordinates<std::size_t, 3> z(r.Lt);
  Coordinates<std::size_t, 3> x(r.Ld);
    
  for(std::size_t  io = 0; io < (std::size_t) r.Ld[2]; io++)
    for(std::size_t i1 = NGHOSTS; i1 < (std::size_t) r.Ld[1] - NGHOSTS ; i1++)
      for(std::size_t i0 = NGHOSTS; i0 < (std::size_t) r.Ld[0] - NGHOSTS ; i0++)
        {
          r.convertCoordinates(z, x.set({i0,i1,io}) );
          v(x.set({i0,i1,io}).index, 0) = aux_wr(z.index);
        }
    
  Exchange_Boundaries();
#pragma omp barrier
    
#pragma omp critical
  {
    for(std::size_t  io = 0; io < (std::size_t) r.Ld[2]*0 + 1; io++)
      for(std::size_t i1 = 0; i1 < (std::size_t) r.Ld[1] ; i1++)
        {
          for(std::size_t i0 = 0; i0 < (std::size_t) r.Ld[0]; i0++)
            {
              r.convertCoordinates(z, x.set({i0,i1,io}) );
              T val = aux_wr(z.index);
		
              if( aux_test(v(x.index , 0), val ) )
                {
                  std::cout << "Problems---->" << v(x.index , 0) << " " << val << " ";
                  std::cout << "\t wrong " << std::real(v(x.index , 0)) << " " << z.index << " " << x.index << "\t\t";
                  x.print();
		    
                }
            };
        }
  }
  std::cout << "Thread : " << r.thread_id << std::endl;
#pragma omp barrier
  exit(1);
}
template <typename T>
void KPM_Vector <T, 2>::empty_ghosts(int mem_index) {
  /* This function takes the kpm vector that's being used, 'v' and sets to zero the part corresponding
   * to the ghosts, that is, the part of the vector that actually belongs to a different thread.
   * This is done so that when we take the dot product 'v' with another vector only terms pertraining 
   * to the current thread are considered.
   * */
  
  Coordinates<long, 3u> x(r.Ld);
  
  
  // There are four sides, so set the ghosts in each side to zero individually.
  // Remember that the size of the ghost boundaries depends on NGHOSTS.
    
  for(long  io = 0; io < (long) r.Ld[2]; io++)
    for(long i0 = 0; i0 < (long) r.Ld[0]; i0++)
      for(int d = 0; d < NGHOSTS; d++)
        v(x.set({i0,(long) d,io}).index, mem_index) *= 0;

  for(long  io = 0; io < (long) r.Ld[2]; io++)
    for(long i0 = 0; i0 < (long) r.Ld[0]; i0++)
      for(int d = 0; d < NGHOSTS; d++)
        v(x.set({i0, (long) (r.Ld[1] - 1 - d),io}).index, mem_index) *= 0;
  
  for(long  io = 0; io < (long) r.Ld[2]; io++)
    for(long i1 = 0; i1 < (long) r.Ld[1]; i1++)
      for(int d = 0; d < NGHOSTS; d++)
        v(x.set({(long) d,i1,io}).index, mem_index) *= 0;

  for(long  io = 0; io < (long) r.Ld[2]; io++)
    for(long i1 = 0; i1 < (long) r.Ld[1]; i1++)
      for(int d = 0; d < NGHOSTS; d++)
        v(x.set({(long) (r.Ld[0] - 1 - d),i1,io}).index, mem_index) *= 0;

}



#define instantiateTYPE(type)               template class KPM_Vector <type,2u>; \
  template void KPM_Vector<type,2u>::template KPM_MOTOR<0u,false>(KPM_Vector<type,2u> * kpm_final, unsigned axis); \
  template void KPM_Vector<type,2u>::template KPM_MOTOR<1u,false>(KPM_Vector<type,2u> * kpm_final, unsigned axis); \
  template void KPM_Vector<type,2u>::template KPM_MOTOR<0u,true>(KPM_Vector<type,2u> * kpm_final,  unsigned axis);

instantiateTYPE(float)
instantiateTYPE(double)
instantiateTYPE(long double)
instantiateTYPE(std::complex<float>)
instantiateTYPE(std::complex<double>)
instantiateTYPE(std::complex<long double>)






