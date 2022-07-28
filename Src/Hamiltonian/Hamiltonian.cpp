/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/



#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Hamiltonian.hpp"
#include "aux.hpp"



extern "C" herr_t getMembers(hid_t loc_id, const char *name, void *opdata);
template <typename T, unsigned D>
Hamiltonian<T,D>::Hamiltonian(char *filename,  LatticeStructure<D> & rr, GLOBAL_VARIABLES <T> & gg) : name(filename), r(rr) , Global(gg),  hr(name, r), cross_mozaic(r.NStr), hV(name, rr, rnd)
{
#pragma omp critical
  {
    H5::H5File *file = new H5::H5File(filename, H5F_ACC_RDONLY);
    get_hdf5<double>(&EnergyScale, file, (char *) "/EnergyScale");
    delete file;
  }
  

  is_custom_local_set = false;
  generate_custom_local();


  /* Anderson disorder */
  build_Anderson_disorder();
  build_vacancies_disorder();
  build_structural_disorder();
}

template <typename T, unsigned D>
void Hamiltonian<T,D>::generate_disorder()
{
  distribute_AndersonDisorder();
  for(std::size_t istr = 0; istr < r.NStr; istr++)
    cross_mozaic[istr] = true;
  cross_mozaic_indexes.clear();

  hV.generate_disorder();
  for(auto id = hd.begin(); id != hd.end(); id++)
    id->generate_disorder();
}

template <typename T, unsigned D>
void Hamiltonian<T,D>::generate_twists()
{
  BoundTwist.setZero(); 
#pragma omp master
  {
    // GlobBTwist is a 3D vector. The last dimension isn't used for 2D
    Global.GlobBTwist.setZero(); // Set the Global Boundary Twists to zero
    for(unsigned i = 0; i < D; i++) {
          //JPPP Fixed Case
          if(r.RandomBoundaries[i] == 0){
	        Global.GlobBTwist(i,0) = r.BdTwist[i];         
          }

            //JPPP Random value between -π and π
          if(r.RandomBoundaries[i] == 1){
	        Global.GlobBTwist(i,0) = 2 * M_PI * (0.5 - rnd.get()); 
          } 

      };
  };
#pragma omp barrier
#pragma omp critical
  {
    for(unsigned i = 0; i < D; i++)
      BoundTwist(i,0) = Global.GlobBTwist(i,0); // Hard Copy Global Boundary Twists to Local Variables
  };
}

template <typename T, unsigned D>
void Hamiltonian<T,D>::build_structural_disorder()
{
  
#pragma omp critical
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    // Test if there is a strutural disorder to build
    H5::Group  grp;
    std::vector<std::string> defects;
    try {
      H5::Exception::dontPrint();
      grp = file->openGroup("/Hamiltonian/StructuralDisorder");
      grp.iterateElems(grp.getObjName(), NULL, getMembers, static_cast<void*>(&defects) );
      
      for(auto id = defects.begin(); id != defects.end(); id++)
        hd.push_back(Defect_Operator<T,D> (*this,  *id, file ));
    }
    catch(H5::Exception& e) {
      // Do nothing
    };
    delete file;
  }
}


template <typename T, unsigned D>
void Hamiltonian<T,D>::build_vacancies_disorder()
{
  r.SizetVacancies = 0;
#pragma omp critical
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    // Test if there is vacancies to build
    H5::Group  grp;
    std::vector<int> tmp;
    double p;
    std::vector<int> orbit;
    std::vector<std::string> vacancies;
    int n;

    try {
      H5::Exception::dontPrint();
      grp = file->openGroup("/Hamiltonian/Vacancy");
      // Get the names of the Vacancies Types
      grp.iterateElems(grp.getObjName(), NULL, getMembers, static_cast<void*>(&vacancies)); 
      for(auto id = vacancies.begin(); id != vacancies.end(); id++)
        {
          std::string field = *id + std::string("/Concentration");
          try {
            H5::Exception::dontPrint();
            get_hdf5<double> ( &p, file, field );
          } catch(H5::Exception& e) {
            // Do nothing
            p = 0.;
          }; 
          

          field = *id + std::string("/FixPosition");  
          try {
            H5::Exception::dontPrint();
            H5::DataSet dataset = H5::DataSet(file->openDataSet(field));
            H5::DataSpace dataspace = H5::DataSpace(dataset.getSpace());  
            std::size_t num = dataspace.getSimpleExtentNpoints ();
            tmp.resize(num);
            dataspace.close();
            dataset.close();
            get_hdf5<int> ( tmp.data(), file, field );
          } catch(H5::Exception& e) {
          };
          
          field = *id + std::string("/NumOrbitals");
          get_hdf5<int> ( &n, file, field );
          orbit.resize(n);
          field = *id + std::string("/Orbitals");
          get_hdf5<int> ( orbit.data(), file, field );
          hV.add_model(p, orbit, tmp);
          tmp.clear();
        }
    }
    catch(H5::Exception& e) {
      // Do nothing
    }
    delete file;
  }
  
}

template <typename T, unsigned D>
void Hamiltonian<T,D>::build_Anderson_disorder() {
  /*
   * Gaussian      : 1
   * Uniform       : 2
   * Deterministic : 3
   */
  hsize_t dim[2];
#pragma omp critical
  {

    H5::H5File    *file      = new H5::H5File(name, H5F_ACC_RDONLY);
    try {
      H5::DataSet   dataset    = H5::DataSet(file->openDataSet("/Hamiltonian/Disorder/OrbitalNum"));
      H5::DataSpace dataspace  = dataset.getSpace();
      dataspace.getSimpleExtentDims(dim, NULL);;
      
      
      orb_num.resize(dim[0]*dim[1]);
      model.resize(dim[1]);
      mu.resize(dim[1]);
      sigma.resize(dim[1]);   
      get_hdf5<int>(orb_num.data(), file, (char *) "/Hamiltonian/Disorder/OrbitalNum");               // read the orbitals that have local disorder
      get_hdf5<int> (model.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderModelType");   // read the the type  of local disorder
      get_hdf5<double> (mu.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderMeanValue");   // read the the mean value
      get_hdf5<double> (sigma.data(), file, (char *) "/Hamiltonian/Disorder/OnsiteDisorderMeanStdv"); // read the the variance

    }
    catch (...){}
    delete file;
  }
    
  Anderson_orb_address.resize(r.Orb);
  U_Orbital.resize(r.Orb);
  Eigen::Array<int,-1,-1> vv = Eigen::Map<Eigen::Array<int,-1,-1>>(orb_num.data(), dim[1], dim[0]);
    
  std::fill_n ( Anderson_orb_address.begin(), r.Orb, -2 );
  std::fill_n ( U_Orbital.begin(), r.Orb,  0 );
  int sum = 0;
  for (unsigned i = 0; i < model.size(); i++){
      if(model.at(i) < 3){
          int count = 0;
          while(unsigned(count) < dim[0] &&  vv(i, count) != -1 )
            {
              int io = vv(i, count);
              Anderson_orb_address.at(io) = sum;
              count++;
            }
          sum++;
        }
      
      if(model.at(i) == 3) // Deterministic
        {
          int count = 0;
          while(count < vv.cols() &&  vv(i, count) != -1 )
            {
              int io = vv(i, count);
              Anderson_orb_address.at(io) = -1;
              U_Orbital.at(io) = mu.at(i);
              count++;
            }
        }
    }
    
  if(sum > 0)
    U_Anderson.resize( sum * r.Nd);    
}

template <typename T, unsigned D>
void Hamiltonian<T,D>::build_velocity(std::vector<unsigned> & components, unsigned n)
{
  hr.build_velocity(components, n);
  for(auto i = hd.begin(); i != hd.end(); i++)
    i->build_velocity(components, n);
}



template <typename T, unsigned D>  
void Hamiltonian<T,D>::distribute_AndersonDisorder()
{
    
  /*
   * Gaussian      : 1
   * Uniform       : 2
   * Deterministic : 3
   */
    
  int sum = 0;
  for (unsigned i = 0; i < model.size(); i++)
    if(model.at(i) == 1 )
      {
        for(unsigned j = 0; j < r.Nd; j++)
          U_Anderson.at(sum + j) = rnd.gaussian(mu.at(i), sigma.at(i)) ;
        sum += r.Nd;
      }
    else if ( model.at(i) == 2 )
      {
        for(unsigned j = 0; j < r.Nd; j++)
          U_Anderson.at(sum + j) = rnd.uniform(mu.at(i), sigma.at(i)) ;
        sum += r.Nd;
      }
}

template <typename T, unsigned D>
void Hamiltonian<T,D>::generate_custom_local(){
    // The user can use two kinds of functions to generate a local energy
    // - type 1: field defined in the whole lattice
    // - type 2: local energy defined around a fixed set of sites
    //
    // If the user chooses type 1, KITE will use a function defined outside of KITE
    // (by the user, inside the folder lib/) to calculate the value of the local energy
    // at each lattice point. --- There can only be one such function ---. 
    // Example: square quantum well in the middle of the lattice
    //
    // If the user chooses type 2, then KITE will use a different function defined
    // outside of KITE in the folder lib/ which defines a local potential around a 
    // set of lattice points. Several different kinds of functions can be used, all 
    // at the same time. These can be used with type 1 functions

    // type 1 details:
    // add_type1()

    // type 2 details: the user specifies the impurity positions and an identifier
    // add_type2(imp_pos1, 1)
    // add_type2(imp_pos2, 2)
    // add_type2(imp_pos3, 3)
    //
    //

    // Check if the user wants to define the local energy. This is an optimization flag
    // KITE does not need to store the local energy if it is not required
    int custom_required = -1;
    int print_flag = false;

#pragma omp barrier
#pragma omp critical
{
    try { 
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5(&custom_required, file, (char*)"Hamiltonian/CustomLocalEnergy");
      get_hdf5(&print_flag, file, (char*)"Hamiltonian/PrintCustomLocalEnergy");
      delete file;
    }
    catch(...) {
        std::cout << "Error: CustomLocalEnergy is not defined in the configuration file. Exiting program.\n";
        exit(1);
    }
}
#pragma omp barrier
    print_custom = print_flag;
    

    if(custom_required != 1) return;

    // Generate the custom local energy from an external library inside /lib
    custom_local = fetch_type1();
    is_custom_local_set = true;

}

template <typename T, unsigned D>
Eigen::Array<T,-1,-1> Hamiltonian<T,D>::fetch_type1(){
//void Hamiltonian<T,D>::fetch_type1(Eigen::Array<T, -1,-1> local){
    // This function could be in LatticeStructure, but that class is not templated
    // with respect to the variable T, so I chose to put it in Hamiltonian

    // The reason for this to be a separate function is that I might want to recycle it

    Eigen::Array<T, -1,-1> vec;
    vec = Eigen::Array<T, -1, -1>::Zero(r.Sized,1); // with ghosts

    Coordinates<std::size_t, D+1> coord_Lt(r.Lt); // total (without ghosts)
    Coordinates<std::size_t, D+1> coord_ld(r.ld); // domain without ghosts
    Coordinates<std::size_t, D+1> coord_Ld(r.Ld); // domain with ghosts

    double Escale, Eshift;


#pragma omp critical
    {
        H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
        get_hdf5(&Escale, file, (char*)"EnergyScale");
        get_hdf5(&Eshift, file, (char*)"EnergyShift");
        delete file;
    }
#pragma omp barrier

    double V;
    unsigned orb;
    T V_converted;

    Eigen::Matrix<double,D,1> vec_real = Eigen::Matrix<double,D,1>::Zero(D);

    std::string filename = std::string("local_potential") + std::to_string(omp_get_thread_num()) + std::string(".dat");
    std::ofstream potential_file;

    // Only runs if the user wants to print the local custom energy
    if(print_custom)
        potential_file = std::ofstream(filename);

    for(unsigned i=0; i<r.N*r.Orb; i++){

        coord_ld.set_coord(i);                          // Convert from index to local coordinates
        r.convertCoordinates(coord_Lt, coord_ld);       // Convert from local to global coordinates, 
        orb = coord_Lt.coord[D];

        for(unsigned j = 0; j < D; j++)
            vec_real(j) = coord_Lt.coord[j];            // Set first coordinates 

        vec_real = r.rLat*vec_real;                     // Convert to real space
        V = uncorrelated_wrapper(vec_real.data(), orb); // and get the value of the potential

        V_converted = T((V - Eshift)/Escale); // type shenanigans

        // Convert to local coordinates with ghosts to get the index of the
        // KPM_Vector in which to store that value of the potential
        r.convertCoordinates(coord_Ld, coord_ld);
        vec(coord_Ld.index) = V_converted;

        // write to a stream
        if(print_custom){
            for(unsigned j = 0; j < D; j++)
                potential_file << vec_real(j) << " ";
            potential_file << orb << " ";
            potential_file << V_converted << "\n";
        }
    }

    if(print_custom)
        potential_file.close();

    return vec;
}







herr_t getMembers(hid_t loc_id, const char *name, void *opdata)
{
  H5::Group  grp(loc_id);
  std::string Disorder = grp.getObjName();
  std::string  sep = "/";
  std::string Defect = name;
  std::string group = Disorder+sep+name; 
  std::vector<std::string> * v = static_cast<std::vector<std::string> *> (opdata);
  
  try {
    H5::Exception::dontPrint();
    grp.openGroup(group);
    v->push_back(group);
  }
  catch(H5::Exception& e) {
    // Don't do nothing 
  }
  return 0;
}

#define instantiate(type, dim)               template class Hamiltonian<type,dim>;
#include "instantiate.hpp"

