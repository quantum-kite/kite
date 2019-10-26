#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Hamiltonian.hpp"


extern "C" herr_t getMembers(hid_t loc_id, const char *name, void *opdata);
template <typename T, unsigned D>
Hamiltonian<T,D>::Hamiltonian(char *filename,  LatticeStructure<D> & rr, GLOBAL_VARIABLES <T> & gg) : name(filename), r(rr) , Global(gg),  hr(name, r), cross_mozaic(r.NStr), hV(name, rr, rnd)
{
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
  for (unsigned i = 0; i < model.size(); i++)
    {
      if(model.at(i) < 3)
        {
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

template class Hamiltonian<float,1u>;
template class Hamiltonian<double,1u>;
template class Hamiltonian<long double,1u>;
template class Hamiltonian<std::complex<float>,1u>;
template class Hamiltonian<std::complex<double>,1u>;
template class Hamiltonian<std::complex<long double>,1u>;

template class Hamiltonian<float,2u>;
template class Hamiltonian<double,2u>;
template class Hamiltonian<long double,2u>;
template class Hamiltonian<std::complex<float>,2u>;
template class Hamiltonian<std::complex<double>,2u>;
template class Hamiltonian<std::complex<long double>,2u>;

template class Hamiltonian<float,3u>;
template class Hamiltonian<double,3u>;
template class Hamiltonian<long double,3u>;
template class Hamiltonian<std::complex<float>,3u>;
template class Hamiltonian<std::complex<double>,3u>;
template class Hamiltonian<std::complex<long double>,3u>;
