#include "Generic.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "ComplexTraits.hpp"
#include "Random.hpp"
#include "HamiltonianVacancies.hpp"
template <typename T,unsigned D>
Vacancy_Operator<T,D>::Vacancy_Operator(char * filename, LatticeStructure <D> & rr, KPMRandom <T> & rndB) : r(rr), rnd(rndB), name(filename), position(r.NStr)
{

}

template <typename T,unsigned D>
void Vacancy_Operator<T,D>::generate_disorder()
{
  Coordinates<std::size_t,D + 1> latt(r.ld), Latt(r.Ld), latStr(r.lStr), x(r.nd);
  // Clean former vacancy distribution
  
  for(unsigned i = 0; i < r.NStr ; i++)
    position.at(i).clear();
  vacancies_with_defects.clear();
  // Distribute Vacancies
  
  for(unsigned k = 0; k < concentration.size(); k++)
    {
      std::size_t vacancies_number  = concentration.at(k) * r.N, count = 0;
      // Test how many vacancies where in this subdomain
      if(vacancies_number < positions_fixed.at(k).size())
        vacancies_number =  positions_fixed.at(k).size();
      
      while(count < vacancies_number)
        {
          std::size_t i;
          if(count <  positions_fixed.at(k).size() )
            i = positions_fixed.at(k).at(count);
          else
            i = rnd.get()*r.N;
          
          latt.set_coord(i + orbitals.at(k).at(0) * r.N );
          r.convertCoordinates(latStr,latt);                  // Get tile position
          r.convertCoordinates(Latt,latt);                    // Get Domain coordinates
          auto & pos = position.at(latStr.index);
          if(!any_of(pos.begin(), pos.end(), std::bind2nd(std::equal_to<std::size_t>(), Latt.index )))
            {
              for(auto o = orbitals.at(k).begin(); o != orbitals.at(k).end(); o++)
                {
                  latt.set_coord(i + std::size_t(*o) * r.N);
                  r.convertCoordinates(Latt,latt);         
                  pos.push_back(Latt.index);
                }
              count++;
            }
        }
    }
  
  for(unsigned i = 0; i < r.NStr ; i++)
    std::sort (position.at(i).begin(), position.at(i).end());
}

template <typename T,unsigned D>
void Vacancy_Operator<T,D>::add_model(double p, std::vector <int> & orb, std::vector<int> & postmp)
{
  Coordinates<std::size_t,D + 1> latt(r.ld),  x(r.nd), LATT(r.Lt);
  std::size_t size_aux  = p * r.Nt * orb.size();
  std::vector<int> tmp;
  size_aux = (size_aux > postmp.size() ? size_aux : postmp.size());
  orbitals.push_back(orb);
  concentration.push_back(p);
  r.SizetVacancies += size_aux;
  
  for(size_t i = 0; i < postmp.size(); i++)
    {
      LATT.set_coord(std::size_t(postmp.at(i)));
      r.convertCoordinates(x, LATT);
      if(x.index == r.thread_id)
        {
          r.convertCoordinates(latt, LATT);
          if(!any_of(tmp.begin(), tmp.end(), std::bind2nd(std::equal_to<std::size_t>(), latt.index )))
            tmp.push_back(latt.index);
        }
    };
  positions_fixed.push_back(tmp);
}

template <typename T,unsigned D>
void Vacancy_Operator<T,D>::add_conflict_with_defect(std::size_t element, unsigned istride)
{
  std::vector<std::size_t> & v = position.at(istride);  
  for(unsigned i = 0; i < v.size(); i++)
    if(element == v.at(i))
      vacancies_with_defects.push_back(element);
}

template <typename T,unsigned D>
bool Vacancy_Operator<T,D>::test_vacancy(Coordinates<std::size_t,D + 1> & Latt)
{
  /*
    1 if is a vacancy
    0 if not
  */
  Coordinates<std::size_t,D + 1> latStr(r.lStr);
  r.convertCoordinates(latStr, Latt);
  auto & vc = position.at(latStr.index);
  if( !any_of(vc.begin(), vc.end(), std::bind2nd(std::equal_to<std::size_t>(), Latt.index)))
    return 0;
  else
    return 1;
}

template <typename T,unsigned D>
void Vacancy_Operator<T,D>::test_field( T * phi0 )
{
  // The field should be zero in the vacancies
  for(unsigned i = 0; i  < r.NStr; i++)
    for(auto vc = position.at(i).begin(); vc != position.at(i).end(); vc++)
      if( abs(phi0[*vc]) > __DBL_EPSILON__)
        {
          std::cout << "Disparate" << std::endl;
          exit(1);
        }
}


template struct Vacancy_Operator<float,1u>;
template struct Vacancy_Operator<double,1u>;
template struct Vacancy_Operator<long double,1u>;
template struct Vacancy_Operator<std::complex<float>,1u>;
template struct Vacancy_Operator<std::complex<double>,1u>;
template struct Vacancy_Operator<std::complex<long double>,1u>;


template struct Vacancy_Operator<float,2u>;
template struct Vacancy_Operator<double,2u>;
template struct Vacancy_Operator<long double,2u>;
template struct Vacancy_Operator<std::complex<float>,2u>;
template struct Vacancy_Operator<std::complex<double>,2u>;
template struct Vacancy_Operator<std::complex<long double>,2u>;

template struct Vacancy_Operator<float,3u>;
template struct Vacancy_Operator<double,3u>;
template struct Vacancy_Operator<long double,3u>;
template struct Vacancy_Operator<std::complex<float>,3u>;
template struct Vacancy_Operator<std::complex<double>,3u>;
template struct Vacancy_Operator<std::complex<long double>,3u>;
