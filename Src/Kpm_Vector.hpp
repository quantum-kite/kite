template <typename T, unsigned D>
class KPM_VectorBasis {
protected:
  int index;
  const int memory;
  const Hamiltonian <T,D> & h;
  LatticeStructure <D> & r;
  KPMRandom <T> & rng;

  std::vector<T> & Border;
  std::vector<T> & GlobalBorder;
  
  
public:
  Eigen::Matrix <T, Eigen::Dynamic,  Eigen::Dynamic > v;
  KPM_VectorBasis(int mem, Hamiltonian <T,D> & h1, LatticeStructure <D>  & r1, std::vector<T> & B, std::vector<T> & GB) :
    memory(mem), h(h1), r(r1), rng(h1.rng), Border(B), GlobalBorder(GB) {
    index  = 0;
    v = Eigen::Matrix <T, Eigen::Dynamic,  Eigen::Dynamic >::Zero(r.Sized, memory);
  };
  
  void set_index(int i) {index = i;};
  void inc_index() {index = (index + 1) % memory;};  
  unsigned get_index(){return index;};
  template <typename U = T>
  typename std::enable_if<is_tt<std::complex, U>::value, U>::type aux_wr(unsigned long x ) {
    typedef typename extract_value_type<U>::value_type value_type;
    return U(value_type(x), value_type(2*x));
  };
  
  template <typename U = T>
  typename std::enable_if<!is_tt<std::complex, U>::value, U>::type aux_wr(unsigned long x ) {
    return U(x);
  };
  

  bool aux_test(T & x, T & y ) {
    return (abs(x - y) > std::numeric_limits<double>::epsilon());
  };
  
  
};


template <typename T, unsigned D>
class KPM_Vector : public KPM_VectorBasis <T,D> {
public:
  KPM_Vector(int mem, Hamiltonian <T,D> & h1, LatticeStructure <D>  & r1, std::vector<T> & B, std::vector<T> & GB ) :
    KPM_VectorBasis<T,D>(mem, h1, r1, B, GB ){};
  
  void initiate_vector() {};
  template <unsigned MULT>
  void Multiply(const int model){};
  void test_boundaries_system() {};
  void Exchange_Boundaries() {};
  inline void  HaIteration( const int model) { Multiply<0>(model); };
  inline void  ChIteration( const int model) { Multiply<1>(model); };
  
};


