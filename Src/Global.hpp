template <typename T>
struct GLOBAL_VARIABLES {
  std::vector<T> ghosts;
  std::vector<std::size_t>    element1;
  std::vector<std::ptrdiff_t> element2_diff;
  std::vector<T> hopping;
  std::vector<std::size_t>    element;
  std::vector<T> U;

  // Averages
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> mu;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> gammaxx;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> gammaxy;
  GLOBAL_VARIABLES() { };
  void addbond( std::size_t  ele1, std::ptrdiff_t ele2, T hop ) {
    element1.push_back(ele1);
    element2_diff.push_back(ele2); 
    hopping.push_back(hop);
  }

  void addlocal( std::size_t  ele,  T u ) {
    element.push_back(ele);
    U.push_back(u);
  }
};
