template <typename T>
struct GLOBAL_VARIABLES {
  std::vector<T> ghosts;
  std::vector<unsigned long> element1;
  std::vector<unsigned long> element2;
  std::vector<T> hopping;
  std::vector<unsigned long> element;
  std::vector<T> U;

  // Averages
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> mu;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> gammaxx;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> gammaxy;
  GLOBAL_VARIABLES() { };
};
