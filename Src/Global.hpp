/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

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
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> gamma;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> lambda;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> singleshot_cond;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> general_gamma;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_x;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_y;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_z;
  double kpm_iteration_time;
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
