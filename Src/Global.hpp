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
  T soma; 

  // Averages
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> mu;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> gamma;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> lambda;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> singleshot_cond;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> general_gamma;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> smaller_gamma;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_x;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_y;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_z;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_ident;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_results;
  double kpm_iteration_time;
  
  bool calculate_arpes;
  bool calculate_ldos;
  bool calculate_wavepacket;
  bool calculate_dos;
  bool calculate_conddc;
  bool calculate_condopt;
  bool calculate_condopt2;
  bool calculate_singleshot;

  GLOBAL_VARIABLES();
  void addbond ( std::size_t, std::ptrdiff_t, T );
  void addlocal( std::size_t,  T);
};
