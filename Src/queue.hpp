/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

class measurement_queue{
public:
  std::string direction_string;
  std::vector<int> NMoments;
  int NDisorder;
  int NRandom;
  std::string label;
  double time_length;
  measurement_queue(std::string, std::vector<int>, int, int, std::string );
  void embed_time(double avg_duration);
};

class singleshot_measurement_queue{
public:
  std::string direction_string;
  int NDisorder;
  int NRandom;
  Eigen::Array<double, -1, 1> singleshot_energies;
  Eigen::Array<double, -1, 1> singleshot_gammas;
  Eigen::Array<double, -1, 1> singleshot_preserve_disorders;
  Eigen::Array<int, -1, 1> NMoments;
  Eigen::Array<double, -1, -1> singleshot_energiesgammas;
  std::string label;
  double time_length;
  singleshot_measurement_queue(std::string, Eigen::Array<int, -1, 1>, int, int, std::string, Eigen::Array<double, -1, 1>, 
                               Eigen::Array<double, -1, 1>, Eigen::Array<double, -1, 1> );
  void  embed_time(double avg_duration);
};

std::string num2str2(int dir_num);
std::string num2str3(int dir_num);
std::string print_time(double);
std::vector<measurement_queue> fill_queue(char *);
std::vector<singleshot_measurement_queue> fill_singleshot_queue(char *name);
