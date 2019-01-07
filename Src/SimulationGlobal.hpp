/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/



template <typename T, unsigned D>
class GlobalSimulation {
private:
  GLOBAL_VARIABLES <T> Global;
  LatticeStructure <D> rglobal;
  // Regular quantities to calculate, such as DOS and CondXX
  Eigen::Array<double, -1, 1> singleshot_energies;
  double EnergyScale;
public:
  GlobalSimulation( char *);
};


