/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T, unsigned DIM>
class system_info{
	H5::H5File file;
	public:
		int dim;
		int num_orbitals;
    int isComplex;
		Eigen::Array<int,1,-1> size;
		Eigen::Array<double,-1,-1> vectors;
		Eigen::Array<double, -1, -1> orbital_positions;
			
		double unit_cell_area;
		double spin_degeneracy;
		double energy_scale;
		double energy_shift;
    std::string filename;
    int NumThreads;
        
    // These two quantities may only be known after the DOS is calculated
    bool EnergyLimitsKnown;
    T minEnergy, maxEnergy;

		system_info(std::string);
    system_info();
    void print_info();
		void read();
	
};
