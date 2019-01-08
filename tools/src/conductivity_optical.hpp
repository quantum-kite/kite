/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T, unsigned DIM>
class conductivity_optical{
	H5::H5File file;
	public:


    bool isRequired; // was this quantity (conductivity_optical) asked for?
    bool isPossible; // do we have all we need to calculate the conductivity?


	// Functions to calculate. They will require the objects present in
    // the configuration file
    int direction;
    int NumDisorder;
    int NumMoments;
    int NumPoints;
    int NumRandoms;
    double temperature;
    double units;
    std::string filename;

    T beta;
    T e_fermi; 
    T scat; 
	int N_energies; 
	int N_omegas; 
    double minFreq; 
    double maxFreq; 
    T lim;


    // information about the Hamiltonian
    system_info<T, DIM> systemInfo;

    // Objects required to successfully calculate the conductivity
    shell_input variables;

    // Objects required to successfully calculate the conductivity
	Eigen::Array<std::complex<T>, -1, -1> Gamma;
	Eigen::Array<std::complex<T>, -1, -1> Lambda;

	std::string dirName;


    conductivity_optical(system_info<T, DIM>&, shell_input &);
	void fetch_parameters();
    void override_parameters();
    void calculate();
    void calculate_efficient();
	
};
