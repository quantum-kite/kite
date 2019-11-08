/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T, unsigned DIM>
class conductivity_optical{
	public:


    bool isRequired; // was this quantity (conductivity_optical) asked for?
    bool isPossible; // do we have all we need to calculate the conductivity?


	// Functions to calculate. They will require the objects present in
    // the configuration file
    int direction;
    int NumDisorder;

    int NumMoments;
    int Moments_D, Moments_G;
    int Convergence_D, Convergence_G;
    bool default_Convergence_D, default_Convergence_G;
    //unsigned Nmax, Mmax;
    //unsigned nmax, mmax;
    bool Moments_divisible;
    int MaxMoments;
    bool default_NumMoments;

    int NumPoints;
    int NumRandoms;
    double units;

    double temperature;
    bool default_temperature;
    T beta;

    bool default_efermi;
    T e_fermi; 

    bool default_scat;
    T scat; 

    bool default_lim;
    T lim;

    int N_energies; 
    bool default_NEnergies;

    unsigned int N_omegas;
    double minFreq; 
    double maxFreq; 
    bool default_Nfreqs;
    bool default_minfreqs;
    bool default_maxfreqs;

    std::string name;
    std::string filename;
    bool default_filename;

    // information about the Hamiltonian
    system_info<T, DIM> systemInfo;

    // Objects required to successfully calculate the conductivity
    shell_input variables;

    // Objects required to successfully calculate the conductivity
    Eigen::Array<std::complex<T>, -1, -1> Gamma;
    Eigen::Array<std::complex<T>, -1, -1> Gamma_Padded;
    Eigen::Array<std::complex<T>, -1, -1> Lambda;
    Eigen::Array<std::complex<T>, -1, -1> Lambda_Padded;

    conductivity_optical(system_info<T, DIM>&, shell_input &);
    bool is_required();
    void set_default_parameters();
    bool fetch_parameters();
    void override_parameters();
    void printOpt();
    void calculate();
    void calculateBlocks();
	
};
