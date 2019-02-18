/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T, unsigned DIM>
class conductivity_nonlinear{
    H5::H5File file;
    public:


    bool isRequired; // was this quantity (conductivity_dc) asked for?
    bool isPossible; // do we have all we need to calculate the conductivity?

    // KPM parameters
    int complex;
    int direction;
    int NumDisorder;
    int NumMoments;
    int NumPoints;
    int NumRandoms;
    int special;
    std::string dirString; 

    // Post-processing parameters
    std::complex<T> imaginary;
    std::string filename;
    double temperature;
    double ratio;
    int print_all;
    T beta;
    T e_fermi;
    T scat;

    bool default_temperature;
    bool default_scat;
    bool default_efermi;
    bool default_minfreq;
    bool default_maxfreq;
    bool default_Nfreqs;
    bool default_ratio;
    bool default_NEnergies;
    bool default_filename;
    
    // Energy parameters needed to run the simulation
    int N_energies;
    double lim;
    Eigen::Matrix<T, -1, 1> energies;

    // Frequency parameters needed to run the simulation
    int N_omegas;
    double minFreq;
    double maxFreq;
    Eigen::Matrix<T, -1, 1> frequencies;
    Eigen::Matrix<T, -1, 2> frequencies2;

    system_info<T, DIM> systemInfo; // information about the Hamiltonian
    shell_input variables;          // Input from the shell to override the configuration file

    // Objects required to successfully calculate the conductivity
    Eigen::Array<std::complex<T>, -1, -1> Gamma0;
    Eigen::Array<std::complex<T>, -1, -1> Gamma1;
    Eigen::Array<std::complex<T>, -1, -1> Gamma2;
    Eigen::Array<std::complex<T>, -1, -1> Gamma3;

    std::string dirName;



    // Functions
    conductivity_nonlinear(system_info<T, DIM>&, shell_input &);
    
    int fetch_gamma0();
    int fetch_gamma1();
    int fetch_gamma2();
    int fetch_gamma3();

    void printOpt2();
    bool fetch_parameters();
    void override_parameters();
    void set_default_parameters();
    void calculate();

    void calculate_photo();
    void calculate_general();

    Eigen::Matrix<std::complex<T>, -1, -1> Gamma0contract();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma1contractAandR();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma2contractAandR();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3Contract_RRandAA();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3Contract_RRandAAblocks();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3Contract_RA();
	
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma1shgcontractAandR();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma2shgcontractAandR();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3shgContract_RA();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3shgContract_RR();
    Eigen::Matrix<std::complex<T>, -1, -1> Gamma3shgContract_AA();

};

