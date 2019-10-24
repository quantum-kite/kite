/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T, unsigned DIM>
class conductivity_dc{
	public:


    bool isRequired; // was this quantity (conductivity_dc) asked for?
    bool isPossible; // do we have all we need to calculate the conductivity?


    double temperature;
    T beta;
    bool default_temp;

		// Functions to calculate. They will require the objects present in
    // the configuration file
    int direction;
    int NumDisorder;
    int NumMoments;
    int MaxMoments;
    int Moments_D;
    int Moments_G;
    bool default_NumMoments;
    int full_range;
    bool default_full_range;
    double units;

    int NumThreads;
    bool default_NumThreads;
    
    std::string filename;
    bool default_filename;

    T scat, deltascat;
    bool default_scat;
    bool default_deltascat;

    int NEnergies; 
    bool default_NEnergies;

    T minEnergy;
    T maxEnergy;
    bool default_energy_limits;
    Eigen::Matrix<T, -1, 1> energies;

    int NFermiEnergies;
    double minFermiEnergy;
    double maxFermiEnergy;
    Eigen::Matrix<T, -1, 1> fermiEnergies;
    bool default_NFermi;
    bool default_mFermi;
    bool default_MFermi;

    // information about the Hamiltonian
    system_info<T, DIM> systemInfo;

    // Input from the shell to override the configuration file
    shell_input variables;

    // Objects required to successfully calculate the conductivity
    //Eigen::Array<std::complex<T>, -1, -1, Eigen::RowMajor> Gamma;
    Eigen::Array<std::complex<T>, -1, -1> Gamma;
    Eigen::Array<std::complex<T>, -1, -1> Gamma_Padded;




    conductivity_dc(system_info<T, DIM>&, shell_input &);
    void printDC();
    void set_energy_limits();
    bool is_required();
    void set_default_parameters();
	  bool fetch_parameters();
	  void override_parameters();
    void calculate();
    void calculate2();
    void calculate_imag();
    Eigen::Matrix<std::complex<T>, -1, -1, Eigen::ColMajor> fill_delta();
    Eigen::Matrix<std::complex<T>, -1, -1, Eigen::RowMajor> fill_dgreenR();
    Eigen::Matrix<std::complex<T>, -1, -1> triple_product(
      Eigen::Matrix<std::complex<T>, -1, -1, Eigen::ColMajor>,
      Eigen::Matrix<std::complex<T>, -1, -1, Eigen::RowMajor>);

    Eigen::Matrix<std::complex<T>, -1, 1> calc_cond(
        Eigen::Matrix<std::complex<T>, -1, -1>);

    void save_to_file(Eigen::Matrix<std::complex<T>, -1, -1>);


};
