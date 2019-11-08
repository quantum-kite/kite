

/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T, unsigned DIM>
class arpes{
	H5::H5File file;
	public:


    // Flags that will be accessed from the other classes
    bool isRequired;                        // was this quantity (ARPES) asked for?
    bool isPossible;                        // do we have all we need to calculate ARPES?
    bool default_filename;

    bool calculate_full_arpes;

    // Variables needed to compute the Local density of states
    int NumMoments;                                    // Number of Chebyshev moments
    bool default_NumMoments;
    unsigned NumVectors;                                  // Number of lattice sites in which to compute the LDoS

    int NumEnergies;                                   // Number of energies
    double Emin, Emax;
    bool default_energies;

    double freq;
    bool default_freq;

    Eigen::Array<double, -1, -1> incident_vector;
    bool default_incident;

    double fermi;
    bool default_fermi;

    double temperature;
    double beta;
    bool default_temp;

    double scale;
    double shift;

    std::string kernel;
    double kernel_parameter;
    bool default_kernel;
    bool default_kernel_parameter;

    std::string filename;                          // Saving results to file with this name
    Eigen::Matrix<double, -1, -1> arpes_k_vectors;     // Position of the lattice sites
    Eigen::Matrix<float, -1, -1> energies;                  // Energies specified to be calculated
    std::string dirName;                                    // Name of the hdf5 dataset where the arpes moments are saved

    // Aditional variables
    system_info<T, DIM> *systemInfo;            // information about the Hamiltonian
    shell_input variables;                      // Input from the shell to override the configuration file

    // Objects required to successfully calculate the conductivity
    Eigen::Matrix<std::complex<T>, -1, -1> kMU;

    std::string name;

    // Class methods
    arpes(system_info<T, DIM>&, shell_input &);  // Constructor
	  bool fetch_parameters();                    // Get the parameters from the hdf5 file
    bool is_required();
    void printARPES();
    void set_default_parameters();
	  void override_parameters();                 // If shell variables were given, this function overrides the current parameters
    void calculate();                           // Compute the local density of states
	
};

