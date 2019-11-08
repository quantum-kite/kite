

/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T, unsigned DIM>
class ldos{
	H5::H5File file;
	public:


    // Flags that will be accessed from the other classes
    bool isRequired;                        // was this quantity local density of states asked for?
    bool isPossible;                        // do we have all we need to calculate the density of states?
    bool default_filename;

    // Variables needed to compute the Local density of states
    int NumMoments;                                    // Number of Chebyshev moments
    int MaxMoments;
    bool default_NumMoments;
    unsigned NumPositions;                                  // Number of lattice sites in which to compute the LDoS
    int NumEnergies;                                   // Number of energies
    std::string filename;                          // Saving results to file with this name
    Eigen::Matrix<unsigned long, -1, -1> global_positions;
    Eigen::Matrix<unsigned long, -1, -1> ldos_Orbitals;     // Position of the lattice sites
    Eigen::Matrix<unsigned long, -1, -1> ldos_Positions;     // Position of the lattice sites
    Eigen::Matrix<float, -1, -1> energies;                  // Energies specified to be calculated
    std::string dirName;                                    // Name of the hdf5 dataset where the ldos moments are saved

    std::string kernel;
    double kernel_parameter;
    bool default_kernel;
    bool default_kernel_parameter;

    // Aditional variables
    system_info<T, DIM> *systemInfo;            // information about the Hamiltonian
    shell_input variables;                      // Input from the shell to override the configuration file

    // Objects required to successfully calculate the conductivity
    Eigen::Matrix<std::complex<T>, -1, -1> lMU;

    std::string name;

    // Class methods
    ldos(system_info<T, DIM>&, shell_input &);  // Constructor
    bool fetch_parameters();                    // Get the parameters from the hdf5 file
    bool is_required();
    void printLDOS();
    void set_default_parameters();
    void override_parameters();                 // If shell variables were given, this function overrides the current parameters
    void calculate();                           // Compute the local density of states
	
};

