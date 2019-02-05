/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/


template <typename T, unsigned DIM>
class dos{
	public:

        system_info<T, DIM> *systemInfo;    // information about the Hamiltonian
        shell_input variables;              // Input from the shell to override the configuration file

        bool isRequired; // was this quantity density of states asked for?
        bool isPossible; // do we have all we need to calculate the density of states?

        double Emax, Emin;
        bool default_Emax, default_Emin;

            // Functions to calculate. They will require the objects present in
        // the configuration file
        int NumMoments;
        int MaxMoments;
        bool default_NumMoments;

        std::string kernel;
        double kernel_parameter;
        bool default_kernel;
        bool default_kernel_parameter;
        
        std::string filename;
        bool default_filename;

        bool default_NEnergies;
        int NEnergies;
        Eigen::Matrix<T, -1, 1> energies;


        Eigen::Array<std::complex<T>, -1, -1> MU;   // Objects required to successfully calculate the conductivity

        Eigen::Array<std::complex<T>, -1, -1> GammaE;
        bool dos_finished;

        dos(system_info<T, DIM>&, shell_input &);
        dos();
        bool is_required();
        void printDOS();
        void find_limits();
        void set_default_parameters();
        bool fetch_parameters();
        void override_parameters();
        void calculate();
	
};


