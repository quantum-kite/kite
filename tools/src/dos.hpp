/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/


template <typename T, unsigned DIM>
class dos{
	H5::H5File file;
	public:


        bool isRequired; // was this quantity density of states asked for?
        bool isPossible; // do we have all we need to calculate the density of states?


            // Functions to calculate. They will require the objects present in
        // the configuration file
        int NumMoments;
        int NumPoints;
        std::string filename;

        int NEnergies;

        // information about the Hamiltonian
        system_info<T, DIM> *systemInfo;

        // Input from the shell to override the configuration file
        shell_input variables;

        // Objects required to successfully calculate the conductivity
        Eigen::Array<std::complex<T>, -1, -1> MU;

        std::string dirName;


        dos(system_info<T, DIM>&, shell_input &);
        dos();
        void fetch_parameters();
        void override_parameters();
        void calculate();
	
};


