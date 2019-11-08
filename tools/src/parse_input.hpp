
class shell_input{
    // Set of parameters that may be passed through the shell

    public:
        // DC conductivity
        double CondDC_Temp; 
        int CondDC_NumEnergies; 
        int CondDC_NumMoments; 
        int CondDC_integrate;
        int CondDC_nthreads;
        double CondDC_Scat; 
        double CondDC_deltaScat;
        double CondDC_FermiMin; 
        double CondDC_FermiMax; 
        int CondDC_NumFermi; 
        std::string CondDC_Name;
        bool CondDC_Exclusive;
        bool CondDC_is_required;

        // Optical Conductivity
        double CondOpt_Temp; 
        int CondOpt_NumEnergies; 
        int CondOpt_NumMoments;
        int CondOpt_Convergence_D;
        int CondOpt_Convergence_G;
        double CondOpt_Scat; 
        double CondOpt_Fermi; 
        double CondOpt_FreqMin; 
        double CondOpt_FreqMax; 
        int CondOpt_NumFreq; 
        std::string CondOpt_Name;
        bool CondOpt_Exclusive;
        bool CondOpt_is_required;


        // Density of states
        int DOS_NumEnergies;
        int DOS_NumMoments;
        std::string DOS_kernel;
        double DOS_kernel_parameter;
        std::string DOS_Name;
        bool DOS_Exclusive;
        bool DOS_is_required;

        // Local density of states
        std::string lDOS_Name;
        bool lDOS_Exclusive;
        bool lDOS_is_required;
        int lDOS_NumMoments;
        std::string lDOS_kernel;
        double lDOS_kernel_parameter;

        // ARPES
        std::string ARPES_Name;
        bool ARPES_Exclusive;
        bool ARPES_is_required;
        bool ARPES_calculate_full_arpes;
        int ARPES_NumMoments;
        std::string ARPES_kernel;
        double ARPES_kernel_parameter;
        double ARPES_Emin;
        double ARPES_Emax;
        double ARPES_NumEnergies;
        double ARPES_Temp;
        double ARPES_Fermi;
        double ARPES_freq;
        Eigen::Array<double, -1, 1> ARPES_vec;


        // 2nd order optical conductivity
        double CondOpt2_Temp; 
        int CondOpt2_NumEnergies; 
        double CondOpt2_Scat; 
        double CondOpt2_Fermi; 
        double CondOpt2_FreqMin; 
        double CondOpt2_FreqMax; 
        int CondOpt2_NumFreq; 
        std::string CondOpt2_Name;
        bool CondOpt2_Exclusive;
        bool CondOpt2_is_required;
        double CondOpt2_ratio;
        int CondOpt2_print_all;

        // Help menu
        bool help;

        // function names
        std::vector<std::string> valid_keys;
        int len;
        std::vector<int> keys_pos;
        std::vector<int> keys_len;

        shell_input(int, char**);
        shell_input();
        void printHelp();
        void printInfo();
        void printDC();
        void printOpt();
        void printOpt2();
        void printDOS();
        void printlDOS();
        void printARPES();

        void parse_CondDC(int argc, char *argv[]);
        void parse_CondOpt(int argc, char *argv[]);
        void parse_CondOpt2(int argc, char *argv[]);
        void parse_DOS(int argc, char *argv[]);
        void parse_lDOS(int argc, char *argv[]);
        void parse_ARPES(int argc, char *argv[]);
        int get_num_exclusives();

};
