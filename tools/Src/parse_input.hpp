/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

#include <string>
#include <vector>

class command{
  public:

    std::vector<std::string> keys;
    std::vector<std::string> values;
    int num_commands;
    std::string func_to_calculate;

    // These are the valid command names
    std::string help_key = "--help";
    std::vector<std::string> valid_keys {"DOS","CondOpt","CondDC"};
    int num_keys = valid_keys.size();


    // Simulation parameters
    double fermi_energy;  
    double         temp;   //temperature
    int      N_energies;
    double     max_freq;
    int         N_freqs;
    double   scattering;

    bool override_config = true;
    bool use_defaults = false;


    // Class methods
    command(std::vector<std::string>, std::vector<std::string>);
    int check_function();
    void print_help();
    std::string find_keys(std::string);
    void fetch_parameters();

    void ifCondOpt();
    void ifDOS();


};

command::command(std::vector<std::string> input_keys, std::vector<std::string> input_values){
  // Class constructor

  // check if the sizes of the input keys and values are the same
  if(input_values.size() != input_keys.size()){
    std::cout << "There was an error processing the input to the program. The length of \
      the keys and values obtained from parsing the input to the program is not the same. \
      Please check that your input is correct. Exiting program. \n";
    exit(1);
  }

  // initialize the class's keys and values to the programs's
  num_commands = input_values.size();

  keys.resize(num_commands);
  values.resize(num_commands);

  keys = input_keys;
  values = input_values;
}

void command::print_help(){
  // This is the text that appears when the user uses --help as a parameter
  std::cout << "Usage:\n";
  std::cout << "--calc \t\tspecify the function to calculate.\n";
  
  // exit the program successfully
  exit(0);
};


int command::check_function(){
  // check if the function being asked to calculate is one of the valid ones
  // or if the user used --help
  // If it's one of the valid ones, then that's the function the user wants to calculate
  // If there's more than one, the user doesn't know what they're doing and the
  // program should exit.
  verbose_message("Checking function to evaluate.\n");


  int num_matches = 0;
  bool need_help = false;
  for(int i = 0; i < num_commands; i++){
    for(int j = 0; j < num_keys; j++)
      if(values.at(i) == valid_keys.at(j)){
        func_to_calculate = values.at(i);
        num_matches++;
      }
    if(keys.at(i) == help_key)
      need_help = true;
  }


  verbose_message("Number of matches: "); verbose_message(num_matches);
  verbose_message("\nUsed --help? "); verbose_message(need_help);
  verbose_message("\nNumber of commands: ");verbose_message(num_commands);
  verbose_message("\n");


  if(need_help and num_commands == 1)
    print_help();

  else if(need_help and num_commands > 1){
    std::cout << "You requested a function to calculate but also the help menu. The";
    std::cout << " program will print the help menu but will not calculate the function";
    std::cout << "requested. To calculate the function, please do not request the help menu.\n";
    print_help();
  }

  else if(num_matches > 1){
    std::cout << "Please enter only one function to calculate. Leaving program.\n";
    exit(1);
  }

  else if(num_matches == 0){
    std::cout << "No valid function to calculate was requested. Leaving program.\n";
    exit(1);
  }

  verbose_message("The function requested to calculate was ");
  verbose_message(func_to_calculate);verbose_message("\n");

  return num_matches;

}

void command::fetch_parameters(){
  // After we know which function to calculate, this method will try to obtian the
  // required parameters to calculate it from the shell input
  verbose_message("Entered fetch_parameters\n");

  if(func_to_calculate == "CondOpt") ifCondOpt();
  if(func_to_calculate == "DOS") ifDOS();

  verbose_message("Left fetch_parameters\n.");

}


void command::ifDOS(){
  // In case the user requested DOS be calculated, some parameters must be present.
  // This function checks for those parameters. Parameters:
  // -NE       Number of energy points to use in the integration
  verbose_message("Entered ifDOS\n");

  // check if indeed the function to calculate is the correct one
  if(func_to_calculate != "DOS"){
    std::cout << "Bad usage of the ifDOS function. ";
    std::cout << "This function is checking for the parameters of 'DOS' but ";
    std::cout << "func_to_calculate is not 'DOS'. Leaving program.\n";
    exit(1);
  }

  N_energies   = stoi(find_keys("-NE"));
}

void command::ifCondOpt(){
  // In case the user requested CondOpt be calculated, some parameters must be present.
  // This function checks for those parameters. Parameters:
  // -EF       Fermi energy (eV)
  // -T        Temperature  (K)
  // -NE       Number of energy points to use in the integration
  // -MF       Maximum frequency
  // -NF       Number of frequency points to calculate
  // -G        Finite scattering parameter (eV)
  verbose_message("Entered ifCondOpt\n");

  // check if indeed the function to calculate is the correct one
  if(func_to_calculate != "CondOpt"){
    std::cout << "Bad usage of the ifCondOpt function. ";
    std::cout << "This function is checking for the parameters of 'CondOpt' but ";
    std::cout << "func_to_calculate is not 'CondOpt'. Leaving program.\n";
    exit(1);
  }
  
  

  fermi_energy = stod(find_keys("-EF"));
  temp         = stod(find_keys("-T"));
  N_energies   = stoi(find_keys("-NE"));
  max_freq     = stod(find_keys("-MF"));
  N_freqs      = stoi(find_keys("-NF"));
  scattering   = stod(find_keys("-G"));
}

std::string command::find_keys(std::string name){
  // Tries to find 'name' among the parameters passed on to the program
  verbose_message("Entered find_keys.\n");



  // Try to find 'name'
  std::string val;
  int num_matches = 0;
  for(int i = 0; i < num_commands; i++){
    if(name == keys.at(i)){
      val = values.at(i);
      num_matches++;
    }
  }

  verbose_message("Number of matches: "); verbose_message(num_matches);
  // Check if there's only one match. Any number different from one means something's wrong
  if(num_matches == 0){
    std::cout << "Unable to find " << name << " in the list of parameters for the ";
    std::cout << "function " << func_to_calculate << ". Leaving program.\n";
    exit(1);
  }
  else if(num_matches > 1){
    std::cout << "More than one occurence of " << name << " was found. Please provide";
    std::cout << " only one of each parameter.\n";
    exit(1);
  }
  else{
    verbose_message("The value found was "); verbose_message(val);verbose_message("\n");
    verbose_message("Left find_keys.\n");
    return val;
  }

  verbose_message("Left find_keys after return.\n");
}

void parser(int argc, char *argv[]){
	// Processes the input that this program recieves from the command line	
	
	std::string *arguments; // Convert it all to strings
	arguments = new std::string[argc];
	for(int i=0; i<argc; i++){
		arguments[i] = argv[i];	
		std::cout << arguments[i] << std::endl;
	}

  std::cout << "Finished converting to strings\n" << std::flush;

  // Determine the number of commands (the number of - )
  int num_commands = 0;
	for(int cursor = 1; cursor < argc; cursor++)
    if(arguments[cursor][0]=='-')
      num_commands++;  






  std::vector<std::string> keys(num_commands);
  std::vector<std::string> values(num_commands);
  //std::string *keys;
  //std::string *values;

  //keys = new std::string[num_commands];
  //values = new std::string[num_commands];

  // organize the entries as a series of key-value pairs
  // the last entry is the name of the configuration file
  // and the first entry is the name of the executable file
  bool prev_is_command = 0;
  int n = -1;
	for(int cursor = 1; cursor < argc-1; cursor++){
    std::cout << arguments[cursor] << "\n" <<std::flush;
    if(arguments[cursor][0]=='-' and !prev_is_command){
      std::cout << "first if";
      prev_is_command = 1;
      n++;
      keys.at(n) = arguments[cursor];
    }

    else if(arguments[cursor][0]=='-' and prev_is_command){
      
      std::cout << "second if";
      prev_is_command = 1;
      values.at(n) = "";
      n++;
      keys.at(n) = arguments[cursor];
    }

    
    else if(arguments[cursor][0]!='-' and prev_is_command){
      std::cout << "third if";
      prev_is_command = 0;
      values.at(n) = arguments[cursor];
    }
    

    else if(arguments[cursor][0]!='-' and !prev_is_command){
      std::cout << "fourth if";

      std::cout << "Wrong string format. Aborting.\n";
      exit(0);
    }
  }

  std::cout << "before print\n" << std::flush;
  for(int i = 0; i < num_commands; i++)
    std::cout << keys.at(i) << " " << values.at(i) << "\n";


  command comm(keys, values);
  comm.check_function();
  comm.fetch_parameters();
  



  // the number of matches has to be exactly one
  //std::cout << "number of matches: " << num_matches << "\n";
  //std::cout << "need help? " << need_help << "\n";

  //if(need_help)
    //std::cout << "Help Menu: \n to be continued..." << std::endl;



  
  //std::cout << "CSTR" << "\n";    
  
	
  
  
  
  //delete[] keys;
  //delete[] values;
  delete[] arguments;
};


