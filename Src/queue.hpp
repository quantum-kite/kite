
std::string num2str3(int dir_num){
  std::string dir;
 
  switch(dir_num){
    case 0:
      dir = "xxx"; break;
    case 1:
      dir = "xxy"; break;
    case 2:
      dir = "xxz"; break;
    case 3:
      dir = "xyx"; break;
    case 4:
      dir = "xyy"; break;
    case 5:
      dir = "xyz"; break;
    case 6:
      dir = "xzx"; break;
    case 7:
      dir = "xzy"; break;
    case 8:
      dir = "xzz"; break;
    case 9:
      dir = "yxx"; break;
    case 10:
      dir = "yxy"; break;
    case 11:
      dir = "yxz"; break;
    case 12:
      dir = "yyx"; break;
    case 13:
      dir = "yyy"; break;
    case 14:
      dir = "yyz"; break;
    case 15:
      dir = "yzx"; break;
    case 16:
      dir = "yzy"; break;
    case 17:
      dir = "yzz"; break;
    case 18:
      dir = "zxx"; break;
    case 19:
      dir = "zxy"; break;
    case 20:
      dir = "zxz"; break;
    case 21:
      dir = "zyx"; break;
    case 22:
      dir = "zyy"; break;
    case 23:
      dir = "zyz"; break;
    case 24:
      dir = "zzx"; break;
    case 25:
      dir = "zzy"; break;
    case 26:
      dir = "zzz"; break;
    default:
      std::cout << "Invalid direction in num2str_dir3.\n"; exit(1);
  }
  return dir;
}

std::string num2str2(int dir_num){
  std::string dir;
 
  switch(dir_num){
    case 0:
      dir = "xx"; break;
    case 1:
      dir = "xy"; break;
    case 2:
      dir = "xz"; break;
    case 3:
      dir = "yx"; break;
    case 4:
      dir = "yy"; break;
    case 5:
      dir = "yz"; break;
    case 6:
      dir = "zx"; break;
    case 7:
      dir = "zy"; break;
    case 8:
      dir = "zz"; break;
    default:
      std::cout << "Invalid direction for the optical conductivity.\n"; exit(1);
  }
  return dir;
}


class measurement_queue{
  public:
    std::string direction_string;
    std::vector<int> NMoments;
    int NDisorder;
    int NRandom;
    std::string label;
    measurement_queue(std::string dir_string, std::vector<int> moments, int disorder, int random, std::string name){
      direction_string = dir_string;
      NMoments = moments;
      NDisorder = disorder;
      NRandom = random;
      label = name;
    };

    double time_required(double avg_duration){
      int prod = 1;
      for(unsigned int i = 0; i < NMoments.size(); i++)
        prod *= NMoments.at(i);

      return prod*avg_duration;
    }
};


std::vector<measurement_queue> fill_queue(char *name){
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
    std::vector<measurement_queue> queue;

    // Now we must check which functions are asked by the configuration file.
    // Checking for the density of states
    try{
      //H5::Exception::dontPrint();
      int NMoments, NRandom, NDisorder;
       debug_message("DOS: checking if we need to calculate DOS.\n");
       get_hdf5<int>(&NMoments,  file, (char *)   "/Calculation/dos/NumMoments");
       //get_hdf5<int>(&NDisorder, file, (char *)   "/Calculation/dos/NumDisorder");
       NDisorder = 1;
       get_hdf5<int>(&NRandom,   file, (char *)   "/Calculation/dos/NumRandoms");
       //dos = true;
       queue.push_back(measurement_queue("", {NMoments}, NDisorder, NRandom, "/Calculation/dos/MU"));
    } catch(H5::Exception& e) {debug_message("DOS: no need to calculate DOS.\n");}

    // Checking for the optical conductivity
    try{
      int direction, NMoments, NRandom, NDisorder;
       debug_message("Optical conductivity: checking if we need to calculate it.\n");
       get_hdf5<int>(&direction, file, (char *) "/Calculation/conductivity_optical/Direction");
       get_hdf5<int>(&NMoments, file, (char *)  "/Calculation/conductivity_optical/NumMoments");
       get_hdf5<int>(&NRandom, file, (char *)   "/Calculation/conductivity_optical/NumRandoms");
       //get_hdf5<int>(&NDisorder, file, (char *)   "/Calculation/conductivity_optical/NumDisorder");
       NDisorder = 1;
       //conductivity_optical = true;
      
       // convert the numerical value for the direction into the string that represents it
       std::string dir(num2str2(direction));
       // same string, but separated by commas. This indicates a different gamma function
       std::string dirc = dir.substr(0,1)+","+dir.substr(1,2);
      
       queue.push_back(measurement_queue(dir,  {NMoments},           NDisorder, NRandom, "/Calculation/conductivity_optical/Lambda"+dir));
       queue.push_back(measurement_queue(dirc, {NMoments, NMoments}, NDisorder, NRandom, "/Calculation/conductivity_optical/Gamma" +dir));
    } catch(H5::Exception& e) {debug_message("Optical conductivity: no need to calculate it.\n");}


    // Checking for the dc conductivity
    try{
      int direction, NMoments, NRandom, NDisorder;
       debug_message("dc conductivity: checking if we need to calculate it.\n");
       get_hdf5<int>(&direction, file, (char *) "/Calculation/conductivity_dc/Direction");
       get_hdf5<int>(&NMoments, file, (char *)  "/Calculation/conductivity_dc/NumMoments");
       get_hdf5<int>(&NRandom, file, (char *)   "/Calculation/conductivity_dc/NumRandoms");
       //get_hdf5<int>(&NDisorder, file, (char *) "/Calculation/conductivity_dc/NumDisorder");
       NDisorder = 1;
       //conductivity_dc = true;

       // convert the numerical value for the direction into the string that represents it
       std::string dir(num2str2(direction));
       // same string, but separated by commas. This indicates a different gamma function
       std::string dirc = dir.substr(0,1)+","+dir.substr(1,2);
      
       queue.push_back(measurement_queue(dirc, {NMoments, NMoments}, NDisorder, NRandom, "/Calculation/conductivity_dc/Gamma"+dir));
    } catch(H5::Exception& e) {debug_message("dc conductivity: no need to calculate it.\n");}


    // nonlinear optical conductivity    
    try{
      int direction, NMoments, NRandom, NDisorder, special;
       debug_message("nonlinear optical cond: checking if we need to calculate it.\n");
       get_hdf5<int>(&direction, file, (char *)   "/Calculation/conductivity_optical_nonlinear/Direction");
       get_hdf5<int>(&NMoments, file, (char *)  "/Calculation/conductivity_optical_nonlinear/NumMoments");
       get_hdf5<int>(&NRandom, file, (char *)   "/Calculation/conductivity_optical_nonlinear/NumRandoms");
       //get_hdf5<int>(&NDisorder, file, (char *)   "/Calculation/conductivity_optical_nonlinear/NumDisorder");
       get_hdf5<int>(&special, file, (char *)   "/Calculation/conductivity_optical_nonlinear/Special");
       //conductivity_optical_nonlinear = true;
       NDisorder = 1;

       // convert the numerical value for the direction into the string that represents it
       std::string dir(num2str3(direction));                                                // xxx Gamma0
       std::string dirc1 = dir.substr(0,1) + "," + dir.substr(1,3);                         // x,xx Gamma1
       std::string dirc2 = dir.substr(0,2) + "," + dir.substr(2,3);                         // xx,x Gamma2
       std::string dirc3 = dir.substr(0,1) + "," + dir.substr(1,2) + "," + dir.substr(2,3); // x,x,x Gamma3

       std::string directory = "/Calculation/conductivity_optical_nonlinear/";
       
      // regular nonlinear calculation
      if(special != 1){
       queue.push_back(measurement_queue(dir,   {NMoments},                     NDisorder, NRandom, directory+"Gamma0"+dir));
       queue.push_back(measurement_queue(dirc1, {NMoments, NMoments},           NDisorder, NRandom, directory+"Gamma1"+dir));
       queue.push_back(measurement_queue(dirc2, {NMoments, NMoments},           NDisorder, NRandom, directory+"Gamma2"+dir));
       queue.push_back(measurement_queue(dirc3, {NMoments, NMoments, NMoments}, NDisorder, NRandom, directory+"Gamma3"+dir));
      }

      // special nonlinear calculation. In this case, it's going to be HBN, which is nonlinear
      // but only has simple objects that need calculating
      if(special == 1){
       queue.push_back(measurement_queue(dirc1, {NMoments, NMoments}, NDisorder, NRandom, directory+"Gamma1"+dir));
       queue.push_back(measurement_queue(dirc2, {NMoments, NMoments}, NDisorder, NRandom, directory+"Gamma2"+dir));
      }

    } catch(H5::Exception& e) {debug_message("nonlinear optical conductivity: no need to calculate it.\n");}

    delete file;
    

  return queue;

};




class singleshot_measurement_queue{
  public:
    std::string direction_string;
    int NMoments;
    int NDisorder;
    int NRandom;
    Eigen::Array<double, -1, 1> singleshot_energies;
    double single_gamma;
    std::string label;
    double time_length;


    double time_required(double avg_duration){
      return NMoments*NMoments*avg_duration;
    }



    singleshot_measurement_queue(std::string dir_string, int moments, int disorder, 
        int random, std::string name, Eigen::Array<double, -1, 1> energies, double gamma){
      direction_string = dir_string;
      NMoments = moments;
      NDisorder = disorder;
      NRandom = random;
      label = name;
      singleshot_energies = energies;
      single_gamma = gamma;
    };

};

std::vector<singleshot_measurement_queue> fill_singleshot_queue(char *name){
    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDONLY);
    
    Eigen::Array<double, -1, 1> energies;
    int NDisorder, NRandom, NMoments, direction;
    std::string direction_string;
    double gamma;
    
    NDisorder = 1;

    std::vector<singleshot_measurement_queue> queue;
    try{
       debug_message("single_shot dc checking if we need to calculate it.\n");
       get_hdf5<int>(&direction, file, (char *)   "/Calculation/singleshot_conductivity_dc/Direction");
       get_hdf5<int>(&NMoments, file, (char *)  "/Calculation/singleshot_conductivity_dc/NumMoments");
       get_hdf5<int>(&NRandom, file, (char *)   "/Calculation/singleshot_conductivity_dc/NumRandoms");
       get_hdf5<double>(&gamma, file, (char *)   "/Calculation/singleshot_conductivity_dc/Gamma");
       
       if(direction == 0)
         direction_string = "x,x";
       else if(direction == 1)
         direction_string = "y,y";
       else{
         std::cout << "Invalid singleshot direction. Exiting.\n";
         exit(1);
       }

      // We also need to determine the number of energies that we need to calculate
      H5::DataSet * dataset_energy     	= new H5::DataSet(file->openDataSet("/Calculation/singleshot_conductivity_dc/Energy"));
      H5::DataSpace * dataspace_energy 	= new H5::DataSpace(dataset_energy->getSpace());
      hsize_t dims_out[2];		
      dataspace_energy->getSimpleExtentDims(dims_out, NULL);	
      energies = Eigen::Array<double, -1, 1>::Zero(dims_out[0]*dims_out[1], 1);	
      delete dataspace_energy;
      delete dataset_energy;
      
      get_hdf5<double>(energies.data(),  	file, (char *)   "/Calculation/singleshot_conductivity_dc/Energy");
      
      queue.push_back(singleshot_measurement_queue(direction_string, NMoments,
            NDisorder, NRandom, "/Calculation/singleshot_conductivity_dc/SingleShot",
            energies, gamma));

    } catch(H5::Exception& e) {debug_message("singleshot dc: no need to calculate it.\n");}

    delete file;
    
  return queue;
}

std::string print_time(double duration){
  if(duration < 500)
    return std::to_string(duration) + " seconds.";
  else if(duration >= 500 and duration < 60*3*60)
    return std::to_string(int(duration/60)) + " minutes.";
  else if(duration >= 60*60*3 and duration < 60*60*50)
    return std::to_string(int(duration/(60*60))) + " hours.";
  else
    return std::to_string(int(duration/(60*60*24))) + " days.";

}
