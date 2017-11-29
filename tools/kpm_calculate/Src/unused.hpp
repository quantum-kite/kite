/*
class info{
	H5std_string FILE_NAME;
	public:
		info(std::string);
		template <typename T> void read(std::string, T);
};

info::info(std::string name){
	FILE_NAME = name;

}

template <typename T>
H5::DataType convert(T irrelevant){
	H5T_class_t type_class = dataset.getTypeClass();
	int size;
	switch(type_class){
		case H5T_INTEGER:
			size = dataset.getIntType().getSize();
			std::cout << "size:" << size << "\n";
			return PredType::NATIVE_INT;
		case H5T_FLOAT:
			size = dataset.getFloatType().getSize();
			std::cout << "size:" << size << "\n";
			
			switch(size){
				case(4):
					return PredType::NATIVE_FLOAT;
				case(8):
					return PredType::NATIVE_DOUBLE;
				default:
					std::cout << "Unknown Data Size. Leaving program.\n";
					exit(0);
				}
		
		default :
			std::cout << "Unknown DataType. Leaving program.\n";
			exit(0);
	}
}

template <typename T>
void info::read(std::string NAME_DATASET_STR, T data_buffer){
	std::cout << "entered\n";
				
	const H5std_string DATASET_NAME1(NAME_DATASET_STR);
	H5File file("SDS.h5", H5F_ACC_RDONLY);
	DataSet dataset = file.openDataSet(DATASET_NAME1);
	H5T_class_t type_class = dataset.getTypeClass();
	DataSpace dataspace = dataset.getSpace();
	
	int rank = dataspace.getSimpleExtentNdims();
	
	
			
	hsize_t dims_out[2];
	int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
	cout << "rank " << rank << ", dimensions " <<
			(unsigned long)(dims_out[0]) << " x " <<
			(unsigned long)(dims_out[1]) << endl; 
			
	std::cout << "EFE\n";
	
	dataset.read(data_buffer, convert(dataset));
	std::cout << "Left\n";
}



template <typename T> void h5_fetch(H5std_string FILE_NAME, H5std_string NAME_DATASET, T *data){
	H5::CompType complex_datatype(sizeof(std::complex<double>));
  complex_datatype.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
  complex_datatype.insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);
  
  H5::H5File file(FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
  
  H5::DataSet dataset;
	dataset = file.openDataSet(NAME_DATASET);
	dataset.read(data->data(), complex_datatype);	
	dataset.close();
	file.close();
}

void single_shot(Eigen::Array<double, -1, 1> energies){
	std::cout << "#########3\n\n\nEntered single_shot\n";fflush(stdout);
	std::cout << energies << "\n\n";
	// Read the data from the HDF5 file
	info config("test_f.h5");
	config.read();
	int n_moments = config.num_moments(0);
	
	std::cout << "here1\n";fflush(stdout);
	double finite_lambda = 0.01;
	std::complex<double> energy; 

	std::cout << "here2\n";fflush(stdout);

	Eigen::Array<double, -1, -1> cond_DC;
	cond_DC = Eigen::Array<double, -1, -1>::Zero(energies.rows(),energies.cols());
	for(int e=0; e < energies.rows()*energies.cols(); e++){
		energy = std::complex<double>(energies(e), finite_lambda);
		std::cout << energies(e) << "\n";
		std::cout << energy << "\n";
		for(int i=0; i < n_moments; i++){
			for(int j=0; j < n_moments; j++){
				cond_DC(e) += config.GammaXX(i*n_moments + j).real()*green(i, 1, energy).imag()*green(j, 1, energy).imag();
			}
		}
	}
	
	std::ofstream myfile;
  myfile.open ("output/single_shot.dat");
	for(int i=0; i < energies.rows()*energies.cols(); i++)
		myfile  << energies(i) << " " << cond_DC(i) << "\n";
	
	myfile.close();
	
	std::cout << "Left single_shot\n";fflush(stdout);
}
 
 
*/


