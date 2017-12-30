#include "headers.hpp"

using std::cout;
using namespace H5;
using std::endl;



int custom_find(int *arr, int arr_size, int value_to_search){
	/* This function searches array 'arr' for any occurence of number 'value to search'.
	 * If it exists, it exits and returns the index where it occurs.
	 * Otherwise, returns -1
	 */
	for(int i = 0; i < arr_size; i++){
		if(*(arr+i) == value_to_search)
			return i;
	}
	
	return -1;	
}


template <typename T, unsigned DIM>
class info{
	H5::H5File file;
	public:
		int dim;
		int num_orbitals;
		Eigen::Array<int,1,-1> size;
		Eigen::Array<double,-1,-1> vectors;
		Eigen::Array<double, -1, -1> orbital_positions;
		Eigen::Array<int,1,-1> num_moments;
		Eigen::Array<int, 1, -1> Quantities;
			
		double unit_cell_area;
		double spin_degeneracy;
		double energy_scale;
		
		int DOS;
		int CondXX;
		int CondXY;
		int OptCondXX;
		int OptCondXY;
		int SpinCond;
		
		bool MU_needed;
		bool GammaXX_needed;
		bool GammaXY_needed;
		bool LambdaXX_needed;
		bool LambdaXY_needed;

		Eigen::Array<std::complex<T>, -1, -1> LambdaXX;
		Eigen::Array<std::complex<T>, -1, -1> LambdaXY;
		Eigen::Array<std::complex<T>, -1, -1> GammaXX;
		Eigen::Array<std::complex<T>, -1, -1> GammaXY;
		Eigen::Array<std::complex<T>, -1, -1> MU;
		
		info(std::string);
		void read();
	
};

template <typename T, unsigned DIM>
info<T, DIM>::info(std::string name){
	file = H5::H5File(name, H5F_ACC_RDONLY);
	
}
	

template <typename T, unsigned DIM>
void info<T, DIM>::read(){
	debug_message("Entered info::read.\n");
	/* This function reads all the data from the hdf5 file that's needed to 
	 * calculate the quantities we want.*/
	 
	
	// Basic information about the lattice 
	verbose_message("Reading basic information about the lattice: Dimension DIM, Length L and primitive lattice vectors LattVectors\n");
	dim = DIM;										// two-dimensional or three-dimensional
	size = Eigen::Array<int,1,-1>::Zero(1,dim);		// size of the sample
	get_hdf5(size.data(), &file, (char*)"L");
	
	vectors = Eigen::Array<double,-1,-1>::Zero(dim,dim);	// Basis of primitive vectors that generate the lattice
	get_hdf5(vectors.data(), &file, (char*)"LattVectors");
	unit_cell_area = fabs(vectors.matrix().determinant());	// Use the basis vectors to determine the area of the unit cell
	
	verbose_message("Reading the energy scale, number of orbitals NOrbitals and their positions OrbPositions\n");
	get_hdf5(&energy_scale, &file, (char*)"EnergyScale");								// energy scale
	get_hdf5(&num_orbitals, &file, (char*)"NOrbitals");									// number of orbitals in each unit cell	
	orbital_positions = Eigen::Array<double,-1,-1>::Zero(num_orbitals, num_orbitals);	// position of each of those orbitals
	get_hdf5(orbital_positions.data(), &file, (char*)"OrbPositions");
	
	spin_degeneracy = 2; // put by hand?
	
	
	// Information about the data types
	verbose_message("Reading data type and checking whether it is complex.\n");
	int precision = 1, complex;
	get_hdf5(&complex, &file, (char *) "/IS_COMPLEX");
	get_hdf5(&precision,  &file, (char *) "/PRECISION");
	
	
	// List of quantities to calculate	
	verbose_message("Reading the list of quantities to calculate, Calculation/FunctionNum.\n");
    H5::DataSet * dataset     = new H5::DataSet(file.openDataSet("/Calculation/FunctionNum"));
    H5::DataSpace * dataspace = new H5::DataSpace(dataset->getSpace());
    size_t NQuantities        = dataspace->getSimpleExtentNpoints();
    Quantities  = Eigen::Array<int, 1, -1>::Zero(1, NQuantities);
    num_moments = Eigen::Array<int, 1, -1>::Zero(1, NQuantities);
    free(dataset);
    free(dataspace);
    get_hdf5<int>(Quantities.data(), &file, (char *)   "/Calculation/FunctionNum");
	get_hdf5(num_moments.data(), &file, (char*)"Calculation/NumMoments");
	
	
	// Check which quantities need to be calculated. If they are not needed, the corresponding number will be -1.
	verbose_message("Searching which quantities need to be calculated.\n");
	DOS = custom_find(Quantities.data(), NQuantities, 1);
	CondXX = custom_find(Quantities.data(), NQuantities, 2);
	CondXY = custom_find(Quantities.data(), NQuantities, 3);
	OptCondXX = custom_find(Quantities.data(), NQuantities, 4); //?? needs to be clarified
	OptCondXY = custom_find(Quantities.data(), NQuantities, 4); //?? needs to be clarified
	//SpinCond = custom_find(Quantities.data(), NQuantities, 5);
	
	// opt cond needs to be clarified
	MU_needed			= DOS >= 0;
	GammaXX_needed 		= CondXX >= 0 or OptCondXX >= 0;
	GammaXY_needed 		= CondXY >= 0 or OptCondXY >= 0;
	LambdaXX_needed 	= CondXX >= 0;
	LambdaXY_needed 	= CondXY >= 0;
	
	
	verbose_message("MU needed? "); verbose_message(MU_needed);
	verbose_message("\nGammaXX needed? "); verbose_message(GammaXX_needed);
	verbose_message("\nGammaXY needed? "); verbose_message(GammaXY_needed);
	verbose_message("\nLambdaXX needed? "); verbose_message(LambdaXX_needed);
	verbose_message("\nLambdaXY needed? "); verbose_message(LambdaXY_needed);
	verbose_message("\n");


	// Fill the MU matrix if it is needed. 
	// For simplicity, let it always be complex	
	if(MU_needed){
		verbose_message("Filling the MU matrix.\n");
		MU = Eigen::Array<std::complex<T>,-1,-1>::Zero(1,num_moments(DOS));
		
		if(complex)
			get_hdf5(MU.data(), &file, (char*)"MU");
		
		if(!complex){
			Eigen::Array<T,-1,-1> MUReal;
			MUReal = Eigen::Array<T,-1,-1>::Zero(1,num_moments(DOS));
			get_hdf5(MUReal.data(), &file, (char*)"MU");
			
			MU = MUReal.template cast<std::complex<T>>();
		}		
	}
	
	// Fill the GammaXX matrix if it is needed. 
	// For simplicity, let it always be complex	
	if(GammaXX_needed){
		verbose_message("Filling the GammaXX matrix.\n");
		GammaXX = Eigen::Array<std::complex<T>,-1,-1>::Zero(1,num_moments(CondXX)*num_moments(CondXX));
		
		if(complex)
			get_hdf5(GammaXX.data(), &file, (char*)"GammaXX");
		
		if(!complex){
			Eigen::Array<T,-1,-1> GammaXXReal;
			GammaXXReal = Eigen::Array<T,-1,-1>::Zero(1,num_moments(CondXX)*num_moments(CondXX));
			get_hdf5(GammaXXReal.data(), &file, (char*)"GammaXX");
			
			GammaXX = GammaXXReal.template cast<std::complex<T>>();
		}				
	}
	
	// Fill the GammaXY matrix if it is needed. 
	// For simplicity, let it always be complex	
	if(GammaXY_needed){
		verbose_message("Filling the GammaXY matrix.\n");
		GammaXY = Eigen::Array<std::complex<T>,-1,-1>::Zero(1,num_moments(CondXY)*num_moments(CondXY));
		
		if(complex)
			get_hdf5(GammaXY.data(), &file, (char*)"GammaXY");
		
		if(!complex){
			Eigen::Array<T,-1,-1> GammaXYReal;
			GammaXYReal = Eigen::Array<T,-1,-1>::Zero(1,num_moments(CondXY)*num_moments(CondXY));
			get_hdf5(GammaXYReal.data(), &file, (char*)"GammaXY");
			
			GammaXY = GammaXYReal.template cast<std::complex<T>>();
		}				
	}
	
	
	// Fill the LambdaXX matrix if it is needed. 
	// For simplicity, let it always be complex	
	if(LambdaXX_needed){
		verbose_message("Filling the LambdaXX matrix.\n");
		LambdaXX = Eigen::Array<std::complex<T>,-1,-1>::Zero(1,num_moments(CondXX));
		
		if(complex)
			get_hdf5(LambdaXX.data(), &file, (char*)"LambdaXX");
		
		if(!complex){
			Eigen::Array<T,-1,-1> LambdaXXReal;
			LambdaXXReal = Eigen::Array<T,-1,-1>::Zero(1,num_moments(CondXX));
			get_hdf5(LambdaXXReal.data(), &file, (char*)"LambdaXX");
			
			LambdaXX = LambdaXXReal.template cast<std::complex<T>>();
		}		
	}
	
	
	// Fill the LambdaXY matrix if it is needed. 
	// For simplicity, let it always be complex	
	if(LambdaXY_needed){
		verbose_message("Filling the LambdaXY matrix.\n");
		LambdaXY = Eigen::Array<std::complex<T>,-1,-1>::Zero(1,num_moments(CondXY));
		
		if(complex)
			get_hdf5(LambdaXY.data(), &file, (char*)"LambdaXY");
		
		if(!complex){
			Eigen::Array<T,-1,-1> LambdaXYReal;
			LambdaXYReal = Eigen::Array<T,-1,-1>::Zero(1,num_moments(CondXY));
			get_hdf5(LambdaXYReal.data(), &file, (char*)"LambdaXY");
			
			LambdaXY = LambdaXYReal.template cast<std::complex<T>>();
		}				
	}
	
	file.close();
	debug_message("Left info::read.\n");
}

