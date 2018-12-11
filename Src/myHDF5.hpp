/****************************************************************/
/*                                                              */
/*  Copyright (C) 2018, M. Andelkovic, L. Covaci, A. Ferreira,  */
/*                    S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                              */
/****************************************************************/

template <typename T>
class DataTypeFor
{
public:
  static  H5::DataType value;
};
template<>
H5::DataType DataTypeFor<int>::value = H5::PredType::NATIVE_INT;
template<>
H5::DataType DataTypeFor<unsigned int>::value = H5::PredType::NATIVE_UINT;
template<>
H5::DataType DataTypeFor<float>::value = H5::PredType::NATIVE_FLOAT;
template<>
H5::DataType DataTypeFor<double>::value = H5::PredType::NATIVE_DOUBLE;
template<>
H5::DataType DataTypeFor<long double>::value = H5::PredType::NATIVE_LDOUBLE;




 
template <typename T>
typename std::enable_if<is_tt<std::complex, T>::value, void>::type get_hdf5(T * l, H5::H5File *  file,  char * name) {
  H5::DataSet dataset = H5::DataSet(file->openDataSet(name));
  H5::CompType complex_data_type(sizeof(l[0]));
  typedef typename extract_value_type<T>::value_type value_type;
  
  complex_data_type.insertMember("r", 0, DataTypeFor<value_type>::value);
  complex_data_type.insertMember( "i", sizeof(value_type), DataTypeFor<value_type>::value);
  dataset.read(l, complex_data_type);
};

template <typename T>
typename std::enable_if<!is_tt<std::complex, T>::value, void>::type get_hdf5(T * l, H5::H5File *  file,  char * name) {
  H5::DataSet dataset = H5::DataSet(file->openDataSet(name));
  dataset.read(l, DataTypeFor<T>::value);
};

template <typename T>
typename std::enable_if<is_tt<std::complex, T>::value, void>::type get_hdf5(T * l, H5::H5File *  file,  std::string & name) {
  H5::DataSet dataset = H5::DataSet(file->openDataSet(name));
  H5::CompType complex_data_type(sizeof(l[0]));
  typedef typename extract_value_type<T>::value_type value_type;
  
  complex_data_type.insertMember("r", 0, DataTypeFor<value_type>::value);
  complex_data_type.insertMember( "i", sizeof(value_type), DataTypeFor<value_type>::value);
  dataset.read(l, complex_data_type);
};

template <typename T>
typename std::enable_if<!is_tt<std::complex, T>::value, void>::type get_hdf5(T * l, H5::H5File *  file,  std::string & name) {
  H5::DataSet dataset = H5::DataSet(file->openDataSet(name));
  dataset.read(l, DataTypeFor<T>::value);
};




template <typename T>
typename std::enable_if<!is_tt<std::complex, T>::value, void>::type write_hdf5(const Eigen::Array<T, -1, -1 > & mu,
									       H5::H5File *  file,
									       const std::string  name) {
  hsize_t    dims[2], chunk_dims[2]; // dataset dimensions
  dims[0] = chunk_dims[0] = mu.cols();
  dims[1] = chunk_dims[1] = mu.rows();      
  H5::DataSet dataset;
  H5::DataSpace dataspace = H5::DataSpace(2, dims );
  H5::DSetCreatPropList plist;
  plist.setChunk(2, chunk_dims);
  plist.setDeflate(6);
  
  try {
    H5::Exception::dontPrint();
    dataset = file->createDataSet(name, DataTypeFor<T>::value, dataspace);
  }
  catch (H5::FileIException & E) { 
    dataset = file->openDataSet(name);
  }
  
  dataset.write(mu.data(), DataTypeFor<T>::value);
};


template <typename T>
typename std::enable_if<is_tt<std::complex, T>::value, void>::type write_hdf5(const Eigen::Array<T, -1, -1 > & mu,
									      H5::H5File * file,
									      const std::string name) {
  hsize_t    dims[2], chunk_dims[2]; // dataset dimensions
  dims[0] = chunk_dims[0] = mu.cols();
  dims[1] = chunk_dims[1] = mu.rows();      
  H5::DataSet dataset;
  H5::DataSpace dataspace = H5::DataSpace(2, dims );
  typedef typename extract_value_type<T>::value_type value_type;
  
  H5::CompType complex_datatype(sizeof(T));
  complex_datatype.insertMember("r", 0, DataTypeFor<value_type>::value);
  complex_datatype.insertMember( "i", sizeof(value_type), DataTypeFor<value_type>::value);
  H5::DSetCreatPropList plist;
  
  plist.setChunk(2, chunk_dims);
  plist.setDeflate(6);
  
  try {
    H5::Exception::dontPrint();
    dataset = file->createDataSet(name, complex_datatype, dataspace);
  }
  catch (H5::FileIException & E) { 
    dataset = file->openDataSet( name);
  }
  
  dataset.write(mu.data(), complex_datatype);
};




