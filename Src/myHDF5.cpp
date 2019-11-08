#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"

template<>
H5::DataType DataTypeFor<int>::value = H5::PredType::NATIVE_INT;
template<>
H5::DataType DataTypeFor<unsigned int>::value = H5::PredType::NATIVE_UINT;
template<>
H5::DataType DataTypeFor<unsigned long>::value = H5::PredType::NATIVE_ULONG;
template<>
H5::DataType DataTypeFor<float>::value = H5::PredType::NATIVE_FLOAT;
template<>
H5::DataType DataTypeFor<double>::value = H5::PredType::NATIVE_DOUBLE;
template<>
H5::DataType DataTypeFor<long double>::value = H5::PredType::NATIVE_LDOUBLE;




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
}


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
}



template <typename T> 
void instantiateHDF<T>:: get_hdf5A(T *l, H5::H5File *file,  std::string & name) {
  get_hdf5<T>(l, file, name);
}

template <typename T> 
void instantiateHDF<T>:: get_hdf5A(T * l, H5::H5File * file,  char *name) {
  get_hdf5<T>(l, file, name);
}

template <typename T> 
void instantiateHDF<T>:: write_hdf5A(const Eigen::Array<T, -1, -1 > & mu, H5::H5File * file, const std::string name) {
  write_hdf5<T>(mu, file, name);
}

template struct instantiateHDF<int>;
template struct instantiateHDF<unsigned>;
template struct instantiateHDF<unsigned long>;
template struct instantiateHDF<float>;
template struct instantiateHDF<double>;
template struct instantiateHDF<long double>;
template struct instantiateHDF<std::complex<float>>;
template struct instantiateHDF<std::complex<double>>;
template struct instantiateHDF<std::complex<long double>>;

template void get_hdf5<int>(int*, H5::H5File*, char*);
template void get_hdf5<unsigned>(unsigned*, H5::H5File*, char*);
template void get_hdf5<unsigned long>(unsigned long*, H5::H5File*, char*);
template void get_hdf5<float>(float*, H5::H5File*, char*);
template void get_hdf5<double>(double*, H5::H5File*, char*);
template void get_hdf5<long double>(long double*, H5::H5File*, char*);
template void get_hdf5<std::complex<float>>(std::complex<float>*, H5::H5File*, char*);
template void get_hdf5<std::complex<double>>(std::complex<double>*, H5::H5File*, char*);
template void get_hdf5<std::complex<long double>>(std::complex<long double>*, H5::H5File*, char*);


template void get_hdf5<int>(int*, H5::H5File*, std::string &);
template void get_hdf5<unsigned>(unsigned*, H5::H5File*, std::string &);
template void get_hdf5<unsigned long>(unsigned long*, H5::H5File*, std::string &);
template void get_hdf5<float>(float*, H5::H5File*, std::string &);
template void get_hdf5<double>(double*, H5::H5File*, std::string &);
template void get_hdf5<long double>(long double*, H5::H5File*, std::string &);
template void get_hdf5<std::complex<float>>(std::complex<float>*, H5::H5File*, std::string &);
template void get_hdf5<std::complex<double>>(std::complex<double>*, H5::H5File*, std::string &);
template void get_hdf5<std::complex<long double>>(std::complex<long double>*, H5::H5File*, std::string &);
