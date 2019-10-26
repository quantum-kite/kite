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





template <typename T>
typename std::enable_if<is_tt<std::complex, T>::value, void>::type get_hdf5(T * l, H5::H5File *  file,  char * name) {
  H5::DataSet dataset = H5::DataSet(file->openDataSet(name));
  H5::CompType complex_data_type(sizeof(l[0]));
  typedef typename extract_value_type<T>::value_type value_type;
  
  complex_data_type.insertMember("r", 0, DataTypeFor<value_type>::value);
  complex_data_type.insertMember( "i", sizeof(value_type), DataTypeFor<value_type>::value);
  dataset.read(l, complex_data_type);
}

template <typename T>
typename std::enable_if<!is_tt<std::complex, T>::value, void>::type get_hdf5(T * l, H5::H5File *  file,  char * name) {
  H5::DataSet dataset = H5::DataSet(file->openDataSet(name));
  dataset.read(l, DataTypeFor<T>::value);
}

template <typename T>
typename std::enable_if<is_tt<std::complex, T>::value, void>::type get_hdf5(T * l, H5::H5File *  file,  std::string & name) {
  H5::DataSet dataset = H5::DataSet(file->openDataSet(name));
  H5::CompType complex_data_type(sizeof(l[0]));
  typedef typename extract_value_type<T>::value_type value_type;
  
  complex_data_type.insertMember("r", 0, DataTypeFor<value_type>::value);
  complex_data_type.insertMember( "i", sizeof(value_type), DataTypeFor<value_type>::value);
  dataset.read(l, complex_data_type);
}

template <typename T>
typename std::enable_if<!is_tt<std::complex, T>::value, void>::type get_hdf5(T * l, H5::H5File *  file,  std::string & name) {
  H5::DataSet dataset = H5::DataSet(file->openDataSet(name));
  dataset.read(l, DataTypeFor<T>::value);
}


/*

template <typename T>
typename std::enable_if< is_tt<std::complex, T>::value, void>::type get_hdf5(T *, H5::H5File *,  char *);

template <typename T>
typename std::enable_if<!is_tt<std::complex, T>::value, void>::type get_hdf5(T *, H5::H5File *,  char *);

template <typename T>
typename std::enable_if< is_tt<std::complex, T>::value, void>::type get_hdf5(T *, H5::H5File *,  std::string &);

template <typename T>
typename std::enable_if<!is_tt<std::complex, T>::value, void>::type get_hdf5(T *, H5::H5File *,  std::string &);
*/
template <typename T>
typename std::enable_if<!is_tt<std::complex, T>::value, void>::type write_hdf5(const Eigen::Array<T, -1, -1 > & ,H5::H5File *, const std::string);

template <typename T>
typename std::enable_if<is_tt<std::complex, T>::value, void>::type write_hdf5(const Eigen::Array<T, -1, -1 > &, H5::H5File *, const std::string);

template <typename T>
struct instantiateHDF {
  void write_hdf5A(const Eigen::Array<T, -1, -1 > &, H5::H5File *, const std::string);
  void get_hdf5A(T *, H5::H5File *,  std::string &);
  void get_hdf5A(T *, H5::H5File *,  char *);
};





 




