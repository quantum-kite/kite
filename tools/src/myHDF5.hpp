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
typename std::enable_if< is_tt<std::complex, T>::value, void>::type get_hdf5(T *, H5::H5File *,  char *);

template <typename T>
typename std::enable_if<!is_tt<std::complex, T>::value, void>::type get_hdf5(T *, H5::H5File *,  char *);

template <typename T>
typename std::enable_if< is_tt<std::complex, T>::value, void>::type get_hdf5(T *, H5::H5File *,  std::string &);

template <typename T>
typename std::enable_if<!is_tt<std::complex, T>::value, void>::type get_hdf5(T *, H5::H5File *,  std::string &);

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
