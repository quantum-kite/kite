#include <H5Cpp.h>
#include <typeinfo>
#include <iostream>
#include <complex>

template <typename T>
void get_hdf5(T * l, H5::H5File *  file,  char * name)
{
  H5::DataSet dataset = H5::DataSet(file->openDataSet(name));
  H5::CompType complex_data_type(sizeof(l[0]));
  
  if(typeid(T) == typeid(int))
    dataset.read(l,H5::PredType::NATIVE_INT);

  if(typeid(T) == typeid(unsigned int))
    dataset.read(l,H5::PredType::NATIVE_UINT);
  
  if(typeid(T) == typeid(float))
    dataset.read(l,H5::PredType::NATIVE_FLOAT);
  
  if(typeid(T) == typeid(double))
    dataset.read(l,H5::PredType::NATIVE_DOUBLE);
  
  if(typeid(T) == typeid(long double))
    dataset.read(l,H5::PredType::NATIVE_LDOUBLE);
  
  if(typeid(T) == typeid(std::complex<float>))
    {
      complex_data_type.insertMember("r", 0, H5::PredType::NATIVE_FLOAT);
      complex_data_type.insertMember( "i", sizeof(float), H5::PredType::NATIVE_FLOAT);
      dataset.read(l, complex_data_type);  
    }
  
  if(typeid(T) == typeid(std::complex<double>))
    {
      complex_data_type.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
      complex_data_type.insertMember( "i", sizeof(double), H5::PredType::NATIVE_DOUBLE);
      dataset.read(l, complex_data_type);  
    }
  
  if(typeid(T) == typeid(std::complex<long double>))
    {
      complex_data_type.insertMember("r", 0, H5::PredType::NATIVE_LDOUBLE);
      complex_data_type.insertMember( "i", sizeof(long double), H5::PredType::NATIVE_LDOUBLE);
      dataset.read(l, complex_data_type);  
    }
}

template void get_hdf5<int>(int * l, H5::H5File *  file,  char * name);
template void get_hdf5<unsigned int>(unsigned int * l, H5::H5File *  file,  char * name);
template void get_hdf5<float>(float * l, H5::H5File *  file,  char * name);
template void get_hdf5<double>(double * l, H5::H5File *  file,  char * name);
template void get_hdf5<long double>(long double * l, H5::H5File *  file,  char * name);
template void get_hdf5<std::complex<float>>(std::complex<float> * l, H5::H5File *  file,  char * name);
template void get_hdf5<std::complex<double>>(std::complex<double> * l, H5::H5File *  file,  char * name);
template void get_hdf5<std::complex<long double>>(std::complex<long double> * l, H5::H5File *  file,  char * name);

