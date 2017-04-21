
#include "H5cpp.h"
#ifndef DTK_HDF5_HPP
#define DTK_HDF5_HPP
namespace dtk{
  template <class T>
  H5::PredType convert_type(T data);
  template <>
  H5::DataType convert_type(float data){
    return H5::PredType::NATIVE_FLOAT;
  }
  H5::DataType convert_type(double data){
    return H5::PredType::NATIVE_DOUBLE;
  }
  H5::DataType convert_type(int data){
    return H5::PredType::NATIVE_INT;
  }
  H5::DataType convert_type(int64_t data){
    return H5::PredType::NATIVE_INT64;
  }

template <class T>
void read_hdf5(std::string file_name,std::string var_name,std::vector<T>& out){
    H5::H5 file(file_name);
    H5::Data dataset = file(var_name);
    out.resize(dataset.getStorageSize()/sizeof(T));
    dataset.read(&out[0],var
  }
}
#endif //DTK_HDF5_HPP
