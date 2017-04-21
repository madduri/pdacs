
#include "H5cpp.h"

namespace dtk{
  template <class T>
  H5::Type convert_type(T data);
  H5::Type convert_type(T data){

  }

template <class T>
void read_hdf5(std::string file_name,std::string var_name,std::vector<T>& out){
    H5::H5 file(file_name);
    H5::Data dataset = file(var_name);
    out.resize(dataset.getStorageSize()/sizeof(T));
    dataset.read(&out[0],var
  }
}
