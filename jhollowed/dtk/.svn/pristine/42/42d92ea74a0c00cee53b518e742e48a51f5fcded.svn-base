

#include <stdint.h>
#define GENERICIO_NO_MPI
#include "GenericIO.h"
#include <sstream>

template <class T>
void read_gio(char* file_name, std::string var_name, T*& data){
  gio::GenericIO reader(file_name);
  reader.openAndReadHeader(gio::GenericIO::MismatchAllowed);
  int num_ranks = reader.readNRanks();
  uint64_t max_size = 0;
  uint64_t rank_size[num_ranks];
  for(int i =0;i<num_ranks;++i){
    rank_size[i] = reader.readNumElems(i);
    if(max_size < rank_size[i])
      max_size = rank_size[i];
  }
  T* rank_data = new T[max_size+reader.requestedExtraSpace()/sizeof(T)];
  int64_t offset =0;
  reader.addVariable(var_name,rank_data,true);
  for(int i=0;i<num_ranks;++i){
    reader.readData(i,false);
    std::copy(rank_data,rank_data+rank_size[i],data+offset);
    offset +=rank_size[i];
  }
  delete [] rank_data;
  reader.close();
}
extern "C" int64_t get_elem_num(char* file_name);

extern "C" void read_gio_float (char* file_name, char* var_name, float* data);
extern "C" void read_gio_double(char* file_name, char* var_name, double* data);
extern "C" void read_gio_int32 (char* file_name, char* var_name, int* data); 
extern "C" void read_gio_int64 (char* file_name, char* var_name, int64_t* data);
enum var_type{
  float_type=0,
  double_type=1,
  int32_type=2,
  int64_type=3,
  type_not_found=9,
  var_not_found=10
};
extern "C" var_type get_variable_type(char* file_name,char* var_name);
extern "C" void inspect_gio(char* file_name);
