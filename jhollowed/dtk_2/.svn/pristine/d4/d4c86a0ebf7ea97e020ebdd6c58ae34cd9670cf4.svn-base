 

#include "gio.hpp"
#include <GenericIO.h>
#include <iostream>
  void read_gio_float(char* file_name, char* var_name, float* data){
    read_gio<float>(file_name,var_name,data);
  }
  void read_gio_double(char* file_name, char* var_name, double* data){
    read_gio<double>(file_name,var_name,data);
  }
  void read_gio_int32(char* file_name, char* var_name, int* data){
    read_gio<int>(file_name,var_name,data);
  }
  void read_gio_int64(char* file_name, char* var_name, int64_t* data){
    read_gio<int64_t>(file_name,var_name,data);
  }
  
  int64_t get_elem_num(char* file_name){
    gio::GenericIO reader(file_name);
    reader.openAndReadHeader(gio::GenericIO::MismatchAllowed);
    int num_ranks = reader.readNRanks();
    uint64_t size = 0;
    for(int i =0;i<num_ranks;++i)
      size +=reader.readNumElems(i);
    reader.close();
    return size;
  }

  var_type get_variable_type(char* file_name,char* var_name){
   gio::GenericIO reader(file_name);
   std::vector<gio::GenericIO::VariableInfo> VI;
   reader.openAndReadHeader(gio::GenericIO::MismatchAllowed);
   reader.getVariableInfo(VI);

   int num =VI.size();
    for(int i =0;i<num;++i){
      gio::GenericIO::VariableInfo vinfo = VI[i];
      if(vinfo.Name == var_name){
	if(vinfo.IsFloat && vinfo.Size == 4)
	  return float_type;
	else if(vinfo.IsFloat && vinfo.Size == 8)
	  return double_type;
	else if(!vinfo.IsFloat && vinfo.Size == 4)
	  return int32_type;
	else if(!vinfo.IsFloat && vinfo.Size == 8)
	  return int64_type;
	else
	  return type_not_found;
      }
    }
    return var_not_found;
      
  }

extern "C" void inspect_gio(char* file_name){
  int64_t size = get_elem_num(file_name);
  gio::GenericIO reader(file_name);
  std::vector<gio::GenericIO::VariableInfo> VI;
  reader.openAndReadHeader(gio::GenericIO::MismatchAllowed);
  reader.getVariableInfo(VI);
  std::cout<<"Number of Elements: "<<size<<std::endl;
  int num =VI.size();
  std::cout<<"[data type] Variable name"<<std::endl;
  std::cout<<"---------------------------------------------"<<std::endl;
  for(int i =0;i<num;++i){
    gio::GenericIO::VariableInfo vinfo = VI[i];

    if(vinfo.IsFloat)
      std::cout<<"[f";
    else
      std::cout<<"[i";
    std::cout<<" "<<vinfo.Size*8<<"] ";
    std::cout<<vinfo.Name<<std::endl;
  }
  std::cout<<"\n(i=integer,f=floating point, number bits size)"<<std::endl;
}

