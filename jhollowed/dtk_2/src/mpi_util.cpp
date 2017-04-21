

#include <mpi.h>
#include <string.h>
#include "dtk/mpi_util.hpp"
namespace dtk{
  void broadcast_param(dtk::Param& param,int root,MPI_Comm comm){
    int myrank;
    MPI_Comm_rank(comm,&myrank);
    std::string param_str;
    char* param_cstr=0;
    int   param_cstr_size;
    if(myrank == root){
      param_str = param.stringify();
      param_cstr_size = param_str.size()+1; //+1 for null char
      param_cstr = new char[param_cstr_size];
      strncpy(param_cstr,param_str.c_str(),param_cstr_size);
    }
    MPI_Bcast(&param_cstr_size,1,MPI_INT,root,comm);
    if(myrank != root)
      param_cstr = new char[param_cstr_size];
    MPI_Bcast(param_cstr,param_cstr_size,MPI_CHAR,root,comm);
    if(myrank != root){
      param_str.assign(param_cstr,param_cstr+param_cstr_size);
      param.parse_string(param_str);
    }
    delete [] param_cstr;
  }
}
