
#include <sstream>

#include "dtk/util.hpp"


namespace dtk{
  std::string cat_strings(std::vector<std::string> strings){
    std::stringstream ss;
    for(uint i =0;i<strings.size();++i){
      ss<<strings[i];
    }
    return ss.str();
  }


  bool ensure_dir(std::string dir){
#ifdef unix
    std::stringstream ss;
    ss<<"mkdir -p "<<dir;
    int ret = system(ss.str().c_str());
    return ret==0;
#else 
    return false;
#endif //unix
  }

  bool ensure_file_path(std::string file){
    size_t end = file.find_last_of("\\/");
    if(end == std::string::npos)
      return false;
    return ensure_dir(file.substr(0,end));
  }
  

}
