#ifndef DTK_MATH_HPP
#define DTK_MATH_HPP

namespace dtk{
  
  template<class T,class U>
  T max(T* data, U size){
    T max_val = data[0];
    for(U i =1;i<size;++i){
      if(data[i]>max_val)
	max_val = data[i];
    }
    return max_val;
  }
  template<class T,class U>
  T min(T* data, U size){
    T min_val = data[0];
    for(U i =1;i<size;++i){
      if(data[i]>min_val)
	min_val = data[i];
    }
    return min_val;
  }
  template<class T,class U>
  T average(T* data, U size){
    T result = 0;
    for(U i =0;i<size;++i){
      result += data[i];
    }
    return result/static_cast<T>(size);
  }
  template<class T>
  T average(std::vector<T> data){
    size_t size = data.size();
    T result = 0;
    for(size_t i =0;i<size;++i){
      result += data[i];
    }
    return result/static_cast<T>(size);

  }
}






#endif
