#ifndef DTK_TIMER_HPP
#define DTK_TIMER_HPP

#include <sys/time.h>
#include <string>

namespace dtk{
  class Timer{
  public:
    Timer();
    Timer(bool autostart);
    virtual ~Timer();
    Timer&    start();
    Timer&    stop();
    double  const get_seconds();
    double  const get_mseconds();
    double  const get_useconds();
    std::string timef() const;
  private:
    struct timeval start_;
    struct timeval end_;

  };
}

std::ostream& operator<<(std::ostream& os, dtk::Timer const &t);

#endif
