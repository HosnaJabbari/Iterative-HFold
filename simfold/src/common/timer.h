#include <sys/times.h>
#include <limits.h>

// timer class for measuring cpu-time of program fragments */
class timer {
  clock_t t_start_u, t_start_s;
  clock_t t_stop_u, t_stop_s;
  enum state {RUNNING,STOPPED};
  state st;
  
public:
  timer(void);
  void start(void);
  void stop(void);
  void reset(void);
  float resolution(void); 
  // returns resolution of timer in CPU-seconds */
  float elapsed(void);  // returns stopped time in CPU-seconds */
};

#ifndef CLK_TCK
#define CLK_TCK 60
#endif
