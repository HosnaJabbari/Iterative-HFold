#include <math.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
//#include <strstream.h>
//#include <fstream.h>
#include <sys/times.h>
#include <limits.h>
#include "timer.h"

// timer class for measuring cpu-time of program fragments */

timer::timer(void) {
  t_start_u=0;
  t_start_s=0;
  t_stop_u=0;
  t_stop_s=0;
  st = STOPPED;
  }

void timer::start(void) {
  struct tms buf;
  clock_t t;
  t = times(&buf);
  t_start_u = buf.tms_utime;
  t_start_s = buf.tms_stime;
  st = RUNNING;
  }

void timer::stop(void) {
  struct tms buf;
  clock_t t;
  t = times(&buf);
  t_stop_u = buf.tms_utime;
  t_stop_s = buf.tms_stime;  
  st = STOPPED;
  }

void timer::reset(void) {
  struct tms buf;
  clock_t t;
  t = times(&buf);
  t_start_u = buf.tms_utime;
  t_start_s = buf.tms_stime;
  }

float timer::elapsed(void) {
  struct tms buf;
  if (st==STOPPED)
    return(1.0*((t_stop_u+t_stop_s)-(t_start_u+t_start_s))/CLK_TCK);
  else {
    clock_t t;
    t = times(&buf);
    return(1.0*((buf.tms_utime+buf.tms_stime)-(t_start_u+t_start_s))/CLK_TCK);
    }
  }

float timer::resolution(void) {
  return(1.0/CLK_TCK);
  }

