#ifndef H_STRUCT_H_
#define H_STRUCT_H_

#include "constants.hh"
#include <vector>

// the data structure stored in the V array
typedef struct minimum_fold
{
    short int pair;
    char type;                  // type can be 'H', 'S', 'I', 'M'
    minimum_fold()
    {
        pair = -1;
        type = NONE;
    }
} minimum_fold;

// This node is used to keep the intervals that need to be further backtracked
struct seq_interval
{
  int i;
  int j;
  int energy;                        // it is used
  char type;
  seq_interval* next;

  void copy (seq_interval *other)
  {
    other->i = i;
    other->j = j;
    other->energy = energy;
    other->type = type;
  }
};

struct free_energy_node
{
    int energy;
    char type;          // type may be: N (NONE), H (HAIRPIN), S (STACKED), I (INTERNAL), M (MULTI)
    free_energy_node()
    {
        energy = 10000; // INF
        type = NONE;
    }
};

#endif /*H_STRUCT_H_*/
