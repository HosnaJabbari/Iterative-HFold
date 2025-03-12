#ifndef H_STRUCT_H_
#define H_STRUCT_H_

#include "structs.h"
// another way to represent a structure. Used to measure free energy of a give structure etc.
typedef struct h_str_features
{
    short int pair;
    char type;                   // type can be 'H', 'S', 'I', 'M' etc
    short int num_branches;
    int bri[MAX_BRANCHES];      // the i of each branch
    int arc;					// keeps the left base pair of the arc
    h_str_features()
    {
        pair = -1;
        type = NONE;
        num_branches = 0;
        arc = -1;
    }
} h_str_features;

#endif /*H_STRUCT_H_*/
