/***************************************************************************
                          structs.h  -  description
                             -------------------
    begin                : Thu Apr 11 2002
    copyright            : (C) 2002 by Mirela Andronescu
    email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef STRUCTS_H
#define STRUCTS_H

#include "constants.h"
#include <vector>

typedef struct
{
    short int i;
    short int j;
} pair;

// not used
/*
struct W_node
{
    int energy;
    short int num_branches;           // 0 if this does not have any branch
    pair rbranch;                     // we need only the right branch, to calculate dangling energies with the next branch
    short int next_back;

    W_node ()
    {
        energy = MAXENERGY;
        num_branches = 0;
        next_back = -1;
    }

};
*/

// info from miscloop.dat
typedef struct
{
    // Extrapolation for large loops based on polymer theory
    // internal, bulge or hairpin loops > 30: dS(T)=dS(30)+param*ln(n/30)
    double param_greater30;    // we keep this fixed for learning
    // helix
    PARAMTYPE terminal_AU_penalty;
    // hairpin data from miscloop file
    PARAMTYPE hairpin_GGG;          // bonus for GGG hairpin
    PARAMTYPE hairpin_c1;           // c hairpin slope
    PARAMTYPE hairpin_c2;           // c hairpin intercept
    PARAMTYPE hairpin_c3;           // c hairpin of 3
    // internal loops
    PARAMTYPE asymmetry_penalty_max_correction;
    PARAMTYPE asymmetry_penalty_array[4];    //for param learning, only [0] and [1] are variable, [2] and [3] are fixed
    int gail_rule;
    // multi-branched loops
    PARAMTYPE multi_offset;
    PARAMTYPE multi_helix_penalty;
    PARAMTYPE multi_free_base_penalty;
    PARAMTYPE intermolecular_initiation;

    // 3 params to replace the tstackh table if MODEL==EXTENDED and parsi_tstackh is 1
    PARAMTYPE hairpin_AU_closure;   // does not need terminal_AU_penalty to be added to it
    PARAMTYPE hairpin_AG_mismatch;
    PARAMTYPE hairpin_GA_mismatch;
    PARAMTYPE hairpin_UU_mismatch;

    // 3 parameters which replace the tstacki table - added on Dec 20, 2006
    PARAMTYPE internal_AU_closure;   // does not need terminal_AU_penalty to be added to it
    PARAMTYPE internal_GA_AG_mismatch;
    // use GA and AG separately, and add GG, according to Schroeder_Turner_2000
    PARAMTYPE internal_AG_mismatch;
    PARAMTYPE internal_GA_mismatch;
    PARAMTYPE internal_GG_mismatch;

    // 6 more parameters for internal loops 3x3 or larger, following Chen_Turner_2006b
    // I only added the ones that would not be included in tstacki
    // These parameters are to be used when !parsi_special
    PARAMTYPE internal_special_3GA;     // 5'-YGGA/GAAR-3' or 5'-GGAR/YGAA-3', in loops 3x3 and larger
    PARAMTYPE internal_special_2GA;     // 5'-GA/GA-3' next to a closing base pair, or 5'-GG/AA-3' next to a closing base pair, for 3x3, 3x4, 4x4, and 4x5 loops;
                                //  ALSO 5'-RGGA/GAAY-3' or 5'-GGAY/YGAA-3' for 3x5, 3x6 and 4x6 loops.
                                // internal_2GA is USED only if internal_3GA was not used
    PARAMTYPE internal_special_2xGA_GC;  // 5'-GANGC/GANGC-3' in 3x3 loops
    PARAMTYPE internal_special_midGA;       // middle GA adjacent to RY in 3x3 loops
                                    // internal_midGA is USED only if none of internal_3GA and internal_2GA is used.
    PARAMTYPE internal_special_UG_AG;       // once or twice, for each 5'-UG/AG-3' at the terminus of loops 3x3 or larger
    PARAMTYPE internal_special_GU_A;        // first mismatch is GA, and U is 3' of G, for loops 3x3

    PARAMTYPE internal_UU_mismatch;
    PARAMTYPE internal22_delta_same_size;
    PARAMTYPE internal22_delta_different_size;
    PARAMTYPE internal22_delta_1stable_1unstable;
    PARAMTYPE internal22_delta_AC;
    PARAMTYPE internal22_match;   // if it has a CG or AU match in the middle

    PARAMTYPE internal11_basic_mismatch;
    PARAMTYPE internal11_GG_mismatch;

    PARAMTYPE internal11_AU_closure;
    PARAMTYPE internal11_GU_closure;
    PARAMTYPE internal11_AG_mismatch;
    PARAMTYPE internal11_UU_mismatch;
    PARAMTYPE internal11_5YRR_5YRR;
    PARAMTYPE internal11_5RYY_5RYY;
    PARAMTYPE internal11_5YYR_5YYR;
    PARAMTYPE internal11_5YRY_5RYR;
    PARAMTYPE internal11_5RRY_5RYY;

    PARAMTYPE internal21_match;   // if it has a CG or AU or GU match in the middle
    PARAMTYPE internal21_AU_closure;   // does not need terminal_AU_penalty to be added to it

    PARAMTYPE internal21_initiation;
    PARAMTYPE internal21_GU_closure;
    PARAMTYPE internal21_AG_mismatch;    // applied once per loop, not applied to 5'RA/3'YG loops
    PARAMTYPE internal21_GG_mismatch;    // applied once per loop
    PARAMTYPE internal21_UU_mismatch;    // applied once per loop
    PARAMTYPE internal22mid_group1;      // group 1 according to Christiansen_Znosko_2008
    // That is: a U · U pair adjacent to an R · R pair, a G · A  or A · G pair adjacent to a Y · Y pair, or  any combination of A · C, U · C, C · U,  C · C, C · A, or A · A pairs
    PARAMTYPE internal22mid_group2;      // group 2 according to Christiansen_Znosko_2008
    // That is: any combination of adjacent G · A and A · G pairs or two U · U pairs
    PARAMTYPE internal22mid_group3;      // group 3 according to Christiansen_Znosko_2008
    // That is: a U · U pair adjacent to a Y · Y (not U · U),  C · A, or A · C pair
    PARAMTYPE internal22mid_group4;      // group 4 according to Christiansen_Znosko_2008
    // That is: a G · G pair not adjacent to a U · U pair

    // 2 more parameters for the case parsi_int22 == 1
    PARAMTYPE internal22_AU_closure;     // as suggested by Christiansen_Znosko_2008
    PARAMTYPE internal22_GU_closure;     // as suggested by Christiansen_Znosko_2008

} miscinfo;

// info from tloop.dat
typedef struct
{
    char seq[10];
    PARAMTYPE energy;
} hairpin_tloop;


// the data structure stored in the V array
typedef struct minimum_fold
{
    short int pair;
    char type;                  // type can be 'H', 'S', 'I', 'M'
    char filled;                // I think this is not used any more
    minimum_fold()
    {
        pair = -1;
        type = NONE;
        filled = 'N';
    }
} minimum_fold;

// another way to represent a structure. Used to measure free energy of a give structure etc.
typedef struct str_features
{
    short int pair;
    char type;                   // type can be 'H', 'S', 'I', 'M' etc
    short int num_branches;
    int bri[MAX_BRANCHES];      // the i of each branch

    // Ian Wark August 8 2017
    // exists restricted is a very common function
    // that is ultimately very time consuming
    // precompute the results and save in this array
    std::vector< std::vector<int> > exists_restricted_arr;

    str_features()
    {
        pair = -1;
        type = NONE;
        num_branches = 0;
    }
} str_features;



typedef struct
{
        int top;
        int elem[MAXSLEN];
} stack_ds;



// This node is used to keep the intervals that need to be further backtracked
struct seq_interval
{
  int i;
  int j;
  PARAMTYPE energy;                        // it is used
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


struct struct_node
{
    minimum_fold* f;                    // an array
    seq_interval* intervals;            // M: a linked list
    PARAMTYPE bot_en;                         // not used?
    PARAMTYPE energy;                         // M: min energy of any structure starting with the partial structure so far
    char* structure;
    struct_node* previous;              // M: made doubly linked list, to be able to keep the size < limit
    struct_node* next;

    struct_node()
    {
        f = NULL;
        intervals = NULL;
        structure = NULL;
        previous = NULL;
        next = NULL;
    }
};


typedef struct seq_node
//class seq_node
{
    char* structure;
    seq_node* next;
}seq_node;


struct free_energy_node
{
    PARAMTYPE energy;
    char type;          // type may be: N (NONE), H (HAIRPIN), S (STACKED), I (INTERNAL), M (MULTI)
    free_energy_node()
    {
        energy = INF;
        type = NONE;
    }
};


// create double parameters, so that we have better precision for Maximum likelihood
// OBSOLETE
/*
typedef struct
{
    // Extrapolation for large loops based on polymer theory
    // internal, bulge or hairpin loops > 30: dS(T)=dS(30)+param*ln(n/30)
    double param_greater30;    // we keep this fixed for learning   // let's not make this complex yet
    // helix
    double terminal_AU_penalty;
    // hairpin data from miscloop file
    double hairpin_GGG;          // bonus for GGG hairpin
    double hairpin_c1;           // c hairpin slope
    double hairpin_c2;           // c hairpin intercept
    double hairpin_c3;           // c hairpin of 3
    // internal loops
    double asymmetry_penalty_max_correction;
    double asymmetry_penalty_array[4];    //for param learning, only [0] and [1] are variable, [2] and [3] are fixed
    int gail_rule;  // this is not a complex parameter
    // multi-branched loops
    double multi_offset;
    double multi_helix_penalty;
    double multi_free_base_penalty;
    double intermolecular_initiation;

    // 3 parameters which replace the tstackh table - added on Jun 16, 2008, for the parsimonious model
    double hairpin_AU_closure;   // does not need terminal_AU_penalty to be added to it
    double hairpin_AG_mismatch;
    double hairpin_UU_mismatch;

    // 3 parameters which replace the tstacki table - added on Dec 20, 2006
    double internal_AU_closure;   // does not need terminal_AU_penalty to be added to it
    double internal_AG_mismatch;
    double internal_UU_mismatch;
    double internal22_delta_same_size;
    double internal22_delta_different_size;
    double internal22_delta_1stable_1unstable;
    double internal22_delta_AC;
    double internal22_match;   // if it has a CG or AU match in the middle
    double internal21_match;   // if it has a CG or AU or GU match in the middle
    double internal21_AU_closure;   // does not need terminal_AU_penalty to be added to it
    double internal11_basic_mismatch;
    double internal11_GG_mismatch;
} miscinfo_double;


// info from tloop.dat
typedef struct
{
    char seq[6];  // pointer to the hairpin_tloop data structure, because no needed is necessary to the seq
    // the seq was commented out for the complex class
    double energy;
} hairpin_tloop_double;
*/

#endif

