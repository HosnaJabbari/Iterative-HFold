/***************************************************************************
                          simfold.h  -  description
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


#ifndef SIMFOLD_H
#define SIMFOLD_H

#include "constants.h"
#include "init.h"

double simfold (char *sequence, char *structure);
// PRE:  the init_data function has been called;
//       the space for structure has been allocated
// POST: fold sequence, return the MFE structure in structure, and return the MFE
// this function is defined in s_specific_functions.cpp

double simfold_loss_augmented (char *sequence, char *known_structure, char *pred_structure);
// like simfold, but loss-augmented
// Returns arg min_y (dG(x,y) - loss(y,y*))
//  where x is sequence and y* is known_structure            

double simfold_restricted (char *sequence, char *restricted, char *structure);
// PRE:  the init_data function has been called;
//       the space for structure has been allocated
//       In restricted, parentheses mean forced base pairs, dot means forced not to pair, 
//            underscore or space means not restricted      
// POST: fold sequence, return the MFE structure in structure, and return the MFE


double free_energy_simfold (char *sequence, char *structure);
// PRE:  sequence and structure are given as input
// !!! structure should not be sth like "...(((....)))", but a variable!!!
// POST: returns the free energy of sequence folded into structure
//       Worst case complexity = n^2 (n = length of sequence)
// (the function is implemented in s_specific_functions.cpp)

double free_energy_simfold_restricted (char *sequence, char *structure, char *restricted);
// PRE:  sequence and structure are given as input
// !!! structure should not be sth like "...(((....)))", but a variable!!!
//       In restricted, parentheses mean forced base pairs, dot means forced not to pair, 
//            underscore or space means not restricted      
// POST: returns the free energy of sequence folded into structure
//       Worst case complexity = n^2 (n = length of sequence)
// (the function is implemented in s_specific_functions.cpp)

double enthalpy_simfold (char *sequence, char *structure);
// PRE:  sequence and structure are given as input
// POST: returns the enthalpy of sequence folded into structure
//       Worst case complexity = n^2 (n = length of sequence)


// in s_specific_functions.cpp
int simfold_restricted_all_mfe_structures (char *sequence, char *restricted, char structures[][MAXSLEN], double energies[]);


int simfold_unordered_suboptimals (char *sequence, int number, char structures[][MAXSLEN], double energies[]);
// compute "number" suboptimal structures, with dangling energies added everywhere
// return the number of suboptimal structures (between 1 and number)

int simfold_ordered_suboptimals (char *sequence, int number, char structures[][MAXSLEN], double energies[]);
// compute "number" suboptimal structures, with dangling energies added everywhere
// in fact, compute 2*number structures, then compute the good free energy, reorder, and get the best "number" structures
// return the number of suboptimal structures (between 1 and number)

int simfold_restricted_unordered_suboptimals (char *sequence, char *restricted, int number, char structures[][MAXSLEN], double energies[]);


int simfold_restricted_all_suboptimals (char *sequence, char *restricted, char structures[][MAXSLEN], double energies[]);
// return MAXSUBSTR suboptimal structures, which have free energy <= 0. Note that the suboptimal structure computation adds both 5' and 3' dangling ends in any case. The suboptimal structures are not reordered in this function.

int simfold_unordered_suboptimals_range (char *sequence, double max_energy, char structures[][MAXSLEN], double energies[], int num_subopt);
// Compute all suboptimal structures (limited by MAXSUBSTR), which have free energy <= max_energy
//     the free energies include the 5' and 3' dangling energies in all cases

int simfold_restricted_unordered_suboptimals_range (char *sequence, char *restricted, double max_energy, char structures[][MAXSLEN], double energies[]);
// Compute all (restricted) suboptimal structures (limited by MAXSUBSTR), which have free energy <= max_energy
//     the free energies include the 5' and 3' dangling energies in all cases


// partition function calculations - in s_specific_functions.cpp
PFTYPE simfold_partition_function_smart (char *sequence, int ignore_dangles=0, int compute_gradient_dangles=1);

PFTYPE simfold_partition_function_smart_numerical (char *sequence, int ignore_dangles=0, int compute_gradient_dangles=1);

PFTYPE simfold_partition_function_smart_double (char *sequence, int ignore_dangles=0, int compute_gradient_dangles=1);

void simfold_gradient_smart (char *sequence, PFTYPE *grad, int ignore_dangles=0, int compute_gradient_dangles=1);

PFTYPE simfold_f_and_gradient_smart (char *sequence, char *restricted, PFTYPE *grad, int ignore_dangles=0, int compute_gradient_dangles=1);

PFTYPE simfold_f_and_gradient_smart_double (char *sequence, PFTYPE *grad, int ignore_dangles=0, int compute_gradient_dangles=1);

void simfold_gradient_smart_numerical (char *sequence, PFTYPE *grad, int ignore_dangles=0, int compute_gradient_dangles=1);

PFTYPE simfold_f_and_gradient_smart_numerical (char *sequence, PFTYPE *grad, int ignore_dangles=0, int compute_gradient_dangles=1, int which_param=-1);

PFTYPE simfold_partition_function_approximately (char *sequence);

void simfold_partition_function_both (char *sequence);

// in params.cpp
double simfold_restricted_partition_function (char *sequence, char *restricted, double &min_energy, double &max_energy, int &actual_num_str);
double compute_probability (double energy, double Z);


#endif
