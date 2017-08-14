/***************************************************************************
                          common.h  -  description
                             -------------------
    begin                : Thu Sep 5 2002
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

#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include "structs.h"
#include "s_partition_function.h"

#define isY(i)  (i==U || i==C)
#define isR(i)  (i==A || i==G)

#define MIN(A, B)      ((A) < (B) ? (A) : (B))
#define MAX(A, B)      ((A) > (B) ? (A) : (B))

// I precomputed these values in s_partition_function.cpp
//#define EXPA   (exp (misc.multi_offset * oneoverRT))
//#define EXPB(X)   (exp (((X)*misc.multi_helix_penalty) * oneoverRT))
//#define EXPC(X)   (exp (((X)*misc.multi_free_base_penalty) * oneoverRT))

#define IFD    if (ignore_dangles)

#define AU_penalty(X,Y)  ((((X) != C || (Y) != G) && ((X) != G || (Y) != C))?misc.terminal_AU_penalty:0)
#define has_AU_penalty(X,Y)  ((((X) != C || (Y) != G) && ((X) != G || (Y) != C))?1:0)
#define AU_penalty_enthalpy(X,Y)  ((((X) != C || (Y) != G) && ((X) != G || (Y) != C))?enthalpy_misc.terminal_AU_penalty:0)

//#define asymmetry_penalty(size1, size2) (MIN (misc.asymmetry_penalty_max_correction, abs (size1-size2) * misc.asymmetry_penalty_array [MIN (2, MIN ((size1), (size2)))-1]))

PARAMTYPE asymmetry_penalty (int size1, int size2);

#define asymmetry_penalty_enthalpy(size1, size2) (MIN (enthalpy_misc.asymmetry_penalty_max_correction, abs (size1-size2) * enthalpy_misc.asymmetry_penalty_array [MIN (2, MIN ((size1), (size2)))-1]))

#define IGINF(x) (((x) == INF)?0:(x))
// ignore infinite values

void get_sorted_positions (int n, double numbers[], int positions[]);
// does not modify numbers


double compute_accuracy (char *ref_structure, char *pred_structure);
double compute_sensitivity (char *ref_structure, char *pred_structure);
double compute_ppv (char *ref_structure, char *pred_structure);

double compute_pf_sensitivity (char *ref_structure, s_partition_function *part, double threshold);
// compute the sensitivity obtained after thresholding the base pair probabilities
// part is the partition function object, which contains base pair probabilities
// returns -1 if undefined (denominator is 0)

double compute_pf_ppv (char *ref_structure, s_partition_function *part, double threshold);
// compute the positive predictive value obtained after thresholding the base pair probabilities
// part is the partition function object, which contains base pair probabilities
// returns -1 if undefined (denominator is 0)

void giveup (const char *string1,const char *string2);
// to add: variable nb of parameters, as in scanf, printf

void giveup2 (const char *string1,const char *string2, FILE *file);
// to add: variable nb of parameters, as in scanf, printf

void create_random_sequence (int length, char *sequence);
// function to create uniformly random sequences - for demonstration purposes

void create_random_restricted (char *sequence, char *restricted);
// sequence is an input argument
// restricted is the output argument

void remove_space (char *structure);
// PRE: none
// POST: remove the space(s) from structure, if any; modifies structure

void empty_string (char * str);

int can_pair (int base1, int base2);
// PRE:  base1 and base2 are nucleotides over the alphabet {A, C, G, T, U}
// POST: return 1 if they can pair, 0 otherwise

int watson_crick (int base1, int base2);
// PRE:  base1 and base2 are nucleotides over the alphabet {A, C, G, T, U}
// POST: return 1 if they are watson crick pair, 0 otherwise

int nuc_to_int (char nucleotide);
// PRE:  nucleotide is 'A', 'C', 'G' or 'T'
// POST: Return 0 for A, 1 for C, 2 for G, 3 for T


char int_to_nuc (int inuc);


int is_nucleotide (char base);
// PRE:  base is a character
// POST: return true if base is a nucleotide (A, C, G or T)
//       return false otherwise


void check_sequence (char *sequence);
// check sequence for length and alphabet


PARAMTYPE penalty_by_size (int size, char type);
// PRE:  size is the size of the loop
//       type is HAIRP or INTER or BULGE
// POST: return the penalty by size of the loop

//PARAMTYPE IL_penalty_by_size_2D (int size1, int size2);

PARAMTYPE penalty_by_size_enthalpy (int size, char type);

void substr (char *source, int begin, int end, char *dest);
// PRE:  begin and end are smaller than strlen(source)
// POST: Put in dest what is in source between position begin and position end

void replace_str_piece (char *sequence, int position, char *seq);
// PRE:  begin + strlen(seq) < strlen (sequence)
// POST: In sequence, at position, replace what is was by seq

void reverse_complement_of_seq (const char *seq, char *complem);
// PRE:  seq is a sequence
// POST: complement and reverse sequence and put the result into compl

void insert_space (char *structure, int place);
// PRE:  None
// POST: insert a space at the specified place, in structure

void detect_original_pairs(char *structure, int *p_table);
// PRE:  structure contains the desired structure
// POST: p_table will contain the index of each base pair
//               or -1 if it does not pair
// Feb 28, 2008: structure can also have:
//  - angles: < or >, which denote the ends of a pseudoknot. In that case, p_table would still be filled in the same way.
//      The assumption is that the <> pairs are always nested within parentheses.
//      That is, a structure like this (<)> is not possible.
//  - x, which denotes that I should ignore that part. p_table would be -3 in that case/

int valid_structure (int i, int j, char *structure);
// returns 1 if this structure is valid (i.e. complete), 0 if it's partial

void detect_structure_features (char *structure, str_features *f);
// PRE:  None
// POST: The variable f is filled with structure features, i.e. what type of elementary structure
//       this base is closing (such as stacked pair, hairpin loop etc.)

int complementary_bases (char b1, char b2);
// returns 1 if b1 and b2 are complementary bases

int self_complementary (char *sequence);
// return 1 if this sequence is self-complementary
// self_complementary means the first half is the reverse complement of the second half
// if length (sequence) is an odd number, the middle base does not matter

int exists_restricted (int i, int j, str_features *fres);
int exists_restricted_ptable (int i, int j, int *ptable);

int is_structured (int i, int j, char *structure);
// return 1 if structure has some parentheses between i and j inclusive
// return 0 otherwise

void print_stacking_energies();
void print_tstacki_dangling_energies();
void print_stack_dangling_energies();
void print_stack_equation_dangling();
void print_int22_tstacki();

void read_parsi_options_from_file (char *filename);
// the file should contain values (0 or 1) for each of the parsi options. For example the following is all-parsimonious options
//     parsi_tstackh = 1;
//     parsi_tstacki = 1;
//     parsi_asymmetry = 1;
//     parsi_int11 = 1;
//     parsi_int21 = 1;
//     parsi_int22 = 1;
//     parsi_bulge1 = 1;
//     parsi_dangles = 1;
//     parsi_others = 1;
//     parsi_length = 1;
//     parsi_special = 1;


int loss (int first, int last);
// known_pairings contains the pairings from 0 to n-1, of the reference structure
// pred_pairings contains the pairings of a potential structure on the region first-last inclusive
//      the other regions don't matter
// Returns the "Hamming" distance between known_pairings and pred_pairings on the region first-last inclusive
// This function is used for the loss-augmented prediction
// Written on August 9, 2008
// Note: Maybe this measure is better than the Hamming distance:
//      (# correctly predicted bp - # incorrectly predicted bp) / # true bp.
//      This will be in (-inf,1], but it only includes the base pairs,
//      whereas the Hamming distance measure also includes the unpaired bases.

// Added on Sep 3, 2008, for loss-augmented prediction
double compute_distance (char *ref_structure, char *pred_structure);
// It has to be the same mathematical function as the one implemented in "loss"


#endif
