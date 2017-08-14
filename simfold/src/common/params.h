/***************************************************************************
                          params.h  -  description
                             -------------------
    begin                : Sun May 8 2005
    copyright            : (C) 2005 by Mirela Andronescu
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

#ifndef PARAMS_H
#define PARAMS_H

#include "constants.h"

void print_stacking_energies();
// prints the stacking energies

void print_tstacki_dangling_energies();
// prints the tstacki and dangling energies

void print_stack_dangling_energies();
// prints the stack and dangling energies

void print_stack_equation_dangling();
// prints the stack and dangling energies

void print_int22_tstacki();
// prints the int22 and tstacki energies

void print_int22();
// prints the int22 energies

void test_int11_differences();
// they are the same as turner's, within 0.05 kcal/mol

void test_int22_differences();
// there are some differences here, because in the 53-param model I didn't include all the measured internal loops

int check_stability_and_size (int k, int l, int o, int p);
// helper function, to detect which delta we need for the int22 parameters

void count_AU_penalty (int base_i, int base_j, double *counter);
// PRE:  base_i and base_j make a pair
// POST: increment counter at the appropriate position
// Mirela: Nov 23, 2003

void count_penalty_by_size (int size, char type, double *counter);
// PRE:  size is the size of the loop
//       type is HAIRP or INTER or BULGE
// POST: return the penalty by size of the loop
// Mirela: Nov 23, 2003

void IL_count_penalty_by_size_2D (int size1, int size2, double *counter);

int structure_type_index (const char input[]);
  // Mirela: Nov 23, 2003
  // Given the type as a string, return the index in string_params

void count_asymmetry_penalty (int size1, int size2, double *counter);
// PRE:  size1 and size2 are the sizes of the two free base groups
// POST: Calculate the asymmetry penalty for internal loops
//       Note that if size1 == size2, pen is 0
// Mirela: Nov 23, 2003

void save_parameters (const char *filename);
  // Mirela: Dec 16, 2003
  // save all parameters in the given file

void save_parameters_in_array (PARAMTYPE *array);
// PRE: parameters have been read
// save all parameters in the given array


int create_building_block_strings ();
// Mirela: Sep 20, 2005
// For most of the building blocks, write which sequence(s) correspond to each building block
// also fill string_params

int create_string_params ();
// Mirela: Nov 23, 2003
// writes each parameter type, excluding duplicates, in a long vector, containing the names of the parameters

int get_num_params ();
// returns the number of parameters in the model
// right now calls create_string_params, which is inefficient, but

void save_paramtypes (const char *filename);
// PRE: call create_string_params ()
// save all parameter types in the given file
void save_paramtypes_machine_readable (const char *filename);

void set_starters ();
// set all the starters to the appropriate values, to be used in structure_type_index

double get_feature_counts (char *sequence, char *structure, char *restricted, double *c, double &f);
// This function computes the c vector and the f value of a linear energy function c' x + f
// a wrapper around count_each_structure_type
// returns the energy function

double get_feature_counts_quadratic (char *sequence, char *structure, char *restricted, double **quadratic, double *linear, double &f);
// This function computes the P matrix (quadratic), the c vector (linear) and the f value
//    of a quadratic energy function x'Px + c'x + f
// just a dummy function (we don't have quadratic energy model in MultiRNAFold, equivalent with get_feature_counts
// Used to test create_structural_constraints_simfold.cpp for the quadratic case

double get_feature_counts_restricted (char *sequence, char *structure, double *c, double &f, int reset_c, int ignore_dangles, int ignore_first_AU_penalty);
// This function computes the c vector and the f value of a linear energy function c' x + f
// structure contains the following characters:
//      - left and right parentheses mean base pairs, to be considered in the energy model
//      - dots mean unpaired bases, to be considered in the energy model
//      - a sequence of x's means ignore that region from the energy model
//      - left and right angles < > mean they enclose a pseudoknot, and ignore that region from the energy model.
// I'm assuming the following are true about structure:
//  - Regions with x are always inside parentheses or angle brackets: (xxxx) or <xxxx>.
//  - (..xx) or <...> or <...((..))> is fine too, I ignore everything that's inside <>
//      and I ignore everything that's inside a hairpin loop if it has at least one x.
//  - For a case like: (..(...)[[)]], this function will be called 3 times,
//      so structure can never be: (..(...)xx)
//  - An angle pair is always nested within parentheses: you can have this: (<xxx>) but not this (<)>
//  - If an angle pair is inside a multi-loop, consider it is a multi-loop branch
//      Example: (((...(((...)))..<xxxxxxx>.)))
//  - If an angle pair is nested inside another base pair, consider it is a multi-loop with 2 branches.
//      Example: ..(((...<xxxxxxx>..)))
//  - A multi-loop branch of type <xxxx> contributes 1 branch penalty
//
// f is the free value of the energy function c' x + f
// If reset_counter is 1, then counter[i] is set to 0 for all i, 0 <= i < num_params.
//      Otherwise, the new values get added to whatever was in counter.
// If ignore_dangles is 1, then dangling ends are not included anywhere
//  Otherwise, dangling ends are included whereever needed only in the parts containing parentheses and dots.
//      For example, in the case (..<xxx>..), dangling ends are included near the parentheses, but not near the angles.
// If ignore_AU is 1, then AU_penalty is not added at the end of the stem.
//      Otherwise, they are.
// Returns the value of the energy function
// Added on Feb 28, 2008
// Modification on May 26, 2008: if c is NULL, then it doesn't compute the c and f, it only returns the energy

int check_counts_linear (int numpars, double *params, double *c, double f, double energy);
// Return 1 if energy ~= c'x + f    (where x is params)
// Return 0 otherwise
// Added on Feb 29, 2008

int check_counts_quadratic (int numpars, double *params, double **P, double *c, double f, double energy);
// Return 1 if energy ~= x'Px + c'x + f    (where x is params)
// Return 0 otherwise
// Added on Feb 29, 2008

double count_each_structure_type (char *sequence, char *structure, char *restricted, double *counter, double &free_value, int reset);
// Mirela: Nov 23, 2003
// PRE: string_params have been filled, i.e. by calling create_string_params() or create_building_block_strings()
// GIven the sequence and the structure, returns the vector counter, with the #of elementary structures
// The sequence and structure may have a space inside, meaning duplex
// The free_value is f, the free value of the energy function c'x + f
// If reset is 1, then counter will be reset to 0 first, otherwise, the values will be added to the current value of counter
// returns the energy value

void print_counter (double *counter, double free_value);
// prints to the screen the counts of each parameter
// Added on Feb 29, 2008

double simfold_restricted_logZ (char *sequence, char *real_structure, char *restricted, double &min_energy, double &max_energy, int &actual_num_str);

void fill_data_structures_with_new_parameters (const char *filename);
  // Mirela: Dec 16, 2003
  // reads parameters from a file, and writes them in the internal data structures
  // PRE: first read the actual standard parameters, to be able to figure out which of them are
  // < INF, and maybe to also keep some old values.

void fill_data_structures_with_new_parameters_from_array (PARAMTYPE *array);
  // Mirela: 25 Aug 2008
  // reads parameters from and array, and writes them in the internal data structures
  // PRE: first read the actual standard parameters, to be able to figure out which of them are
  // < INF, and maybe to also keep some old values.


void fill_data_structures_with_new_parameters_double (const char *filename);

void fill_data_structures_with_new_parameters_fixed_dangles (const char *filename, char *dangfilename);
// reads all params from filename, except the dangling parameters, which reads from a different file
// assumes the fm363 model

int get_info_from_file (FILE *file, char *sequence, char *real_structure, char *restricted);
// return 1 if everything was ok
// return 0 instead of continue;

PFTYPE compute_f (char *input_file);
// PRE:  the parameters have been read
// POST: computes the f function, which is sum_{i=1}^N (1/RT G_{i,nat,theta} + logZ)

PFTYPE compute_likelihood_exactly (char *input_file);
// created on Oct 7, 2006
// PRE:  the parameters have been read
// POST: computes the likelihood function, which is \prod_{i=1}^N (exp (1/RT G_{i,nat,theta}) / Z)
// the partition function is computed exactly.
// the restriction string is not included so far

PFTYPE compute_log_likelihood_smart (char *input_file);

void compute_gradient_f (char *input_file, PFTYPE *f_gradient);
// the exhaustive way
// PRE:  the parameters have been read
// POST: compute the gradient of f, which is 1/RT sum_{i=1}^N (c_{nat}^i - sum_k c_k^i exp (-1/RT G_{i,k,theta}) / sum_k exp (-1/RT G_{i,k,theta}))

void compute_gradient_f_smart (char *input_file, PFTYPE *f_gradient);
// using dynamic programming
// PRE:  the parameters have been read
// POST: compute the gradient of f, which is 1/RT sum_{i=1}^N (c_{nat}^i - sum_k c_k^i exp (-1/RT G_{i,k,theta}) / sum_k exp (-1/RT G_{i,k,theta}))

PFTYPE compute_f_and_gradient_f_smart (char *input_file, PFTYPE *f_gradient);
// PRE:  the parameters have been read
// POST: compute the gradient of f, which is 1/RT sum_{i=1}^N (c_{nat}^i - sum_k c_k^i exp (-1/RT G_{i,k,theta}) / sum_k exp (-1/RT G_{i,k,theta}))
// returns f
// unrestricted for now
// DEPRICATED, another function exists in evaluate_function_smart_map

PFTYPE compute_f_and_gradient_f (char *input_file, PFTYPE *f_gradient);
// PRE:  the parameters have been read
// POST: computes the f function, which is sum_{i=1}^N (1/RT G_{i,nat,theta} + logZ)
// POST: compute the gradient of f, which is 1/RT sum_{i=1}^N (c_{nat}^i - sum_k c_k^i exp (-1/RT G_{i,k,theta}) / sum_k exp (-1/RT G_{i,k,theta}))

void compute_counts_vector_LP (char *input_file, double *total_counter);
// PRE:  the parameters have been read
// POST: given p input sequence/structure pairs, compute sum_{j=1}^p (c_j, cbar_j)
//          where c_j is the counts for the known structures
//          and cbar_j is the counts for the predicted mfe structure, with the given set of parameters

void compute_counts_matrix_LP (char *input_file, int train_samples);
// PRE:  the parameters have been read
//       train_samples = # of training instances. -1 if all.
// POST: given p input sequence/structure pairs, compute and display the counts c_j - cbar_j
//          where c_j is the counts for the known structures
//          and cbar_j is the counts for the predicted mfe structure, with the given set of parameters

int restricted_compatible (char *given_restricted, char *restricted);
// assume given_restricted doesn't have dots
// return 1 if the two are compatible
// given_restricted is input variable
// restricted is input-output, and it will contain the joint restriction

int generate_structure_withbb (char *sequence, char *known_structure, char *restricted, char *turner_structure, char structures[][MAXSLEN], double *old_counts, int threshold);
// written by Mirela in Sep, 2005
// PRE: the bbseq_left, bbstr_left, bbseq_right and bbstr_right are filled
//        i.e. call num_params = create_building_block_strings();
// First try to find if the building block bb_index can exist in any structure.
//   If found, restrict that part, and call simfold_restricted
//   Return the number of structures found
//   Look for at most threshold structures
//   num_params = total number of parameters
//   true_fe = the estimated free energy of the true structure.
//        We consider only those structures whose energies are <= true_fe
//   Update old_counts every time I consider a structure


int generate_structure_withbb_many_thresholds (char *sequence, char structures[][MAXSLEN], double energies[], double *old_counts, int many_thresholds[], double true_fe, int num_params);
// many_thresholds is a threshold per parameter (INF of no threshols wanted)


void search_bb (char *sequence, double *old_counts, int threshold, int num_params);

void print_parameters_in_almost_mfold_format ();
// save the parameters in a format that is easy to parse by an external perl script, and which puts the parameters in mfold (or RNAstructure) format.
// written on Oct 23, 2007

void print_parameters_in_ViennaRNA_format ();
// writes the parameters in Vienna RNA format, on the screen
// this function is somewhat similar to the function write_parameter_file from the Vienna RNA library
// see format documentation at http://www.tbi.univie.ac.at/~ivo/RNA/RNAlib/Param-Files.html
// DO NOT ADD ENTHALPIES, because my method doesn't produce enthalpies, only free energies
// THE VALUES for non-standard bases or base pairs are the same as in the original Vienna RNA file
// IT LOOKS LIKE THERE IS A MISTAKE int the Vienna RNA original file in the dangle sections!!!
// written on Nov 19, 2007

void fill_similarity_rule_with_optical_melting_reference (char *xml_filename);
// PRE: create_string_params must be called
// reads data from xml_filename, and fills up the array similarity_rule with the experiment id.
// started on Mar 18, 2008.

void fill_similarity_rules ();

int get_data_from_buffer (char *buffer, const char *header, const char last_char, char *output);
// function to get the sequence, structure etc data from the XML lines

int apply_rule_1 (int i, int j, int &i_rule1, int &j_rule1);
// check if i-j form a base pair: A-U, C-G or G-U.
// If they do, return 1, and write the replacement in i_rule1 and j_rule1
// If they don't, return 0, and i_rule1=i, j_rule1=j.

int is_int22_group_1 (int i, int j, int k, int l);
int is_int22_group_2 (int i, int j, int k, int l);
int is_int22_group_3 (int i, int j, int k, int l);
int is_int22_group_4 (int i, int j, int k, int l);

void extrapolate_parameters ();
// Start from the basic set of parameters, according to the model complexity,
//      and fill up the remaining structures.

int is_special_internal_1 (int *sequence, int i, int j, int ip, int jp);
int is_special_internal_2 (int *sequence, int i, int j, int ip, int jp);
int is_special_internal_3 (int *sequence, int i, int j, int ip, int jp);
int is_special_internal_4 (int *sequence, int i, int j, int ip, int jp);
int is_special_internal_5 (int *sequence, int i, int j, int ip, int jp);
int is_special_internal_6 (int *sequence, int i, int j, int ip, int jp);

PARAMTYPE special_energy_internal (int *sequence, int i, int j, int ip, int jp);
// Return the energy obtained when we consider 6 additional parameters for internal loop 3x3 and larger,
//  as described in Chen_Turner_2006b.
// the arguments are positions in sequence

PARAMTYPE count_special_internal (double *counter, int *sequence, int i, int j, int ip, int jp);
// Return the energy and counts obtained when we consider 6 additional parameters for internal loop 3x3 and larger,
//  as described in Chen_Turner_2006b.
// the arguments are positions in sequence

void check_int11_parameters (int i, int j, int k, int l, int m, int n);
        // check if int11 is the sum up of experimental addition etc.


#endif

