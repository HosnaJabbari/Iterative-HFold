#ifndef H_COMMON_H_
#define H_COMMON_H_

#include "h_struct.h"
#include "common.h"

//Hosna: June 22, 2007

// #define PS_penalty = 960 		//exterior pseudoloop initiation penalty (9.6 Kcal/mol)
// #define PSM_penalty 		1500		//penalty for introducing pseudoknot inside a multiloop (15 Kcal/mol)
// #define PSP_penalty 		1500		//penalty for introducing pseudoknot inside a pseudoloop (15 Kcal/mol)
// #define PB_penalty 			20			//band penalty (0.2 Kcal/mol)
// #define PUP_penalty			10			//penalty for an un-paired base in a pseudoloop or a band (0.1 Kcal/mol)
// #define PPS_penalty 		10			//penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band(0.1 Kcal/mol)

// #define a_penalty			340			//penalty for introducing a multiloop (3.4 Kcal/mol)
// #define b_penalty			40			//penalty for base pair in a multiloop (0.4 Kcal/mol)
// #define c_penalty			0			//penalty for un-paired base in a multi-loop

// #define e_stP_penalty		0.83		// e_stP = 0.83 * e_s
// #define e_intP_penalty		0.83		// e_intP = 0.83 * e_int

// #define ap_penalty			340			//penalty for introducing a multiloop that spans a band (3.4 Kcal/mol)
// #define bp_penalty			40			//base pair penalty for a multiloop that spans a band (0.4 Kcal/mol)
// #define cp_penalty          0			//penalty for unpaired base in a multiloop that spans a band



#define P_WMB			'R'
#define P_VP			'D'
#define P_VPP			'E'
#define P_WI			'G'
#define P_BE			'J'
#define P_WIP			'L'
#define P_WMBP			'T'
#define P_V				'A' // This case is only for the cases that we have some pairings in input structure that are not valid in terms of simfold restrictions


#define NOT_COVERED		-1
#define STACK_EMPTY -1 // originally this value is 0, which I think is wrong! Hosna, March 8, 2012
#define RESTRICTED_UNPAIR -1
#define FREE_TO_PAIR	-2

// Hosna, March 19, 2012
// the original value of the matrices should be set to -INF and then be changed to their correct value
#define MINUS_INF             -1600000      // a very small value (minus infinity)

void detect_original_pairs_arcs(char *structure, int *p_table, int *arc_table);
void detect_original_PKed_pairs(char *structure, int *p_table);
double compute_h_sensitivity (char *ref_structure, char *pred_structure);
double compute_h_ppv (char *ref_structure, char *pred_structure);
void h_fill_data_structures_with_new_parameters (char *filename);
void h_fill_data_structures_with_new_parameters (double *param);
int h_create_string_params();

// Hosna: helper function to fill in the weakly_closed array
void detect_weakly_closed(h_str_features *fres, int *weakly_closed, int nb_nucleotides, int *index);

// Hosna: helper function to fill in the not_paired_all array
void detect_not_paired_all(h_str_features *fres, int *not_paired_all, int nb_nucleotides, int *index);

// Hosna: this function fills the bs table which keeps track of
// bs and Bs for each i and l
void detect_border_bs(h_str_features *fres, int** border_bs, int nb_nucleotides);

// Hosna: this function filld the bps table which keeps track of
// b' and B' for each l and j
void detect_border_bps(h_str_features *fres, int** border_bps, int nb_nucleotides);

void h_init (stack_ds *st);
void h_push (stack_ds *st, int el);
int h_pop (stack_ds *st);

// Hosna June 26, 2007
// I need functions to convert str_features to h_str_features and the other way around

h_str_features *convert_str_features_to_h_str_features(str_features *f);
str_features *convert_h_str_features_to_str_features(h_str_features *f);

void detect_h_structure_features (char *structure, h_str_features *f);

#endif /*H_COMMON_H_*/
