/***************************************************************************
                          globals.h  -  description
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

// this file contains mainly the global variables. It's the mirror of externs.h.
// Include globals.h in the driver, and externs.h in the library files. If you include globals.h
// twice, you get linking errors.

#ifndef GLOBALS_H
#define GLOBALS_H

#include "structs.h"
#include "constants.h"

int debug = 0;        // set debug 1 and recompile if you want to get a more verbose folding
int ignore_internal = 0;      // for debugging purposes, ignore internal loops
int ignore_multi = 0;   // for debugging purposes, ignore multi-loops
int max_internal_loop = MAXLOOP;
int *constraints;
int fix_dangles = 0;

int *known_pairings = NULL;    // used for the loss-augmented prediction
int *pred_pairings = NULL;

char similarity_rule[MAXNUMPARAMS][5000];

// these variables are used for all models: turner99, parsimonious and lavish
//  For each of the variables, LAVISH (0) means lavish, PARSI (1) means parsimonious, T99 (2) means turner99 (sometimes they overlap)
// In addition, I incorporated the Zhang_Liang_2008 model for bulge loop and internal loop length.
//      In this case, parsi_length is ZL (3). The only modifications are in params.cpp:extrapolate_parameters.
// For parsi_special, for now, by default it's only the new ones. To add the old tetraloops, do this:
//  cp params/turner-special-hairpins-rna.dat params/turner-special-hairpins-rna-only-new.dat
//  cp params/turner-special-hairpins-rna-with-T99.dat params/turner-special-hairpins-rna.dat
int parsi_tstackh = T99;   // 0: lavish and turner99 model for tstackh; 1: parsimonious model
int parsi_tstacki = T99;
int parsi_asymmetry = T99;    // category 3
int parsi_int11 = T99;      // category 4
int parsi_int21 = T99;      // category 5
int parsi_int22 = T99;      // category 6
int parsi_bulge1 = T99;    // 0: lavish for bulge1; 1: parsi for bulge1, i.e. 4 features; 2: turner99, i.e. no bulge features
int parsi_dangles = T99;      // category 8
int parsi_length = T99;      // category L
int parsi_special = T99;      // category S
int use_similarity_rules = 0;
int creating_model = 0;     // this is 1 only when I use the thermo xml set to figure out the experimental additions

// it's at least the specified number, updated after it goes once

//#if (MODEL == EXTENDED)   // I'm not using SIMPLE and EXTENDED model any more
int start_bulge = 21;
int start_internal = 25;
int start_internal11 = 25;
int start_internal11_C = 25;    //130;
int start_internal11_G = 25;    //204;
int start_internal11_U = 25;    //304;

int start_internal21 = 40;
int start_internal21_AUA = 40;      //46;   // first, second and third
int start_internal21_AUC = 40;      //142;
int start_internal21_AUG = 40;      //238;
int start_internal21_AUU = 40;      //334;
int start_internal21_CGA = 40;      //430;
int start_internal21_CGC = 40;      //526;
int start_internal21_CGG = 40;      //622;
int start_internal21_CGU = 40;      //718;
int start_internal21_GCA = 40;      //814;
int start_internal21_GCC = 40;      //910;
int start_internal21_GCG = 40;      //1006;
int start_internal21_GCU = 40;      //1102;
int start_internal21_GUA = 40;      //1198;
int start_internal21_GUC = 40;      //1198;
int start_internal21_GUG = 40;      //1198;
int start_internal21_GUU = 40;      //1198;
int start_internal21_UAA = 40;      //1582;
int start_internal21_UAC = 40;      //1582;
int start_internal21_UAG = 40;      //1582;
int start_internal21_UAU = 40;      //1582;
int start_internal21_UGA = 40;      //1966;
int start_internal21_UGC = 40;      //1966;
int start_internal21_UGG = 40;      //1966;
int start_internal21_UGU = 40;      //1966;

int start_internal22 = 46;
int start_internal22_AUA = 46;      //52;
int start_internal22_AUC = 46;      //52;
int start_internal22_AUG = 46;      //52;
int start_internal22_AUU = 46;      //52;
int start_internal22_CGA = 46;      //1468;
int start_internal22_CGC = 46;      //1468;
int start_internal22_CGG = 46;      //1468;
int start_internal22_CGU = 46;      //1468;
int start_internal22_GCA = 46;      //2628;
int start_internal22_GCC = 46;      //2628;
int start_internal22_GCG = 46;      //2628;
int start_internal22_GCU = 46;      //2628;
int start_internal22_GUA = 46;      //3532;
int start_internal22_GUC = 46;      //3532;
int start_internal22_GUG = 46;      //3532;
int start_internal22_GUU = 46;      //3532;
int start_internal22_UAA = 46;      //4180;
int start_internal22_UAC = 46;      //4180;
int start_internal22_UAG = 46;      //4180;
int start_internal22_UAU = 46;      //4180;
int start_internal22_UGA = 46;      //4572;
int start_internal22_UGC = 46;      //4572;
int start_internal22_UGG = 46;      //4572;
int start_internal22_UGU = 46;      //4572;

int start_dangle = 52;
int start_internal_size = 52;
int start_bulge_size = 61;
int start_hairpin_size = 63;
int start_misc_last = 70;       // the last misc: from misc.terminal_AU_penalty on
int start_special_hl = 75;

//**********************************************************************//
// Energies information for RNA
PARAMTYPE stack [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE tstackh [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE tstacki [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE int11   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE int21   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE int22   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE dangle_top  [NUCL] [NUCL] [NUCL];
PARAMTYPE dangle_bot  [NUCL] [NUCL] [NUCL];
PARAMTYPE internal_penalty_by_size [MAXLOOP+1];
PARAMTYPE bulge_penalty_by_size [MAXLOOP+1];
PARAMTYPE hairpin_penalty_by_size [MAXLOOP+1];
miscinfo misc;

//#if (MODEL == SIMPLE)
hairpin_tloop triloop[MAXTRILOOPNO];
hairpin_tloop tloop[MAXTLOOPNO];
int nb_triloops = 0;
int nb_tloops = 0;

//#elif (MODEL == EXTENDED)
hairpin_tloop special_hl[MAX_SPECIAL_LOOP_NO];
int nb_special_hl = 0;

// middle of asymmetric internal loops 2x2
//PARAMTYPE int22mid[NUCL] [NUCL] [NUCL] [NUCL];    // I'm not using this any more, I'm using int22_experimental_addition instead

PARAMTYPE int11_experimental_addition   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];       // values to be added to the simple 10-parameter model proposed by Davis_Znosko_2007, so that we use the experimental values for these parameters
PARAMTYPE int21_experimental_addition   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];       // values to be added to the simple 6-parameter model proposed by Badhwar_Znosko_2007, so that we use the experimental values for these parameters
PARAMTYPE int22_experimental_addition   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];       // values to be added to the simple 6-parameter model proposed by Christiansen_Znosko_2008, so that we use the experimental values for these parameters

PARAMTYPE internal_asymmetry_initiation;
PARAMTYPE internal_asymmetry_slope;
PARAMTYPE internal_asymmetry [MAXLOOP+1];
//PARAMTYPE internal_symmetry [MAXLOOP+1];
//PARAMTYPE internal_penalty_by_size_2D[MAXLOOP_I][MAXLOOP_I];

PARAMTYPE bulgeA;
PARAMTYPE bulgeC;
PARAMTYPE bulgeG;
PARAMTYPE bulgeU;
PARAMTYPE bulge1[NUCL] [NUCL] [NUCL] [NUCL] [NUCL];     // bulge of size 1     [i][j][k][ip][jp], where k=i+1, ip=k+1, j=jp+1
//#endif
//**********************************************************************//

//**********************************************************************//
// Enthalpies information for RNA
PARAMTYPE enthalpy_stack [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE enthalpy_tstackh [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE enthalpy_tstacki [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE enthalpy_int11   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE enthalpy_int21   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE enthalpy_int22   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
PARAMTYPE enthalpy_dangle_top  [NUCL] [NUCL] [NUCL];
PARAMTYPE enthalpy_dangle_bot  [NUCL] [NUCL] [NUCL];
PARAMTYPE enthalpy_internal_penalty_by_size [MAXLOOP+1];
PARAMTYPE enthalpy_bulge_penalty_by_size [MAXLOOP+1];
PARAMTYPE enthalpy_hairpin_penalty_by_size [MAXLOOP+1];
miscinfo enthalpy_misc;
hairpin_tloop enthalpy_triloop[MAXTRILOOPNO];
hairpin_tloop enthalpy_tloop[MAXTLOOPNO];
int enthalpy_nb_triloops = 0;
int enthalpy_nb_tloops = 0;
//**********************************************************************//

// parameters from the configuration file
int nb_params = 58;
char par_name [58] [100] = {
            "STD_DIR",
            "RNA_STACK_ENERGY37_V31",
            "RNA_STACK_ENERGY37_V23",
            "RNA_STACK_ENTHALPY",
            "RNA_LOOP_ENERGY37_V31",
            "RNA_LOOP_ENERGY37_V23",
            "RNA_LOOP_ENTHALPY",
            "RNA_TRILOOP_ENERGY37_V31",
            "RNA_TRILOOP_ENERGY37_V23",
            "RNA_TRILOOP_ENTHALPY",
            "RNA_TLOOP_ENERGY37_V31",
            "RNA_TLOOP_ENERGY37_V23",
            "RNA_TLOOP_ENTHALPY",
            "RNA_TSTACKH_ENERGY37_V31",
            "RNA_TSTACKH_ENERGY37_V23",
            "RNA_TSTACKH_ENTHALPY",
            "RNA_TSTACKI_ENERGY37_V31",
            "RNA_TSTACKI_ENERGY37_V23",
            "RNA_TSTACKI_ENTHALPY",
            "RNA_INT11_ENERGY37_V31",
            "RNA_INT11_ENERGY37_V23",
            "RNA_INT11_ENTHALPY",
            "RNA_INT21_ENERGY37_V31",
            "RNA_INT21_ENERGY37_V23",
            "RNA_INT21_ENTHALPY",
            "RNA_INT22_ENERGY37_V31",
            "RNA_INT22_ENERGY37_V23",
            "RNA_INT22_ENTHALPY",
            "RNA_MISCLOOP_ENERGY37_V31",
            "RNA_MISCLOOP_ENERGY37_V23",
            "RNA_MISCLOOP_ENTHALPY",
            "RNA_DANGLE_ENERGY37_V31",
            "RNA_DANGLE_ENERGY37_V23",
            "RNA_DANGLE_ENTHALPY",

            "DNA_STACK_ENERGY37",
            "DNA_STACK_ENTHALPY",
            "DNA_LOOP_ENERGY37",
            "DNA_LOOP_ENTHALPY",
            "DNA_TRILOOP_ENERGY37",
            "DNA_TRILOOP_ENTHALPY",
            "DNA_TLOOP_ENERGY37",


            "DNA_TLOOP_ENTHALPY",
            "DNA_TSTACKH_ENERGY37",
            "DNA_TSTACKH_ENTHALPY",
            "DNA_TSTACKI_ENERGY37",
            "DNA_TSTACKI_ENTHALPY",
            "DNA_INT11_ENERGY37",
            "DNA_INT11_ENTHALPY",
            "DNA_INT21_ENERGY37",
            "DNA_INT21_ENTHALPY",
            "DNA_INT22_ENERGY37",
            "DNA_INT22_ENTHALPY",
            "DNA_MISCLOOP_ENERGY37",
            "DNA_MISCLOOP_ENTHALPY",
            "DNA_DANGLE_ENERGY37",
            "DNA_DANGLE_ENTHALPY",
            "RNA_SPECIAL_HL_ENERGY37_V31",
            "RNA_SPECIAL_T99_L_ENERGY37_V31"
            };

char par_value [58] [1000];
char * std_dir_par                  = par_value[0];
char * rna_stack_energy37_v31_par   = par_value[1];
char * rna_stack_energy37_v23_par   = par_value[2];
char * rna_stack_enthalpy_par       = par_value[3];
char * rna_loop_energy37_v31_par    = par_value[4];
char * rna_loop_energy37_v23_par    = par_value[5];
char * rna_loop_enthalpy_par        = par_value[6];
char * rna_triloop_energy37_v31_par = par_value[7];
char * rna_triloop_energy37_v23_par = par_value[8];
char * rna_triloop_enthalpy_par     = par_value[9];
char * rna_tloop_energy37_v31_par   = par_value[10];
char * rna_tloop_energy37_v23_par   = par_value[11];
char * rna_tloop_enthalpy_par       = par_value[12];
char * rna_tstackh_energy37_v31_par = par_value[13];
char * rna_tstackh_energy37_v23_par = par_value[14];
char * rna_tstackh_enthalpy_par     = par_value[15];
char * rna_tstacki_energy37_v31_par = par_value[16];
char * rna_tstacki_energy37_v23_par = par_value[17];
char * rna_tstacki_enthalpy_par     = par_value[18];
char * rna_int11_energy37_v31_par   = par_value[19];
char * rna_int11_energy37_v23_par   = par_value[20];
char * rna_int11_enthalpy_par       = par_value[21];
char * rna_int21_energy37_v31_par   = par_value[22];
char * rna_int21_energy37_v23_par   = par_value[23];
char * rna_int21_enthalpy_par       = par_value[24];
char * rna_int22_energy37_v31_par   = par_value[25];
char * rna_int22_energy37_v23_par   = par_value[26];
char * rna_int22_enthalpy_par       = par_value[27];
char * rna_miscloop_energy37_v31_par= par_value[28];
char * rna_miscloop_energy37_v23_par= par_value[29];
char * rna_miscloop_enthalpy_par    = par_value[30];
char * rna_dangle_energy37_v31_par  = par_value[31];
char * rna_dangle_energy37_v23_par  = par_value[32];
char * rna_dangle_enthalpy_par      = par_value[33];

char * dna_stack_energy37_par       = par_value[34];
char * dna_stack_enthalpy_par       = par_value[35];
char * dna_loop_energy37_par        = par_value[36];
char * dna_loop_enthalpy_par        = par_value[37];
char * dna_triloop_energy37_par     = par_value[38];
char * dna_triloop_enthalpy_par     = par_value[39];
char * dna_tloop_energy37_par       = par_value[40];
char * dna_tloop_enthalpy_par       = par_value[41];
char * dna_tstackh_energy37_par     = par_value[42];
char * dna_tstackh_enthalpy_par     = par_value[43];
char * dna_tstacki_energy37_par     = par_value[44];
char * dna_tstacki_enthalpy_par     = par_value[45];
char * dna_int11_energy37_par       = par_value[46];
char * dna_int11_enthalpy_par       = par_value[47];
char * dna_int21_energy37_par       = par_value[48];
char * dna_int21_enthalpy_par       = par_value[49];
char * dna_int22_energy37_par       = par_value[50];
char * dna_int22_enthalpy_par       = par_value[51];
char * dna_miscloop_energy37_par    = par_value[52];
char * dna_miscloop_enthalpy_par    = par_value[53];
char * dna_dangle_energy37_par      = par_value[54];
char * dna_dangle_enthalpy_par      = par_value[55];
char * rna_special_hl_energy_par    = par_value[56];
char * rna_special_t99_l_energy_par    = par_value[57];


// The following are used for parameter learning only - you don't need them for folding

int simple_internal_energy = 0;
int simple_dangling_ends = 1;
// Until Sep 2, 2008, simple_dangling_ends was 0. -- This might cause problems in CG!!
// Sep 2, 2008: Changed this to 1. This should be 1 by default, otherwise the counts get screwed up
//      if a 3' and a 5' dangle are the same. Then one can be picked once and another can be picked another time.
// if 1, don't do minimization in multi-loops and exterior loops, just add the one at the 3' end
int no_dangling_ends = 0;    // if 1, don't add dangling ends at all
char string_params[MAXNUMPARAMS][MAXPNAME];    // for playing with the parameters
char string_params_human_readable[MAXNUMPARAMS][MAXPNAME];    // another version of string_params, more human readable
int num_params;

char bbseq_left[MAXNUMPARAMS][10];   // the 5' sequence of the corresponding building block
char bbseq_right[MAXNUMPARAMS][10];   // the 3' sequence of the corresponding building block
char bbstr_left[MAXNUMPARAMS][10];   // the 5' structure of the corresponding building block
char bbstr_right[MAXNUMPARAMS][10];   // the 3' structure of the corresponding building block

int counter_min_dangle[NUM_DANG][NUM_DANG];    // cell[0][1] says how many times min(dangle0,dangle1) appear, where dangle0 is x213 in the 318 nomenclature
// there are 48 dangling ends: 24 top and 24 bottom.

// OBSOLETE
// energies information
// double stack_double [NUCL] [NUCL] [NUCL] [NUCL];
// double tstackh_double [NUCL] [NUCL] [NUCL] [NUCL];
// double tstacki_double [NUCL] [NUCL] [NUCL] [NUCL];
// double int11_double   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
// double int21_double   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
// double int22_double   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
// double dangle_top_double  [NUCL] [NUCL] [NUCL];
// double dangle_bot_double  [NUCL] [NUCL] [NUCL];
// double internal_penalty_by_size_double [MAXLOOP+1];
// double bulge_penalty_by_size_double [MAXLOOP+1];
// double hairpin_penalty_by_size_double [MAXLOOP+1];
// miscinfo_double misc_double;
// hairpin_tloop_double triloop_double[MAXTRILOOPNO];
// hairpin_tloop_double tloop_double[MAXTLOOPNO];

char linker[6] ="XXXXX";
int linker_pos=0;
int linker_length = 5;
char structure_one_type[3];
char structure_two_type[3];
bool is_hybrid = false;

// AP
double PMO_RNA_penalty = 0 ;

#endif
