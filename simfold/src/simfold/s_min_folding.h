/***************************************************************************
                          s_min_folding.cpp  -  description
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

// This class is the main class to compute the MFE prediction


#ifndef S_MIN_FOLDING_H
#define S_MIN_FOLDING_H


#include <stdio.h>
#include <string.h>
#include "structs.h"
#include "s_energy_matrix.h"
#include "s_multi_loop.h"


class s_min_folding
{
    public:
  
        s_min_folding (char *seq);
        // constructor for the unrestricted mfe case 
        
        s_min_folding (char *seq, char* restricted);
        // constructor for the restricted mfe case
                      
        ~s_min_folding ();
        // The destructor
        
        double s_simfold ();
        // PRE:  the init_data function has been called;
        //       the space for structure has been allocated
        // POST: fold sequence, return the MFE structure in structure, and return the MFE        
        
        double s_simfold_restricted ();
        // PRE:  the init_data function has been called;
        //       the space for structure has been allocate
        // POST: fold sequence, return the MFE structure in structure, and return the MFE
        
        void return_structure (char *structure) { strcpy (structure, this->structure); }        
        // writes the predicted MFE structure into structure       

    // better to have protected variable rather than private, it's necessary for Hfold
    protected:
    //private:
    
        char* structure;        // MFE structure
        s_hairpin_loop *H;      // hairpin loop object  
        s_stacked_pair *S;      // stack pair object 
        s_internal_loop *VBI;   // internal loop object
        s_multi_loop *VM;       // multi loop object
        s_energy_matrix *V;     // the V object
        int* int_sequence;      // sequence in integer representation (faster)  

        minimum_fold *f;        // the minimum folding, see structs.h
        PARAMTYPE *W;                 // the W exterior loop array
        int nb_nucleotides;     // sequence length (number of nucleotides)
        char* sequence;
        seq_interval *stack_interval;  // used for backtracking
        char *restricted;    // restricted structure given as input - restricts base pairs eg (________)    
        char *known_structure;  // known structure, used for the loss-augmented prediction

        void allocate_space();
        // allocate the necessary memory
        
        double fold_sequence ();
        double fold_sequence_restricted ();
        void insert_node (int i, int j, char type);
        
        void backtrack (seq_interval *cur_interval);
        // backtrack to retreive the MFE structure
        // PRE:  All matrixes V, VM, WM and W have been filled
        // POST: Discover the MFE path
        
        void backtrack_restricted (seq_interval *cur_interval, str_features *fres);
        // backtrack, the restricted case
        
        void compute_W (int j);
        // fill the W array
        
        PARAMTYPE compute_W_br2 (int j);
        // fill the second branch of the W array
        
        void compute_W_restricted (int j, str_features *fres);
        // fill the W array, the restricted case
        
        PARAMTYPE compute_W_br2_restricted (int j, str_features *fres, int &must_choose_this_branch);
        // fill the second branch of the W array, the restricted case
        
        void print_result ();
        // PRE:  The matrix V has been calculated and the results written in f
        // POST: Prints details of each elementary structure        

};

#endif //S_MIN_FOLDING_H
