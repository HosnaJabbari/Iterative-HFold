/***************************************************************************
                          s_sub_folding.h  -  description
                             -------------------
    begin                : Thu May 06 2005
    copyright            : (C) 2003 orginally started by Zhi Chuan Zhang, 
                           continued by Mirela Andronescu
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

// this class represents free energy minimization and backtracking for pairfold
// WITH suboptimal structures. Recurrences are from Wutchy complete suboptimal 
// folding paper. The dangling free energies for suboptimal structures are not like
// in the MFE folding, but the 5' and 3' dangling ends are added all in any case.
 
 
#ifndef S_SUB_FOLDING_H
#define S_SUB_FOLDING_H


#include <stdio.h>
#include <string.h>
#include "structs.h"
#include "s_energy_matrix.h"
#include "s_multi_loop_sub.h"


class s_sub_folding
{
    public:
 
        s_sub_folding (char *seq, PARAMTYPE en_var);
        // CONSTRUCTOR
        
        s_sub_folding (char *seq, char *restricted, PARAMTYPE en_var);
        // CONSTRUCTOR, the restricted case
        
        ~s_sub_folding ();
        // The destructor

        double s_simfold (double &enthalpy);
        // fold sequence
        
        double s_simfold_restricted (double &enthalpy);
        // fold sequence, the restricted case
        
        long get_num_partial_structures_thrown_away()  { return num_partial_structures_thrown_away; }; 
        PARAMTYPE compute_W_br2 (int j);

        void compute_W (int j);
        void compute_W_restricted (int j, str_features *fres);
        
        PARAMTYPE calculate_enthalpy (int i);
        void backtrack_hairpin(int i, int j);
        void backtrack_VBI(int i, int j);        
        void backtrack_stack(int i, int j);
        void backtrack_multi(int i, int j);
        void backtrack_freebases(int i, int j);
        void backtrack_MFM1(int i, int j);
        void backtrack_MFM(int i, int j);
        void backtrack();
        
        void backtrack_restricted (str_features *fres);              
        void backtrack_hairpin_restricted (int i, int j, str_features *fres);
        void backtrack_VBI_restricted (int i, int j, str_features *fres);
        void backtrack_freebases_restricted (int i, int j, str_features *fres);
        void backtrack_MFM1_restricted (int i, int j, str_features *fres);
        void backtrack_MFM_restricted (int i, int j, str_features *fres);

                
        struct_node* copy_struct();
        void release_struct(struct_node* sn);
        void print_result(int flag);
        void insert_node(struct_node* sn);
        void set_limit(int limit);
        void adjust();
        
        //int return_structures (char **structures, double energies[]);
        int return_structures (char structures[][MAXSLEN], double energies[]);

    // better to have protected variable rather than private, it's necessary for Hfold
    protected:
    //private:
    
        //char* structure;
        s_hairpin_loop *H;
        s_stacked_pair *S;
        s_internal_loop *VBI;
        s_multi_loop_sub *VM_sub;
        s_energy_matrix *V;  
        int* int_sequence;

        //minimum_fold *f;// the minimum folding; 
        PARAMTYPE *W;
        int nb_nucleotides;
        char* sequence;
        char *restricted;
  
        // A stack of valid structures, each represented by a list
        struct_node* folding_list;

        // M: added
        struct_node *tail_folding_list;      // the last node of the folding_list

        // Use two temporary variable cur_folding for the struct poped out
        struct_node* cur_folding;
        // Use two temporary variable cur_intevral for the interval that is being backtracked 
        seq_interval* cur_interval;

        // Keep all the valid structure so that it can be printed out
        struct_node *result_list;
        struct_node *last_list;

        // Min_energy 
        PARAMTYPE min_energy;

        // how much energy from the MFE we allow
        PARAMTYPE en_var;

        // Number of valid structures needed
        int limit;
    
        // M: added
        int num_partial_structures;

        // M: added
        int num_complete_structures;

        // M: added
        PARAMTYPE max_energy;
        // having "limit" partial structures, the energy of the last suboptimal structure (i.e. max_energy) 
        // will be at least as low as the energy of the last partial structure.

        // M: added on July 5th
        long num_partial_structures_thrown_away;

        /* Mybe need to keep pointers to V, FM, FM1, W */
        /* Best way is to put the code in energy.cpp and energy.h in this class */

        double fold_sequence (double &enthalpy);        
        
        double fold_sequence_restricted (double &enthalpy);        
        
        void allocate_space (char *seq, PARAMTYPE var);
        
        void copy_list (seq_interval* from, seq_interval *& to);
        // PRE: from is a linked list
        // POST: copy all linked list from "from" to "to"
        // M: added   


};

#endif //S_SUB_FOLDING_H
