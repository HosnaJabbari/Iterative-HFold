/***************************************************************************
                          s_energy_matrix.h  -  description
                             -------------------
    begin                : Fri Apr 12 2002
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

// the V matrix

#ifndef ENERGY_MATRIX_H
#define ENERGY_MATRIX_H

#include "s_stacked_pair.h"
#include "s_hairpin_loop.h"
#include "s_internal_loop.h"
#include "s_multi_loop.h"
#include "s_multi_loop_sub.h"


class s_energy_matrix
{
    public:

        friend class s_stacked_pair;
        friend class s_internal_loop;
        friend class s_multi_loop;

        s_energy_matrix (int *seq, int length);
        // The constructor

        ~s_energy_matrix ();
        // The destructor

        void set_loops (s_hairpin_loop *H, s_stacked_pair *S,
                        s_internal_loop *VBI, s_multi_loop *VM, s_multi_loop_sub *VM_sub)
        // Set the local loops to the given values
        {
            this->H = H;
            this->S = S;
            this->VBI = VBI;
            this->VM = VM;
            this->VM_sub = VM_sub;
        }
        // VM_sub should be NULL if you don't want suboptimals

        void compute_energy (int i, int j);
        // compute the V(i,j) value

        void compute_energy_restricted (int i, int j, str_features *fres);

		void compute_energy_restricted_pkonly (int i, int j, str_features *fres);
		// April 18, 2012, Hosna:
		// this is the necessary change for having the pkonly variation of HFold

        void compute_energy_sub (int i, int j);
        // compute the V(i,j) value

        void compute_energy_sub_restricted (int i, int j, str_features *fres);

        free_energy_node* get_node (int i, int j) { int ij = index[i]+j-i; return &nodes[ij]; }
        // return the node at (i,j)

        // May 15, 2007. Added "if (i>=j) return INF;"  below. It was miscalculating the backtracked structure.
        PARAMTYPE get_energy (int i, int j) { if (i>=j) return INF; int ij = index[i]+j-i; return nodes[ij].energy; }
        // return the value at V(i,j)

        char get_type (int i, int j) { int ij = index[i]+j-i; return nodes[ij].type; }
        // return the type at V(i,j)


    // better to have protected variable rather than private, it's necessary for Hfold
    protected:
    //private:
        s_hairpin_loop *H;
        s_stacked_pair *S;
        s_internal_loop *VBI;
        s_multi_loop *VM;
        s_multi_loop_sub *VM_sub;

        int *sequence;             // the entire sequence for which we compute the energy.
                                   //     Each base is converted into integer, because it's faster.
        int seqlen;                // sequence length
        int *index;                // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
        free_energy_node *nodes;   // the free energy and type (i.e. base pair closing a hairpin loops, stacked pair etc), for each i and j
};



#endif
