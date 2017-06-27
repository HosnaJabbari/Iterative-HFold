/***************************************************************************
                          s_multi_loop_sub.h  -  description
                             -------------------
    begin                : Mon Apr 15 2002
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

// this class represents multi-loops for suboptimal pairfold, it implements the Wuchty complete suboptimal paper 
 
#ifndef S_MULTI_LOOP_SUB_H
#define S_MULTI_LOOP_SUB_H

#include "structs.h"

class s_energy_matrix;

class s_multi_loop_sub
{
    public:

        friend class s_energy_matrix;

        s_multi_loop_sub (int *seq, int length);
        // The constructor
        
        ~s_multi_loop_sub ();
        // The destructor
        
        void set_energy_matrix (s_energy_matrix *V) { this->V = V; }
        // Set local energy matrix to V

        PARAMTYPE compute_energy (int i, int j);
        // computes the free energy of a multi-loop closed by i and j
      
        void compute_energy_FM (int j);
        // computes a partial multi-loop

        void compute_energy_FM1 (int j);
        // compute the multi-loop branch which is right at the left of j
      
        void compute_energy_FM1_restricted (int j, str_features *fres);
        // compute the multi-loop branch which is right at the left of j, the restricted case
        
        void compute_energy_FM_restricted (int j, str_features *fres);
        // computes a partial multi-loop, the restricted case        
        
        PARAMTYPE get_FM_energy (int i, int j) { if (i>=j) return INF; return FM[index[i]-i+j]; }
        // return the energy of FM(i,j)
        
        PARAMTYPE get_FM1_energy (int i, int j) { if (i>=j) return INF; return FM1[index[i]-i+j]; }
        // return the energy of FM1(i,j)
        

    // better to have protected variable rather than private, it's necessary for Hfold
    protected:
    //private:
    
        int *sequence;                 // the entire sequence for which we compute the energy. 
                                       //     Each base is converted into integer, because it's faster.
        int seqlen;                     // sequence length

        s_energy_matrix *V;               // a pointer to the free energy matrix V

        int *index;    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2

        PARAMTYPE* FM;  // 2D array to keep the energy for FM (actually 1D of length n*(n-1)/2)
        PARAMTYPE* FM1;  // 2d array to keep the energy for FM1 (actually 1D of length n*(n-1)/2)
                    // FM1 contains the rightmost multi-loop branch at this point

};

#include "s_energy_matrix.h"
#endif

