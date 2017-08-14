/***************************************************************************
                          s_multi_loop.h  -  description
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

 // The main class for multi-loop related functions, the mfe case

#ifndef S_MULTI_LOOP_H
#define S_MULTI_LOOP_H

#include "structs.h"

class s_energy_matrix;

class s_multi_loop
{
    public:

        friend class s_energy_matrix;

        s_multi_loop (int *seq, int length);
        // The constructor

        ~s_multi_loop ();
        // The destructor

        void set_energy_matrix (s_energy_matrix *V) { this->V = V; }
        // Set local energy matrix to V

        PARAMTYPE compute_energy (int i, int j);
        // compute the MFE of a multi-loop closed at (i,j)

        PARAMTYPE compute_energy_restricted (int i, int j, str_features *fres);
        // compute the MFE of a multi-loop closed at (i,j), the restricted case

        void compute_energy_WM (int j);
        // compute de MFE of a partial multi-loop closed at (i,j)

        void compute_energy_WM_restricted (int j, str_features *fres);
        // compute de MFE of a partial multi-loop closed at (i,j), the restricted case

		// Hosna, April 18, 2012
		void compute_energy_WM_restricted_pkonly (int j, str_features *fres);
		// compute de MFE of a partial multi-loop closed at (i,j), the restricted case

        // May 15, 2007. Added "if (i>=j) return INF;"  below. It was miscalculating the backtracked structure.
        PARAMTYPE get_energy_WM (int i, int j) { if (i>=j) return INF; int ij = index[i]+j-i; return WM[ij]; }
        // returns the previously computed free energy of WM(i,j)


    // better to have protected variable rather than private, it's necessary for Hfold
    protected:
    //private:

        int *sequence;                 // the entire sequence for which we compute the energy.
                                       //     Each base is converted into integer, because it's faster.
        int seqlen;                    // sequence length

        s_energy_matrix *V;            // a pointer to the free energy matrix V

        int *index;    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
        PARAMTYPE *WM;      // WM - 2D array (actually n*(n-1)/2 long 1D array)

};

#include "s_energy_matrix.h"
#endif

