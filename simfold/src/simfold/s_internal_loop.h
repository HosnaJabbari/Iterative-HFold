/***************************************************************************
                          s_internal_loop.h  -  description
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

#ifndef S_INTERNAL_LOOP_H
#define S_INTERNAL_LOOP_H

// a class for internal loop related functions

#include "constants.h"
#include "structs.h"

class s_energy_matrix;

class s_internal_loop
{
    public:

        friend class s_energy_matrix;

        s_internal_loop (int *seq, int length);
        // The constructor

        ~s_internal_loop ();
        // The destructor

        void set_energy_matrix (s_energy_matrix *V) { this->V = V; }
        // PRE:  None
        // POST: Set local energy matrix to V

        PARAMTYPE compute_energy (int i, int j);
        // computes the MFE of the structure closed by an internal loop closed at (i,j)

        PARAMTYPE compute_energy_restricted (int i, int j, str_features *fres);
        // computes the MFE of the structure closed by a restricted internal loop closed by (i,j)

		// Hosna, April 18, 2012
		PARAMTYPE compute_energy_restricted_pkonly (int i, int j, str_features *fres);
		// computes the MFE of the structure closed by a restricted internal loop closed by (i,j)

        PARAMTYPE get_energy_str (int i, int j, int ip, int jp);
        // returns the free energy of the structure closed by the internal loop (i,j,ip,jp)


		// Hosna, March 26, 2012
		// This function is added for non-cannonical base pairs in the restricted struture
		PARAMTYPE get_energy_str_restricted (int i, int j, int ip, int jp, str_features *fres);
		// returns the free energy of the structure closed by the internal loop (i,j,ip,jp)

        static PARAMTYPE get_energy (int i, int j, int ip, int jp, int *sequence, int *ptable=NULL);
        // returns the free energy of the internal loop closed at (i,j,ip,jp)

        static PARAMTYPE get_energy_00 (int i, int j, int ip, int jp, int *sequence);

        static PARAMTYPE get_enthalpy (int i, int j, int ip, int jp, int *sequence);
        // returns the enthalpy of the internal loop closed at (i,j,ip,jp)

        static void count_get_energy (int i, int j, int ip, int jp, int *sequence, double *counter);
        // this function is needed for parameter learning, not for folding
        // fill the counter vectro accordingly

    // better to have protected variable rather than private, it's necessary for Hfold
    protected:
    //private:
        int *sequence;                  // the entire sequence for which we compute the energy.
                                        //     Each base is converted into integer, because it's faster.
        int seqlen;                     // sequence length
        s_energy_matrix *V;             // a pointer to the free energy matrix V

        // we don't need to store the energy value(i,j), we just compute and return it
};

#include "s_energy_matrix.h"
#endif
