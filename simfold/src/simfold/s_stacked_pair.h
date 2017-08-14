/***************************************************************************
                          s_stacked_pair.h  -  description
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

// class representing a stacked pair in simfold

#ifndef S_STACKED_PAIR_H
#define S_STACKED_PAIR_H


class s_energy_matrix;

class s_stacked_pair
{
    public:

        friend class s_energy_matrix;

        s_stacked_pair (int *seq, int length);
        // The constructor

        ~s_stacked_pair ();
        // The destructor

        void set_energy_matrix (s_energy_matrix *V) { this->V = V; }
        // Set local energy matrix to V

        PARAMTYPE compute_energy (int i, int j);
        // compute the free energy of the structure closed by this stacked pair

		// Hosna, March 26, 2012
		PARAMTYPE compute_energy_restricted (int i, int j, str_features *fres);
		//  compute the free energy if this is restricted to be a stacked pair closed at (i,j)

		// Hosna, April 18, 2012
		PARAMTYPE compute_energy_restricted_pkonly (int i, int j, str_features *fres);
		//  compute the free energy if this is restricted to be a stacked pair closed at (i,j)

        static PARAMTYPE get_energy (int i, int j, int *sequence);
        // returns the free energy of the stacked pair closed at (i,j)

        static PARAMTYPE get_enthalpy (int i, int j, int *sequence);
        // returns the enthalpy of the stacked pair closed at (i,j)

        static void count_get_energy (int i, int j, int *sequence, double *counter);
        // used for parameter learning, not for folding
        // update the counter vector accordingly

    // better to have protected variable rather than private, it's necessary for Hfold
    protected:
    //private:

        int *sequence;             // the entire sequence for which we compute the energy.
                                   //     Each base is converted into integer, because it's faster.
        int seqlen;                 // sequence length
        s_energy_matrix *V;           // a pointer to the free energy matrix V

        // we don't need to store the energy value(i,j), we just compute and return it
};

#include "s_energy_matrix.h"
#endif
