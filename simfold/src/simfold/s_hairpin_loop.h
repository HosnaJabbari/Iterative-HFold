/***************************************************************************
                          s_hairpin_loop.h  -  description
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

#ifndef S_HAIRPIN_LOOP_H
#define S_HAIRPIN_LOOP_H


// a class for hairpin loop related functions

class s_hairpin_loop
{
    public:
        s_hairpin_loop (char * char_seq, int *seq, int length);
        // The constructor

        ~s_hairpin_loop ();
        // The destructor

        PARAMTYPE compute_energy (int i, int j);
        // compute the free energy if this hairpin loop, closed at (i,j)

        PARAMTYPE compute_energy_restricted (int i, int j, str_features *fres);
        // compute the free energy if this is restricted to be a hairpin loop, closed at (i,j)

        static PARAMTYPE get_energy (int i, int j, int *sequence, char *csequence, int *ptable=NULL);
        // returns the free energy of the hairpin loop closed at (i,j)

        static PARAMTYPE get_enthalpy (int i, int j, int *sequence, char *csequence);
        // returns the enthalpy of the hairpin loop closed at (i,j)

        static void count_get_energy (int i, int j, int* sequence, char *csequence, double *counter);
        // PRE:  csequence is the sequence of the loop; important for tetraloops
        //       I assume i-j can pair
        // POST: Increment the counter vector accordingly


    // better to have protected variable rather than private, it's necessary for Hfold
    protected:
    //private:

        int *sequence;             // the entire sequence for which we compute the energy.
                                   //     Each base is converted into integer, because it's faster.

        char *csequence;           // the sequence a array of chars. We need it here for triloops and tetraloops

        int seqlen;                // sequence length

        // we don't need to store the energy value(i,j), we just compute and return it

};


#endif
