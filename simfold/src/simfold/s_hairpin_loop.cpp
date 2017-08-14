 /***************************************************************************
                          s_hairpin_loop.cpp  -  description
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

// a class for hairpin loop related functions

#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "structs.h"
#include "externs.h"
#include "common.h"
#include "simfold.h"
#include "s_hairpin_loop.h"
#include "params.h"

// TODO: Not sure what is correct. If size < 3  and is restricted, if I add AU_penalty, it gives inconsistencies with free_energy_simfold_restricted.
// But I think I had a reason to add AU_penalty - maybe for ML restricted?
// Right now it does not add them.

s_hairpin_loop::s_hairpin_loop (char * char_seq, int *seq, int length)
// The constructor
{
    csequence = char_seq;
    // make sure all bases are uppercase, otherwise it's not good for triloops and tetraloops
    for (int i=0; i < length; i++)
    {
        toupper(csequence[i]);
    }
    sequence = seq;     // just refer it from where it is in memory
    seqlen = length;
}



s_hairpin_loop::~s_hairpin_loop ()
// The destructor
{
}


PARAMTYPE s_hairpin_loop::compute_energy (int i, int j)
// compute the free energy if this is a hairpin loop, closed by (i,j)
{
    PARAMTYPE energy=0;
    PARAMTYPE terminal_mismatch_energy = 0, bonus = 0, special_bonus = 0, AU_pen = 0;
    int k, is_poly_C;
    int size;
    char seq[10] = "";

    size = j-i-1;

     // TODO
     //if (size < 3)
     //    return INF;
     //return 0;

    if (size < 3)
        return INF;

    else if (size == 3)
    {
        terminal_mismatch_energy = 0;
        AU_pen = AU_penalty (sequence[i], sequence[j]);
    }
    else
    {
        terminal_mismatch_energy =
             tstackh[sequence[i]]
                    [sequence[j]]
                    [sequence[i+1]]
                    [sequence[j-1]];
    }

    if (parsi_special == T99)
    {
        // check if it is a triloop
        if (size == 3)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_triloops; k++)
            {
                if (strcmp (seq, triloop[k].seq) == 0)
                    bonus = triloop[k].energy;
            }
        }

        // check to see it is a tetraloop in tloop
        else if (size == 4)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_tloops; k++)
            {
                if (strcmp (seq, tloop[k].seq) == 0)
                    bonus = tloop[k].energy;
            }
        }
    }
    else if (parsi_special == LAVISH || parsi_special == T99_LAVISH)
    {
        // nb_special_hl is (should be) 0 if parsi_special
        if (size <= 6)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_special_hl; k++)
            {
                if (strcmp (seq, special_hl[k].seq) == 0)
                    bonus = special_hl[k].energy;
            }
        }
    }

    // special_bonus from miscloop file
    // check if we have to apply "GGG" loop special bonus
    // to come back - Vienna doesn't have it

    // if MODEL is EXTENDED and parsi_special, then misc.hairpin_GGG, misc.hairpin_c1-3 are 0
    if (i > 1)
    {
        if (sequence[i-2] == G && sequence[i-1] == G &&
            sequence[i] == G && sequence[j] == U)
            special_bonus += misc.hairpin_GGG;
    }


    // check for the special case of "poly-C" hairpin loop
    is_poly_C = 1;
    for (k=i+1; k<j; k++)
    {
        if (sequence[k] != C)
        {
            is_poly_C = 0;
            break;
        }
    }
    if (is_poly_C)
    {
        if (size == 3)
            special_bonus += misc.hairpin_c3;
        else
            special_bonus += misc.hairpin_c2 + misc.hairpin_c1 * size;
    }

    energy = penalty_by_size (size, 'H') +
             terminal_mismatch_energy + bonus + special_bonus + AU_pen;

    // add the loss
    if (pred_pairings != NULL)
    {
        pred_pairings[i] = j;   pred_pairings[j] = i;
        for (k=i+1; k < j; k++) pred_pairings[k] = -1;
        energy = energy - loss (i,j);
    }
    return energy;
}


PARAMTYPE s_hairpin_loop::compute_energy_restricted (int i, int j, str_features *fres)
// compute the free energy if this is restricted to be a hairpin loop, closed at (i,j)
// Hosna, March 27, 2012
// in restricted cases, we might have hairpins with non-canonical closing base pairs
// This method needs to take care of those cases as well
{
    PARAMTYPE energy=0;
    PARAMTYPE terminal_mismatch_energy = 0, bonus = 0, special_bonus = 0, AU_pen = 0;
    int k, is_poly_C;
    int size;
    char seq[10] = "";

    // don't allow the formation of a hairpin if there are restricted base pairs inside
    if (exists_restricted (i, j, fres))
        return INF;

    size = j-i-1;
    if (size < 3)
    {
        if (fres[i].pair == j)
            return 0;
            //return (AU_penalty (sequence[i], sequence[j]));   // not sure if we should return AU_penalty - I think we should, since we do it in get_energy and count_energy (Apr 18, 2008)

        else
            return INF;
    }
    else if (size == 3)
    {
        terminal_mismatch_energy = 0;
        AU_pen = AU_penalty (sequence[i], sequence[j]);

    }
    else
    {
        terminal_mismatch_energy =
             tstackh[sequence[i]]
                    [sequence[j]]
                    [sequence[i+1]]
                    [sequence[j-1]];

    }

    if (parsi_special == T99)
    {
        // check if it is a triloop
        if (size == 3)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_triloops; k++)
            {
                if (strcmp (seq, triloop[k].seq) == 0)
                    bonus = triloop[k].energy;
            }
        }

        // check to see it is a tetraloop in tloop
        else if (size == 4)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_tloops; k++)
            {
                if (strcmp (seq, tloop[k].seq) == 0)
                    bonus = tloop[k].energy;
            }
        }
    }
    else if (parsi_special == LAVISH || parsi_special == T99_LAVISH)
    {
        if (size <= 6)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_special_hl; k++)
            {
                if (strcmp (seq, special_hl[k].seq) == 0)
                    bonus = special_hl[k].energy;
            }
        }
    }

    // special_bonus from miscloop file
    // check if we have to apply "GGG" loop special bonus
    // to come back - Vienna doesn't have it


    if (i > 1)
    {
        if (sequence[i-2] == G && sequence[i-1] == G &&
            sequence[i] == G && sequence[j] == U)
            special_bonus += misc.hairpin_GGG;
    }


    // check for the special case of "poly-C" hairpin loop
    is_poly_C = 1;
    for (k=i+1; k<j; k++)
    {
        if (sequence[k] != C)
        {
            is_poly_C = 0;
            break;
        }
    }
    if (is_poly_C)
    {
        if (size == 3)
            special_bonus += misc.hairpin_c3;
        else
            special_bonus += misc.hairpin_c2 + misc.hairpin_c1 * size;
    }
	// Hosna, March 27, 2012
	if (fres[i].pair == j && fres[j].pair ==i && !can_pair(sequence[i],sequence[j])){
		//printf("in H(%d,%d) non-canonical: terminal_mistmatch = %d \n", i,j,terminal_mismatch_energy);
		terminal_mismatch_energy = MIN(0,terminal_mismatch_energy);

	}
	/*if (fres[i].pair == j && fres[j].pair ==i && !can_pair(sequence[i],sequence[j])){
		printf("in H(%d,%d) non-canonical: penalty by size = %d, bonus = %d, special_bonus=%d  and AU_pen\n", i,j,penalty_by_size(size, 'H'),bonus, special_bonus, AU_pen);

	}*/

    energy = penalty_by_size (size, 'H') +
             terminal_mismatch_energy + bonus + special_bonus + AU_pen;


    return energy;
}

PARAMTYPE s_hairpin_loop::get_energy (int i, int j, int* sequence, char *csequence, int *ptable)
// returns the free energy of the hairpin loop closed at (i,j)
// this function is redundant with compute_energy
// includes restricted cases
// This is also used by the partition function
{
    PARAMTYPE energy=0;
    PARAMTYPE terminal_mismatch_energy = 0, bonus = 0, special_bonus = 0, AU_pen = 0;
    int k, is_poly_C;
    int size;
    char seq[10] = "";

    // don't allow the formation of a hairpin if there are restricted base pairs inside
    if (exists_restricted_ptable (i, j, ptable))
        return INF;

    size = j-i-1;

    // TODO
//      if (size < 3)
//          return INF;
//      return 500;

    if (size < 3)
    {
        // Mar 6, 2008: I changed "return INF" into "return 0", to match count_get_energy - these are called by params.cpp:count_types, and by get_features_counts
        // Mar 21, 2008: Not sure why I did that, but it doesn't work when I compute the partition function, because if considers very small hairpin loops. They should only be considered in the restricted case.
        //return INF;
        if (ptable != NULL)     // meaning we have the restricted case
        {
            if (ptable[i] == j)
                //return (AU_penalty (sequence[i], sequence[j]));
                return 0;
            return INF;
        }
        return INF;
    }
    else if (size == 3)
    {
        terminal_mismatch_energy = 0;
        AU_pen = AU_penalty (sequence[i], sequence[j]);
    }
    else
    {
        terminal_mismatch_energy =
             tstackh[sequence[i]]
                    [sequence[j]]
                    [sequence[i+1]]
                    [sequence[j-1]];
        //printf ("IN GET: tstackh[%d][%d][%d][%d] = %Lf\n", sequence[i], sequence[j], sequence[i+1], sequence[j-1], terminal_mismatch_energy);
    }

    if (parsi_special == T99)
    {
        // check if it is a triloop
        if (size == 3)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_triloops; k++)
            {
                if (strcmp (seq, triloop[k].seq) == 0)
                    bonus = triloop[k].energy;
            }
        }

        // check to see it is a tetraloop in tloop
        else if (size == 4)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_tloops; k++)   // nb_tloops should be the same as for energy
            {
                if (strcmp (seq, tloop[k].seq) == 0)
                    bonus = tloop[k].energy;
            }
        }
    }
    else if (parsi_special == LAVISH || parsi_special == T99_LAVISH)
    {
        if (size <= 6)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_special_hl; k++)
            {
                if (strcmp (seq, special_hl[k].seq) == 0)
                    bonus = special_hl[k].energy;
            }
        }
    }

    // special_bonus from miscloop file
    // check if we have to apply "GGG" loop special bonus
    // Vienna package doesn't have it

    if (i > 1)
    {
        if (sequence[i-2] == G && sequence[i-1] == G &&
            sequence[i] == G && sequence[j] == U)
            special_bonus += misc.hairpin_GGG;
    }

    // check for the special case of "poly-C" hairpin loop
    is_poly_C = 1;
    for (k=i+1; k<j; k++)
    {
        if (sequence[k] != C)
        {
            is_poly_C = 0;
            break;
        }
    }
    if (is_poly_C)
    {
        if (size == 3)
            special_bonus += misc.hairpin_c3;
        else
            special_bonus += misc.hairpin_c2 + misc.hairpin_c1 * size;
    }
    //printf ("In GETEN, penalty=%Lf, terminal=%Lf, bonus=%Lf, specialb=%Lf, AUpen=%Lf\n", penalty_by_size (size, 'H'),
    //        terminal_mismatch_energy, bonus, special_bonus, AU_pen);

    energy = penalty_by_size (size, 'H') +
             terminal_mismatch_energy + bonus + special_bonus + AU_pen;
    //printf ("term_me = %lf\n", terminal_mismatch_energy);
    return energy;
}



// this is needed for parameter learning
void s_hairpin_loop::count_get_energy (int i, int j, int* sequence, char *csequence, double *counter)
// PRE:  csequence is the sequence of the loop; important for tetraloops
//       I assume i-j can pair
// POST: Increment the counter vector accordingly
// Mirela: Nov 23, 2003
{
    char type[100];
    int index;
    PARAMTYPE energy=0;
    PARAMTYPE terminal_mismatch_energy = 0, bonus = 0, special_bonus = 0, AU_pen = 0;
    int k, is_poly_C;
    int size;
    char seq[10] = "";

    size = j-i-1;

    // TODO: not sure this is correct
    if (size < 3)    // so if the size was < 3 (i.e. restricted hairpin loop), then just don't count anything for hairpin loop
    // actually do count AU_penalty
    {
        //count_AU_penalty (sequence[i], sequence[j], counter);
        return;
    }
    else if (size == 3)
    {
        terminal_mismatch_energy = 0;
        AU_pen = AU_penalty (sequence[i], sequence[j]);
        count_AU_penalty (sequence[i], sequence[j], counter);
        energy += AU_pen;
    }
    else
    {
        terminal_mismatch_energy =
             IGINF(tstackh[sequence[i]]
                    [sequence[j]]
                    [sequence[i+1]]
                    [sequence[j-1]]);
        if (parsi_tstackh == T99 || parsi_tstackh == LAVISH)
        {
            // We don't apply rule 1 here: if i+1 and j-1 pair, the tstackh for that one is used
            sprintf (type, "tstackh[%d][%d][%d][%d]", sequence[i], sequence[j], sequence[i+1], sequence[j-1]);
            index = structure_type_index(type);
            counter[index]++;
            //printf ("In COUNT: tstackh[%d][%d][%d][%d] = %Lf\n", sequence[i], sequence[j], sequence[i+1], sequence[j-1], tstackh[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]]);
            energy += tstackh[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]];
        }
        else if (parsi_tstackh == PARSI)
        {
            // there are only 4 parameters
            if (((sequence[i]==A || sequence[i]==G) && sequence[j]==U) || (sequence[i]==U && (sequence[j]==A || sequence[j]==G)))
            {
                index = structure_type_index ("misc.hairpin_AU_closure");
                counter[index]++;
                energy += misc.hairpin_AU_closure;
            }
            if (sequence[i+1] == A && sequence[j-1] == G)
            {
                index = structure_type_index ("misc.hairpin_AG_mismatch");
                counter[index]++;
                energy += misc.hairpin_AG_mismatch;
            }
            if (sequence[i+1] == G && sequence[j-1] == A)
            {
                index = structure_type_index ("misc.hairpin_GA_mismatch");
                counter[index]++;
                energy += misc.hairpin_GA_mismatch;
            }
            if (sequence[i+1] == U && sequence[j-1] == U)
            {
                index = structure_type_index ("misc.hairpin_UU_mismatch");
                counter[index]++;
                energy += misc.hairpin_UU_mismatch;
            }
            //printf ("Energy tstackh equiv %lf\n", energy);
        }
    }

    if (parsi_special == T99)
    {
        // check if it is a triloop
        if (size == 3)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_triloops; k++)
            {
                if (strcmp (seq, triloop[k].seq) == 0)
                {
                    bonus = triloop[k].energy;
                    sprintf (type, "triloop[%d].energy", k);
                    index = structure_type_index(type);
                    counter[index]++;
                    energy += triloop[k].energy;
                    break;
                }
            }
        }

        // check to see it is a tetraloop in tloop
        if (size == 4)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_tloops; k++)   // nb_tloops should be the same as for energy
            {
                if (strcmp (seq, tloop[k].seq) == 0)
                {
                    bonus = tloop[k].energy;
                    sprintf (type, "tloop[%d].energy", k);
                    index = structure_type_index(type);
                    counter[index]++;
                    energy += tloop[k].energy;
                    break;
                }
            }
        }
    }
    else if (parsi_special == LAVISH || parsi_special == T99_LAVISH)
    {
        if (size <= 6)
        {
            substr (csequence, i, j, seq);
            for (k=0; k < nb_special_hl; k++)
            {
                if (strcmp (seq, special_hl[k].seq) == 0)
                {
                    bonus = special_hl[k].energy;
                    sprintf (type, "special_hl[%d].energy", k);
                    index = structure_type_index(type);
                    counter[index]++;
                    energy += special_hl[k].energy;
                    break;

                }
            }
        }
    }

  // special_bonus from miscloop file
    // check if we have to apply "GGG" loop special bonus
    // to come back - Vienna doesn't have it

    if (parsi_special == T99 || parsi_special == LAVISH || parsi_special == T99_LAVISH)
    {
        if (i > 1)
        {
            if (sequence[i-2] == G && sequence[i-1] == G &&
                sequence[i] == G && sequence[j] == U)
            {
                special_bonus += misc.hairpin_GGG;
                index = structure_type_index("misc.hairpin_GGG");
                counter[index]++;
                energy += misc.hairpin_GGG;
            }
        }


        // check for the special case of "poly-C" hairpin loop
        is_poly_C = 1;
        for (k=i+1; k<j; k++)
        {

            if (sequence[k] != C)
            {
                is_poly_C = 0;
                break;
            }
        }
        if (is_poly_C)
        {
            if (size == 3)
            {
                special_bonus += misc.hairpin_c3;
                index = structure_type_index("misc.hairpin_c3");
                counter[index]++;
                energy += misc.hairpin_c3;
            }
            else
            {
                special_bonus += misc.hairpin_c2 + misc.hairpin_c1 * size;
                index = structure_type_index("misc.hairpin_c2");
                counter[index]++;
                energy += misc.hairpin_c2;
                index = structure_type_index("misc.hairpin_c1");
                counter[index]+= size;
                energy += misc.hairpin_c1 * size;
            }
        }
    }

    //printf ("In COUNT, penalty size is %Lf\n",  penalty_by_size (size, 'H'));
    energy += penalty_by_size (size, 'H');
    count_penalty_by_size (size, 'H', counter);

    PARAMTYPE energy2 = get_energy (i, j, sequence, csequence, NULL);
    if (fabs (energy/100.0-energy2/100.0) > 0.01)
    {
        printf ("ERROR! The way I compute get_energy and the way I count in s_hairpin_loop.cpp is different!\n");
#ifdef DOUBLEPARAMS
        printf ("By counts energy is %.2lf, by get_energy is %.2lf\n", energy/100.0, energy2/100.0);
#elif LDOUBLEPARAMS
        printf ("By counts energy is %.2Lf, by get_energy is %.2Lf\n", energy/100.0, energy2/100.0);
#endif
        fprintf (stderr, "Size is %d ", j-i-1);
        for (int myi = i; myi <= j ; myi++)
            fprintf (stderr, "%c", csequence[myi]);
        fprintf (stderr, "\n");
        exit(1);
    }
    else
    {
        //printf ("Counts and energy in s_hairpin_loop.cpp are equal! %.2lf\n", energy/100.0);
    }
}






PARAMTYPE s_hairpin_loop::get_enthalpy (int i, int j, int* sequence, char *csequence)
// returns the enthalpy of the hairpin loop closed at (i,j)
{
    PARAMTYPE energy=0;
    PARAMTYPE terminal_mismatch_energy = 0, bonus = 0, special_bonus = 0, AU_pen = 0;
    int k, is_poly_C;
    int size;
    char seq[10] = "";


    size = j-i-1;

    if (size < 3)
        return INF;
    else if (size == 3)
    {
        terminal_mismatch_energy = 0;
        AU_pen = AU_penalty_enthalpy (sequence[i], sequence[j]);
    }
    else
    {
        terminal_mismatch_energy =
                enthalpy_tstackh[sequence[i]]
                [sequence[j]]
                [sequence[i+1]]
                [sequence[j-1]];
    }

    // TODO MODEL
    // check if it is a triloop
    if (size == 3)
    {
        substr (csequence, i, j, seq);
        for (k=0; k < enthalpy_nb_triloops; k++)
        {
            if (strcmp (seq, enthalpy_triloop[k].seq) == 0)
                bonus = enthalpy_triloop[k].energy;
        }
    }

    // check to see it is a tetraloop in tloop
    if (size == 4)
    {
        substr (csequence, i, j, seq);
        for (k=0; k < enthalpy_nb_tloops; k++)   // nb_tloops should be the same as for energy
        {
            if (strcmp (seq, enthalpy_tloop[k].seq) == 0)
                bonus = enthalpy_tloop[k].energy;
        }
    }

    // special_bonus from miscloop file
    // check if we have to apply "GGG" loop special bonus
    // Vienna package doesn't have it

    if (i > 1)
    {
        if (sequence[i-2] == G && sequence[i-1] == G &&
            sequence[i] == G && sequence[j] == U)
            special_bonus += enthalpy_misc.hairpin_GGG;
    }


    // check for the special case of "poly-C" hairpin loop
    is_poly_C = 1;
    for (k=i+1; k<j; k++)
    {

        if (sequence[k] != C)
        {
            is_poly_C = 0;
            break;
        }
    }
    if (is_poly_C)
    {
        if (size == 3)
            special_bonus += enthalpy_misc.hairpin_c3;
        else
            special_bonus += enthalpy_misc.hairpin_c2 + enthalpy_misc.hairpin_c1 * size;
    }

    energy = penalty_by_size_enthalpy (size, 'H') +
            terminal_mismatch_energy + bonus + special_bonus + AU_pen;
    return energy;
}

