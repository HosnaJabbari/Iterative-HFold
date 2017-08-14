/***************************************************************************
                          s_energy_matrix.cpp  -  description
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

 // This is the V matrix

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
// Hosna, March 5, 2012
// malloc.h is not needed in my mac as stdlib.h does the same
//#include <malloc.h>

#include "constants.h"
#include "structs.h"
#include "externs.h"
#include "common.h"
#include "simfold.h"
#include "s_energy_matrix.h"
#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"


s_energy_matrix::s_energy_matrix (int *seq, int length)
// The constructor
{
    this->H = NULL;
    this->S = NULL;
    this->VBI = NULL;
    this->VM = NULL;

    sequence = seq;     // just refer it from where it is in memory
    seqlen = length;

    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    index = new int [length];
    int total_length = (length *(length+1))/2;
    index[0] = 0;
    for (int i=1; i < length; i++)
        index[i] = index[i-1]+length-i+1;

    // this array holds V(i,j), and what (i,j) encloses: hairpin loop, stack pair, internal loop or multi-loop
    nodes = new free_energy_node [total_length];
    if (nodes == NULL) giveup ("Cannot allocate memory", "s_energy_matrix");
}


s_energy_matrix::~s_energy_matrix ()
// The destructor
{
    delete [] index;
    delete [] nodes;
}



void s_energy_matrix::compute_energy (int i, int j)
// compute the V(i,j) value
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;

    if (can_pair (sequence[i], sequence[j]))    // if i and j can pair
    {
        // compute free energy of hairpin loop, stack pair, internal loop and multi-loop
        min_en[0] = H->compute_energy (i, j);
        if (i<=j-TURN-1)
        {
            min_en[1] = S->compute_energy (i, j);

            // TODO: uncomment
            if (!ignore_internal)
                min_en[2] = VBI->compute_energy (i, j);
            if (!ignore_multi)
                min_en[3] = VM->compute_energy (i, j);
        }
    }

    // see which of them is the minimum
    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = STACK; break;
        case  2: type = INTER; break;
        case  3: type = MULTI; break;
        default: type = NONE;
    }

    if (min_rank > -1)
    {
        if (debug)
            printf ("V(%d,%d) type %c energy %d\n", i, j, type, min);
    }

    if (min < INF/2)
    {
        int ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }
}


void s_energy_matrix::compute_energy_restricted (int i, int j, str_features *fres)
// compute the V(i,j) value, if the structure must be restricted
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;

	// Hosna, March 26, 2012
	// if the restricted base pairs are non-canonical then checking for can_pair only will cause missing those base pairs
	if (can_pair (sequence[i], sequence[j]) || (fres[i].pair == j && fres[j].pair ==i))
    {
        if (fres[i].pair == i+1)
            min_en[0] = 0;
        else
        {
            if (!exists_restricted (i, j, fres))
                min_en[0] = H->compute_energy_restricted (i, j, fres); // there was a stupid bug here, I was calling H->compute_energy instead of the restricted version. Fixed on June 30, 2007.

            min_en[1] = S->compute_energy_restricted (i, j,fres);//S->compute_energy (i, j); Hosna, March 26, 2012
            min_en[2] = VBI->compute_energy_restricted (i, j, fres);
            min_en[3] = VM->compute_energy_restricted (i, j, fres);
        }
    }

    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = STACK; break;
        case  2: type = INTER; break;
        case  3: type = MULTI; break;
        default: type = NONE;
    }

    if (min_rank > -1 && debug) {
    	printf ("V(%d,%d) type %c energy %d\n", i, j, type, min);
   	}

    if (min < INF/2) {
        int ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }
}

// Hosna, April 18, 2012
// pkonly version
void s_energy_matrix::compute_energy_restricted_pkonly (int i, int j, str_features *fres)
// compute the V(i,j) value, if the structure must be restricted
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;

	if (fres[i].pair == j && fres[j].pair ==i)
    {
        if (fres[i].pair == i+1)
            min_en[0] = 0;
        else
        {
            if (!exists_restricted (i, j, fres))
                min_en[0] = H->compute_energy_restricted (i, j, fres);

            min_en[1] = S->compute_energy_restricted_pkonly (i, j,fres);
            min_en[2] = VBI->compute_energy_restricted_pkonly (i, j, fres);
            min_en[3] = VM->compute_energy_restricted (i, j, fres); // should be left as is, Hosna April 18, 2012
        }
    }

    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = STACK; break;
        case  2: type = INTER; break;
        case  3: type = MULTI; break;
        default: type = NONE;
    }

    if (min_rank > -1)
    {
        if (debug)
            printf ("V(%d,%d) type %c energy %d\n", i, j, type, min);
    }

    if (min < INF/2)
    {
        int ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }
}


void s_energy_matrix::compute_energy_sub (int i, int j)
// suboptimals computation for V(i,j)
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;
    int ij;

    min = INF;
    min_rank = -1;


    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;

    if (can_pair (sequence[i], sequence[j]))
    {
        min_en[0] = H->compute_energy (i, j);
        if (i<=j-TURN-1)
        {
            min_en[1] = S->compute_energy (i, j);
            // TODO: uncomment
            if (!ignore_internal)
                min_en[2] = VBI->compute_energy (i, j);
            if (!ignore_multi)
                min_en[3] = VM_sub->compute_energy (i, j);
        }
    }

    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
      {
      case  0: type = HAIRP; break;
      case  1: type = STACK; break;
      case  2: type = INTER; break;
      case  3: type = MULTI; break;
      default: type = NONE;
    }

    if (min < INF/2)
    {
        ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
        //    printf ("V(%d,%d) = %d, type=%c\n", i,j, nodes[ij].energy, nodes[ij].type);
    }
}




void s_energy_matrix::compute_energy_sub_restricted (int i, int j, str_features *fres)
// compute the V(i,j) value - suboptimals and restricted
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    min_en[3] = INF;

    if (can_pair (sequence[i], sequence[j]))
    {
        min_en[0] = H->compute_energy_restricted (i, j, fres);
        min_en[1] = S->compute_energy (i, j);
        min_en[2] = VBI->compute_energy_restricted (i, j, fres);
        // I don't need restricted for VM_sub because I include dangling ends all the time
        min_en[3] = VM_sub->compute_energy (i, j);
    }

    for (k=0; k<4; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = STACK; break;
        case  2: type = INTER; break;
        case  3: type = MULTI; break;
        default: type = NONE;
    }

    if (min_rank > -1)
    {
        if (debug)
            printf ("V(%d,%d) type %c energy %d\n", i, j, type, min);
    }

    if (min < INF/2)
    {
        int ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }
}
