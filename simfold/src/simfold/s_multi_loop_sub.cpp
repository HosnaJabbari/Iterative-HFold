/***************************************************************************
                          s_multi_loop_sub.cpp  -  description
                             -------------------
    begin                : Mon Apr 15 2002
    copyright            : (C) 2002 by Zhi-Chuan Zhang and Mirela Andronescu
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

#include <stdio.h>
//Hosna, March 5, 2012
// malloc.h is not needed on my mac as stdlib.h does the same
//#include <malloc.h>

#include "externs.h"
#include "common.h"
#include "s_multi_loop_sub.h"
#include "s_energy_matrix.h"

s_multi_loop_sub::s_multi_loop_sub (int *seq, int length)
// The constructor
{
    int i;
    sequence = seq;
    seqlen = length;
    this->V = NULL;

    index = new int[length];
    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    int total_length;
    total_length = (length *(length+1))/2;
    index[0] = 0;
    for (i=1; i < length; i++)
        index[i] = index[i-1]+length-i+1;

    FM = new PARAMTYPE [total_length];
    if (FM == NULL) giveup ("Cannot allocate memory", "s_multi_loop_sub");
    for (i=0; i < total_length; i++) FM[i] = INF;
    FM1 = new PARAMTYPE [total_length];
    if (FM1 == NULL) giveup ("Cannot allocate memory", "s_multi_loop_sub");
    for (i=0; i < total_length; i++) FM1[i] = INF;

}

s_multi_loop_sub::~s_multi_loop_sub ()
// The destructor
{
    delete [] index;
    delete [] FM;
    delete [] FM1;
}

void s_multi_loop_sub::compute_energy_FM1 (int j)
// compute the multi-loop branch which is right at the left of j
{
    int i, tmp_k, ij, k;
    PARAMTYPE tmp, min_e;

    tmp=INF;
    for (i=j-1; i>=0; i--)
    {
        ij = index[i]+j-i;
        tmp_k = -1;
        tmp = INF;
        min_e = INF;
        for (k=i+1; k <= j; k++) // RANGE OF K<=j-TURN-2
          {
            tmp = V->get_energy(i, k) + (j-k)*misc.multi_free_base_penalty+
              misc.multi_helix_penalty+
              AU_penalty(sequence[i], sequence[k]);

            if(i>0)
              {
                tmp += dangle_bot [sequence[k]]  // M: changed: it was i,k,i-1
                  [sequence[i]]
                  [sequence[i-1]];
              }

            if(k<seqlen-1)
              {
                tmp += dangle_top [sequence [k]]  // M: changed: it was k,i,k+1
                  [sequence [i]]
                  [sequence [k+1]];
              }
            if (tmp < min_e)
              {
                min_e = tmp;
              }
          }// k for loop
        if(min_e < INF)
          {
            FM1[ij] = min_e;
          }
    }//i for loop
}



void s_multi_loop_sub::compute_energy_FM (int j)
// computes a partial multi-loop
{
    int i, tmp_k, k, ij;
    PARAMTYPE tmp, min_e;

    for (i=j-1; i>=0; i--)
      {
        //    int ii = index[i]-i;
        ij = index[i]+j-i;
        tmp_k = -1;
        min_e = INF;

        // first branch
        for (k=i+1; k < j; k++)
          {
            tmp = get_FM_energy (i,k) + get_FM1_energy (k+1,j);
            if (tmp < min_e)
                min_e = tmp;
          }// k for loop

        // second branch
        for (k=i; k < j; k++)
          {
            tmp = get_FM1_energy (k,j) + (k-i)*misc.multi_free_base_penalty;
            if (tmp < min_e)
                min_e = tmp;
          }// k for loop

        if (min_e < INF)
          FM[ij] = min_e;
      }
}



PARAMTYPE s_multi_loop_sub::compute_energy (int i, int j)
// computes the free energy of a multi-loop closed by i and j
{
    int k;
    PARAMTYPE min_en = INF;
    PARAMTYPE tmp;
    for(k=i+1; k<j; k++)
      {
        tmp = get_FM_energy (i+1,k) + get_FM1_energy (k+1,j-1);
        if(tmp < min_en)
          {
            min_en = tmp;
          }
      }
    if(min_en < INF)
      {
        return min_en+misc.multi_offset+misc.multi_helix_penalty+
          AU_penalty (sequence[i], sequence[j])+
          dangle_top[sequence[i]][sequence[j]][sequence[i+1]]+
          dangle_bot[sequence[i]][sequence[j]][sequence[j-1]];
      }
    return INF;
}




void s_multi_loop_sub::compute_energy_FM1_restricted (int j, str_features *fres)
// compute the multi-loop branch which is right at the left of j, the restricted case
{
    int i, tmp_k, ij, k;
    PARAMTYPE tmp, min_e;

    tmp=INF;
    for (i=j-1; i>=0; i--)
    {
        ij = index[i]+j-i;
        tmp_k = -1;
        tmp = INF;
        min_e = INF;
        for (k=i+1; k <= j; k++) // RANGE OF K<=j-TURN-2
          {
            if (exists_restricted (k, j+1, fres))
                continue;
            tmp = V->get_energy(i, k) + (j-k)*misc.multi_free_base_penalty+
              misc.multi_helix_penalty+
              AU_penalty(sequence[i], sequence[k]);

            if(i>0)
              {
                tmp += dangle_bot [sequence[k]]  // M: changed: it was i,k,i-1
                  [sequence[i]]
                  [sequence[i-1]];
              }

            if(k<seqlen-1)
              {
                tmp += dangle_top [sequence [k]]  // M: changed: it was k,i,k+1
                  [sequence [i]]
                  [sequence [k+1]];
              }
            if (tmp < min_e)
              {
                min_e = tmp;
              }
          }// k for loop
        if(min_e < INF)
          {
            FM1[ij] = min_e;
          }
    }//i for loop
}


void s_multi_loop_sub::compute_energy_FM_restricted (int j, str_features *fres)
// computes a partial multi-loop, the restricted case
{
    int i, tmp_k, k, ij;
    PARAMTYPE tmp, min_e;

    for (i=j-1; i>=0; i--)
      {
        //    int ii = index[i]-i;
        ij = index[i]+j-i;
        tmp_k = -1;
        min_e = INF;

        // first branch
        for (k=i+1; k < j; k++)
          {
            tmp = get_FM_energy (i,k) + get_FM1_energy (k+1,j);
            if (tmp < min_e)
                min_e = tmp;
          }// k for loop

        // second branch
        for (k=i; k < j; k++)
          {
            if (i > 0 && exists_restricted (i-1, k, fres))
                continue;
            tmp = get_FM1_energy (k,j) + (k-i)*misc.multi_free_base_penalty;
            if (tmp < min_e)
                min_e = tmp;
          }// k for loop

        if (min_e < INF)
          FM[ij] = min_e;
      }
}

