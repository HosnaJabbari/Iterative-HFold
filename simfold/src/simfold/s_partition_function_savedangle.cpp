/***************************************************************************
                          s_partition_function.cpp  -  description
                             -------------------
    begin                : Wed Mar 15 2006
    copyright            : (C) 2006 by Mirela Andronescu
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <ctype.h>

#include "constants.h"
#include "structs.h"
#include "externs.h"
#include "common.h"
#include "simfold.h"
#include "s_energy_matrix.h"
#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"
#include "s_partition_function.h"
#include "s_min_folding.h"
#include "s_sub_folding.h"
#include "params.h"

#define PREC 1.0e-13

s_partition_function::s_partition_function (char *cseq, int ignore)
// The constructor
{
    int i;
    seqlen = strlen (cseq);
    ignore_dangles = ignore;    
    
    // for the exhaustive calculations
    IFD
        no_dangling_ends = 1;
    else
        simple_dangling_ends = 1;
            
    // NOTE: valgrind doesn't like if I just do: csequence = cseq
    csequence = cseq;
    //strcpy (csequence, cseq);    // just refer it from where it is in memory
    
    for (i=0; i < seqlen; i++)
    {
        toupper(csequence[i]);
    }          
    sequence = new int[seqlen];
    if (sequence == NULL) giveup ("Cannot allocate memory", "s_partition_function");
    for (i=0; i < seqlen; i++) 
    {
        sequence[i] = nuc_to_int(csequence[i]);
    }
    
    oneoverRT = -10.0/(1.98717*310.15);
    //oneoverRT = 1000.0/(1.98717*310.15);

    // faster to compute eAU
    eAU = exp (AU_penalty(A,U)* oneoverRT);
    EXPA = exp (misc.multi_offset * oneoverRT);
    EXPB = new double [seqlen/(TURN+2)];    // there could be at most seqlen/(TURN+2) branches
    if (EXPB == NULL) giveup ("Cannot allocate memory", "s_partition_function");
    double scaled_helix_penalty = misc.multi_helix_penalty * oneoverRT; //
    for (int i=0; i < seqlen/(TURN+2); i++)
    {
        EXPB[i] = exp (i * scaled_helix_penalty);
    }
    EXPC = new double [seqlen];  
    if (EXPC == NULL) giveup ("Cannot allocate memory", "s_partition_function");
    double scaled_free_base_penalty = misc.multi_free_base_penalty * oneoverRT; //
    for (int i=0; i < seqlen; i++)
    {
        EXPC[i] = exp (i * scaled_free_base_penalty);
    }    
    
    index = new int [seqlen];
    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    int total_length = (seqlen *(seqlen+1))/2;
    index[0] = 0;
    for (int i=1; i < seqlen; i++)
        index[i] = index[i-1]+seqlen-i+1;

    up = new double [total_length];
    if (up == NULL) giveup ("Cannot allocate memory", "s_partition_function");

    upm = new double [total_length];
    if (upm == NULL) giveup ("Cannot allocate memory", "s_partition_function");   

    p = new double [total_length];
    if (p == NULL) giveup ("Cannot allocate memory", "s_partition_function");

    IFD
    {
        u = new double[total_length];
        if (u == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        u1 = new double[total_length];
        if (u1 == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        s1 = new double[total_length];
        if (s1 == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        s2 = new double[total_length];
        if (s2 == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        s3 = new double[total_length];
        if (s3 == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        pm = new double[total_length];
        if (pm == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        pm1 = new double[total_length];
        if (pm1 == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        pm2 = new double[total_length];
        if (pm2 == NULL) giveup ("Cannot allocate memory", "s_partition_function");
    }
    else
    {
        u_ip_jp = new double [total_length];
        if (u_ip_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function");
    
        u_ip_ju = new double [total_length];
        if (u_ip_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        
        u_iu_jp = new double [total_length];
        if (u_iu_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        
        u_iu_ju = new double [total_length];
        if (u_iu_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function");            
        
        s1_jp = new double [total_length];
        if (s1_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function");    
    
        s1_ju = new double [total_length];
        if (s1_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function");        
                                    
        u1_ip_jp = new double [total_length];
        if (u1_ip_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function");
    
        u1_ip_ju_jm1p = new double [total_length];
        if (u1_ip_ju_jm1p == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        
        u1_ip_ju = new double [total_length];
        if (u1_ip_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        
        u1_iu_jp = new double [total_length];
        if (u1_iu_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function");
    
        u1_iu_ju_jm1p = new double [total_length];
        if (u1_iu_ju_jm1p == NULL) giveup ("Cannot allocate memory", "s_partition_function");    
            
        u1_iu_ju = new double [total_length];
        if (u1_iu_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function");            
    
        s2_jp = new double [total_length];
        if (s2_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        
        s2_ju = new double [total_length];
        if (s2_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function");    
        
        s3_jp = new double [total_length];
        if (s3_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        
        s3_ju_jm1p = new double [total_length];
        if (s3_ju_jm1p == NULL) giveup ("Cannot allocate memory", "s_partition_function");
            
        s3_ju = new double [total_length];
        if (s3_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function");
                               
        pmnod3_needmidd3 = new double[total_length];
        if (pmnod3_needmidd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function");
    
        pmnod3_noneedmidd3 = new double[total_length];
        if (pmnod3_noneedmidd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function");    
            
        pmd3_needmidd3 = new double[total_length];
        if (pmd3_needmidd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function");
    
        pmd3_noneedmidd3 = new double[total_length];
        if (pmd3_noneedmidd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function");    
            
        pm1nod3_needendd3 = new double[total_length];
        if (pm1nod3_needendd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function");        
        
        pm1d3_needendd3 = new double[total_length];
        if (pm1d3_needendd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function");            
                
        pm2d5_needmidd5 = new double[total_length];
        if (pm2d5_needmidd5 == NULL) giveup ("Cannot allocate memory", "s_partition_function");            
    
        pm2d5_noneedmidd5 = new double[total_length];
        if (pm2d5_noneedmidd5 == NULL) giveup ("Cannot allocate memory", "s_partition_function");            
    
        pm2nod5_needmidd5 = new double[total_length];
        if (pm2nod5_needmidd5 == NULL) giveup ("Cannot allocate memory", "s_partition_function");            
    
        pm2nod5_noneedmidd5 = new double[total_length];
        if (pm2nod5_noneedmidd5 == NULL) giveup ("Cannot allocate memory", "s_partition_function");            
    }                        
    
    // if we want to go the exhaustive way, just for verification
    pexhaustive = new double[total_length];
    if (pexhaustive == NULL) giveup ("Cannot allocate memory", "s_partition_function");        

    uexhaustive = new double[total_length];
    if (uexhaustive == NULL) giveup ("Cannot allocate memory", "s_partition_function");                        
        
    upexhaustive = new double[total_length];
    if (upexhaustive == NULL) giveup ("Cannot allocate memory", "s_partition_function");                    
    
    // TODO: to erase in the final version
    num_params = create_string_params();
    
    //printf ("num_params=%d\n", num_params);
    GlogZ = new double[num_params];
    if (GlogZ == NULL) giveup ("Cannot allocate memory", "s_partition_function");                        
    
    GlogZexhaustive = new double[num_params];
    if (GlogZexhaustive == NULL) giveup ("Cannot allocate memory", "s_partition_function");
        
    initialize_arrays ();    
}


s_partition_function::~s_partition_function ()
// The destructor
{
    //delete [] csequence;
    delete [] sequence;
    delete [] index;     
    
    delete [] up;
    delete [] upm;
    delete [] p;
        
    IFD
    {
        delete [] u;
        delete [] u1;
        delete [] s1;
        delete [] s2;
        delete [] s3;
        delete [] pm;
        delete [] pm1;
        delete [] pm2;
    }
    else
    {        
        delete [] u_ip_jp;
        delete [] u_iu_jp;
        delete [] u_ip_ju;
        delete [] u_iu_ju;        
        delete [] s1_jp;
        delete [] s1_ju;    
        delete [] u1_ip_jp;
        delete [] u1_ip_ju_jm1p;
        delete [] u1_iu_jp;
        delete [] u1_iu_ju_jm1p;
        delete [] u1_ip_ju;
        delete [] u1_iu_ju;
        delete [] s2_jp;
        delete [] s2_ju;    
        delete [] s3_jp;
        delete [] s3_ju_jm1p;
        delete [] s3_ju;
        
        delete [] pmnod3_needmidd3;
        delete [] pmnod3_noneedmidd3;
        delete [] pmd3_needmidd3;
        delete [] pmd3_noneedmidd3;
        delete [] pm1nod3_needendd3;
        delete [] pm1d3_needendd3;

        delete [] pm2d5_needmidd5;
        delete [] pm2d5_noneedmidd5;
        delete [] pm2nod5_needmidd5;
        delete [] pm2nod5_noneedmidd5;        
    }
        
    // if we go the exhaustive way
    delete [] pexhaustive;
    delete [] uexhaustive;
    delete [] upexhaustive;
    
    delete [] GlogZ;
        
    delete [] GlogZexhaustive;        
}


void s_partition_function::initialize_arrays ()
{
    int i, j, ij;
    for (i = 0; i < seqlen; i++)
    {
        for (j = i; j < seqlen; j++)
        //for (j = i; j <= MIN(i + TURN, seqlen-1); j++)
        {
            ij = index[i] + j - i;
            up[ij] = 0;
            upm[ij] = 0;
            p[ij] = 0;
            
            IFD
            {
                u[ij] = 1;
                u1[ij] = 0;
                s1[ij] = 0;
                s2[ij] = 0;
                s3[ij] = 0;
                pm[ij] = 0;
                pm1[ij] = 0;
                pm2[ij] = 0;
            }
            else
            {
                u_ip_jp[ij] = 0;
                u_ip_ju[ij] = 0;
                u_iu_jp[ij] = 0;
                u_iu_ju[ij] = 1;
                s1_jp[ij] = 0;
                s1_ju[ij] = 0;
                            
                u1_ip_jp[ij] = 0;
                u1_ip_ju_jm1p[ij] = 0;
                u1_ip_ju[ij] = 0;
                u1_iu_jp[ij] = 0;
                u1_iu_ju_jm1p[ij] = 0;
                u1_iu_ju[ij] = 0;
                
                s2_jp[ij] = 0;
                s2_ju[ij] = 0;
                
                s3_jp[ij] = 0;
                s3_ju_jm1p[ij] = 0;
                s3_ju[ij] = 0;
                
                pmnod3_needmidd3[ij] = 0;
                pmnod3_noneedmidd3[ij] = 0;
                pmd3_needmidd3[ij] = 0;
                pmd3_noneedmidd3[ij] = 0;
                pm1nod3_needendd3[ij] = 0;
                pm1d3_needendd3[ij] = 0;
    
                pm2d5_needmidd5[ij] = 0;
                pm2d5_noneedmidd5[ij] = 0;
                pm2nod5_needmidd5[ij] = 0;
                pm2nod5_noneedmidd5[ij] = 0;        
            }    
            // if we go the exhaustive way
            pexhaustive[ij] = 0;
            uexhaustive[ij] = 1;
            upexhaustive[ij] = 0;
        }
    }
}


double s_partition_function::compute_partition_function ()
// the recursions are taken from Ding and Lawrence, "A statistical sampling algorithm for RNA secondary structure prediction", NAR 2003
{ 
    int i, j, ij;
    for (j=TURN+1; j < seqlen; j++)
    {
        for (i=j-TURN-1; i>=0; i--)
        {
            ij = index[i] + j - i;

            IFD
            {
                compute_upm_nodangles (i, j);    // doesn't matter where it is, all dependencies have been computed at previous steps
                compute_up (i, j);    // must be after upm            
                compute_s1 (i,j);   
                compute_u (i,j);    // must be after s1            
                compute_u1 (i,j);
                compute_s2 (i,j);
                compute_s3 (i,j);
            }
            else
            {            
                compute_upm(i, j);    // doesn't matter where it is, all dependencies have been computed at previous steps
                compute_up (i, j);    // must be after upm                
                compute_u_ip_jp (i, j);    // must be after up
                compute_u_ip_ju (i, j);    // must be after up
                compute_u_iu_jp (i, j);    // must be after up
                compute_u_iu_ju (i, j);    // must be after up
                compute_s1_jp (i, j);    // must be after up
                compute_s1_ju (i, j);    // must be after up
                            
                compute_u1_ip_jp (i, j);    // must be after up
                compute_u1_ip_ju_jm1p (i, j);    // must be after up
                compute_u1_ip_ju (i, j);    // must be after up
                compute_u1_iu_jp (i, j);    // must be after up
                compute_u1_iu_ju_jm1p (i, j);    // must be after up
                compute_u1_iu_ju (i, j);    // must be after up
                
                compute_s2_jp (i, j);    // must be after up
                compute_s2_ju (i, j);    // must be after up
                
                compute_s3_jp (i, j);    // must be after up
                compute_s3_ju_jm1p (i, j);    // must be after up
                compute_s3_ju (i, j);    // must be after up
                
                //printf ("u[%d,%d] = %g\n", i, j, u_ip_jp[ij] + u_ip_ju[ij] + u_iu_jp[ij] + u_iu_ju[ij]);
    //                 printf ("u_ip_jp[%d,%d] = %g\n", i, j, u_ip_jp[ij]);
    //                 printf ("u_ip_ju[%d,%d] = %g\n", i, j, u_ip_ju[ij]);
    //                 printf ("u_iu_jp[%d,%d] = %g\n", i, j, u_iu_jp[ij]);
    //                 printf ("u_iu_ju[%d,%d] = %g\n", i, j, u_iu_ju[ij]);
            }
        }
    }    

    int firstlast = index[0]+seqlen-1;                    
//     Z = u[firstlast];
//     return u[firstlast];
    IFD
    {
        Z = u[firstlast];
    }
    else
    {
        Z = u_ip_jp[firstlast] + u_ip_ju[firstlast] + u_iu_jp[firstlast] + u_iu_ju[firstlast];
    }
    return Z;
}


double s_partition_function::exp_AUpenalty (int i, int j)
{    
    //double AUpen = AU_penalty (sequence[i], sequence[j]);
    //return (double)(exp (AUpen* oneoverRT));
    // Note: doing exp(sum log) is much slower than doing exp(product)

    // this way is faster than computing exp(AUpen*oneoverRT) every time
    if (has_AU_penalty(sequence[i], sequence[j]))
        return eAU;
    return 1.0;
}


double s_partition_function::exp_dangle5 (int i, int j, int k)
// dangle_bot
{    
    // added the next lines on Oct 7, 2006
    // make sure we are not out of bounds
    if (i < 0 || i >= seqlen || j < 0 || j >= seqlen || k < 0 || k >= seqlen)
        return 1.0;    
    double dang = dangle_bot[sequence[i]][sequence[j]][sequence[k]];
    return (double) (exp (dang * oneoverRT));
}


double s_partition_function::exp_dangle3 (int i, int j, int k)
// dangle_top
{
    if (i < 0 || i >= seqlen || j < 0 || j >= seqlen || k < 0 || k >= seqlen)
        return 1.0;        
    double dang = dangle_top[sequence[i]][sequence[j]][sequence[k]];
    return (double) (exp (dang * oneoverRT));
}


void s_partition_function::compute_u (int i, int j)
// called only when no dangling ends
{
    int ij = index[i]+j-i;
    int hj, il, lp2j, lp1j;
    int h, l;

    // we don't need separate cases to compute u, because we don't care about dangling ends

    for (h = i; h < j; h++)                 // case ...(...)...---
    {
        hj = index[h]+j-h;
        u[ij] += s1[hj];
    }
}


void s_partition_function::compute_s1 (int h, int j)
// called only when no dangling ends
{
    // if (h < 1)  return;  // we don't want this here, because we don't care about dangling ends
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s1[hj] = 0;

    // I split it into 2 for loops to avoid the if inside. Maybe this is faster.
    for (l = h+1; l < j-2; l++)        // (...)---
    {
        hl = index[h] + l - h;
        lp1j = index[l+1]+j-l-1;
        s1[hj] += up[hl] * exp_AUpenalty (h, l) * u[lp1j];
    }
    
    for (l = j-2; l <= j; l++)        // (...)---
    {
        hl = index[h] + l - h;
        s1[hj] += up[hl] * exp_AUpenalty (h, l);
    }


//     for (l = h+1; l < j; l++)        // .(...)---
//     {
//         hl = index[h] + l - h;
//         if (l+2 < j)
//         {
//             lp1j = index[l+1]+j-l-1;
//             s1[hj] += up[hl] * exp_AUpenalty (h, l) * u[lp1j];
//         }
//         else
//         {
//             s1[hj] += up[hl] * exp_AUpenalty (h, l);
//         }
//     }
    
}


void s_partition_function::compute_u1 (int i, int j)
// called only when no dangling ends
// contains at least one branch of a multi-loop
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1[ij] = 0;

    for (h=i; h <= j-1; h++)    // ...(...)---
    {
        hj = index[h]+j-h;
        u1[ij] += EXPC[h-i] * s3[hj];
    }
    u1[ij] *= EXPB[1];
}

void s_partition_function::compute_s3 (int h, int j)
// called only when no dangling ends
// must contain at least one branch
// s3 doesn't contain the helix penalty
{
    if (h < 1) return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s3[hj] = 0;   
    
    for (l = h+1; l <= j; l++)
    {
        hl = index[h] + l - h;        
        
        if (l+2 < j)                  
        {
            lp1j = index[l+1]+j-l-1;
            s3[hj] += up[hl] * exp_AUpenalty (h, l) * (u1[lp1j] + EXPC[j-l]);
        }
        else
        {
            s3[hj] += up[hl] * exp_AUpenalty (h, l) * EXPC[j-l];
        }
    }
}


void s_partition_function::compute_upm_nodangles (int i, int j)
// i and j close a multi-loop
{
    int ij = index[i]+j-i;
    int l, h;
    int ip1l, ip2l, lp2jm1, lp1jm1, lp2jm2, lp1jm2, lp2jm3, lp1jm3, hj;
    double upm_temp;
    
    if (!can_pair(sequence[i], sequence[j]))
        return;
    upm[ij] = 0;
    // took this out for speed
    double upm_common = exp_AUpenalty (i,j) * EXPA * EXPB[2];

    /*
    // I don't think we need all these cases, since we don't care about dangling ends
    // there must be at least one more branch (i.e. 5 nucleotides) from l+1 to j-1: l+1+TURN <= j-1
    for (l=i+2; l < j-TURN-2; l++)    // case ((...)--(--)-)
    {
        ip1l = index[i+1]+l-i-1;
        lp1jm1 = index[l+1]+j-1-l-1;

        upm[ij] += up[ip1l] * exp_AUpenalty (i+1,l) * u1[lp1jm1];
    }
    for (l=i+3; l < j-TURN-2; l++)    // case [.(...)--(--)-]
    {
        ip2l = index[i+2]+l-i-2;
        lp1jm1 = index[l+1]+j-1-l-1;
        upm[ij] += up[ip2l] * EXPC[1] * exp_AUpenalty (i+2, l) * u1[lp1jm1];
    }
    */
    
    for (h=i+1; h < j-TURN-2; h++)
    //for (h=i+3; h < j-TURN-2; h++)    // case (....(...)--(--)-)
    {
        hj = index[h]+j-h;
        int hjm1 = index[h] + j-1 -h;
        upm[ij] += EXPC[h-i-1] * s2[hjm1];
    }
    upm[ij] *= upm_common;
    //printf ("i=%d, j=%d, upm=%g\n", i, j, upm[ij]);
}

void s_partition_function::compute_s2 (int h, int j)
// called only when no dangling ends
// must contain at least 2 branches
// helper for computing upm, which closes a multi-loop
// modified from Ding and Lawrence, to be able to add the right most d5 dangling end: (--(...)..)
// has at least 2 branches, and h is paired with some l in between h and j
{
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s2[hj] = 0;
    
    for (l = h+1; l < j-3; l++)    // .(...)
    {
        hl = index[h] + l - h;
        lp1j = index[l+1]+j-l-1;
        s2[hj] += up[hl] * exp_AUpenalty (h, l) * u1[lp1j];
    }
}



void s_partition_function::compute_u_ip_jp (int i, int j)
// i paired
// j paired or j-1 paired
{
    int ij = index[i]+j-i;
    int hj, il, lp2j, lp1j;
    int h, l;    
    int ijm1 = index[i]+j-1-i;

    u_ip_jp[ij] = up[ij] * exp_AUpenalty (i,j);    // case (...)
    u_ip_jp[ij] += up[ijm1] * exp_AUpenalty (i, j-1) * exp_dangle3 (j-1, i, j);    // case (...).
    
    for (l = i+1; l < j-2; l++)                    // case (...)-(--)
    {
        il = index[i]+l-i;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        u_ip_jp[ij] += up[il] * exp_AUpenalty (i, l) * (u_ip_jp[lp1j] + exp_dangle3 (l, i, l+1) * (u_ip_jp[lp2j] + u_iu_jp[lp2j]));
        //if (u_ip_jp[lp2j] != 1)            
        //    u_ip_jp[ij] += up[il] * exp_AUpenalty (i, l) * exp_dangle3 (l, i, l+1) * u_ip_jp[lp2j];
        //if (u_iu_jp[lp2j] != 1)
        //    u_ip_jp[ij] += up[il] * exp_AUpenalty (i, l) * exp_dangle3 (l, i, l+1) * u_iu_jp[lp2j];
    }                                   
}


void s_partition_function::compute_u_ip_ju (int i, int j)
// i paired
// j must be unpaired and j-1 must be unpaired
{
    int ij = index[i]+j-i;
    int hj, il, lp2j, lp1j;
    int h, l;    
    
    u_ip_ju[ij] = 0;
    
    
    l=j-2;
    il = index[i]+l-i;
    u_ip_ju[ij] += up[il] * exp_AUpenalty (i, l) * exp_dangle3 (l, i, l+1);
    
    for (l = i+1; l < j-2; l++)                    // case (...)...---
    {
        il = index[i]+l-i;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;    
        u_ip_ju[ij] += up[il] * exp_AUpenalty (i, l) *
            (u_ip_ju[lp1j] + exp_dangle3 (l, i, l+1) *(u_ip_ju[lp2j]+u_iu_ju[lp2j]));   // put them together to be faster        
    }                                   
    
}


void s_partition_function::compute_u_iu_jp (int i, int j)
// i unpaired
// j paired or j-1 paired
{
    int ij = index[i]+j-i;
    int hj, il, lp2j, lp1j, hjm1;
    int h, l;    
    
    u_iu_jp[ij] = 0;
    //for (h = i+1; h < j-1; h++)                    // case ...(...) or ...(...).
    //{
    //    hj = index[h]+j-h;
    //    hjm1 = index[h]+j-1-h;
        //u_iu_jp[ij] += up[hj] * exp_dangle5 (j, h, h-1) * exp_AUpenalty (h, j);
        //u_iu_jp[ij] += up[hjm1] * exp_dangle5 (j-1, h, h-1) * exp_AUpenalty (h, j-1) * exp_dangle3(j-1,h,j);
    //}                    
    for (h = i+1; h < j-1; h++)                 // case ...(...)...---
    {
        hj = index[h]+j-h;
        u_iu_jp[ij] += s1_jp[hj];
    }
}

void s_partition_function::compute_u_iu_ju (int i, int j)
// i unpaired
// j unpaired and j-1 unpaired
{
    int ij = index[i]+j-i;
    int hj, il, lp2j, lp1j;
    int h, l;    
    
    u_iu_ju[ij] = 1;
    for (h = i+1; h < j-1; h++)                 // case ...(...)...---
    {
        hj = index[h]+j-h;
        u_iu_ju[ij] += s1_ju[hj];                      
    }
}


void s_partition_function::compute_s1_jp (int h, int j)
// j paired or j-1 paired
{
    if (h < 1)  return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    int hjm1 = index[h]+j-1-h;
    s1_jp[hj] = 0;
    
    // case .(...)
    s1_jp[hj] += up[hj] * exp_dangle5 (j, h, h-1) * exp_AUpenalty (h, j);        
    // case .(...).
    s1_jp[hj] += up[hjm1] * exp_dangle5 (j-1, h, h-1) * exp_AUpenalty (h, j-1) * exp_dangle3 (j-1, h, j);        
    
    for (l = h+1; l < j-2; l++)        // .(...)---
    {
        hl = index[h] + l - h;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        // put them together to be faster
        s1_jp[hj] += up[hl] * exp_dangle5(l,h,h-1) * exp_AUpenalty(h,l) *
            (u_ip_jp[lp1j] + exp_dangle3(l,h,l+1) * (u_ip_jp[lp2j] + u_iu_jp[lp2j]));
    }
}


void s_partition_function::compute_s1_ju (int h, int j)
// j unpaired and j-1 unpaired
{
    if (h < 1)  return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s1_ju[hj] = 0;
    
    l=j-2;
    hl = index[h] + l - h;
    s1_ju[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) * exp_dangle3 (l, h, l+1);
    
    for (l = h+1; l < j-2; l++)        // .(...)---
    {
        hl = index[h] + l - h;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        // put them together to be faster
        s1_ju[hj] += up[hl] * exp_dangle5(l,h,h-1) * exp_AUpenalty(h,l) *
            (u_ip_ju[lp1j] + exp_dangle3(l,h,l+1) * (u_ip_ju[lp2j] + u_iu_ju[lp2j]));      
    }
}



void s_partition_function::compute_up (int i, int j)
{
    int ij = index[i]+j-i;
    int ip1jm1 = index[i+1]+j-1-i-1;
    int ip, jp, minq;
    int en_hairpin, en_stack, en_internal;

    up[ij] = 0;
    
    if (! can_pair (sequence[i], sequence[j]))
    {
        return;
    }    
    
    // hairpin loop
    en_hairpin = s_hairpin_loop::get_energy (i, j, sequence, csequence);
    /*
    // weirdly, if this has the GGG hairpin, it's added here, although it shouldn't
    // so try to remove it in case it's added
    //    that would only give me equality between the exhaustive and the true u's and up's
    if (j-i-1 == 3 && i > 1)
    {
        if (sequence[i-2] == G && sequence[i-1] == G &&
            sequence[i] == G && sequence[j] == U)
            en_hairpin -= misc.hairpin_GGG;
    }
    */
    
    if (en_hairpin >= INF/2)
    {
        //printf ("** Infinite hairpin (%d, %d) !\n", i, j);
    }
    else
        up[ij] += exp (en_hairpin * oneoverRT);                   
        
//     if (i==0 && j==8)    printf ("r1 up[0,8] = %g\n", up[ij]);
    // stack pair
    if (can_pair (sequence[i+1], sequence[j-1]))    
    {
        en_stack = s_stacked_pair::get_energy (i, j, sequence);    
        if (en_stack >= INF/2)
        {
            //printf ("** Infinite stack   (%d, %d) !\n", i, j);
        }
        else            
            up[ij] += exp (en_stack * oneoverRT) * up [ip1jm1];
    }             
//     if (i==0 && j==8)    printf ("r2 up[0,8] = %g\n", up[ij]);
             
    // internal loop/bulge
    for (ip = i+1; ip <= MIN(j-2-TURN,i+MAXLOOP+1) ; ip++)  // j-2-TURN
    {
        minq = MAX (j-i+ip-MAXLOOP-2, ip+1+TURN);    // ip+1+TURN);
        for (jp = minq; jp < j; jp++)
        {        
            if (sequence[ip]+sequence[jp] == 3 ||
                sequence[ip]+sequence[jp] == 5)        
            {
                if (ip == i+1 && jp == j-1) continue;    // we don't want stacked pairs here
                int ipjp = index[ip] + jp - ip;
                en_internal = s_internal_loop::get_energy (i, j, ip, jp, sequence);
                if (en_internal >= INF/2)
                {
                    //printf ("** Infinite internal(%d, %d) !\n", i, j);
                }
                else
                {
                    up[ij] += exp (en_internal * oneoverRT) * up[ipjp];
//                     if (i==0 && j==8)    printf ("r3, ip=%d, jp=%d, up[0,8] = %g, en_internal=%d\n", ip, jp, up[ij], en_internal);
                }
            }
        }        
    }    
    
    // multi-loop
    up[ij] += upm[ij];
}


void s_partition_function::compute_upm (int i, int j)
// this assumes we include dangling ends. For the no dangling end version, see compute_upm_nodangles
// i and j close a multi-loop
{
    int ij = index[i]+j-i;
    int l, h;
    int ip1l, ip2l, lp2jm1, lp1jm1, lp2jm2, lp1jm2, lp2jm3, lp1jm3, hj;
    double upm_temp;
    
    if (!can_pair(sequence[i], sequence[j]))
        return;
    upm[ij] = 0;
    // took this out for speed
    double upm_common = exp_AUpenalty (i,j) * EXPA * EXPB[2];
    // there must be at least one more branch (i.e. 5 nucleotides) from l+1 to j-1: l+1+TURN <= j-1
    for (l=i+2; l < j-TURN-2; l++)    // case ((...)--(--)-)
    {
        ip1l = index[i+1]+l-i-1;
        lp2jm1 = index[l+2]+j-1-l-2;
        lp1jm1 = index[l+1]+j-1-l-1;

        upm[ij] += up[ip1l] * exp_AUpenalty (i+1,l) *
                    (  u1_ip_jp[lp1jm1]      // [(...)(.-..)] or [(...)(.-..).]                    
                       + exp_dangle3 (l, i+1, l+1) * EXPC[1] *
                            ( u1_ip_jp[lp2jm1] + u1_iu_jp[lp2jm1] +     // [(...).-(...)] or [(...).-(...).] 
                                exp_dangle5(i,j,j-1) * (u1_ip_ju[lp2jm1] + u1_iu_ju[lp2jm1])) +    // [(...).-(...)-..]
                       + exp_dangle5(i,j,j-1) * u1_ip_ju[lp1jm1]);    // [(...)(...)-..]
    }
    for (l=i+3; l < j-TURN-2; l++)    // case [.(...)--(--)-]
    {
        ip2l = index[i+2]+l-i-2;
        lp2jm1 = index[l+2]+j-1-l-2;
        lp1jm1 = index[l+1]+j-1-l-1;
        lp2jm2 = index[l+2]+j-2-l-2;
        lp1jm2 = index[l+1]+j-2-l-1;        
        lp2jm3 = index[l+2]+j-3-l-2;
        lp1jm3 = index[l+1]+j-3-l-1;      
        
        upm[ij] += up[ip2l] * EXPC[1] *
                    exp_dangle3 (i, j, i+1) * exp_AUpenalty (i+2, l) * 
                    (  u1_ip_jp[lp1jm1]      // [.(...)(.-..)] or [.(...)(.-..).]                    
                       + exp_dangle3 (l, i+2, l+1) * EXPC[1] *
                            ( u1_ip_jp[lp2jm1] + u1_iu_jp[lp2jm1] +     // [.(...).-(...)] or [.(...).-(...).] 
                                exp_dangle5(i,j,j-1) * (u1_ip_ju[lp2jm1] + u1_iu_ju[lp2jm1])) +    // [.(...).-(...)-..]
                       + exp_dangle5(i,j,j-1) * u1_ip_ju[lp1jm1]);    // [.(...)(...)-..] 
    }
    upm_temp = 0;
    for (h=i+3; h < j-TURN-2; h++)    // case (....(...)--(--)-)
    {
        hj = index[h]+j-h;
        int hjm1 = index[h] + j-1 -h;
        int hjm2 = index[h] + j-2 -h;
        int hjm3 = index[h] + j-3 -h;
            
        upm_temp += EXPC[h-i-1] *
                ( s2_jp[hjm1]    // --(...)) or --(...).) 
                  + s2_ju[hjm1]* exp_dangle5(i,j,j-1));    // --(...)..)
    }   
    upm[ij] += exp_dangle3 (i, j, i+1) * upm_temp;
    upm[ij] *= upm_common;
    //printf ("i=%d, j=%d, upm=%g\n", i, j, upm[ij]);
}


void s_partition_function::compute_s2_jp (int h, int j)
// helper for computing upm, which closes a multi-loop
// modified from Ding and Lawrence, to be able to add the right most d5 dangling end: (--(...)..)
// has at least 2 branches, and h is paired with some l in between h and j
{
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s2_jp[hj] = 0;
    
    for (l = h+1; l < j-3; l++)    // .(...)
    {
        hl = index[h] + l - h;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        s2_jp[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) * 
                ( u1_ip_jp[lp1j] + exp_dangle3 (l, h, l+1) * EXPC[1]*(u1_ip_jp[lp2j] + u1_iu_jp[lp2j]));
    }
}



void s_partition_function::compute_s2_ju (int h, int j)
// helper for computing upm, which closes a multi-loop
// modified from Ding and Lawrence, to be able to add the right most d5 dangling end: (--(...)..)
// has at least 2 branches, and h is paired with some l in between h and j
{
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s2_ju[hj] = 0;
    
    for (l = h+1; l < j-3; l++)    // .(...)
    {
        hl = index[h] + l - h;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        s2_ju[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) * 
                ( u1_ip_ju[lp1j] + exp_dangle3 (l, h, l+1) * EXPC[1]*(u1_ip_ju[lp2j] + u1_iu_ju[lp2j]));
    }
}


void s_partition_function::compute_u1_ip_jp (int i, int j)
// contains at least one branch of a multi-loop
// i must be paired
// j must be paired, or j-1 must be paired
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_ip_jp[ij] = 0;


            
    for (l=j-1; l <= j; l++)    // (...) or (...).
    {
        il = index[i]+l-i;        
        u1_ip_jp[ij] += up[il] * EXPB[1] * exp_AUpenalty(i,l) * ( fd3(j+1,i,l) * EXPC[j-l] );
    }
                
    for (l=i+1; l < j-2; l++)    // (...)-(---) or (...)-(---).
    {
        lp1j = index[l+1]+j-l-1;
        lp2j = index[l+2]+j-l-2;                      
        il = index[i]+l-i;        
        u1_ip_jp[ij] += up[il] * EXPB[1] * exp_AUpenalty(i,l) *
                    (u1_ip_jp[lp1j] + exp_dangle3 (l,i,l+1)*EXPC[1]*(u1_ip_jp[lp2j] + u1_iu_jp[lp2j]));
    }
}


void s_partition_function::compute_u1_ip_ju_jm1p (int i, int j)
// contains at least one branch of a multi-loop
// i must be paired
// j must be paired, or j-1 must be paired
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_ip_ju_jm1p[ij] = 0;
    int ijm1 = index[i]+j-1-i;       
    u1_ip_ju_jm1p[ij] += up[ijm1] * EXPB[1] * exp_AUpenalty(i,j-1)* fd3(j+1,i,j-1) * EXPC[1];
                
    for (l=i+1; l < j-2; l++)    // (...)-(---).
    {
        lp1j = index[l+1]+j-l-1;
        lp2j = index[l+2]+j-l-2;                      
        il = index[i]+l-i;        
        u1_ip_ju_jm1p[ij] += up[il] * EXPB[1] * exp_AUpenalty(i,l) *
                    (u1_ip_ju_jm1p[lp1j] + exp_dangle3 (l,i,l+1)*EXPC[1]*(u1_ip_ju_jm1p[lp2j] + u1_iu_ju_jm1p[lp2j]));
    }
}


void s_partition_function::compute_u1_ip_ju (int i, int j)
// contains at least one branch of a multi-loop
// i must be paired
// j must be unpaired, and j-1 must be unpaired
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_ip_ju[ij] = 0;

    for (l=i+1; l <= j-2; l++)    // (...)----
    {
        il = index[i]+l-i;        // (...)....
        u1_ip_ju[ij] += up[il] * EXPB[1] * exp_AUpenalty(i,l) *
                ( fd3(j+1,i,l) * EXPC[j-l] );
                  
        if (l+2 < j)            // (...)-(--)-
        {          
            lp1j = index[l+1]+j-l-1;
            lp2j = index[l+2]+j-l-2;                  
            u1_ip_ju[ij] += up[il] * EXPB[1] * exp_AUpenalty(i,l) *
                        (u1_ip_ju[lp1j] + exp_dangle3 (l,i,l+1)*EXPC[1]*(u1_ip_ju[lp2j] + u1_iu_ju[lp2j]));
        }                        
    }
}

    
void s_partition_function::compute_u1_iu_jp (int i, int j)
// contains at least one branch of a multi-loop
// i must be unpaired
// j must be paired, or j-1 must be paired
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_iu_jp[ij] = 0;

    for (l=j-1; l <= j; l++)    // .(...) or .(...).
    {
        ip1l = index[i+1]+l-i-1;        
        u1_iu_jp[ij] += up[ip1l] * EXPB[1] * exp_AUpenalty(i+1,l) * exp_dangle5 (l, i+1, i) * EXPC[1] *
                ( fd3(j+1,i+1,l) * EXPC[j-l] );
    }
    
    for (l=i+2; l < j-2; l++)    // .(...)-(---) or .(...)-(---).
    {
        ip1l = index[i+1]+l-i-1;
        lp1j = index[l+1]+j-l-1;
        lp2j = index[l+2]+j-l-2;  
        u1_iu_jp[ij] += up[ip1l] * EXPB[1] * exp_AUpenalty(i+1,l) * exp_dangle5 (l, i+1, i) * EXPC[1] *
                    (u1_ip_jp[lp1j] + exp_dangle3 (l,i+1,l+1)*EXPC[1]*(u1_ip_jp[lp2j] + u1_iu_jp[lp2j]));
    }

    for (h=i+2; h <= j-1; h++)    // ..-(...)-(--) or ..-(...)-(--).
    {
        hj = index[h]+j-h;
        // d5 is included in s3, but helix penalty is not 
        u1_iu_jp[ij] += EXPB[1] * EXPC[h-i] * s3_jp[hj];
    }
}


void s_partition_function::compute_u1_iu_ju_jm1p (int i, int j)
// contains at least one branch of a multi-loop
// i must be unpaired
// j must be paired, or j-1 must be paired
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_iu_ju_jm1p[ij] = 0;

    int ip1jm1 = index[i+1]+j-1-i-1;        // .(...).
    u1_iu_ju_jm1p[ij] += up[ip1jm1] * EXPB[1] * exp_AUpenalty(i+1,j-1) * exp_dangle5 (j-1, i+1, i) * fd3(j+1,i+1,j-1) * EXPC[2];
    
    for (l=i+2; l < j-2; l++)    // .(...)-(---) or .(...)-(---).
    {
        ip1l = index[i+1]+l-i-1;
        lp1j = index[l+1]+j-l-1;
        lp2j = index[l+2]+j-l-2;  
        u1_iu_ju_jm1p[ij] += up[ip1l] * EXPB[1] * exp_AUpenalty(i+1,l) * exp_dangle5 (l, i+1, i) * EXPC[1] *
                    (u1_ip_ju_jm1p[lp1j] + exp_dangle3 (l,i+1,l+1)*EXPC[1]*(u1_ip_ju_jm1p[lp2j] + u1_iu_ju_jm1p[lp2j]));
    }

    for (h=i+2; h <= j-1; h++)    // ..-(...)-(--) or ..-(...)-(--).
    {
        hj = index[h]+j-h;
        // d5 is included in s3, but helix penalty is not 
        u1_iu_ju_jm1p[ij] += EXPB[1] * EXPC[h-i] * s3_ju_jm1p[hj];
    }
}



void s_partition_function::compute_u1_iu_ju (int i, int j)
// contains at least one branch of a multi-loop
// i must be unpaired
// j must be unpaired, and j-1 must be unpaired
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_iu_ju[ij] = 0;

    for (l=i+2; l < j-1; l++)    // .(...)-..
    {
        ip1l = index[i+1]+l-i-1;        
        u1_iu_ju[ij] += up[ip1l] * EXPB[1] * exp_AUpenalty(i+1,l) * exp_dangle5 (l, i+1, i) * EXPC[1] *
                ( fd3(j+1,i+1,l) * EXPC[j-l] );
    }
        
    for (l=i+2; l < j-2; l++)    // .(...)-(--)-..
    {
        ip1l = index[i+1]+l-i-1;
        lp1j = index[l+1]+j-l-1;
        lp2j = index[l+2]+j-l-2;                 
        u1_iu_ju[ij] += up[ip1l] * EXPB[1] * exp_AUpenalty(i+1,l) * exp_dangle5 (l, i+1, i) * EXPC[1] *
                    (u1_ip_ju[lp1j] + exp_dangle3 (l,i+1,l+1)*EXPC[1]*(u1_ip_ju[lp2j] + u1_iu_ju[lp2j]));
    }

    for (h=i+2; h <= j-1; h++)    // ..-(...)-(--)-..
    {
        hj = index[h]+j-h;
        // d5 is included in s3, but helix penalty is not 
        u1_iu_ju[ij] += EXPB[1] * EXPC[h-i] * s3_ju[hj];
    }
}



double s_partition_function::fd3 (int jplus1, int h, int l)
{
    if (l > jplus1-1)    giveup ("Error, l > j", "f function, partition_function");
    if (l == jplus1-1)
        return 1.0;
    return exp_dangle3 (l, h, l+1);    
}


void s_partition_function::compute_s3_jp (int h, int j)
// s3 doesn't contain the helix penalty
// j paired or j-1 paired
{
    if (h < 1) return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s3_jp[hj] = 0;   

    for (l = j-1; l <= j; l++)
    {
        hl = index[h] + l - h;
        s3_jp[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) * 
                    ( fd3 (j+1, h, l) * EXPC[j-l] );
    }
    for (l = h+1; l < j-2; l++)
    {
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        hl = index[h] + l - h;
        s3_jp[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) *
                (u1_ip_jp[lp1j] + exp_dangle3 (l,h,l+1)*EXPC[1]*(u1_ip_jp[lp2j] + u1_iu_jp[lp2j]));
    }
}


void s_partition_function::compute_s3_ju_jm1p (int h, int j)
// s3 doesn't contain the helix penalty
// j paired or j-1 paired
{
    if (h < 1) return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s3_ju_jm1p[hj] = 0;   

    int hjm1 = index[h]+j-1-h;
    s3_ju_jm1p[hj] += up[hjm1] * exp_dangle5 (j-1, h, h-1) * exp_AUpenalty (h, j-1) * fd3 (j+1, h, j-1) * EXPC[1];

    for (l = h+1; l < j-2; l++)
    {
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        hl = index[h] + l - h;
        s3_ju_jm1p[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) *
                (u1_ip_ju_jm1p[lp1j] + exp_dangle3 (l,h,l+1)*EXPC[1]*(u1_ip_ju_jm1p[lp2j] + u1_iu_ju_jm1p[lp2j]));
    }
}


void s_partition_function::compute_s3_ju (int h, int j)
// s3 doesn't contain the helix penalty
// j unpaired and j-1 unpaired
{
    if (h < 1) return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s3_ju[hj] = 0;   
    
    for (l = h+1; l < j-1; l++)
    {
        hl = index[h] + l - h;
        s3_ju[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) * 
                    ( fd3 (j+1, h, l) * EXPC[j-l] );
    }
    
    for (l = h+1; l < j-2; l++)
    {
        hl = index[h] + l - h;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        s3_ju[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) *
                (u1_ip_ju[lp1j] + exp_dangle3 (l,h,l+1)*EXPC[1]*(u1_ip_ju[lp2j] + u1_iu_ju[lp2j]));
    }
}


void s_partition_function::compute_base_pair_probabilities ()
// Nov 9, 2006. Computes base pair probabilities
// PRE: the partition function arrays have been filled
{
    int h, l;
    for (h=0; h < seqlen; h++)   
    {
        for (l=seqlen-1; l > h+TURN; l--)
        {
            if (can_pair (sequence[h], sequence[l]))    
            {
                compute_p (h, l);
            }
            IFD
            {
                compute_pm (h, l);
                compute_pm1 (h, l);
                compute_pm2 (h, l);
            }
            else
             {
                compute_pmnod3_needmidd3 (h, l);
                compute_pmnod3_noneedmidd3 (h, l);
                compute_pmd3_needmidd3 (h, l);
                compute_pmd3_noneedmidd3 (h, l);
                compute_pm1nod3_needendd3 (h, l);
                compute_pm1d3_needendd3 (h, l);            
                
                // the next 4 are only needed to compute the partial derivative wrt multi_free_base_penalty
                compute_pm2d5_needmidd5 (h, l);
                compute_pm2d5_noneedmidd5 (h, l);
                compute_pm2nod5_needmidd5 (h, l);
                compute_pm2nod5_noneedmidd5 (h, l);
            }
        }
    }
}


void s_partition_function::print_base_pair_probabilities (double threshold)
//prints all the base pair probabilities above the given threshold
{
    int i, j;
    printf ("Indeces start from 0.\n");
    printf ("i\tj\tprobability\n");
    for (i=0; i < seqlen-TURN-1; i++)
    {
        for (j=i+TURN+1; j < seqlen; j++)
        {
            int ij = index[i]+j-i;
            if (p[ij] > threshold)
            {
                printf ("%d\t%d\t%g\n", i, j, p[ij]);
            }
        }
    }
}



void s_partition_function::compute_p (int h, int l)
{
    int hl, zeronminus1, hm1lp1;
    int en_stack;
    int i,j;
    int en_internal;
    double term1, term2;
    
    hl = index[h]+l-h;

    //exterior base pair
    // first, 5' end
    term1 = 0.0;
    
    IFD
    {
        if (h > 0)      term1 = u[h-1];
        else            term1 = 1;
    }
    else
    {
        if (h > TURN)
        {
            term1 +=  (u_ip_jp[h-1] + u_iu_jp[h-1] + exp_dangle5(l,h,h-1)*(u_ip_ju[h-1] + u_iu_ju[h-1]));
        }
        else if (h > 0)       
        {
            term1 += exp_dangle5(l,h,h-1);  //-.[...]
        }
        else                  term1 = 1;  
    }
    
    // then, 3' end
    term2 = 0.0;
    IFD
    {
        if (l < seqlen-1)   term2 = u[index[l+1]+seqlen-1-(l+1)];
        else                term2 = 1;
    }
    else
    {
        if (l < seqlen-3)
        {
            int lp1n = index[l+1]+seqlen-1-(l+1);
            int lp2n = index[l+2]+seqlen-1-(l+2);            
            term2 += (u_ip_jp[lp1n] + u_ip_ju[lp1n] + 
                        exp_dangle3(l,h,l+1)*(u_ip_jp[lp2n] + u_ip_ju[lp2n] + u_iu_jp[lp2n] + u_iu_ju[lp2n]));
        }
        else if (l < seqlen-1) 
        {    
            term2 += exp_dangle3(l,h,l+1);
            //printf ("Add4 d3(%d,%d,%d)\n", l, h, l+1);
        }
        else                  term2 = 1;
    }
    
    zeronminus1 = index[0] + seqlen-1;
    IFD
    {
        p[hl] = term1 * up[hl] * term2 * exp_AUpenalty (h, l) / u[zeronminus1];
    }
    else
    {
        p[hl] = term1 * up[hl] * term2 * exp_AUpenalty (h, l) /
            (u_ip_jp[zeronminus1] + u_ip_ju[zeronminus1] + u_iu_jp[zeronminus1] + u_iu_ju[zeronminus1]);
    }
    
    // case stacking energies
    if (h > 0 && l < seqlen-1)
    {
        hm1lp1 = index[h-1] + l+1 - h+1;
        if (can_pair (sequence[h-1], sequence[l+1]))    
        {
            en_stack = s_stacked_pair::get_energy (h-1, l+1, sequence);    
            if (en_stack >= INF/2)
            {
                //printf ("** Infinite stack   (%d, %d) !\n", h-1, l+1);
            }
            else            
                p[hl] += p[hm1lp1] * up[hl] / up[hm1lp1] * exp (en_stack * oneoverRT);
        }
    }
    
    // case internal loop
    for (i = MAX(0, h-MAXLOOP-1); i < h; i++)
    {
        for (j = l+1; j <= MIN(MAXLOOP+l-h+i-2, seqlen-1); j++)
        {        
            if (sequence[i]+sequence[j] == 3 ||
                sequence[i]+sequence[j] == 5)        
            {
                if (i == h-1 && j == l+1) continue;    // we don't want stacked pairs here
                int ij = index[i] + j - i;
                en_internal = s_internal_loop::get_energy (i, j, h, l, sequence);
                if (en_internal >= INF/2)
                {
                    //printf ("** Infinite internal(%d, %d) !\n", i, j);
                }
                else
                    p[hl] += p[ij] * up[hl] /up[ij] * exp (en_internal * oneoverRT); 
            }
        }        
    }

    // case multi-loop
    //     u1 has multi_helix_penalty                    
        

    double pml = 0;
    for (i=0; i < h; i++)
    {
        int il = index[i] + l - i;
                
        // the case when h-l is the first branch of the multi-loop    //  (-[...] ..      
        IFD
            pml += EXPC[h-i-1] * pm[il];
        else
        {
            if (i < h-1)  // we must add d3 (which is added in pmd3); also add d5 here, if i < h-2
            {
                 pml += EXPC[h-i-1] * (i<h-2?exp_dangle5(l,h,h-1):1) *
                    (pmd3_noneedmidd3[il] + pmd3_needmidd3[il]  * exp_dangle3(l,h,l+1) * EXPC[1]);
                    // (.-[...](---)-)    i.-h...l(---)-j 
                    // (.-[...].-(---)-)    i.-h...l.-(---)-j
            }
            else    // case ((... , no dangling end
            {
                pml += pmnod3_noneedmidd3[il] + pmnod3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1];
                    // ([...](---)-)    ih...l(---)-j 
                    // ([...].-(---)-)    ih...l.-(---)-j
            }
        }
                
        if (i < h - TURN -2)    // only now we can have a branch to the left
        {      
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            IFD
            {
                pml += (pm1[il] + pm[il]) * u1[ip1hm1];
                // no branch to the right of h-l
                // branch to the left and to the right of h-l
            }
            else
            {   
                int ip2hm1 = index[i+2] + h-1 - (i+2);
                int ilp1 = index[i] + l+1 - i;
                double term1 = 0.0;
                if (up[ilp1] != 0)    // j is l+1
                {
                    term1 = p[ilp1] / up[ilp1] * exp_AUpenalty (i,l+1) *
                                (  // first, the case ((..-)-[...])
                                  ( u1_ip_jp[ip1hm1]    // ((..-)[...]) or ((..-).[...])
                                    + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1] )   // ((..-).-[...])  
                                        // no need to add EXPC[1], it's in u1_ip_ju
                                        
                                    // next, the case (.-(...)-[...])      i..-(...)h...lj
                                  + exp_dangle3(i,l+1,i+1) * EXPC[1] *
                                    ( u1_ip_jp[ip2hm1] + u1_iu_jp[ip2hm1] 
                                            + exp_dangle5(l,h,h-1)*(u1_ip_ju[ip2hm1]+u1_iu_ju[ip2hm1]) ) 
                                        // (..-(..-).-[...])     
                                );
                }
                
                pml +=
                
                // first, when h-l is the last branch to the right
                ( term1 + 
                
                // case ((..-)-[...].-)  
                pm1nod3_needendd3[il] * exp_dangle3(l,h,l+1) *    // EXPC added in pm1nod3
                    (u1_ip_jp[ip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1])
                            // c is added in pm1nod3_needendd3
                            // ((...).-[...].-)         i(...)h...l.-j
                // case (.(..-)-[...].-)
                + pm1d3_needendd3[il] * exp_dangle3(l,h,l+1) *     // don't need EXPC[1]
                    (u1_ip_jp[ip2hm1] + exp_dangle5(l,h,h-1)* u1_ip_ju[ip2hm1])
                            // (.(...).-[...].-)        i.(...)h...l.-j
                            
                // case (..-(..-)-[...].-)            
                + pm1d3_needendd3[il]  * exp_dangle3(l,h,l+1) *     // don't need EXPC[1]
                    (u1_iu_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_iu_ju[ip2hm1])
                            // (..-(...)..-[...].-)      i..-(...)h...l.-j

                            
                // we have branches to the left and to the right of h-l
                // let's do same as above, but with pm instead of pm1
                                // case ((..-)-[...].-(--)-) or  ((..-)-[...](--)-)
                + (pmnod3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1] + pmnod3_noneedmidd3[il]) *
                    (u1_ip_jp[ip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1])
                            // c is added in pm1nod3_needendd3
                            // ((...).-[...].-(--)-)         i(...)h...l.-j
                // case (.(..-)-[...].-(--)-) or (.(..-)-[...](--)-)
                + EXPC[1] * (pmd3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1] + pmd3_noneedmidd3[il]) *
                        (u1_ip_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip2hm1])
                            // (.(...).-[...].-)        i.(...)h...l.-j
                            
                // case (..-(..-)-[...].-(--)-) or (..-(..-)-[...](--)-)
                + EXPC[1] * (pmd3_needmidd3[il]  * exp_dangle3(l,h,l+1) * EXPC[1] + pmd3_noneedmidd3[il])*
                    (u1_iu_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_iu_ju[ip2hm1]) ); 
                            
            }    // end if-else ignore_dangles
        }      // end if (i < h - TURN -2)  
    }
    p[hl] +=  up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) * pml;
     
}


void s_partition_function::compute_pm (int i, int l)
// called only when no dangling ends
//  region l+1  - j-1 has at least one branch
//  i.-h    l-(---)-j
//  (.-(    )-(---)-)
{
    int j;
    int ij, lp1jm1, il;
    il = index[i]+l-i;   
    for (j = l+TURN+3; j < seqlen; j++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            lp1jm1 = index[l+1] + j-1 - l-1;
            pm[il] += exp_AUpenalty (i,j) * p[ij] / up[ij] * u1[lp1jm1];            
        }
    }    
}


void s_partition_function::compute_pmd3_needmidd3 (int i, int l)
// add the dangling end d3(i,j,i+1)
//  region l+1  - j-1 has at least one branch
//  i.-h    l-(---)-j
//  (.-(    )-(---)-)
{
    int j;
    int ij, lp2jm1, lp2jm2, lp2jm3, il;
    il = index[i]+l-i;   
    for (j = l+TURN+3; j < seqlen; j++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            lp2jm1 = index[l+2] + j-1 - l-2;
            lp2jm2 = index[l+2] + j-2 - l-2;
            lp2jm3 = index[l+2] + j-3 - l-2;            
            pmd3_needmidd3[il] += exp_AUpenalty (i,j) * exp_dangle3(i,j,i+1) * p[ij] / up[ij] * // don't need to add EXPC[1] here, it's added in the calling function
                            ( u1_ip_jp[lp2jm1] + u1_iu_jp[lp2jm1]      // ..))  or ..).)
                                + exp_dangle5(i,j,j-1) * (u1_ip_ju[lp2jm1] + u1_iu_ju[lp2jm1])  );    //)-..)            
        }
    }
}


void s_partition_function::compute_pmd3_noneedmidd3 (int i, int l)
// add the dangling end d3(i,j,i+1)
//  region l+1  - j-1 has at least one branch
//  i.-h    l-(---)-j
//  (.-(    )-(---)-)
{
    int j;
    int ij, lp1jm1, lp1jm2, lp1jm3, lp2jm1, lp2jm2, lp2jm3, il;
    il = index[i]+l-i;
    for (j = l+TURN+3; j < seqlen; j++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            lp1jm1 = index[l+1] + j-1 - l-1;
            lp1jm2 = index[l+1] + j-2 - l-1;
            lp1jm3 = index[l+1] + j-3 - l-1;            
            lp2jm1 = index[l+2] + j-1 - l-2;
            lp2jm2 = index[l+2] + j-2 - l-2;
            lp2jm3 = index[l+2] + j-3 - l-2;

            pmd3_noneedmidd3[il] += exp_AUpenalty (i,j) * exp_dangle3(i,j,i+1) * p[ij] / up[ij] * // shouldn't add EXPC[1] here because it's in the calling function
                            (u1_ip_jp[lp1jm1] + exp_dangle5(i,j,j-1) * u1_ip_ju[lp1jm1]);  
                                // ..))  or ..).),  and then )-..)            
        }
    }
}

void s_partition_function::compute_pmnod3_needmidd3 (int i, int l)
// the pm from McCaskill, which does not have d3(i,j,i+1)
//  region l+1  - j-1 has at least one branch
//  ih    l-(---)-j
//  ((    )-(---)-)
{
    int j;
    int ij, lp2jm1, lp2jm2, lp2jm3, il;
    il = index[i]+l-i;   
    for (j = l+TURN+3; j < seqlen; j++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            lp2jm1 = index[l+2] + j-1 - l-2;
            lp2jm2 = index[l+2] + j-2 - l-2;
            lp2jm3 = index[l+2] + j-3 - l-2;
            pmnod3_needmidd3[il] += exp_AUpenalty (i,j) * p[ij] / up[ij] * // no need to add EXPC[1] here, it's added in the calling function
                            ( u1_ip_jp[lp2jm1] + u1_iu_jp[lp2jm1]      // ..))  or ..).)
                                + exp_dangle5(i,j,j-1) * (u1_ip_ju[lp2jm1] + u1_iu_ju[lp2jm1]));    //)-..)
        }
    }
}


void s_partition_function::compute_pmnod3_noneedmidd3 (int i, int l)
// the pm from McCaskill, which does not have d3(i,j,i+1)
//  region l+1  - j-1 has at least one branch
//  ih    l(---)-j
//  ((    )(---)-)    // right branch follows right after, no free base in between
{
    int j;
    int ij, lp1jm1, lp1jm2, lp1jm3, lp2jm1, lp2jm2, lp2jm3, il;
    il = index[i]+l-i;   
    for (j = l+TURN+3; j < seqlen; j++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            lp1jm1 = index[l+1] + j-1 - l-1;
            lp1jm2 = index[l+1] + j-2 - l-1;
            lp1jm3 = index[l+1] + j-3 - l-1;
            lp2jm1 = index[l+2] + j-1 - l-2;
            lp2jm2 = index[l+2] + j-2 - l-2;
            lp2jm3 = index[l+2] + j-3 - l-2;            

            pmnod3_noneedmidd3[il] += p[ij] / up[ij] * exp_AUpenalty (i,j) * 
                            (u1_ip_jp[lp1jm1] + exp_dangle5(i,j,j-1) * u1_ip_ju[lp1jm1]);                               
        }
    }
}


void s_partition_function::compute_pm1 (int i, int l)
// called only when no dangling ends
// the pm1 from McCaskill, no dangling ends
// region l+1 ... j-1 is unpaired
//  i h    l.-j   
//  ( (    ).-)
{
    int j;
    int ij, il;
    il = index[i]+l-i;   
    // the case j=l+1 is dealt with separately, directly in compute_p
    for (j = l+1; j < seqlen; j++)    
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;           
            pm1[il] += p[ij] / up[ij] * exp_AUpenalty (i,j) * EXPC[j-l-1];
        }
    }
}


void s_partition_function::compute_pm1nod3_needendd3 (int i, int l)
// all free bases between l and j
// the pm1 from McCaskill, which does not have d3(i,j,i+1)
// region l+1 ... j-1 is unpaired
//  i h    l.-j    - we need to add d3(l,h,l+1)
//  ( (    ).-)
{
    int j;
    int ij, il;
    il = index[i]+l-i;   
    // the case j=l+1 is dealt with separately, directly in compute_p
    for (j = l+2; j < seqlen; j++)    
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;           
            pm1nod3_needendd3[il] += p[ij] / up[ij] * exp_AUpenalty (i,j) * 
                                     EXPC[j-l-1] * (l+2<j?exp_dangle5(i,j,j-1):1);
        }
    }
}


void s_partition_function::compute_pm1d3_needendd3 (int i, int l)
// all free bases between l and j
// add the dangling end d3(i,j,i+1)
// region l+1 ... j-1 is unpaired
//  i h    l........j
//  ( (    )........)
{
    int j;
    int ij, il;
    il = index[i]+l-i;   
    // the case j=l+1 is dealt with separately, directly in compute_p
    for (j = l+2; j < seqlen; j++)    
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;           
            pm1d3_needendd3[il] += exp_dangle3(i,j,i+1) * EXPC[1] * p[ij] / up[ij] * exp_AUpenalty (i,j) *
                                   EXPC[j-l-1] * (l+2<j?exp_dangle5(i,j,j-1):1);
        }
    }
}


void s_partition_function::compute_pm2 (int h, int j)
// called only when no dangling ends
//  region i+1  - h-3 has at least one branch, and there must at least 2 free bases before h
// i---------h...l.-.j
// (-(---)-..[...].-.)
{
    int i, ij;
    int hj, ip1;
    hj = index[h]+j-h;   
    
    for (i=0; i < h-TURN-2; i++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            pm2[hj] += p[ij] / up[ij] * exp_AUpenalty (i,j) * u1[ip1hm1];
        }
    }
}


void s_partition_function::compute_pm2d5_needmidd5 (int h, int j)
//  region i+1  - h-3 has at least one branch, and there must at least 2 free bases before h
// also add the dangling end d5(i,j,j-1)
// i---------h...l.-.j
// (-(---)-..[...].-.)
{
    int i, ij;
    int hj, ip1;
    hj = index[h]+j-h;   
    
    for (i=0; i < h-TURN-3; i++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            int ip2hm1 = index[i+2] + h-1 - (i+2);
            pm2d5_needmidd5[hj] += p[ij] / up[ij] * exp_AUpenalty (i,j) * exp_dangle5(i,j,j-1) *
                // the following is the same in compute_pm2nod5_needmidd5
                // we have to make sure the middle d3 is included
                // first, ((---)-..[...].-.)
                ( u1_ip_ju[ip1hm1] + exp_dangle3(i,j,i+1) * EXPC[1] * (u1_ip_ju[ip2hm1] + u1_iu_ju[ip2hm1]));
        }
    }
}


void s_partition_function::compute_pm2d5_noneedmidd5 (int h, int j)
//  region i+1  - h-1 has at least one branch, and there must be 0 or 1 free bases left of h
// also add the dangling end d5(i,j,j-1)
// i------.h...l.-.j
// (-(---).[...].-.) or (-(---)[...].-.)
{
    int i, ij;
    int hj, ip1;
    hj = index[h]+j-h;   
    
    for (i=0; i < h-TURN-2; i++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            int ip2hm1 = index[i+2] + h-1 - (i+2);
            pm2d5_noneedmidd5[hj] += p[ij] / up[ij] * exp_AUpenalty (i,j) * exp_dangle5(i,j,j-1) * 
                // first, ((---)[...].-.) or ((---).[...].-.)
                ( u1_ip_jp[ip1hm1] + exp_dangle3(i,j,i+1) * EXPC[1] * (u1_ip_jp[ip2hm1] + u1_iu_jp[ip2hm1]));
        }
    }
}


void s_partition_function::compute_pm2nod5_needmidd5 (int h, int j)
//  region i+1  - h-3 has at least one branch, and there must at least 2 free bases before h
// DO NOT add the dangling end d5(i,j,j-1)
// i---------h...l.j
// (-(---)-..[...]) or (-(---)-..[...].)
// exactly the same as compute_pm2d5_needmidd5, but withour exp_dangle5
{
    int i, ij;
    int hj, ip1;
    hj = index[h]+j-h;   
    
    for (i=0; i < h-TURN-2; i++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            int ip2hm1 = index[i+2] + h-1 - (i+2);
            pm2nod5_needmidd5[hj] += p[ij] / up[ij] * exp_AUpenalty (i,j) * 
                // we have to make sure the middle d3 is included
                // first, ((---)-..[...]) or ((---)-..[...].)
                ( u1_ip_ju[ip1hm1] + exp_dangle3(i,j,i+1) * EXPC[1] * (u1_ip_ju[ip2hm1] + u1_iu_ju[ip2hm1]));
        }
    }
}


void s_partition_function::compute_pm2nod5_noneedmidd5 (int h, int j)
//  region i+1  - h-1 has at least one branch, and there must be 0 or 1 free bases left of h
// also add the dangling end d5(i,j,j-1)
// i------.h...l.-.j
// (-(---).[...]) or (-(---)[...])         or   (-(---).[...].) or (-(---)[...].)
// same as compute_pm2d5_noneedmidd5, just remove the exp_dangle5
{
    int i;
    int hj, ip1, ij;
    hj = index[h]+j-h;   
    
    for (i=0; i < h-TURN-2; i++)
    {
    
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            int ip2hm1 = index[i+2] + h-1 - (i+2);
            pm2nod5_noneedmidd5[hj] += p[ij] / up[ij] * exp_AUpenalty (i,j) * 
                // first, ((---)[...].-.) or ((---).[...].-.)
                ( u1_ip_jp[ip1hm1] + exp_dangle3(i,j,i+1) * EXPC[1] * (u1_ip_jp[ip2hm1] + u1_iu_jp[ip2hm1]));
        }
    }
}




/*---------------------------------------------------------------------------*/
#define PMIN 0.00001
void s_partition_function::PS_dot_plot(char *wastlfile)
// taken from Vienna package 1.4 
// plots a matrix with the square root of probabilities (for better visibility)
{
  /* produce PostScript dot plot from probabilities in p[] array */
   
  FILE *wastl;
  char name[31], *c;
  int i, j, length;
  double tmp;
   
  length= strlen(csequence);
  wastl = fopen(wastlfile,"w");
  if (wastl==NULL) {
    fprintf(stderr, "can't open %s for dot plot\n", wastlfile);
    return; /* return 0 for failure */
  }
  strncpy(name, wastlfile, 30);
  if ((c=strrchr(name, '_'))!=0) *c='\0';
  fprintf(wastl,"%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(wastl,"%%%%Title: RNA DotPlot\n");
  fprintf(wastl,"%%%%Creator: PS_dot.c, ViennaRNA Package\n");
  //fprintf(wastl,"%%%%CreationDate: %s", time_stamp());
  fprintf(wastl,"%%%%BoundingBox: 66 211 518 662\n");
  fprintf(wastl,"%%%%DocumentFonts: Helvetica\n");
  fprintf(wastl,"%%%%Pages: 1\n");
  fprintf(wastl,"%%%%EndComments\n\n");
//  fprintf(wastl,"%%Options: %s\n", option_string());

  fprintf(wastl,"%%This file contains the square roots "
      "of the base pair probabilities in the form\n");
  fprintf(wastl,"%% i  j  sqrt(p(i,j)) ubox\n");
   
  fprintf(wastl,"100 dict begin\n");  /* DSC says EPS should create a dict */
  fprintf(wastl,"\n/logscale false def\n\n");
  fprintf(wastl,"%%delete next line to get rid of title\n"
      "270 665 moveto /Helvetica findfont 14 scalefont setfont "
      "(%s) show\n\n", name);
  fprintf(wastl,"/lpmin {\n"
      "   %g log  %% log(pmin) only probs>pmin will be shown\n"
      "} bind def\n\n",PMIN);

  /* EPS should not contain lines >255 characters */
  fprintf(wastl,"/sequence { (\\\n");
  i=0;
  while (i<length) {
    fprintf(wastl, "%.255s\\\n", csequence+i);
    i+=255;
  }
  fprintf(wastl,") } def\n");
  fprintf(wastl,"/len { sequence length } def\n\n");
   
  fprintf(wastl,"/ubox {\n"     /* upper triangle matrix */
      "   logscale {\n"
      "      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n"
      "   } if\n"
      "   3 1 roll\n"
      "   exch len exch sub 1 add box\n"
      "} bind def\n\n");

  fprintf(wastl,"/lbox {\n"     /* lower triangle matrix */
      "   3 1 roll\n"
      "   len exch sub 1 add box\n"
      "} bind def\n\n");

  fprintf(wastl,"/box { %%size x y box - draws box centered on x,y\n"
      "   2 index 0.5 mul add            %% x += 0.5\n"
      "   exch 2 index 0.5 mul add exch  %% x += 0.5\n"
      "   newpath\n"
      "   moveto\n"
      "   dup neg   0 rlineto\n"
      "   dup neg   0 exch rlineto\n"
      "             0 rlineto\n"
      "   closepath\n"
      "   fill\n"
      "} def\n\n");

  fprintf(wastl,"72 216 translate\n");
  fprintf(wastl,"72 6 mul len 1 add div dup scale\n");
  fprintf(wastl,"/Helvetica findfont 0.95 scalefont setfont\n\n");

   /* print sequence along all 4 sides */
  fprintf(wastl,"%% print sequence along all 4 sides\n");
  fprintf(wastl,"0 1 len 1 sub {\n");
  fprintf(wastl,"    dup\n");
  fprintf(wastl,"    0.7 add -0.3 moveto\n");
  fprintf(wastl,"    sequence exch 1 getinterval\n");
  fprintf(wastl,"    show\n");
  fprintf(wastl,"} for\n");
  fprintf(wastl,"\n");
  fprintf(wastl,"0 1 len 1 sub {\n");
  fprintf(wastl,"    dup\n");
  fprintf(wastl,"    0.7 add 0.7 len add moveto\n");
  fprintf(wastl,"    sequence exch 1 getinterval\n");
  fprintf(wastl,"    show\n");
  fprintf(wastl,"} for\n\n");

  fprintf(wastl,"90  rotate\n");
  fprintf(wastl,"0 1 len 1 sub {\n");
  fprintf(wastl,"    dup\n");
  fprintf(wastl,"    0.7 add -0.2 moveto\n");
  fprintf(wastl,"    len 1 sub exch sub\n");
  fprintf(wastl,"    sequence exch 1 getinterval\n");
  fprintf(wastl,"    show\n");
  fprintf(wastl,"} for\n");
  fprintf(wastl,"270 rotate\n\n");

  fprintf(wastl,"270 rotate\n");
  fprintf(wastl,"0 1 len 1 sub {\n");
  fprintf(wastl,"    dup\n");
  fprintf(wastl,"    -0.3 add len sub  0.7 len add  moveto\n");
  fprintf(wastl,"    sequence exch 1 getinterval\n");
  fprintf(wastl,"    show\n");
  fprintf(wastl,"} for\n");
  fprintf(wastl,"90 rotate\n\n");

  /* do grid */
  fprintf(wastl,"0.5 dup translate\n"
      "%% draw diagonal\n"
      "0.04 setlinewidth\n"
      "0 len moveto len 0 lineto stroke \n\n");
  fprintf(wastl,"%%draw grid\n"
      "0.01 setlinewidth\n"
      "len log 0.9 sub cvi 10 exch exp  %% grid spacing\n"
      "dup 1 gt {\n"
      "   dup dup 20 div dup 2 array astore exch 40 div setdash\n"
      "} { [0.3 0.7] 0.1 setdash } ifelse\n"
      "0 exch len {\n"      /* for (i=0; i<=len; i++) */
      "   dup dup\n"        
      "   0 moveto\n"                     /* i 0 moveto   */
      "   len lineto \n"                  /* i len lineto */
      "   dup\n"
      "   len exch sub 0 exch moveto\n"   /* 0 i moveto   */
      "   len exch len exch sub lineto\n" /* len i lineto */
      "   stroke\n"
      "} for\n"
      "0.5 neg dup translate\n\n");

  /* print boxes */
  for (i=0; i<length; i++)
    for (j=i+1; j<=length; j++) {
      if (p[index[i]+j-i]<PMIN) continue;
      tmp = sqrt(p[index[i]+j-i]);
      fprintf(wastl,"%d %d %1.5f ubox\n", i+1, j+1, tmp);
    }
    // come back to this later    
  /* do mfe */
  /*
  if (base_pair)
    for(i=1; i<=base_pair[0].i; i++) 
      fprintf(wastl,"%d %d 0.95 lbox\n",
          base_pair[i].i, base_pair[i].j); 
  */  
  fprintf(wastl,"showpage\n");
  fprintf(wastl,"end\n");
  fprintf(wastl,"%%%%EOF\n");
  fclose(wastl);
  //return 1; /* success */
}


int s_partition_function::has_base_pair (int i, int j, char *structure)
// returns 1 if structure contains base pair (i,j)
{
    int ptable[MAXSLEN];
    detect_original_pairs (structure, ptable); 
    if (ptable[i] == j)
    {
        //printf ("Structure %s has base pair %d-%d\n", structure, i, j);    
        return 1;
    }
    return 0;
}


int s_partition_function::identical_structure (int i, int j, char *structure1, char *structure2)
// return 1 if structure 1 and structure2 are identical between i and j
{
    int k;
    for (k=i; k <= j; k++)
    {
        if (structure1[k] != structure2[k])
            return 0;
    }
    return 1;
}

double s_partition_function::exp_free_energy_partial (int i, int j, char *seq, char *structure, int removeAU)
// PRE: structure is valid from i to j
// return exp of free energy of sequence and structure, for region i-j
{
    char partial_sequence[MAXSLEN];
    char partial_structure[MAXSLEN];
    int k;
    double energy, beta;
    beta = 1000.0/(1.98717*310.15);
    for (k=i; k <= j; k++)
    {
        partial_sequence[k-i] = seq[k];
        partial_structure[k-i] = structure[k];
    }
    partial_sequence[k-i] = '\0';
    partial_structure[k-i] = '\0';
    //printf ("Partial: %s %s\n", partial_sequence, partial_structure);
    energy = free_energy_simfold (partial_sequence, partial_structure);
    if (removeAU)
    {
        if (has_AU_penalty(nuc_to_int(seq[i]), nuc_to_int(seq[j])))    // remove the AU penalty, because we don't know if this is going to be an exterior loop or not
        { 
            energy -= AU_penalty (nuc_to_int(seq[i]), nuc_to_int(seq[j]))/100.0;
        }
    }
    return exp(-1.0 * energy * beta);
}



double s_partition_function::compute_partition_function_exhaustively ()
{
    char structure[MAXSLEN];
    double enthalpy, energy;
    char tmp_structures[MAXSUBSTR][MAXSLEN];
    double tmp_energies[MAXSUBSTR];
    double min_energy, max_energy;
    int actual_num_str;
    int i,ii,j, k;
    double R, temp, beta;
    //R = 0.00198717;
    //temp = 310.15;
    beta = (double)1000.0/(1.98717*310.15);
    //beta = 1; 


    int ptable[MAXSLEN];
    double strprob;
        
    s_min_folding *min_fold = new s_min_folding (csequence);
    min_energy = min_fold->s_simfold ();
    min_fold->return_structure (structure);
    delete min_fold;      
    
    //s_sub_folding* sub_fold = new s_sub_folding(sequence, -(int)(min_energy*100.0));
    s_sub_folding* sub_fold = new s_sub_folding(csequence, 10000);
    sub_fold->set_limit (MAXSUBSTR);
    sub_fold->s_simfold (enthalpy);
    actual_num_str = sub_fold->return_structures(tmp_structures, tmp_energies);
    delete sub_fold;
    
    max_energy = tmp_energies[actual_num_str-1];
    Zexhaustive = 0.0;
    for (i=0; i < actual_num_str; i++)
    {
        //if (!(i==0 || i==2 || i==4 || i==7)) continue;
        //if (!(i==54)) continue;
        //if (!(i==49 || i==50 || i==51 || i==53 || i==54 || i==55 || i==57)) continue;
        //if (i==57 || i==52 || i==55 || i==54 || i==53 || i==51 || i==50) continue;
        //if (i==56 || i==40 || i==39 || i==38 || i==31 || i==25 || i==15) continue;
        // recompute the free energy, i.e. with the correct dangling ends
        energy = free_energy_simfold (csequence, tmp_structures[i]);
        strprob = (double)(exp(-1.0 * energy * beta)/Z);
        detect_original_pairs (tmp_structures[i], ptable); 
        for (ii=0; ii < seqlen; ii++)
        {
            if (ptable[ii] > ii)
            {
                //printf ("ptable[%d]=%d\n", ii, ptable[ii]);
                int ind = index[ii] + ptable[ii] - ii;
                pexhaustive[ind] += strprob;
            }
        }   
        
        printf ("Substr %d: %s\ten=%.2lf \tprob=%e\n", i, tmp_structures[i], energy, strprob);
        Zexhaustive += exp (-1.0 * energy * beta);
    }
    
    printf ("\n*** Number of substructures: %d\n\n", actual_num_str);
    // try to get the exhaustive u[i,j] and exhaustive up[i,j]
    int l, identical;
    for (i=0; i < seqlen; i++)
    {
        for (j=i+TURN+1; j < seqlen; j++)
        {    
            int ij = index[i] + j - i;
            for (k=0; k < actual_num_str; k++)                
            {
                // compute exhaustive u[i,j]
                // ignore "...." because every u[i,j] 
                if (valid_structure (i, j, tmp_structures[k]) && is_structured (i, j, tmp_structures[k]))
                {
                    identical = 0;
                    for (l=0; l < k; l++)
                    {
                        if (identical_structure (i, j, tmp_structures[l], tmp_structures[k]))
                            identical = 1;
                    }
                    if (!identical)
                    {
                        //printf ("Adding: %g\n", exp_free_energy_partial (i, j, csequence, tmp_structures[k], 0));
                        uexhaustive[ij] += exp_free_energy_partial (i, j, csequence, tmp_structures[k], 0);
                    }
                }
                //else
                //    printf ("%s not valid from %d to %d\n", tmp_structures[k], i, j);
                                
                // compute exhaustive up[i,j]
                if (has_base_pair(i,j, tmp_structures[k]))
                {
                    identical = 0;
                    for (l=0; l < k; l++)
                    {
                        if (identical_structure (i, j, tmp_structures[l], tmp_structures[k]))
                            identical = 1;
                    }
                    if (!identical)
                    {
                        upexhaustive[ij] += exp_free_energy_partial (i, j, csequence, tmp_structures[k], 1);
                    }
                    //else
                    //    printf ("Str %s identical with other for i=%d, j=%d\n", tmp_structures[k], i, j);
                }                        
            }
        }
    }
    return Zexhaustive;        
}


void s_partition_function::verify_partition_function()
{
    //#define PREC 0.0000000001    
    int i,j,ij;    
    printf ("Checking partition function... ");
    if (Z - Zexhaustive < PREC && Zexhaustive - Z < PREC)
        printf ("\t++++++++++++++++++++++\n\tZ = %g, Ze = %g, diff = %g\n", Z, Zexhaustive, Z-Zexhaustive);
    else
    {
        printf ("\t----------------------\n\tZ = %g, Ze = %g, diff = %g\n", Z, Zexhaustive, Z-Zexhaustive);
        exit(1);
    }
        
        
//     printf ("Checking u  .................. ");
//     int ufine = 1;
//     for (i=0; i < seqlen; i++)
//         for (j=i+1; j < seqlen; j++)
//         {
//             ij = index[i]+j-i;
//             if (u[ij] - uexhaustive[ij] > PREC || uexhaustive[ij] - u[ij] > PREC)
//             {
//                 if (ufine)
//                     printf ("\t----------------------\n");
//                 printf ("\tu[%d,%d] = %g, uex = %g, diff = %g\n", i,j,u[ij],uexhaustive[ij], u[ij]-uexhaustive[ij]);
//                 ufine = 0;
//             }
//             //printf ("\tu[%d,%d] = %g, uex = %g, diff = %g\n", i,j,u[ij],uexhaustive[ij],u[ij]-uexhaustive[ij]);
//         }
//     if (ufine) 
//         printf ("\t++++++++++++++++++++++\n");
        
    printf ("Checking up .................. ");
    int upfine = 1;
    for (i=0; i < seqlen; i++)
        for (j=i+1; j < seqlen; j++)
        {
            ij = index[i]+j-i;
            if (up[ij] - upexhaustive[ij] > PREC || upexhaustive[ij] - up[ij] > PREC)
            {
                if (upfine)
                    printf ("\t----------------------\n");
                printf ("\tup[%d,%d] = %g, upex = %g, diff = %g\n", i,j,up[ij],upexhaustive[ij], up[ij]-upexhaustive[ij]);
                upfine = 0;                
            }
            //printf ("\tup[%d,%d] = %g, upex = %g, diff = %g\n", i,j,up[ij],upexhaustive[ij],up[ij]-upexhaustive[ij]);
        }
    if (upfine) 
        printf ("\t++++++++++++++++++++++\n");
    // NOTE: this will be different when we have the special GGG-U hairpin loop. It's okay to be different, it doesn't affect the partition function.

                
    printf ("Checking p ................... ");
    int pfine = 1;
    for (i=0; i < seqlen; i++)
        for (j=i+1; j < seqlen; j++)
        {
            ij = index[i]+j-i;
            if (p[ij] - pexhaustive[ij] > PREC || pexhaustive[ij] - p[ij] > PREC)
            {
                if (pfine)
                    printf ("\t----------------------\n");
                printf ("\tp[%d,%d] = %g, pex = %g, diff = %g\n", i,j,p[ij],pexhaustive[ij], p[ij]-pexhaustive[ij]);
                pfine = 0;
            }
            //printf ("\tp[%d,%d] = %g, pex = %g, diff = %g\n", i,j,p[ij],pexhaustive[ij], p[ij]-pexhaustive[ij]);
        }
    if (pfine) 
        printf ("\t++++++++++++++++++++++\n");
    else
        exit(1);


    printf ("Checking u ................... ");
    int ufine = 1;
    for (i=0; i < seqlen; i++)
        for (j=i+1; j < seqlen; j++)
        {
            ij = index[i]+j-i;
            //u1 should be the same as u1_ip_jp + u1_ip_ju + u1_iu_jp + u1_iu_ju
            double unew;
            IFD
                unew = u[ij];
            else
                unew = u_ip_jp[ij] + u_ip_ju[ij] + u_iu_jp[ij] + u_iu_ju[ij];
            if (uexhaustive[ij] - unew > PREC || unew - uexhaustive[ij] > PREC)
            {
                if (ufine)
                    printf ("\t----------------------\n");
                printf ("\tunew[%d,%d] = %g, u = %g, diff = %g\n", i,j,unew,uexhaustive[ij], unew-uexhaustive[ij]);
                //printf ("\t  u_ip_jp=%g, u_ip_ju=%g, u_iu_jp=%g, u_iu_ju=%g\n", u_ip_jp[ij], u_ip_ju[ij], u_iu_jp[ij], u_iu_ju[ij]);
                ufine = 0;
            }
        }
    if (ufine) 
        printf ("\t++++++++++++++++++++++\n");
    else
        exit(1);
        
    /*
    printf ("Checking s1 ................... ");
    int s1fine = 1;
    for (i=0; i < seqlen; i++)
        for (j=i+1; j < seqlen; j++)
        {
            ij = index[i]+j-i;
            //s3 should be the same as s3_jp + s3_ju
            double s1new = s1_jp[ij] + s1_ju[ij];
            if (s1[ij] - s1new > PREC || s1new - s1[ij] > PREC)
            {
                if (s1fine)
                    printf ("\t----------------------\n");
                printf ("\ts1new[%d,%d] = %g, s1 = %g, diff = %g\n", i,j,s1new,s1[ij], s1new-s1[ij]);
                s1fine = 0;
            }
        }
    if (s1fine) 
        printf ("\t++++++++++++++++++++++\n");
    */
        
//     printf ("Checking u1 ................... ");
//     int u1fine = 1;
//     for (i=0; i < seqlen; i++)
//         for (j=i+1; j < seqlen; j++)
//         {
//             ij = index[i]+j-i;
//             //u1 should be the same as u1_ip_jp + u1_ip_ju + u1_iu_jp + u1_iu_ju
//             double u1new = u1_ip_jp[ij] + u1_ip_ju[ij] + u1_iu_jp[ij] + u1_iu_ju[ij];
//             if (u1[ij] - u1new > PREC || u1new - u1[ij] > PREC)
//             {
//                 if (u1fine)
//                     printf ("\t----------------------\n");
//                 printf ("\tu1new[%d,%d] = %g, u1 = %g, diff = %g\n", i,j,u1new,u1[ij], u1new-u1[ij]);
//                 printf ("\t  u1_ip_jp=%g, u1_ip_ju=%g, u1_iu_jp=%g, u1_iu_ju=%g\n", u1_ip_jp[ij], u1_ip_ju[ij], u1_iu_jp[ij], u1_iu_ju[ij]);
//                 u1fine = 0;
//             }
//         }
//     if (u1fine) 
//         printf ("\t++++++++++++++++++++++\n");
        
//     printf ("Checking s3 ................... ");
//     int s3fine = 1;
//     for (i=0; i < seqlen; i++)
//         for (j=i+1; j < seqlen; j++)
//         {
//             ij = index[i]+j-i;
//             //s3 should be the same as s3_jp + s3_ju
//             double s3new = s3_jp[ij] + s3_ju[ij];
//             if (s3[ij] - s3new > PREC || s3new - s3[ij] > PREC)
//             {
//                 if (s3fine)
//                     printf ("\t----------------------\n");
//                 printf ("\ts3new[%d,%d] = %g, s3 = %g, diff = %g\n", i,j,s3new,s3[ij], s3new-s3[ij]);
//                 s3fine = 0;
//             }
//         }
//     if (s3fine) 
//         printf ("\t++++++++++++++++++++++\n");


}



void s_partition_function::compute_logZ_gradient_exhaustively ()
{        
    char structure[MAXSLEN];
    double enthalpy, energy;
    char tmp_structures[MAXSUBSTR][MAXSLEN];
    double tmp_energies[MAXSUBSTR];
    double min_energy, max_energy;
    int i, k;
    double R, temp, beta;
    int real_str_found, actual_num_str;
    R = 0.00198717;
    temp = 310.15;
    beta = 1.0/(R*temp);
    real_str_found = 0;
    double numerator [MAXNUMPARAMS];
    double counter [MAXNUMPARAMS];
    double denominator;      
    
    s_sub_folding* sub_fold = new s_sub_folding(csequence, 10000);
    sub_fold->set_limit (MAXSUBSTR);
    sub_fold->s_simfold (enthalpy);
    actual_num_str = sub_fold->return_structures(tmp_structures, tmp_energies);
    delete sub_fold;
    
    for (i=0; i < num_params; i++)
    {
        numerator[i] = 0;
    }    
        
    // first compute the denominator, which is the same for all parameters
    denominator = 0;
    for (k=0; k < actual_num_str; k++)
    {        
        // recompute the free energy, i.e. with the correct dangling ends
        energy = free_energy_simfold (csequence, tmp_structures[k]);        
        denominator += exp ((-1) * energy * beta);
        count_each_structure_type (csequence, tmp_structures[k], counter, 1);        
        //printf ("%s, en=%g\n", tmp_structures[k], energy);
        for (i=0; i < num_params; i++)
        {
            numerator[i] += counter[i] * exp ((-1) * energy * beta);
            //if (counter[i] > 0)
            //    printf ("\t%s\t%d times, value=%g\n", string_params[i], counter[i], numerator[i]);
        }
        //printf ("Sub str %d, counter[317] = %d, numerator[317] = %e\n", k, counter[317], numerator[317]);
    }
    // now the denominator and nominator are computed
           
    for (i=0; i < num_params; i++)
    {
        GlogZexhaustive[i] = numerator[i] / denominator;
    }
    //printf ("Denominator = Z = %g\n", denominator);
    //printf ("logZ_gradient[317] = %e\n", GlogZexhaustive[317]);    
}


void s_partition_function::compute_logZ_gradient ()
// u, up and p arrays are filled
{ 

    int index_param;
    int i, j, k, l, m, n, o, r, h;
    int ii, jj, iip, jjp, iijj, iipjjp, iip1jjm1;
    int en_hairpin, en_stack, en_internal;
    int AUpen, dangle, dangle5, dangle3;
    char type[100];
    int tindex;    // type index
    int il_i, il_j, il_ip1, il_jm1;
    int index_should_be;
    
    num_params = create_string_params();

    int AUpen_index = structure_type_index ("misc.terminal_AU_penalty");
   
    for (i=0; i < num_params; i++)    GlogZ[i] = 0;
    
    // stacking energies
    index_param = 0;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (stack[i][j][k][l] < INF)
                    {
                        // exclude duplicates
                        // stack[i][j][k][l] is the same as stack[l][k][j][i]
                        if (i*1000 + j*100 + k*10 + l <= l*1000 + k*100 + j*10 + i)
                        {
                            // traverse the whole sequence and look for this building block
                            // first look from the left
                            for (ii = 0; ii < seqlen; ii++)
                            {
                                for (jj = ii+TURN+3; jj < seqlen; jj++)
                                {
                                    if (sequence[ii] == i && sequence[jj] == j)
                                    {
                                        iijj = index[ii] + jj -ii;
                                        // stacking energies are also involved in bulges of size 1
                                        for (iip = ii + 1; iip <= ii+2; iip++)
                                        {
                                            for (jjp = jj-1; jjp >= jj-2; jjp--)
                                            {
                                                if (iip == ii+2 && jjp == jj-2)     // this is not stacking energy, is internal loop 1x1
                                                    continue;                                                
                                                if (sequence[iip] == k && sequence[jjp] == l) 
                                                {                                                        
                                                    iipjjp = index[iip]+jjp-iip;
                                                    if (iip == ii+1 && jjp == jj-1)    // stack pair
                                                        en_stack = s_stacked_pair::get_energy (ii, jj, sequence);    
                                                    else    // bulge
                                                        en_stack = s_internal_loop::get_energy (ii, jj, iip, jjp,     sequence); 
                                                    GlogZ[index_param] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_stack * oneoverRT);     
                                                }
                                            }                                    
                                        }
                                    }
                                }
                            }
                            // check if it's symmetric. If so, it was already considered when we looked from the left
                            if (i!=l || j!=k)    // then it's not symmetric
                            {
                                // second, look from the right
                                for (ii = seqlen-1; ii > TURN; ii--)
                                {
                                    for (jj = ii-TURN-1; jj > 0; jj--)
                                    {
                                        if (sequence[ii] == i && sequence[jj] == j)
                                        {
                                            iijj = index[jj] + ii -jj;
                                            
                                            // stacking energies are also involved in bulges of size 1
                                            for (iip = ii + 1; iip <= MIN(ii+2,seqlen-1); iip++)
                                            {
                                                for (jjp = jj-1; jjp >= MAX(jj-2,0); jjp--)
                                                {
                                                    if (iip == ii+2 && jjp == jj-2)     // this is not stacking energy, is internal loop 1x1
                                                        continue;                                            
                                                     
                                                    if(sequence[iip] == k && sequence[jjp] == l)       
                                                    {                                                        
                                                        iipjjp = index[jjp]+iip-jjp;
                                                        if (iip == ii+1 && jjp == jj-1)    // stack pair
                                                            en_stack = s_stacked_pair::get_energy (jjp, iip, sequence);    
                                                        else    // bulge
                                                            en_stack = s_internal_loop::get_energy (jjp, iip, jj, ii,     sequence); 
                                                        GlogZ[index_param] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_stack * oneoverRT);     
                                                    }
                                                }
                                            
                                            }       
                                        }
                                    }
                                }                                
                            }    
                            index_param++;                                                        
                        }                        
                    }                    
                }
                 
    index_should_be = structure_type_index("tstackh[0][3][0][0]");                
    if (index_param != index_should_be)
    {
        printf ("Index param after stack = %d, should be %d\n", index_param, index_should_be);
        exit(1);
    }                
    
    // tstackh and hairpin energies
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (tstackh[i][j][k][l] < INF)
                    {
                        // no duplicates here
                        for (ii = 0; ii < seqlen; ii++)
                        {                            
                            for (jj = ii+TURN+1; jj < seqlen; jj++)
                            {
                                if (sequence[ii] == i && sequence[jj] == j && 
                                    sequence[ii+1] == k && sequence[jj-1] == l)
                                {
                                
                                    // compute the derivative for hairpin penalty by size
                                    char s[100];
                                    sprintf (s, "hairpin_penalty_by_size[%d]", jj-ii-1);
                                    // this is very slow, it's linear in the number of parameters
                                    int sizeindex = structure_type_index (s);
                                    iijj = index[ii] + jj -ii;
                                    en_hairpin = s_hairpin_loop::get_energy (ii, jj, sequence, csequence);
                                    GlogZ[sizeindex] += p[iijj] / up[iijj] * exp (en_hairpin * oneoverRT);  
                                    
                                    // add the AU_penalty
                                    
                                    if (jj-ii-1 == 3)
                                    {
                                        if (has_AU_penalty (sequence[ii], sequence[jj]))
                                        {
                                            GlogZ[AUpen_index] += p[iijj] / up[iijj] * exp (en_hairpin * oneoverRT);
                                        }                                        
                                    }
                                    
                                    // check to see if it's a special triloop or tetraloop
                                        // check if it is a triloop
                                    if (jj-ii-1 == 3)
                                    {
                                        char seq[10];
                                        substr (csequence, ii, jj, seq);
                                        for (int kk=0; kk < nb_triloops; kk++)
                                        {
                                            if (strcmp (seq, triloop[kk].seq) == 0)
                                            {
                                                sprintf (s, "triloop[%d].energy", kk);
                                                sizeindex = structure_type_index (s);
                                                GlogZ[sizeindex] += p[iijj] / up[iijj] * exp (en_hairpin * oneoverRT);  
                                            }
                                        }
                                    }
                                    
                                    // check to see it is a tetraloop in tloop
                                    else if (jj-ii-1 == 4)
                                    {
                                        char seq[10];
                                        substr (csequence, ii, jj, seq);                                        
                                        for (int kk=0; kk < nb_tloops; kk++)
                                        {
                                            if (strcmp (seq, tloop[kk].seq) == 0)
                                            {
                                                sprintf (s, "tloop[%d].energy", kk);
                                                sizeindex = structure_type_index (s);
                                                GlogZ[sizeindex] += p[iijj] / up[iijj] * exp (en_hairpin * oneoverRT);  
                                            }
                                        }
                                    }

                                    // compute the derivatives for misc.hairpin_GGG
                                    if (ii > 1)
                                    {
                                        if (sequence[ii]==G && sequence[ii-1]==G && sequence[ii-2]==G and sequence[jj]==U)
                                        {
                                            sizeindex = structure_type_index ("misc.hairpin_GGG");
                                            GlogZ[sizeindex] += p[iijj] / up[iijj] * exp (en_hairpin * oneoverRT);  
                                        }
                                    }
                                    
                                    // check for the special case of "poly-C" hairpin loop
                                    int is_poly_C = 1;
                                    for (int kk=ii+1; kk<jj; kk++)
                                    {
                                        if (sequence[kk] != C)
                                        {
                                            is_poly_C = 0;
                                            break;
                                        }
                                    }
                                    if (is_poly_C)
                                    {
                                        if (jj-ii-1 == 3)
                                        {
                                            sizeindex = structure_type_index ("misc.hairpin_c3");
                                            GlogZ[sizeindex] += p[iijj] / up[iijj] * exp (en_hairpin * oneoverRT);  
                                        }
                                        else
                                        {
                                            sizeindex = structure_type_index ("misc.hairpin_c1");
                                            // we have to multiply by size because this param is multiplied by size in the hairpin loop calculation
                                            GlogZ[sizeindex] += (jj-ii-1) * p[iijj] / up[iijj] * exp (en_hairpin * oneoverRT);
                                            sizeindex = structure_type_index ("misc.hairpin_c2");
                                            GlogZ[sizeindex] += p[iijj] / up[iijj] * exp (en_hairpin * oneoverRT);
                                        }
                                    }
                                    
                                    // compute the derivatives of tstackh
                                    // tstackh is not added to hairpin loops of size 3
                                    if (jj-ii-1 > 3)
                                        GlogZ[index_param] += p[iijj] / up[iijj] * exp (en_hairpin * oneoverRT);
                                }
                            }
                        }
                        index_param++;
                    }
                }
    index_should_be = structure_type_index("misc.internal_AU_closure");                
    if (index_param != index_should_be)
    {
        printf ("Index param after tstackh = %d, should be %d\n", index_param, index_should_be);
        exit(1);
    }                
                                
    // tstacki energies                
    // in fact, we only have 3 parameters                            
    
    char s[100];
    //strcpy (s, "misc.terminal_AU_penalty");
    // this is very slow, it's linear in the number of parameters
    //int index_terminal_AU_penalty = structure_type_index (s);                        
    strcpy (s, "misc.internal_AU_closure");
    int index_internal_AU_closure = structure_type_index (s);                        
    strcpy (s, "misc.internal_AG_mismatch");
    int index_internal_AG_mismatch = structure_type_index (s);
    strcpy (s, "misc.internal_UU_mismatch");
    int index_internal_UU_mismatch = structure_type_index (s);    
                                
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (tstacki[i][j][k][l] < INF)
                    {
                        // no duplicates here
                        for (ii = 0; ii < seqlen; ii++)
                        {
                            for (jj = ii+TURN+4; jj < seqlen; jj++)
                            {
                                // first check the weird case of k==0 and l==0
                                if (sequence[ii] == i && sequence[jj] == j && k==0 && l==0)
                                {
                                    iijj = index[ii] + jj -ii;
                                    
                                    for (iip = ii+1; iip <= MIN(jj-2-TURN,ii+MAXLOOP+1) ; iip++)  // j-2-TURN
                                    {
                                        int minq = MAX (jj-ii+iip-MAXLOOP-2, iip+1+TURN);    // ip+1+TURN);
                                        for (jjp = minq; jjp < jj; jjp++)
                                        {        
                                            if (sequence[iip]+sequence[jjp] == 3 ||
                                                sequence[iip]+sequence[jjp] == 5)        
                                            {
                                                int branch1 = iip-ii-1;
                                                int branch2 = jj-jjp-1;
                                                // compute the derivative for internal penalty by size
                                                    
                                                iipjjp = index[iip] + jjp - iip;  
                                                if (((branch1 == 1 && branch2 > 2) || (branch1 > 2 && branch2 == 1)) && misc.gail_rule)
                                                {
                                                    //char s[100];
                                                    sprintf (s, "internal_penalty_by_size[%d]", branch1+branch2);
                                                    // this is very slow, it's linear in the number of parameters
                                                    int sizeindex = structure_type_index (s);
                                                    //printf ("sizeindex=%d\n", sizeindex);
                                                
                                                    en_internal = s_internal_loop::get_energy (ii, jj, iip, jjp, sequence);
                                                    GlogZ[sizeindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT); 
                                                    
                                                    // now check if any of the 3 internal params is involved
                                                    if (((i == A || i == G) && j == U) ||
                                                        ((j == A || j == G) && i == U))
                                                    {
                                                        //GlogZ[index_terminal_AU_penalty] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT); 
                                                        GlogZ[index_internal_AU_closure] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                    }
                                                    // AG and UU mismatch can't happen here, because k is A and l is A
                                                }
                                            }
                                        }
                                    }                                
                                }

                                // left tstacki, the right one to be added
                                if (sequence[ii] == i && sequence[jj] == j && 
                                    sequence[ii+1] == k && sequence[jj-1] == l)
                                {
                                    iijj = index[ii] + jj -ii;
                                    
                                    for (iip = ii+1; iip <= MIN(jj-2-TURN,ii+MAXLOOP+1) ; iip++)  // j-2-TURN
                                    {
                                        int minq = MAX (jj-ii+iip-MAXLOOP-2, iip+1+TURN);    // ip+1+TURN);
                                        for (jjp = minq; jjp < jj; jjp++)
                                        {        
                                            if (sequence[iip]+sequence[jjp] == 3 ||
                                                sequence[iip]+sequence[jjp] == 5)        
                                            {
                                                int branch1 = iip-ii-1;
                                                int branch2 = jj-jjp-1;
                                                if (branch1 < 0 || branch2 < 0)    // not valid
                                                    continue;
                                                if (branch1 == 0 && branch2 == 0)    // stack pair, not good here
                                                    continue;
                                                if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)
                                                    continue;                                                    
                                                if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
                                                    continue;
                                                if (branch1 == 2 && branch2 == 1 && !simple_internal_energy)
                                                    continue;
                                                if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
                                                    continue;                                                

                                                if (branch1 == 0 || branch2 == 0)    // bulge
                                                {
                                                    // compute the derivative for bulge penalty by size
                                                    //char s[100];
                                                    sprintf (s, "bulge_penalty_by_size[%d]", branch1+branch2);
                                                    // this is very slow, it's linear in the number of parameters
                                                    int sizeindex = structure_type_index (s);
                                                    en_internal = s_internal_loop::get_energy (ii, jj, iip, jjp, sequence);
                                                    iipjjp = index[iip] + jjp - iip;  
                                                    GlogZ[sizeindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT); 
                                                    
                                                    if (branch1 + branch2 > 1)    // add AU_penalty
                                                    {                                                        
                                                        if (has_AU_penalty (sequence[ii], sequence[jj]))
                                                        {
                                                            GlogZ[AUpen_index] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                        }
                                                        if (has_AU_penalty (sequence[iip], sequence[jjp]))
                                                        {
                                                            GlogZ[AUpen_index] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                        }
                                                    }
                                                    continue;
                                                }
                                                    
                                                // compute the derivative for internal penalty by size
                                                //char s[100];
                                                sprintf (s, "internal_penalty_by_size[%d]", branch1+branch2);
                                                // this is very slow, it's linear in the number of parameters
                                                int sizeindex = structure_type_index (s);
                                                //printf ("sizeindex=%d\n", sizeindex);
                                                    
                                                iipjjp = index[iip] + jjp - iip;  
                                                if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
                                                    continue;
                                                    
                                                en_internal = s_internal_loop::get_energy (ii, jj, iip, jjp, sequence);                                                    
                                                GlogZ[sizeindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT); 
                                                
                                                // now check if any of the 3 internal params is involved
                                                if (((i == A || i == G) && j == U) ||
                                                    ((j == A || j == G) && i == U))
                                                {
                                                    //GlogZ[index_terminal_AU_penalty] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT); 
                                                    GlogZ[index_internal_AU_closure] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                }
                                                if ((k == A && l == G) ||
                                                    (l == A && k == G))
                                                {
                                                    GlogZ[index_internal_AG_mismatch] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT); 
                                                }
                                                if (k == U && l == U)
                                                {
                                                    GlogZ[index_internal_UU_mismatch] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT); 
                                                }  
                                            }
                                        }
                                    }                                    
                                }
                            }
                        }           
                        // the upside down orientation
                        // no duplicates here
                        for (ii = seqlen-2; ii > TURN; ii--)
                        {
                            for (jj = ii-TURN-1; jj > 0; jj--)
                            {
                                // first check the weird case of k==0 si l==0
                                if (sequence[ii] == i && sequence[jj] == j && k==0 && l==0)
                                {
                                    // jj < ii
                                    iijj = index[jj] + ii -jj;
                                    
                                    for (iip = ii+1; iip <= MIN(seqlen-1,ii+MAXLOOP+1) ; iip++)  // j-2-TURN
                                    {
                                        int minq = MAX (jj-ii+iip-MAXLOOP-2, 0);    // ip+1+TURN);
                                        for (jjp = minq; jjp < jj; jjp++)
                                        {        
                                            if (sequence[iip]+sequence[jjp] == 3 ||
                                                sequence[iip]+sequence[jjp] == 5)        
                                            {
                                                int branch1 = iip-ii-1;
                                                int branch2 = jj-jjp-1;
                                                // compute the derivative for internal penalty by size
                                                    
                                                iipjjp = index[jjp] + iip - jjp; 
                                                if (((branch1 == 1 && branch2 > 2) || (branch1 > 2 && branch2 == 1)) && misc.gail_rule)
                                                {                                                
                                                    en_internal = s_internal_loop::get_energy (jjp, iip, jj, ii, sequence);               
                                                    // now check if any of the 3 internal params is involved
                                                    if (((i == A || i == G) && j == U) ||
                                                        ((j == A || j == G) && i == U))
                                                    {
                                                        //GlogZ[index_terminal_AU_penalty] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT); 
                                                        GlogZ[index_internal_AU_closure] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT); 
                                                    }
                                                    // AG and UU mismatch can't happen here, because k is A and l is A                                               
                                                }
                                            }
                                        }
                                    }                                
                                }                            
                            
                                // the general case
                                // left tstacki, the right one to be added
                                if (sequence[ii] == i && sequence[jj] == j && 
                                    sequence[ii+1] == k && sequence[jj-1] == l)
                                {
                                    // jj < ii
                                    iijj = index[jj] + ii -jj;
                                    for (iip = ii+1; iip <= MIN(seqlen-1,ii+MAXLOOP+1) ; iip++)  // j-2-TURN
                                    {
                                        int minq = MAX (jj-ii+iip-MAXLOOP-2, 0);    // ip+1+TURN);
                                        for (jjp = minq; jjp < jj; jjp++)
                                        {        
                                            if (sequence[iip]+sequence[jjp] == 3 ||
                                                sequence[iip]+sequence[jjp] == 5)        
                                            {
                                                // same as before
                                                int branch1 = iip-ii-1;
                                                int branch2 = jj-jjp-1;
                                                if (branch1 < 0 || branch2 < 0)    // not valid
                                                    continue;
                                                if (branch1 == 0 && branch2 == 0)    // stack pair, not good here
                                                    continue;
                                                if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
                                                    continue;
                                                if (branch1 == 2 && branch2 == 1 && !simple_internal_energy)
                                                    continue;
                                                if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
                                                    continue;
                                                if (branch1 == 0 || branch2 == 0)    // bulge
                                                    continue;
                                                    
                                                //jjp < iip
                                                iipjjp = index[jjp] + iip - jjp;  
                                                if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
                                                    continue;
                                                en_internal = s_internal_loop::get_energy (jjp, iip, jj, ii, sequence);                                                    
                                                
                                                // now check if any of the 3 internal params is involved
                                                if (((i == A || i == G) && j == U) ||
                                                    ((j == A || j == G) && i == U))
                                                {
                                                    //GlogZ[index_terminal_AU_penalty] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT); 
                                                    GlogZ[index_internal_AU_closure] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT); 
                                                }
                                                if ((k == A && l == G) ||
                                                    (l == A && k == G))
                                                {
                                                    GlogZ[index_internal_AG_mismatch] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT); 
                                                }
                                                if (k == U && l == U)
                                                {
                                                    GlogZ[index_internal_UU_mismatch] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT); 
                                                }                                                 
                                            }
                                        }
                                    }                                    
                                }
                            }
                        }
                        //index_param++;   // there are only 3 tstacki parameters                              
                    }    // end if (tstacki[i][j][k][l] < INF)                                            
                }
                
    index_param += 3;           
    index_should_be = structure_type_index("int11[0][3][3][3][0][3]");     
    if (index_param != index_should_be)
    {
        printf ("Index param after tstacki = %d, should be %d\n", index_param, index_should_be);
        exit(1);
    }                
                    
    // internal loops 1x1                        
    if (!simple_internal_energy)          
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                            {                            
                                if (int11[i][j][k][l][m][n] < INF)
                                {
                                    // exclude duplicates
                                    // int11[i][j][k][l][m][n] is the same as int11[n][m][l][k][j][i]
                                    if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n <= n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
                                    {             
                                       
                                        // look from the left
                                        for (ii = 0; ii < seqlen-TURN-5; ii++)
                                        {
                                            if (sequence[ii] == i && sequence[ii+1]==k && sequence[ii+2] == m)
                                            {
                                                for (jj = ii + TURN+5; jj < seqlen; jj++)
                                                {
                                                    if (sequence[jj] == j && sequence[jj-1] == l && sequence[jj-2] == n)
                                                    {                      
                                                        iijj = index[ii] + jj -ii;    
                                                        iipjjp = index[ii+2] + jj-2 - (ii+2); 
                                                        en_internal = s_internal_loop::get_energy (ii, jj, ii+2, jj-2, sequence);                                             
                              
                                                        int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc;
                                                        iinuc = sequence[ii]; 
                                                        jjnuc = sequence[jj];
                                                        kknuc = sequence[ii+1];
                                                        llnuc = sequence[jj-1];
                                                        mmnuc = sequence[ii+2];
                                                        nnnuc = sequence[jj-2];
                                        
                                                        if ( ((iinuc==C && jjnuc==G) || (iinuc==G && jjnuc==C)) && ((mmnuc==C && nnnuc==G) || (mmnuc==G && nnnuc==C))) 
                                                        {
                                                            if (!can_pair(kknuc,llnuc))
                                                                sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc);
                                                            else
                                                                sprintf (type, "misc.internal11_basic_mismatch");
                                                            tindex = structure_type_index (type);                                                            
                                                            GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                        }        
                                                        else if (watson_crick(iinuc,jjnuc) && watson_crick(mmnuc,nnnuc) && kknuc==U && llnuc==U)
                                                        {
                                                            sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc);
                                                            tindex = structure_type_index (type);
                                                            GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                        }
                                                        else
                                                        {
                                                            if (kknuc==G && llnuc==G)
                                                                sprintf (type, "misc.internal11_GG_mismatch");
                                                            else
                                                                sprintf (type, "misc.internal11_basic_mismatch");
                                                            tindex = structure_type_index (type);
                                                            GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                                                        
                                                            if (has_AU_penalty(iinuc,jjnuc))
                                                            {
                                                                sprintf (type, "misc.internal_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);                                           
                                                            }
                                                            if (has_AU_penalty(mmnuc,nnnuc))
                                                            {
                                                                sprintf (type, "misc.internal_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                            }
                                                        }                                                                                                       
                                                    }    
                                                }
                                            }                   
                                        }
                        
                                        // look from the right, if it's not symmetric
                                        if (!(i==n && k==l && m==j))
                                        {
                                            for (ii = seqlen-3; ii >= TURN+3; ii--)
                                            {
                                                if (sequence[ii] == i && sequence[ii+1]==k && sequence[ii+2] == m)
                                                {
                                                    for (jj = ii-TURN-1; jj >= 2; jj--)
                                                    {
                                                        if (sequence[jj] == j && sequence[jj-1] == l && sequence[jj-2] == n)
                                                        {
                                                            //jj < ii
                                                            iijj = index[jj] + ii -jj;    
                                                            iipjjp = index[jj-2] + ii+2 - (jj-2); 
                                                            en_internal = s_internal_loop::get_energy (jj-2, ii+2, jj, ii, sequence); 
                                                            
                                                            int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc;
                                                            iinuc = sequence[ii]; 
                                                            jjnuc = sequence[jj];
                                                            kknuc = sequence[ii+1];
                                                            llnuc = sequence[jj-1];
                                                            mmnuc = sequence[ii+2];
                                                            nnnuc = sequence[jj-2];
                                            
                                                            if ( ((iinuc==C && jjnuc==G) || (iinuc==G && jjnuc==C)) && ((mmnuc==C && nnnuc==G) || (mmnuc==G && nnnuc==C))) 
                                                            {
                                                                if (!can_pair(kknuc,llnuc))
                                                                    sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc);
                                                                else
                                                                    sprintf (type, "misc.internal11_basic_mismatch");
                                                                tindex = structure_type_index (type);                                                            
                                                                GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                            }        
                                                            else if (watson_crick(iinuc,jjnuc) && watson_crick(mmnuc,nnnuc) && kknuc==U && llnuc==U)
                                                            {
                                                                sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc);
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                            }
                                                            else
                                                            {
                                                                if (kknuc==G && llnuc==G)
                                                                    sprintf (type, "misc.internal11_GG_mismatch");
                                                                else
                                                                    sprintf (type, "misc.internal11_basic_mismatch");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                                                            
                                                                if (has_AU_penalty(iinuc,jjnuc))
                                                                {
                                                                    sprintf (type, "misc.internal_AU_closure");
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);  
                                                                }
                                                                if (has_AU_penalty(mmnuc,nnnuc))
                                                                {
                                                                    sprintf (type, "misc.internal_AU_closure");
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                                }
                                                            }                                                                                                       
                                                                
                                                        }    
                                                    }
                                                }                   
                                            }
                                        }                         
                                    }                                            
                                } // end if (int11[i][j][k][l][m][n] < INF)
                            }    // end int11
        index_param += 33;
        index_should_be = structure_type_index("int21[1][2][0][0][1][2][0]");
        if (index_param != index_should_be)
        {
            printf ("Index param after int11 = %d, should be %d\n", index_param, index_should_be);
            exit(1);
        }                 
        
        
        // int12 energies
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                                for (o=0; o < NUCL; o++)
                                {
                                    if (int21[i][j][k][l][m][n][o] < INF)
                                    {
                                        // no duplicates here
                                        // look from the left
                                        for (ii = 0; ii < seqlen-TURN-6; ii++)
                                        {
                                            if (sequence[ii] == i && sequence[ii+1]==k && sequence[ii+2] == m)
                                            {
                                                for (jj = ii + TURN+6; jj < seqlen; jj++)
                                                {
                                                    if (sequence[jj] == j && sequence[jj-1] == l && sequence[jj-2] == o && sequence[jj-3] == n)
                                                    {
                                                        iijj = index[ii] + jj -ii;    
                                                        iipjjp = index[ii+2] + jj-3 - (ii+2); 
                                                        en_internal = s_internal_loop::get_energy (ii, jj, ii+2, jj-3, sequence);  
                                                        
                                                        int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc;
                                                        iinuc = sequence[ii];
                                                        jjnuc = sequence[jj];
                                                        kknuc = sequence[ii+1];
                                                        llnuc = sequence[jj-1];
                                                        mmnuc = sequence[ii+2];
                                                        nnnuc = sequence[jj-3];
                                                        oonuc = sequence[jj-2];            
                                            
                                                        if ((iinuc==C && jjnuc==G && mmnuc==C && nnnuc==G) ||  // these are already filled above, except what can pair inside
                                                            (iinuc==G && jjnuc==C && mmnuc==G && nnnuc==C))
                                                        {
                                                            if (can_pair(kknuc,llnuc) || can_pair(kknuc,oonuc))
                                                                sprintf (type, "misc.internal21_match");
                                                            else
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc);
                                                            tindex = structure_type_index (type);
                                                            GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                        }
                                                        else
                                                        {
                                                            if (can_pair(kknuc,llnuc) || can_pair(kknuc,oonuc))
                                                            {
                                                                sprintf (type, "misc.internal21_match");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                            }
                                                            else
                                                            {
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", C, G, kknuc, llnuc, C, G, oonuc);
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += 0.5 * p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);             
                                                                
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", G, C, kknuc, llnuc, G, C, oonuc);
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += 0.5 * p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                            }    
                                                            if (has_AU_penalty(iinuc,jjnuc))
                                                            {
                                                                sprintf (type, "misc.internal21_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);             
                                                            }    
                                                            if (has_AU_penalty(mmnuc,nnnuc))    
                                                            {
                                                                sprintf (type, "misc.internal21_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);                 
                                                            }    
                                                        }                                                                                                                                                    
                                                    }    
                                                }
                                            }                   
                                        }
                        
                                        // int21 cannot be symmetric
                                        for (ii = seqlen-3; ii >= TURN+4; ii--)
                                        {
                                            if (sequence[ii] == i && sequence[ii+1]==k && sequence[ii+2] == m)
                                            {
                                                for (jj = ii-TURN-1; jj >= 3; jj--)
                                                {
                                                    if (sequence[jj] == j && sequence[jj-1] == l && sequence[jj-2] == o && sequence[jj-3] == n)
                                                    {
                                                        //jj < ii
                                                        iijj = index[jj] + ii -jj;    
                                                        iipjjp = index[jj-3] + ii+2 - (jj-3); 
                                                        en_internal = s_internal_loop::get_energy (jj-3, ii+2, jj, ii, sequence);
                                                        
                                                        int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc;
                                                        iinuc = sequence[ii];
                                                        jjnuc = sequence[jj];
                                                        kknuc = sequence[ii+1];
                                                        llnuc = sequence[jj-1];
                                                        mmnuc = sequence[ii+2];
                                                        nnnuc = sequence[jj-3];
                                                        oonuc = sequence[jj-2];            
                                            
                                                        if ((iinuc==C && jjnuc==G && mmnuc==C && nnnuc==G) ||  // these are already filled above, except what can pair inside
                                                            (iinuc==G && jjnuc==C && mmnuc==G && nnnuc==C))
                                                        {
                                                            if (can_pair(kknuc,llnuc) || can_pair(kknuc,oonuc))
                                                                sprintf (type, "misc.internal21_match");
                                                            else
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc);
                                                            tindex = structure_type_index (type);
                                                            GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                        }
                                                        else
                                                        {
                                                            if (can_pair(kknuc,llnuc) || can_pair(kknuc,oonuc))
                                                            {
                                                                sprintf (type, "misc.internal21_match");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                            }
                                                            else
                                                            {
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", C, G, kknuc, llnuc, C, G, oonuc);
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += 0.5 * p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                                
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", G, C, kknuc, llnuc, G, C, oonuc);
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += 0.5 * p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                            }
                                                            if (has_AU_penalty(iinuc,jjnuc))
                                                            {
                                                                sprintf (type, "misc.internal21_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                            }    
                                                            if (has_AU_penalty(mmnuc,nnnuc))    
                                                            {
                                                                sprintf (type, "misc.internal21_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                            }    
                                                        }                
                                                    }    
                                                }
                                            }                   
                                        }
                                    //index_param++;                            
                                    }                                            
                                } // end int21
                                
        index_param += 54;                                
        index_should_be = structure_type_index("int22[0][3][0][0][3][0][0][0]");                                
        if (index_param != index_should_be)
        {
            printf ("Index param after int21 = %d, should be %d\n", index_param, index_should_be);
            exit(1);
        }                                     
        // energies int22                                
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                                for(o=0; o < NUCL; o++)
                                    for (r=0; r < NUCL; r++)
                                    {
                                        if (int22[i][j][k][l][m][n][o][r] < INF)
                                        {
                                            // exclude duplicates
                                            // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
                                            if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + r <= n*10000000 + m*1000000 + r*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                                            {
                                                // look from the left
                                                for (ii = 0; ii < seqlen-TURN-7; ii++)
                                                {
                                                    if (sequence[ii]==i && sequence[ii+1]==k && sequence[ii+2]==o && sequence[ii+3]==m)
                                                    {
                                                        for (jj = ii + TURN+7; jj < seqlen; jj++)
                                                        {
                                                            if (sequence[jj]==j && sequence[jj-1]==l && sequence[jj-2]==r && sequence[jj-3]==n)
                                                            {
                                                                iijj = index[ii] + jj -ii;    
                                                                iipjjp = index[ii+3] + jj-3 - (ii+3); 
                                                                en_internal = s_internal_loop::get_energy (ii, jj, ii+3, jj-3, sequence);                                                                   
                                                                
                                                                int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc, ppnuc;
                                                                iinuc = sequence[ii];
                                                                jjnuc = sequence[jj];
                                                                kknuc = sequence[ii+1];
                                                                llnuc = sequence[jj-1];
                                                                mmnuc = sequence[ii+3];
                                                                nnnuc = sequence[jj-3];
                                                                oonuc = sequence[ii+2];            
                                                                ppnuc = sequence[jj-2];
                                                                
                                                                if (nnnuc==iinuc && mmnuc==jjnuc && ppnuc==kknuc && oonuc==llnuc & watson_crick(iinuc,jjnuc) && !watson_crick(kknuc,llnuc))
                                                                {
                                                                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc, ppnuc);
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                                }
                                                            
                                                                int iinuc2, jjnuc2, mmnuc2, nnnuc2;                
                                                                if (iinuc==G && jjnuc==U)   iinuc2 = A;     else iinuc2 = iinuc;
                                                                if (iinuc==U && jjnuc==G)   jjnuc2 = A;     else jjnuc2 = jjnuc;
                                                                if (mmnuc==G && nnnuc==U)   mmnuc2 = A;     else mmnuc2 = mmnuc;
                                                                if (mmnuc==U && nnnuc==G)   nnnuc2 = A;     else nnnuc2 = nnnuc;
                                                                
                                                                if (watson_crick(kknuc,llnuc) || watson_crick(oonuc,ppnuc))
                                                                {
                                                                    sprintf (type, "misc.internal22_match");
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                                }
                                                                else if ( ((iinuc==G && jjnuc==U) || (iinuc==U && jjnuc==G) || (mmnuc==G && nnnuc==U) || (mmnuc==U && nnnuc==G)) &&  
                                                                        (nnnuc2==iinuc2 && mmnuc2==jjnuc2 && ppnuc==kknuc && oonuc==llnuc))  // the UG closing pairs are the same as UA
                                                                {
                                                                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc2, jjnuc2, kknuc, llnuc, mmnuc2, nnnuc2, oonuc, ppnuc);
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                                }                
                                                                else if (!(nnnuc==iinuc && mmnuc==jjnuc && ppnuc==kknuc && oonuc==llnuc))   // was already filled above
                                                                {
                                                                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc2, jjnuc2, kknuc, llnuc, jjnuc2, iinuc2, llnuc, kknuc);
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += 0.5 * p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", nnnuc2, mmnuc2, ppnuc, oonuc, mmnuc2, nnnuc2, oonuc, ppnuc);
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += 0.5 * p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                                    
                                                                    int result = check_stability_and_size (kknuc, llnuc, oonuc, ppnuc);                    
                                                                    switch (result)
                                                                    {
                                                                        case 1: 
                                                                            sprintf (type, "misc.internal22_delta_same_size");
                                                                            break;
                                                                        case 2: 
                                                                            sprintf (type, "misc.internal22_delta_different_size");
                                                                            break;
                                                                        case 3: 
                                                                            sprintf (type, "misc.internal22_delta_1stable_1unstable");
                                                                            break;
                                                                        case 4: 
                                                                            sprintf (type, "misc.internal22_delta_AC");
                                                                            break;
                                                                        default: 
                                                                            printf ("ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, kknuc, llnuc, oonuc, ppnuc);
                                                                            exit(1);
                                                                    }
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                                }                                                                
                                                            }    
                                                        }
                                                    }                   
                                                }
                                
                                                // look from the right, if it's not symmetric
                                                if (!(i==n && k==r && o==l && m==j))
                                                {
                                                    for (ii = seqlen-4; ii >= TURN+4; ii--)
                                                    {
                                                        if (sequence[ii]==i && sequence[ii+1]==k && sequence[ii+2]==o && sequence[ii+3]==m)
                                                        {
                                                            for (jj = ii-TURN-1; jj >= 3; jj--)
                                                            {
                                                                if (sequence[jj]==j && sequence[jj-1]==l && sequence[jj-2]==r && sequence[jj-3]==n)
                                                                {
                                                                    //jj < ii
                                                                    iijj = index[jj] + ii -jj;    
                                                                    iipjjp = index[jj-3] + ii+3 - (jj-3); 
                                                                    en_internal = s_internal_loop::get_energy (jj-3, ii+3, jj, ii, sequence);
                                                                    
                                                                    int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc, ppnuc;
                                                                    iinuc = sequence[ii];
                                                                    jjnuc = sequence[jj];
                                                                    kknuc = sequence[ii+1];
                                                                    llnuc = sequence[jj-1];
                                                                    mmnuc = sequence[ii+3];
                                                                    nnnuc = sequence[jj-3];
                                                                    oonuc = sequence[ii+2];            
                                                                    ppnuc = sequence[jj-2];
                                                                    
                                                                    if (nnnuc==iinuc && mmnuc==jjnuc && ppnuc==kknuc && oonuc==llnuc & watson_crick(iinuc,jjnuc) && !watson_crick(kknuc,llnuc))
                                                                    {
                                                                        sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc, ppnuc);
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                                    }
                                                                
                                                                    int iinuc2, jjnuc2, mmnuc2, nnnuc2;                
                                                                    if (iinuc==G && jjnuc==U)   iinuc2 = A;     else iinuc2 = iinuc;
                                                                    if (iinuc==U && jjnuc==G)   jjnuc2 = A;     else jjnuc2 = jjnuc;
                                                                    if (mmnuc==G && nnnuc==U)   mmnuc2 = A;     else mmnuc2 = mmnuc;
                                                                    if (mmnuc==U && nnnuc==G)   nnnuc2 = A;     else nnnuc2 = nnnuc;
                                                                    
                                                                    if (watson_crick(kknuc,llnuc) || watson_crick(oonuc,ppnuc))
                                                                    {
                                                                        sprintf (type, "misc.internal22_match");
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                                    }
                                                                    else if ( ((iinuc==G && jjnuc==U) || (iinuc==U && jjnuc==G) || (mmnuc==G && nnnuc==U) || (mmnuc==U && nnnuc==G)) &&  
                                                                            (nnnuc2==iinuc2 && mmnuc2==jjnuc2 && ppnuc==kknuc && oonuc==llnuc))  // the UG closing pairs are the same as UA
                                                                    {
                                                                        sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc2, jjnuc2, kknuc, llnuc, mmnuc2, nnnuc2, oonuc, ppnuc);
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                                    }                
                                                                    else if (!(nnnuc==iinuc && mmnuc==jjnuc && ppnuc==kknuc && oonuc==llnuc))   // was already filled above
                                                                    {
                                                                        sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc2, jjnuc2, kknuc, llnuc, jjnuc2, iinuc2, llnuc, kknuc);
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += 0.5 * p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                                        sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", nnnuc2, mmnuc2, ppnuc, oonuc, mmnuc2, nnnuc2, oonuc, ppnuc);
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += 0.5 * p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                                        
                                                                        int result = check_stability_and_size (kknuc, llnuc, oonuc, ppnuc);                    
                                                                        switch (result)
                                                                        {
                                                                            case 1: 
                                                                                sprintf (type, "misc.internal22_delta_same_size");
                                                                                break;
                                                                            case 2: 
                                                                                sprintf (type, "misc.internal22_delta_different_size");
                                                                                break;
                                                                            case 3: 
                                                                                sprintf (type, "misc.internal22_delta_1stable_1unstable");
                                                                                break;
                                                                            case 4: 
                                                                                sprintf (type, "misc.internal22_delta_AC");
                                                                                break;
                                                                            default: 
                                                                                printf ("ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, kknuc, llnuc, oonuc, ppnuc);
                                                                                exit(1);
                                                                        }
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                                    }
                                                                }    
                                                            }
                                                        }                   
                                                    }
                                                }    
                                            //index_param++;                            
                                            }                                            
                                        } // end if (int11[i][j][k][l][m][n] < INF)
                                    }    // end int11
                                        
        index_param += 53;                                    
        index_should_be = structure_type_index("dangle_top[0][3][0]");                                                
        if (index_param != index_should_be)
        {
            printf ("Index param after int22 = %d, should be %d\n", index_param, index_should_be);
            exit(1);
        }              
                                    
    }    // end int11, int12, int22,     // end if (!simple_internal_energy)     
    

    if (!ignore_dangles)    // we don't need to compute dangling ends
    {
    
    // dangling 3', or dangle_top    
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
            {
                if (dangle_top[i][j][k] < INF)
                {
                    // no duplicates here                    
                    // traverse the whole sequence and look for this building block                    

                    // it never gets here but anyway
                    IFD
                    {
                        GlogZ[index_param] = 0.0;
                        index_param++;
                        continue;
                    }
                    
                    // first add the dangling energies that participate to the external loop                     
                    for (ii = 0; ii < seqlen; ii++)
                    {
                        for (jj = ii+TURN+1; jj < seqlen-1; jj++)
                        {
                            if (sequence[jj] == i && sequence[ii] == j && sequence[jj+1] == k)
                            {            
                                double term1 = 0.0;
                                double term2 = 0.0;
                                                            
                                iijj = index[ii]+jj-ii;
                                
                                if (ii > TURN)        
                                    term1 += u_ip_jp[ii-1] + u_iu_jp[ii-1] + exp_dangle5(jj,ii,ii-1)*(u_ip_ju[ii-1] + u_iu_ju[ii-1]);                                  
                                else if (ii > 0)      term1 += exp_dangle5(jj,ii,ii-1);
                                else                  term1 = 1;  
    
                                if (jj < seqlen-3)
                                {
                                    int jjp2n = index[jj+2]+seqlen-1-(jj+2);
                                    term2 += 
                                            exp_dangle3(jj,ii,jj+1) * (u_ip_jp[jjp2n] + u_ip_ju[jjp2n] + u_iu_jp[jjp2n] + u_iu_ju[jjp2n]);
                                }
                                else if (jj < seqlen-1) term2 += exp_dangle3(jj,ii,jj+1);
                                // it will never get here
                                //else                  term2 = 1;
                                
                                GlogZ[index_param] += term1*up[iijj]*term2* exp_AUpenalty (ii,jj) /
                                    (u_ip_jp[seqlen-1] + u_ip_ju[seqlen-1] + u_iu_jp[seqlen-1] + u_iu_ju[seqlen-1]) ;
                            }
                        }
                    }
                    // contribution from multi-loops
                    // the first free bases
                    // taken from compute_upm  
                                      
                    double temp = 0;                    
                    for (jj=TURN+1; jj < seqlen; jj++)
                    {
                        for (ii=jj-TURN-1; ii>=0; ii--)
                        {                     
                            iijj = index[ii] + jj -ii;  
                            if (sequence[ii] == i && sequence[jj] == j && sequence[ii+1] == k)
                            {                 
                                temp = 0;
                                if (upm[iijj] == 0)
                                    continue;                                                              
                                for (l=ii+3; l < jj-TURN-2; l++)    // case (.(...)--(--)-)
                                {
                                    int iip2l = index[ii+2]+l-ii-2;
                                    int lp2jjm1 = index[l+2]+jj-1-l-2;
                                    int lp1jjm1 = index[l+1]+jj-1-l-1;
                                    int lp2jjm2 = index[l+2]+jj-2-l-2;
                                    int lp1jjm2 = index[l+1]+jj-2-l-1;        
                                    int lp2jjm3 = index[l+2]+jj-3-l-2;
                                    int lp1jjm3 = index[l+1]+jj-3-l-1;                
                                    temp += exp_AUpenalty (ii,jj) * up[iip2l] * 
                                            EXPA * EXPB[2] * EXPC[1] *
                                            exp_dangle3 (ii, jj, ii+1) * exp_AUpenalty (ii+2, l) *   
                                            ( u1_ip_jp[lp1jjm1]      // [.(...)(.-..)] or [.(...)(.-..).]                    
                                                + exp_dangle3 (l, ii+2, l+1) * EXPC[1] *
                                                    ( u1_ip_jp[lp2jjm1] + u1_iu_jp[lp2jjm1] +     // [.(...).-(...)] or [.(...).-(...).] 
                                                       exp_dangle5(ii,jj,jj-1) * (u1_ip_ju[lp2jjm1] + u1_iu_ju[lp2jjm1])) +    // [.(...).-(...)-..]
                                                + exp_dangle5(ii,jj,jj-1) * u1_ip_ju[lp1jjm1]);    // [.(...)(...)-..]                                             
                                }
                                double upm_temp = 0;
                                for (h=ii+3; h < jj-TURN-2; h++)    // case (....(...)--(--)-)
                                {
                                    int hjj = index[h]+jj-h;
                                    int hjjm1 = index[h] + jj-1 -h;
                                    int hjjm2 = index[h] + jj-2 -h;
                                    int hjjm3 = index[h] + jj-3 -h;
                                    upm_temp += EXPA * EXPB[2] * EXPC[h-ii-1] *
                                        ( s2_jp[hjjm1]    // --(...)) or --(...).) 
                                          + s2_ju[hjjm1]* exp_dangle5(ii,jj,jj-1));    // --(...)..)
                                    //s2[hj];
                                }    
                                temp += exp_AUpenalty (ii,jj) * exp_dangle3 (ii, jj, ii+1) * upm_temp;
                                //printf ("i=%d, j=%d, upm=%g\n", i, j, upm[ij]);

                                GlogZ[index_param] += temp * p[iijj] / up[iijj];                                
                            }
                        }
                    }
                    

                    // the multi-loop 3' dangling end which can appear to the right of branches
                    // taken from compute_p
                    temp = 0;
                    for (h=0; h < seqlen; h++)   
                    {
                        for (l=seqlen-2; l > h+TURN; l--)    // it's seqlen-2 because l+1 must be valid
                        {
                            //printf ("h=%d, l=%d, seqlen=%d\n", h, l, seqlen);
                            if (can_pair (sequence[h], sequence[l]) && 
                                sequence[l] == i && sequence[h] == j && sequence[l+1] == k)    
                            {                    
                                for (ii=0; ii < h; ii++)
                                {
                                    int iil = index[ii] + l - ii;
                                    int hl = index[h] + l - h;
                                    // the case when h-l is the first branch of the multi-loop
                                    if (ii < h-1)  // we must add d3 (which is added in pmd3); also add d5 here, if i < h-2
                                    {
                                        // removed one line, no dangling end of interest
                                        GlogZ[index_param] += up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                                            EXPC[h-ii-1] * pmd3_needmidd3[iil]  * exp_dangle3(l,h,l+1) * EXPC[1] *
                                                    (ii<h-2?exp_dangle5(l,h,h-1):1);    
                                            // (.-[...].-(---)-)    i.-h...l.-(---)-j                 
                                    }
                                    else    // case ((... , no dangling end
                                    {
                                        // removed one line, no dangling end of interest
                                        GlogZ[index_param] += up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                                            ( pmnod3_needmidd3[iil] * exp_dangle3(l,h,l+1) * EXPC[1] * EXPC[h-ii-1]);
                                            // ([...].-(---)-)    ih...l.-(---)-j                 
                                    }
                                    
                                    if (ii < h - TURN -2)    // only now we can have a branch to the left
                                    {        
                                        int iip1hm1 = index[ii+1] + h-1 - (ii+1);
                                        int iip1hm2 = index[ii+1] + h-2 - (ii+1);
                                        int iip1hm3 = index[ii+1] + h-3 - (ii+1);
                                        int iip2hm1 = index[ii+2] + h-1 - (ii+2);
                                        int iip2hm2 = index[ii+2] + h-2 - (ii+2);
                                        int iip2hm3 = index[ii+2] + h-3 - (ii+2);
                                        int iip3hm1 = index[ii+3] + h-1 - (ii+3);
                                        int iip3hm2 = index[ii+3] + h-2 - (ii+3);
                                        int iip3hm3 = index[ii+3] + h-3 - (ii+3);
                                        int iilp1 = index[ii] + l+1 - ii;
                                        // removed term1, no dangling end of interest
                                        
                                        GlogZ[index_param] += up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                                        
                                        // first, when h-l is the last branch to the right
                                        (                                                    
                                        // case ((..-)-[...].-)  
                                        pm1nod3_needendd3[iil] * exp_dangle3(l,h,l+1) *     // EXPC added in pm1nod3
                                            (u1_ip_jp[iip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[iip1hm1])
                                                    // c is added in pm1nod3_needendd3
                                                    // ((...).-[...].-)         i(...)h...l.-j
                                        // case (.(..-)-[...].-)
                                        + pm1d3_needendd3[iil] * exp_dangle3(l,h,l+1) * //  don't need EXPC[1]
                                            (u1_ip_jp[iip2hm1] + u1_iu_jp[iip2hm1] +
                                                exp_dangle5(l,h,h-1)*(u1_ip_ju[iip2hm1] + u1_iu_ju[iip2hm1]))
                                                    // (.(...).-[...].-)        i.(...)h...l.-j                                                    
                                        // case (..-(..-)-[...].-)            
                                                    // (..-(...)[...].-)      i..-(...)h...l.-j
                                                    
                                        // we have branches to the left and to the right of h-l
                                        // if this is changed in compute_p, it should be changed here too
                                        //+ pmnod3_noneedmidd3[iil]  * u1[iip1hm1]     
                                                    // ((...)[...]-(---)...)    i(...)h...l-(---)...j 
                                        //+ pmnod3_needmidd3[iil] * exp_dangle3(l,h,l+1) * u1[iip1hm1]);        
                                                    // ((...).-[...]-(---)...)    i(...).-h...l-(---)...j 
                                                    
                                        + pmnod3_needmidd3[iil] * exp_dangle3(l,h,l+1) * EXPC[1]*    //+ pmnod3_noneedmidd3[il]) *
                                            (u1_ip_jp[iip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[iip1hm1])
                                                    // c is added in pm1nod3_needendd3
                                                    // ((...).-[...].-(--)-)         i(...)h...l.-j
                                        // case (.(..-)-[...].-(--)-) or (.(..-)-[...](--)-)
                                        + EXPC[1]*(pmd3_needmidd3[iil] * exp_dangle3(l,h,l+1) * EXPC[1])*     //+ pmd3_noneedmidd3[il]) *
                                                (u1_ip_jp[iip2hm1] + u1_iu_jp[iip2hm1] +
                                                    exp_dangle5(l,h,h-1)*(u1_ip_ju[iip2hm1]+u1_iu_ju[iip2hm1])));
                                                    // (.(...).-[...].-)        i.(...)h...l.-j
                                                    
                                        // case (..-(..-)-[...].-(--)-) or (..-(..-)-[...](--)-)
                                    }
                                }     
                            }
                        }
                    }
                    
                    index_param++;
                }
            }
    index_should_be = structure_type_index("dangle_bot[0][3][0]");            
    if (index_param != index_should_be)
    {
        printf ("Index param after dangle_top = %d, should be %d\n", index_param, index_should_be);
        exit(1);
    }      

    // dangling 5', or dangle_bot
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
            {
                if (dangle_bot[i][j][k] < INF)
                {
                    // it never gets here
                    IFD
                    {
                        GlogZ[index_param] = 0.0;
                        index_param++;
                        continue;
                    }
                    // no duplicates here
                    // first add the dangling energies that participate to the external loop
                    // traverse the whole sequence and look for this building block
                    
                    
                    for (ii = 1; ii < seqlen; ii++)
                    {
                        for (jj = ii+TURN+1; jj < seqlen; jj++)
                        {
                            if (sequence[jj] == i && sequence[ii] == j && sequence[ii-1] == k)
                            {                            
                                double term1 = 0.0;
                                double term2 = 0.0;                            
                                iijj = index[ii]+jj-ii;                                
                                if (ii > TURN)        
                                        // case i-1 paired
                                    term1 += exp_dangle5(jj,ii,ii-1)*(u_ip_ju[ii-1] + u_iu_ju[ii-1]);
                                else if (ii > 0)      term1 += exp_dangle5(jj,ii,ii-1);
                                // this case doesn't occur anyway
                                //else                  term1 = 1;  
                                if (jj < seqlen-3)
                                {
                                    int jjp1n = index[jj+1]+seqlen-1-(jj+1);
                                    int jjp2n = index[jj+2]+seqlen-1-(jj+2);
                                    term2 += (u_ip_jp[jjp1n] + u_ip_ju[jjp1n] + 
                                        exp_dangle3(jj,ii,jj+1)*(u_ip_jp[jjp2n] + u_ip_ju[jjp2n] + u_iu_jp[jjp2n] + u_iu_ju[jjp2n]));
                                }    
                                else if (jj < seqlen-1) term2 += exp_dangle3(jj,ii,jj+1);
                                else                  term2 = 1;
                                
                                GlogZ[index_param] += term1 * up[iijj] * term2 * exp_AUpenalty (ii,jj) /
                                    (u_ip_jp[seqlen-1] + u_ip_ju[seqlen-1] + u_iu_jp[seqlen-1] + u_iu_ju[seqlen-1]);
                                //printf ("GlogZ[%d]=%g\n", index_param, GlogZ[index_param]);
                            }
                        }
                    }                    
                    
                    
                    // contribution from multi-loops
                    // the last dangling 5'
                    // taken from compute_upm                    
                    double temp = 0;                    
                    
                    for (jj=TURN+1; jj < seqlen; jj++)
                    {
                        for (ii=jj-TURN-1; ii>=0; ii--)
                        {                     
                            iijj = index[ii] + jj -ii;  
                            if (upm[iijj] == 0)
                                continue;
                            if (sequence[ii] == i && sequence[jj] == j && sequence[jj-1] == k)
                            {       
                                temp = 0;
                                // there must be at least one more branch (i.e. 5 nucleotides) from l+1 to j-1, plus at least 2 free bases: l+1+TURN <= j-3
                                // removed the cases where there isn't a 5' dangling end on the last free base
                                for (l=ii+2; l < jj-TURN-4; l++)    // case ((...)--(--)-..)
                                {
                                    int iip1l = index[ii+1]+l-ii-1;
                                    int lp1jjm1 = index[l+1]+jj-1-l-1;        
                                    int lp2jjm1 = index[l+2]+jj-1-l-2;
                                    temp += up[iip1l] * exp_AUpenalty (ii+1,l) *
                                            (  //u1_ip_jp[lp1jjm1] +     // no dangling end of interest [(...)(.-..)] or [(...)(.-..).]                    
                                             exp_dangle3 (l, ii+1, l+1) * EXPC[1] *
                                                    ( //u1_ip_jp[lp2jjm1] + u1_iu_jp[lp2jjm1] +     // [(...).-(...)] or [(...).-(...).] 
                                                        exp_dangle5(ii,jj,jj-1) * (u1_ip_ju[lp2jjm1] + u1_iu_ju[lp2jjm1])) +    // [(...).-(...)-..]
                                            + exp_dangle5(ii,jj,jj-1) * u1_ip_ju[lp1jjm1]);    // [(...)(...)-..]
                                            // TODO shouldn't there be an EXPC[1] here in the last line?
                                }   
                                                             
                                for (l=ii+3; l < jj-TURN-4; l++)    // case (.(...)--(--)-..)
                                {
                                    int iip2l = index[ii+2]+l-ii-2;
                                    int lp1jjm1 = index[l+1]+jj-1-l-1;        
                                    int lp2jjm1 = index[l+2]+jj-1-l-2;
                                    temp += up[iip2l] * EXPC[1] *
                                            exp_dangle3 (ii, jj, ii+1) * exp_AUpenalty (ii+2, l) * 
                                            (  //u1_ip_jp[lp1jjm1] +     // no dangling end of interest [.(...)(.-..)] or [.(...)(.-..).]                    
                                              exp_dangle3 (l, ii+2, l+1) * EXPC[1] *
                                                    ( //u1_ip_jp[lp2jjm1] + u1_iu_jp[lp2jjm1] +     // no dangling end [.(...).-(...)] or [.(...).-(...).] 
                                                        exp_dangle5(ii,jj,jj-1) * (u1_ip_ju[lp2jjm1] + u1_iu_ju[lp2jjm1]))   // [.(...).-(...)-..]
                                              + exp_dangle5(ii,jj,jj-1) * u1_ip_ju[lp1jjm1]);    // [.(...)(...)-..]                                                                                         
                                }
                                
                                double upm_temp = 0;
                                for (h=ii+3; h < jj-TURN-4; h++)    // case (....(...)--(--)-..)
                                {
                                    int hjjm1 = index[h] + jj-1 -h;
                                    upm_temp +=  EXPC[h-ii-1] *
                                        ( //s2_jp[hjjm1] +   // no dangling end --(...)) or --(...).) 
                                            s2_ju[hjjm1]* exp_dangle5(ii,jj,jj-1));    // --(...)..)                                    
                                }    
                                temp += exp_dangle3 (ii, jj, ii+1) * upm_temp;                                
                                //printf ("i=%d, j=%d, upm=%g\n", i, j, upm[ij]);
                            
                                GlogZ[index_param] += temp * EXPA * EXPB[2] * exp_AUpenalty (ii,jj) * p[iijj] / up[iijj];
                            }
                        }
                    }                                                                   

                    
                    // the multi-loop 5' dangling end which can appear to the left of branches
                    // taken from compute_p
                    temp = 0;
                    
                    for (h=1; h < seqlen; h++)   // starts from 1 because h-1 must be valid
                    {
                        for (l=seqlen-1; l > h+TURN; l--)
                        {
                            if (can_pair (sequence[h], sequence[l]) && 
                                sequence[l] == i && sequence[h] == j && sequence[h-1] == k)    
                            {                    
                                // from now on, taken from compute_p      
                                int hl = index[h] + l - h;
                                
                                for (ii=0; ii < h; ii++)
                                {
                                    int iil = index[ii] + l - ii;
                                    // the case when h-l is the first branch of the multi-loop
                                    if (ii < h-2)  // we must add d3 (which is added in pmd3); also add d5 here, if i < h-2
                                    {
                                        temp += EXPC[h-ii-1] * exp_dangle5(l,h,h-1) *
                                            (pmd3_noneedmidd3[iil] + pmd3_needmidd3[iil]  * exp_dangle3(l,h,l+1) * EXPC[1]);
                                            // (.-[...](---)-)    i.-h...l(---)-j 
                                            // (.-[...].-(---)-)    i.-h...l.-(---)-j
                                    }
                                    
                                    if (ii < h - TURN -2)    // only now we can have a branch to the left
                                    {        
                                        int iip1hm1 = index[ii+1] + h-1 - (ii+1);
                                        int iip1hm2 = index[ii+1] + h-2 - (ii+1);
                                        int iip1hm3 = index[ii+1] + h-3 - (ii+1);
                                        int iip2hm1 = index[ii+2] + h-1 - (ii+2);
                                        int iip2hm2 = index[ii+2] + h-2 - (ii+2);
                                        int iip2hm3 = index[ii+2] + h-3 - (ii+2);
                                        int iip3hm1 = index[ii+3] + h-1 - (ii+3);
                                        int iip3hm2 = index[ii+3] + h-2 - (ii+3);
                                        int iip3hm3 = index[ii+3] + h-3 - (ii+3);
                                        int iilp1 = index[ii] + l+1 - ii;
                                        double term1 = 0.0;
                                        
                                        if (up[iilp1] != 0)
                                        {
                                            //printf ("EXPC=%g\n", exp ((misc.multi_free_base_penalty) * oneoverRT) );
                                            term1 = p[iilp1] / up[iilp1] *  exp_AUpenalty (ii,l+1) *
                                            
                                                        (  // first, the case ((..-)-[...])
                                                        ( //u1_ip_jp[iip1hm1] +  // there's no 5' dangling end here   // ((..-)[...]) or ((..-).[...])
                                                             exp_dangle5(l,h,h-1)*u1_ip_ju[iip1hm1] )   // ((..-).-[...])  
                                                                
                                                            // next, the case (.-(...)-[...])      i..-(...)h...lj
                                                        + exp_dangle3(ii,l+1,ii+1) * EXPC[1] *
                                                            ( //u1_ip_jp[iip2hm1] + u1_iu_jp[iip2hm1] +     // no dangling end
                                                                    exp_dangle5(l,h,h-1)*(u1_ip_ju[iip2hm1]+u1_iu_ju[iip2hm1]) ) 
                                                                // (..-(..-).-[...])     
                                                        );                                                                                                                                                                                                                   
                                        }
                                        
                                        temp += 
                                        
                                        // first, when h-l is the last branch to the right
                                        ( term1 + 
                                        
                                        
                                        // case ((..-)-[...].-)  
                                        pm1nod3_needendd3[iil] * exp_dangle3(l,h,l+1) *   // EXPC added in pm1nod3
                                            (exp_dangle5(l,h,h-1)*u1_ip_ju[iip1hm1])
                                                    // c is added in pm1nod3_needendd3
                                                    // ((...).-[...].-)         i(...)h...l.-j
                                        // case (.(..-)-[...].-)
                                        + pm1d3_needendd3[iil] * exp_dangle3(l,h,l+1) *     // don't need EXPC[1]
                                            (exp_dangle5(l,h,h-1)*(u1_ip_ju[iip2hm1]+u1_iu_ju[iip2hm1]))
                                                    // (.(...).-[...].-)        i.(...)h...l.-j                                                    
                                        // case (..-(..-)-[...].-)            
                                                    // (..-(...)..-[...].-)      i..-(...)h...l.-j
                                                    
                                        // we have branches to the left and to the right of h-l
                                        //+ pmnod3_needmidd3[iil] * exp_dangle3(l,h,l+1) * u1[iip1hm1]);
                                                    // ((...).-[...]-(---)...)    i(...).-h...l-(---)...j 
                                                                                                        
                                        + (pmnod3_needmidd3[iil] * exp_dangle3(l,h,l+1) * EXPC[1]+ pmnod3_noneedmidd3[iil]) *
                                            (exp_dangle5(l,h,h-1)*u1_ip_ju[iip1hm1])
                                                    // c is added in pm1nod3_needendd3
                                                    // ((...).-[...].-(--)-)         i(...)h...l.-j
                                        // case (.(..-)-[...].-(--)-) or (.(..-)-[...](--)-)
                                        + EXPC[1]*(pmd3_needmidd3[iil] * exp_dangle3(l,h,l+1) * EXPC[1] + pmd3_noneedmidd3[iil]) *
                                            (exp_dangle5(l,h,h-1)*(u1_ip_ju[iip2hm1]+u1_iu_ju[iip2hm1])));
                                                    // (.(...).-[...].-)        i.(...)h...l.-j                                                    
                                        // case (..-(..-)-[...].-(--)-) or (..-(..-)-[...](--)-)
                                                    // (..-(...)..-[...].-)      i..-(...)h...l.-j
                                                    
                                    }        
                                }
                                GlogZ[index_param] += up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) * temp;
                            }
                        }
                    }             
                           
                    index_param++;
                }
            }
    index_should_be = structure_type_index("internal_penalty_by_size[4]");            
    if (index_param != index_should_be)
    {
        printf ("Index param after dangle_bot = %d, should be %d\n", index_param, index_should_be);
        exit(1);
    }
          

    }       // end if (fix_dangles)
    
    // compute derivatives for AU_penalties
    // the contribution from hairpin loop of size 3 was computed above
    // now consider the constributions from the exterior loop

    double term1 = 0.0;
    double term2 = 0.0;
        
    
    for (ii = 0; ii < seqlen; ii++)
    {
        for (jj = ii+TURN+1; jj < seqlen; jj++)
        {
            if (can_pair (sequence[ii], sequence[jj]) && has_AU_penalty (sequence[ii], sequence[jj]))
            {   
                iijj = index[ii]+jj-ii;
                AUpen = AU_penalty (sequence[ii], sequence[jj]);
    
                
                // first, exterior loop - same as in compute_p            
                
                term1 = 0.0;
                term2 = 0.0;
                    
                IFD
                {
                    if (ii > 0)     term1 = u[ii-1];
                    else            term1 = 1.0;
                }
                else
                {
                    if (ii > TURN)        
                            // case i-1 paired
                        term1 += (u_ip_jp[ii-1] + u_iu_jp[ii-1] + exp_dangle5(jj,ii,ii-1)*(u_ip_ju[ii-1] + u_iu_ju[ii-1]));
                    else if (ii > 0)      term1 += exp_dangle5(jj,ii,ii-1);
                    else                  term1 = 1;  
                }
                                
                IFD
                {
                    if (jj < seqlen-1)    term2 = u[index[jj+1]+seqlen-1-(jj+1)];
                    else                  term2 = 1.0;
                }
                else
                {
                    int jjp1n = index[jj+1] + seqlen-1 - jj-1;
                    int jjp2n = index[jj+2] + seqlen-1 - jj-2;
                    if (jj < seqlen-3)    
                        term2 += u_ip_jp[jjp1n] + u_ip_ju[jjp1n] + 
                                exp_dangle3(jj,ii,jj+1)*(u_ip_jp[jjp2n] + u_ip_ju[jjp2n] + u_iu_jp[jjp2n] + u_iu_ju[jjp2n]);
                    else if (jj < seqlen-1) term2 += exp_dangle3(jj,ii,jj+1);
                    else                  term2 = 1;
                }

                IFD
                {
                    GlogZ[AUpen_index] += term1 *up[iijj] * term2 * exp_AUpenalty(ii,jj) / u[seqlen-1];                        
                }
                else
                {
                    GlogZ[AUpen_index] += term1 *up[iijj] * term2 * exp_AUpenalty(ii,jj) /
                        (u_ip_jp[seqlen-1] + u_ip_ju[seqlen-1] + u_iu_jp[seqlen-1] + u_iu_ju[seqlen-1]);
                }
                // if it closes a multi loop, or closes a helix of a multi-loop - it's considered at helix below
            }                                        
        }
    }               
    
    int sizeindex = structure_type_index ("misc.multi_offset");
    for (ii = 0; ii < seqlen; ii++)
    {
        for (jj = ii+TURN+1; jj < seqlen; jj++)
        {
            iijj = index[ii]+jj-ii;
            if (upm[iijj] > 0)
            {   
                GlogZ[sizeindex] += upm[iijj] * p[iijj] / up[iijj];
            }
        }
    }
    
    sizeindex = structure_type_index ("misc.multi_helix_penalty");
    // took the multi-loop contribution of p, and upm, like above
    // taken from compute_upm and compute_p
    double temp;
    for (h = 0; h < seqlen; h++)
    {
        for (l = h+TURN+1; l < seqlen; l++)
        {
            if (!can_pair (sequence[h], sequence[l]))
                continue;
            int hl = index[h]+l-h;
            // the exterior multi-loop pair      
            temp = upm[hl] * p[hl] / up[hl];      
            GlogZ[sizeindex] += temp;
            if (has_AU_penalty (sequence[h],sequence[l]))    GlogZ[AUpen_index] += temp;
            
            for (i=0; i < h; i++)
            {
                int il = index[i] + l - i;  
                IFD
                {
                    temp = up[hl] * exp_AUpenalty (h,l) * EXPA * EXPB[2] * EXPC[h-i-1] * pm[il];
                    GlogZ[sizeindex] += temp;
                    if (has_AU_penalty (sequence[h],sequence[l]))    GlogZ[AUpen_index] += temp;
                }
                else
                {                          
                    if (i < h-1)  // we must add d3 (which is added in pmd3); also add d5 here, if i < h-2
                    {
                        temp = up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                            EXPC[h-i-1] * pmd3_noneedmidd3[il] * (i<h-2?exp_dangle5(l,h,h-1):1);
                        GlogZ[sizeindex] += temp; 
                        if (has_AU_penalty (sequence[h],sequence[l]))    GlogZ[AUpen_index] += temp;
                            // (.-[...](---)-)    i.-h...l(---)-j 
                                
                        temp = up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                            EXPC[h-i-1] * pmd3_needmidd3[il]  * exp_dangle3(l,h,l+1) * EXPC[1] *
                                    (i<h-2?exp_dangle5(l,h,h-1):1);    
                        GlogZ[sizeindex] += temp;
                        if (has_AU_penalty (sequence[h],sequence[l]))    GlogZ[AUpen_index] += temp;
                            // (.-[...].-(---)-)    i.-h...l.-(---)-j  
                    }
                    else    // case ((... , no dangling end
                    {
                        temp = up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                            ( pmnod3_noneedmidd3[il]);  
                        GlogZ[sizeindex] += temp;            
                        if (has_AU_penalty (sequence[h],sequence[l]))    GlogZ[AUpen_index] += temp;
                            // ([...](---)-)    ih...l(---)-j 
                            
                        temp = up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                            ( pmnod3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1]);
                        GlogZ[sizeindex] += temp;
                        if (has_AU_penalty (sequence[h],sequence[l]))    GlogZ[AUpen_index] += temp;
                            // ([...].-(---)-)    ih...l.-(---)-j                 
                    }
                }                    
                if (i < h - TURN -2)    // only now we can have a branch to the left
                {        
                    int ip1hm1 = index[i+1] + h-1 - (i+1);                
                    
                    IFD
                    {
                        temp = up[hl] * exp_AUpenalty (h,l) *  EXPA * EXPB[2] *
                                            u1[ip1hm1] * (pm1[il] + pm[il]);
                        GlogZ[sizeindex] += temp;
                        if (has_AU_penalty (sequence[h],sequence[l]))    GlogZ[AUpen_index] += temp;
                    }
                    else
                    {
                        int ip1hm2 = index[i+1] + h-2 - (i+1);
                        int ip1hm3 = index[i+1] + h-3 - (i+1);
                        int ip2hm1 = index[i+2] + h-1 - (i+2);
                        int ip2hm2 = index[i+2] + h-2 - (i+2);
                        int ip2hm3 = index[i+2] + h-3 - (i+2);
                        int ip3hm1 = index[i+3] + h-1 - (i+3);
                        int ip3hm2 = index[i+3] + h-2 - (i+3);
                        int ip3hm3 = index[i+3] + h-3 - (i+3);
                        int ilp1 = index[i] + l+1 - i;
                        double term1 = 0.0;
                        if (up[ilp1] != 0)
                        {
                            //printf ("EXPC=%g\n", exp ((misc.multi_free_base_penalty) * oneoverRT) );
                            term1 = p[ilp1] / up[ilp1] *  exp_AUpenalty (i,l+1) *
                                        (  // first, the case ((..-)-[...])
                                          ( u1_ip_jp[ip1hm1]    // ((..-)[...]) or ((..-).[...])
                                            + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1] )   // ((..-).-[...])  
                                                
                                            // next, the case (.-(...)-[...])      i..-(...)h...lj
                                          + exp_dangle3(i,l+1,i+1) * EXPC[1] *
                                            ( u1_ip_jp[ip2hm1] + u1_iu_jp[ip2hm1] 
                                                    + exp_dangle5(l,h,h-1)*(u1_ip_ju[ip2hm1]+u1_iu_ju[ip2hm1]) ) 
                                                // (..-(..-).-[...])     
                                        ); 
                        }
                        
                        temp = up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                        
                        // first, when h-l is the last branch to the right
                        ( term1 + 
                        
                        // case ((..-)-[...].-)  
                        pm1nod3_needendd3[il] * exp_dangle3(l,h,l+1) * // we don't need EXPC[1] here, because it's in pm1nod3_needendd3
                            (u1_ip_jp[ip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1])
                                    // c is added in pm1nod3_needendd3
                                    // ((...).-[...].-)         i(...)h...l.-j
                        // case (.(..-)-[...].-)
                        + pm1d3_needendd3[il] * exp_dangle3(l,h,l+1) * // we don't need EXPC[1] here, because it's in pm1d3_needendd3
                            (u1_ip_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip2hm1])
                                    // (.(...).-[...].-)        i.(...)h...l.-j
                                    
                        // case (..-(..-)-[...].-)            
                        + pm1d3_needendd3[il]  * exp_dangle3(l,h,l+1) * // don't need EXPC[1]
                            (u1_iu_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_iu_ju[ip2hm1])
                                    // (..-(...)..-[...].-)      i..-(...)h...l.-j
                                    
                        // we have branches to the left and to the right of h-l
                        //+ pmnod3_noneedmidd3[il]  * u1[ip1hm1] +    
                                    // ((...)[...]-(---)...)    i(...)h...l-(---)...j 
                        //pmnod3_needmidd3[il] * exp_dangle3(l,h,l+1) * u1[ip1hm1]);        
                                    // ((...).-[...]-(---)...)    i(...).-h...l-(---)...j 
                          
                        + (pmnod3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1]+ pmnod3_noneedmidd3[il]) *
                            (u1_ip_jp[ip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1])
                                    // c is added in pm1nod3_needendd3
                                    // ((...).-[...].-(--)-)         i(...)h...l.-j
                        // case (.(..-)-[...].-(--)-) or (.(..-)-[...](--)-)
                        + EXPC[1]*(pmd3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1] + pmd3_noneedmidd3[il]) *
                                (u1_ip_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip2hm1])
                                    
                        // case (..-(..-)-[...].-(--)-) or (..-(..-)-[...](--)-)
                        + EXPC[1]*(pmd3_needmidd3[il]  * exp_dangle3(l,h,l+1) * EXPC[1] + pmd3_noneedmidd3[il])*
                            (u1_iu_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_iu_ju[ip2hm1]) ); 
                                    
                        GlogZ[sizeindex] += temp;
                        if (has_AU_penalty (sequence[h],sequence[l]))    GlogZ[AUpen_index] += temp;
                    }                
                }        
            }
        }
    }
    
                
    sizeindex = structure_type_index ("misc.multi_free_base_penalty");
    // traverse the sequence from left to right         
    int ij, hl;
    int ip1l, ip2l, lp2jm1, lp1jm1, lp2jm2, lp1jm2, lp2jm3, lp1jm3, hj;
    double upm_temp;

    // the first free bases
    for (h=0; h < seqlen; h++)   
    {
        for (l=seqlen-1; l > h+TURN; l--)
        {    
            if (!can_pair(sequence[h], sequence[l]))
                continue;
            hl = index[h]+l-h;        
                
            for (i=0; i < h-1; i++)
            {
                int il = index[i] + l - i;                
                IFD
                {
                    GlogZ[sizeindex] += (h-i-1)* up[hl] * exp_AUpenalty (h,l) * EXPA * EXPB[2] * EXPC[h-i-1] * pm[il];
                }
                else
                {
                    // the case when h-l is the first branch of the multi-loop
                    GlogZ[sizeindex] += (h-i-1)*up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                        EXPC[h-i-1] * pmd3_noneedmidd3[il] * (i<h-2?exp_dangle5(l,h,h-1):1);
                        // (.-[...](---)-)    i.-h...l(---)-j
                    GlogZ[sizeindex]  += (h-i-1)*up[hl] * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                        EXPC[h-i-1] * pmd3_needmidd3[il]  * exp_dangle3(l,h,l+1) * EXPC[1] *
                                (i<h-2?exp_dangle5(l,h,h-1):1);
                        // (.-[...].-(---)-)    i.-h...l.-(---)-j
                }
            }
        }
    }
    
    //printf ("GlogZ[%d] = %g\n", sizeindex, GlogZ[sizeindex]);
       
    // second, there is at least a branch to the left of k, and a branch to the right
    // Now, k is the free base inside a multi-loop, closed be i,j     
    for (i=0; i < seqlen-2*TURN-5; i++)
    {
        for (j=i+2*TURN+5; j < seqlen; j++)
        {
            ij = index[i]+j-i;
            if (upm[ij] > 0)
            {
                for (k=i+TURN+3; k < j-TURN-2; k++)
                {
                    int ip1km1 = index[i+1] + k-1 - (i+1);
                    int kp1jm1 = index[k+1] + j-1 - (k+1);
                    
                    IFD
                    {
                        GlogZ[sizeindex] += p[ij] / up[ij] * EXPA * EXPB[1] * EXPC[1] * exp_AUpenalty(i,j)*
                                            u1[ip1km1] * u1[kp1jm1];
                    }
                    else
                    {
                        int ip1k = index[i+1] + k - (i+1);                    
                        int ip1km2 = index[i+1] + k-2 - (i+1);
                        int ip2k = index[i+2] + k - (i+2);
                        int ip2km1 = index[i+2] + k-1 - (i+2);
                        int ip2km2 = index[i+2] + k-2 - (i+2);
                        int kjm1 = index[k] + j-1 - (k);
                        int kjm2 = index[k] + j-2 - (k);
                        int kjm3 = index[k] + j-3 - (k);                    
                        int kp1jm2 = index[k+1] + j-2 - (k+1);
                        int kp1jm3 = index[k+1] + j-3 - (k+1);
                        // first situation (---)-.k---
                                          
                        GlogZ[sizeindex] += p[ij] / up[ij] * EXPA * EXPB[1] * exp_AUpenalty(i,j)*
                            
                            // the branch(es) to the left of k
                            
                            ( exp_dangle3 (i, j, i+1) * EXPC[1]* (u1_ip_ju[ip2k] + u1_iu_ju[ip2k])    // (.-(---)-.k(
                              + u1_ip_ju[ip1k] ) *   // ((---).k(
                            // the branch(es) to the right of k
                            ( u1_iu_jp[kjm1] / EXPC[1] + u1_iu_ju[kjm1]/EXPC[1] *exp_dangle5(i,j,j-1));
                                           
                        // second situation (---)k---
                        GlogZ[sizeindex] += p[ij] / up[ij] * EXPA * EXPB[1] * exp_AUpenalty(i,j) *
                            // the branch(es) to the left of k
                            (exp_dangle3 (i, j, i+1) * EXPC[1]*(u1_ip_ju_jm1p[ip2k] + u1_iu_ju_jm1p[ip2k]) + u1_ip_ju_jm1p[ip1k]) *    // only j should be paired, not j-1
                            // the branch(es) to the right of k
                            (u1_ip_jp[kp1jm1] + u1_iu_jp[kp1jm1] + (u1_ip_ju[kp1jm1] + u1_iu_ju[kp1jm1])*exp_dangle5(i,j,j-1));                                                        
                    }
                }
            }
        }
    }



    // last, no branch to the right of k
    // similar to the first free bases 
    for (h=0; h < seqlen; h++)   
    {
        int hj;
        for (l=seqlen-3; l > h+TURN; l--)
        {    
            if (!can_pair(sequence[h], sequence[l]))
                continue;
            hl = index[h]+l-h;        
                
            j = l+2;
            hj = index[h] + j - h;
            
            IFD
            {
                //GlogZ[sizeindex] += (j-l-1)*up[hl] * exp_AUpenalty (h,l) * EXPA * EXPB[2] * EXPC[j-l-1] * pm2[hj];
                for (j=l+2; j < seqlen; j++)
                {
                    hj = index[h] + j - h;
                    GlogZ[sizeindex] += (j-l-1)*up[hl] * exp_AUpenalty (h,l) * EXPA * EXPB[2] * EXPC[j-l-1] * pm2[hj];
                }
            }
            else
            {
                GlogZ[sizeindex] += (j-l-1)*up[hl] * exp_dangle3(l,h,l+1) * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                    EXPC[j-l-1] * pm2nod5_noneedmidd5[hj];    // (-(---)[...].) or (-(---).[...].)
                    
                GlogZ[sizeindex] += (j-l-1)*up[hl] * exp_dangle3(l,h,l+1) * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                    EXPC[j-l-1] * pm2nod5_needmidd5[hj] * exp_dangle5(l,h,h-1);    // (-(---)-..[...].)
                
                for (j=l+3; j < seqlen; j++)
                {
                    // needs the last 5' dangling end            
                    hj = index[h] + j - h;
                    GlogZ[sizeindex] += (j-l-1)*up[hl] * exp_dangle3(l,h,l+1) * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                        EXPC[j-l-1] * pm2d5_noneedmidd5[hj];    // (-(---)[...].-.) or (-(---).[...].-.)
                    
                    GlogZ[sizeindex] += (j-l-1)*up[hl] * exp_dangle3(l,h,l+1) * EXPA * EXPB[2] * exp_AUpenalty (h,l) *
                        EXPC[j-l-1] * pm2d5_needmidd5[hj] * exp_dangle5(l,h,h-1);//*EXPC[1];    // (-(---)-..[...].-.)
                }
            }
        }
    }    
      
    
    /*
    // make sure all gradients are <=1            // actually it can be > 1, but not by tooo much
    for (i=0; i < num_params; i++)
    {
        if (GlogZ[i] > 1)
        {
            printf ("ERROR: Gradient GlogZ[%d] too large!! It is %g\n", i, GlogZ[i]);
            //exit(1);
        }
    }
    */
}

                
int s_partition_function::correct_gradient ()
// return 1 if it is correct, 0 otherwise
{
    int index_param;
    int i, j, k, l, m, n, o, p;
    index_param = 0;
    int correct = 1;
 
    double diff; 
    for (i=0; i < num_params; i++)
    {
        diff = GlogZexhaustive[i] - GlogZ[i];
        if (diff < 0)     diff *= -1.0;        
        if ((GlogZexhaustive[i] != 0 || GlogZ[i] != 0 ) && diff > PREC)
        {
            printf ("%d %s\texhaust=%g\tsmart=%g\tdiff=%g\n", i, string_params[i], GlogZexhaustive[i], GlogZ[i], diff);
            correct = 0;
        }
        else if (GlogZexhaustive[i] != 0) 
            printf ("%d %s\texhaust=%g\tsmart=%g\t+++++++++\n", i, string_params[i], GlogZexhaustive[i], GlogZ[i]);
    }
    return correct;
}

void s_partition_function::copy_gradient (double *grad)
// write the gradient values into grad
{
    int i;
    for (i=0; i < num_params; i++)
    {
        grad[i] = GlogZ[i];
    }
}

int s_partition_function::correct_gradient_nan ()
// return 1 if no derivative is nan
{
    int i;
    int correct = 1;
    for (i=0; i < num_params; i++)
    {
        if (isnan(GlogZ[i]))
        {
            printf ("Glog[%d]=%g\n", i, GlogZ[i]);
            correct = 0;
        }
    }
    return correct;
}



void s_partition_function::print_gradient ()
// print the gradient
{
    int i;
    for (i=0; i < num_params; i++)
    {
        printf ("Glog[%d]=%g\n", i, GlogZ[i]);
    }
}

/*
// compute_upm which has both IFD and else, just in case I messed up something when I split them

void s_partition_function::compute_upm (int i, int j)
// i and j close a multi-loop
{
    int ij = index[i]+j-i;
    int l, h;
    int ip1l, ip2l, lp2jm1, lp1jm1, lp2jm2, lp1jm2, lp2jm3, lp1jm3, hj;
    double upm_temp;
    
    if (!can_pair(sequence[i], sequence[j]))
        return;
    upm[ij] = 0;
    // took this out for speed
    double upm_common = exp_AUpenalty (i,j) * EXPA * EXPB[2];
    // there must be at least one more branch (i.e. 5 nucleotides) from l+1 to j-1: l+1+TURN <= j-1
    for (l=i+2; l < j-TURN-2; l++)    // case ((...)--(--)-)
    {
        ip1l = index[i+1]+l-i-1;
        lp2jm1 = index[l+2]+j-1-l-2;
        lp1jm1 = index[l+1]+j-1-l-1;

        IFD
            upm[ij] += up[ip1l] * exp_AUpenalty (i+1,l) * u1[lp1jm1];
        else
            upm[ij] += up[ip1l] * exp_AUpenalty (i+1,l) *
                    (  u1_ip_jp[lp1jm1]      // [(...)(.-..)] or [(...)(.-..).]                    
                       + exp_dangle3 (l, i+1, l+1) * EXPC[1] *
                            ( u1_ip_jp[lp2jm1] + u1_iu_jp[lp2jm1] +     // [(...).-(...)] or [(...).-(...).] 
                                exp_dangle5(i,j,j-1) * (u1_ip_ju[lp2jm1] + u1_iu_ju[lp2jm1])) +    // [(...).-(...)-..]
                       + exp_dangle5(i,j,j-1) * u1_ip_ju[lp1jm1]);    // [(...)(...)-..]
    }
    for (l=i+3; l < j-TURN-2; l++)    // case [.(...)--(--)-]
    {
        ip2l = index[i+2]+l-i-2;
        lp2jm1 = index[l+2]+j-1-l-2;
        lp1jm1 = index[l+1]+j-1-l-1;
        lp2jm2 = index[l+2]+j-2-l-2;
        lp1jm2 = index[l+1]+j-2-l-1;        
        lp2jm3 = index[l+2]+j-3-l-2;
        lp1jm3 = index[l+1]+j-3-l-1;      
        IFD
            upm[ij] += up[ip2l] * EXPC[1] * exp_AUpenalty (i+2, l) * u1[lp1jm1];
        else
            upm[ij] += up[ip2l] * EXPC[1] *
                    exp_dangle3 (i, j, i+1) * exp_AUpenalty (i+2, l) * 
                    (  u1_ip_jp[lp1jm1]      // [.(...)(.-..)] or [.(...)(.-..).]                    
                       + exp_dangle3 (l, i+2, l+1) * EXPC[1] *
                            ( u1_ip_jp[lp2jm1] + u1_iu_jp[lp2jm1] +     // [.(...).-(...)] or [.(...).-(...).] 
                                exp_dangle5(i,j,j-1) * (u1_ip_ju[lp2jm1] + u1_iu_ju[lp2jm1])) +    // [.(...).-(...)-..]
                       + exp_dangle5(i,j,j-1) * u1_ip_ju[lp1jm1]);    // [.(...)(...)-..] 
    }
    upm_temp = 0;
    for (h=i+3; h < j-TURN-2; h++)    // case (....(...)--(--)-)
    {
        hj = index[h]+j-h;
        int hjm1 = index[h] + j-1 -h;
        int hjm2 = index[h] + j-2 -h;
        int hjm3 = index[h] + j-3 -h;
        IFD
            upm_temp += EXPC[h-i-1] * s2[hjm1];
        else
            upm_temp += EXPC[h-i-1] *
                ( s2_jp[hjm1]    // --(...)) or --(...).) 
                  + s2_ju[hjm1]* exp_dangle5(i,j,j-1));    // --(...)..)
    }   
    IFD
        upm[ij] += upm_temp;
    else
    {
        upm[ij] += exp_dangle3 (i, j, i+1) * upm_temp;        
    }
    upm[ij] *= upm_common;
    //printf ("i=%d, j=%d, upm=%g\n", i, j, upm[ij]);
}
*/


