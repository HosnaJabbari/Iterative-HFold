/***************************************************************************
                          s_min_folding.cpp  -  description
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

// This class is the main class to compute the MFE prediction

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "constants.h"
#include "structs.h"
#include "externs.h"
#include "common.h"
#include "simfold.h"
#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"
#include "s_multi_loop.h"
#include "s_energy_matrix.h"
#include "s_specific_functions.h"
#include "s_min_folding.h"


s_min_folding::s_min_folding (char *sequence)
// constructor for the unrestricted mfe case
{
    //check_sequence (sequence);
    this->sequence = sequence;
    this->known_structure = NULL;
    allocate_space();
}

s_min_folding::s_min_folding (char *sequence, char *restricted)
// constructor for the restricted mfe case
{
    check_sequence (sequence);
    this->sequence = sequence;
    this->restricted = restricted;
    allocate_space();
}

void s_min_folding::allocate_space()
// allocate the necessary memory
{
    int i;
    nb_nucleotides = strlen(sequence);
    
    f = new minimum_fold [nb_nucleotides];
    if (f == NULL) giveup ("Cannot allocate memory", "energy");
    W = new PARAMTYPE [nb_nucleotides];
    if (W == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) W[i] = 0;

    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);

    H = new s_hairpin_loop (sequence, int_sequence, nb_nucleotides);
    if (H == NULL) giveup ("Cannot allocate memory", "energy");
    S = new s_stacked_pair (int_sequence, nb_nucleotides);
    if (S == NULL) giveup ("Cannot allocate memory", "energy");
    VBI = new s_internal_loop (int_sequence, nb_nucleotides);
    if (VBI == NULL) giveup ("Cannot allocate memory", "energy");
    VM = new s_multi_loop (int_sequence, nb_nucleotides);
    if (VM == NULL) giveup ("Cannot allocate memory", "energy");
    V = new s_energy_matrix (int_sequence, nb_nucleotides);
    if (V == NULL) giveup ("Cannot allocate memory", "energy");

    // initialize the structure

    structure = new char[nb_nucleotides+2];
    for (i=0; i<nb_nucleotides; i++)  structure[i] = '.';
    structure[nb_nucleotides] = '\0';

    S->set_energy_matrix (V);
    VBI->set_energy_matrix (V);
    VM->set_energy_matrix (V);
    V->set_loops (H, S, VBI, VM, NULL);
}


s_min_folding::~s_min_folding()
// release space
{
    delete V;
    delete VM;
    delete VBI;
    delete S;
    delete H;
    delete [] f;
    delete [] W;
    delete [] int_sequence;
    delete [] structure;
}


double s_min_folding::s_simfold ()
// PRE:  the init_data function has been called;
//       the space for structure has been allocated
// POST: fold sequence, return the MFE structure in structure, and return the MFE
{
    double energy;
    energy = fold_sequence ();
    return energy;
}


double s_min_folding::s_simfold_restricted ()
// PRE:  the init_data function has been called;
//       the space for structure has been allocate
// POST: fold sequence, return the MFE structure in structure, and return the MFE
{
    double energy;
    energy = fold_sequence_restricted ();
    return energy;    
}


double s_min_folding::fold_sequence ()
// helper function, it folds sequence, returns structure and free energy
{
    double energy;
    int i, j;

    for (j=0; j < nb_nucleotides; j++)
    {
        // if (constraints[j]) continue;
        for (i=0; i<j; i++)
        {
            // if (constraints[i]) continue;
            V->compute_energy (i,j);
        }
        // if I put this before V calculation, WM(i,j) cannot be calculated, because it returns infinity
        VM->compute_energy_WM (j);
    }
    for (j=1; j < nb_nucleotides; j++)
    {
        compute_W (j);
    }
    energy = W[nb_nucleotides-1]/100.0;        
    
    if (debug)
    {
        for (j=1; j < nb_nucleotides; j++)
        {
            printf ("W(%d) = %d\n", j, W[j]);
        }
    }    
    
    
    // backtrack
    // first add (0,n-1) on the stack
    stack_interval = new seq_interval;
    stack_interval->i = 0;
    stack_interval->j = nb_nucleotides - 1;
    stack_interval->energy = W[nb_nucleotides-1]; 
    stack_interval->type = FREE;
    stack_interval->next = NULL;

    seq_interval *cur_interval = stack_interval;

    while ( cur_interval != NULL)
    {         
        stack_interval = stack_interval->next;
        backtrack (cur_interval);
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }
    if (debug)
    {
        print_result ();
    }    
    //delete stack_interval;
    return energy;
}


double s_min_folding::fold_sequence_restricted ()
// helper function, it folds sequence, returns structure and free energy
{
    double energy;
    int i, j;

    str_features *fres;
    if ((fres = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");   
    // detect the structure features  
    detect_structure_features (restricted, fres);
    
    /*
    for (i=0; i < nb_nucleotides; i++)
        if (fres[i].pair != -1)
            printf ("%d pairs %d, type %c\n", i, fres[i].pair, fres[i].type);
    */
    
    for (j=0; j < nb_nucleotides; j++)
    {
        for (i=0; i<j; i++)
        {            
            // V(i,j) = infinity if i restricted or j restricted and pair of i is not j
            if ((fres[i].pair > -1 && fres[i].pair !=j) || (fres[j].pair > -1 && fres[j].pair != i)) 
                continue;              
            if (fres[i].pair == -1 || fres[j].pair == -1)   // i or j MUST be unpaired
                continue;
            V->compute_energy_restricted (i, j, fres);
        }
        // if I put this before V calculation, WM(i,j) cannot be calculated, because it returns infinity
        VM->compute_energy_WM_restricted (j, fres);
    }
    for (j=1; j < nb_nucleotides; j++)
    {
        compute_W_restricted (j, fres);
    }
    energy = W[nb_nucleotides-1]/100.0;        
    
    if (debug)
    {
        for (j=1; j < nb_nucleotides; j++)
        {
            printf ("W(%d) = %d\n", j, W[j]);
        }
    }    
    
    
    // backtrack
    // first add (0,n-1) on the stack
    stack_interval = new seq_interval;
    stack_interval->i = 0;
    stack_interval->j = nb_nucleotides - 1;
    stack_interval->energy = W[nb_nucleotides-1]; 
    stack_interval->type = FREE;
    stack_interval->next = NULL;

    seq_interval *cur_interval = stack_interval;

    while ( cur_interval != NULL)
    {         
        stack_interval = stack_interval->next;
        backtrack_restricted (cur_interval, fres);
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }
    if (debug)
    {
        print_result ();
    }    
    delete [] fres;
    //delete stack_interval;
    return energy;
}



void s_min_folding::insert_node (int i, int j, char type)
  // insert at the beginning
{      
    seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
}



void s_min_folding::backtrack (seq_interval *cur_interval)
// PRE:  All matrixes V, VM, WM and W have been filled
// POST: Discover the MFE path   
{
    char type;

    if(cur_interval->type == LOOP)
    {
        int i = cur_interval->i;
        int j = cur_interval->j;
        f[i].pair = j;
        f[j].pair = i;
        structure[i] = '(';
        structure[j] = ')';
      
        type = V->get_type (i,j);
        if (debug) 
            printf ("\t(%d,%d) LOOP - type %c\n", i,j,type);
        if (type == STACK)
        {
            f[i].type = STACK;
            f[j].type = STACK;
            if (i+1 < j-1)
                insert_node (i+1, j-1, LOOP);
            else
            {
                printf ("NOT GOOD STACK, i=%d, j=%d\n", i, j);
                exit (0);
            }                    
        }
        else if (type == HAIRP)
        {
            f[i].type = HAIRP;
            f[j].type = HAIRP;
        }
        else if (type == INTER) 
        {
            f[i].type = INTER;
            f[j].type = INTER;
            // detect the other closing pair
            int ip, jp, best_ip, best_jp, minq;
            PARAMTYPE tmp, min = INF;
            for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1); ip++) 
            {
                minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    
                for (jp = minq; jp < j; jp++)
                {
                    tmp = VBI->get_energy_str (i,j,ip,jp);
                    if (tmp < min)
                    {
                        min = tmp;
                        best_ip = ip;
                        best_jp = jp;
                    }
                }
            }
            if (best_ip < best_jp)
                insert_node (best_ip, best_jp, LOOP);
            else
            {
                printf ("NOT GOOD INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
                exit (0);
            }                    
        }
        else if (type == MULTI)
        {
            f[i].type = MULTI;
            f[j].type = MULTI;
            int k, best_k, best_row;
            PARAMTYPE tmp, min = INF;
            for (k = i+TURN+1; k <= j-TURN-2; k++)
            {
                tmp = VM->get_energy_WM (i+1,k) + VM->get_energy_WM (k+1, j-1);
                if (tmp < min)
                {
                    min = tmp;
                    best_k = k;
                    best_row = 1;
                }
                tmp = VM->get_energy_WM (i+2,k) + VM->get_energy_WM (k+1, j-1) + 
                        dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] + 
                        misc.multi_free_base_penalty;
                // add the loss
                if (pred_pairings != NULL)
                {
                    pred_pairings[i+1] = -1;
                    tmp = tmp - loss (i+1,i+1);
                }                        
                if (tmp < min)
                {
                    min = tmp;
                    best_k = k;
                    best_row = 2;
                }
                tmp = VM->get_energy_WM (i+1,k) + VM->get_energy_WM (k+1, j-2) + 
                        dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] + 
                        misc.multi_free_base_penalty;
                // add the loss
                if (pred_pairings != NULL)
                {
                    pred_pairings[j-1] = -1;
                    tmp = tmp - loss (j-1,j-1);
                }                                      
                if (tmp < min)
                {
                    min = tmp;
                    best_k = k;
                    best_row = 3;
                }
                tmp = VM->get_energy_WM (i+2,k) + VM->get_energy_WM (k+1, j-2) + 
                        dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] + 
                        dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] + 
                        2*misc.multi_free_base_penalty;
                // add the loss
                if (pred_pairings != NULL)
                {
                    pred_pairings[i+1] = -1;
                    pred_pairings[j-1] = -1;
                    tmp = tmp - loss (i+1,i+1) - loss (j-1,j-1);
                }                                                    
                if (tmp < min)
                {
                    min = tmp;
                    best_k = k;
                    best_row = 4;
                }
            }
            switch (best_row)
              {
              case 1:                
                  insert_node (i+1, best_k, M_WM);
                  insert_node (best_k+1, j-1, M_WM); break;
              case 2: insert_node (i+2, best_k, M_WM);
                insert_node (best_k+1, j-1, M_WM); break;
              case 3: insert_node (i+1, best_k, M_WM);
                insert_node (best_k+1, j-2, M_WM); break;
              case 4: insert_node (i+2, best_k, M_WM);
                insert_node (best_k+1, j-2, M_WM); break;
              }
        }
    }
    else if(cur_interval->type == FREE)
    {
        int j = cur_interval->j;
    
        if (j <= TURN) return; 
        
        PARAMTYPE min = INF, tmp, acc, energy_ij;
        int best_row, i, best_i;
    
        if (debug)
            printf ("\t(0,%d) FREE\n", j);
        tmp = W[j-1];      
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[j] = -1;
            tmp = tmp - loss (j,j);
        }                            
      
        if (tmp < min)
        {
            min = tmp;          
            best_row = 0;
        }
        for (i=0; i<=j-TURN-1; i++)
        {
            acc = (i-1>0) ? W[i-1] : 0;          
            energy_ij = V->get_energy(i,j);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
                if (tmp < min)
                {
                    min = tmp;
                    best_i = i;
                    best_row = 1;
                }
            }        
            energy_ij = V->get_energy(i+1,j);            
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
                tmp += dangle_bot [int_sequence[j]]
                    [int_sequence[i+1]]
                    [int_sequence[i]];
                // add the loss
                if (pred_pairings != NULL)
                {
                    pred_pairings[i] = -1;
                    tmp = tmp - loss (i,i);
                }                                
                if (tmp < min)
                {
                    min = tmp;
                    best_i = i;
                    best_row = 2;
                }
            }          
            energy_ij = V->get_energy(i,j-1);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j-1]) + acc;
                tmp += dangle_top [int_sequence[j-1]]
                    [int_sequence[i]]
                    [int_sequence[j]];
                // add the loss
                if (pred_pairings != NULL)
                {
                    pred_pairings[j] = -1;
                    tmp = tmp - loss (j,j);
                }                                
                if (tmp < min)
                {
                    min = tmp;
                    best_i = i;
                    best_row = 3;
                }
            }                  
            energy_ij = V->get_energy(i+1,j-1);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
                tmp += dangle_bot [int_sequence[j-1]]
                    [int_sequence[i+1]]
                    [int_sequence[i]];
                tmp += dangle_top [int_sequence[j-1]]
                    [int_sequence[i+1]]
                    [int_sequence[j]];
                // add the loss
                if (pred_pairings != NULL)
                {
                    pred_pairings[i] = -1;
                    pred_pairings[j] = -1;
                    tmp = tmp - loss (i,i) - loss (j,j);
                }            
                if (tmp < min)
                {
                    min = tmp;
                    best_i = i;
                    best_row = 4;
                }
            }
        }
        switch (best_row)
            {
            case 0: insert_node (0, j-1, FREE); break;
            case 1: insert_node (best_i, j, LOOP); 
            if (best_i-1 > TURN) 
                insert_node (0, best_i-1, FREE); 
            break;
            case 2: insert_node (best_i+1, j, LOOP); 
            if (best_i-1 > TURN) 
                insert_node (0, best_i-1, FREE); 
            break;
            case 3: insert_node (best_i, j-1, LOOP); 
            if (best_i-1 > TURN) 
                insert_node (0, best_i-1, FREE); 
            break;
            case 4: insert_node (best_i+1, j-1, LOOP); 
            if (best_i-1 > TURN) 
                insert_node (0, best_i-1, FREE); 
            break;
            }
    }
    else if(cur_interval->type == M_WM)
    {
        int i = cur_interval->i;
        int j = cur_interval->j;
        PARAMTYPE tmp, min = INF;
        int best_k, best_row;
    
        if (debug)
            printf ("\t (%d,%d) M_WM\n", i,j);
    
        tmp = V->get_energy(i,j) +
            AU_penalty (int_sequence[i], int_sequence[j]) +
            misc.multi_helix_penalty;
        if (tmp < min)
        {
            min = tmp;
            best_row = 1;
        }
    
        tmp = V->get_energy(i+1,j) +
                AU_penalty (int_sequence[i+1], int_sequence[j]) +
                dangle_bot [int_sequence[j]]
                [int_sequence[i+1]]
                [int_sequence[i]] +
                misc.multi_helix_penalty +
                misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i] = -1;
            tmp = tmp - loss (i,i);
        }                    
        if (tmp < min)
        {
            min = tmp;
            best_row = 2;
        }
            
        tmp = V->get_energy(i,j-1) +
                    AU_penalty (int_sequence[i], int_sequence[j-1]) +
                    dangle_top [int_sequence[j-1]]
                                [int_sequence[i]]
                                [int_sequence[j]] +
                    misc.multi_helix_penalty +
                    misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[j] = -1;
            tmp = tmp - loss (j,j);
        }                                
        if (tmp < min)
        {
            min = tmp;
            best_row = 3;
        }
    
        tmp = V->get_energy(i+1,j-1) +
                    AU_penalty (int_sequence[i+1], int_sequence[j-1]) +
                    dangle_bot [int_sequence[j-1]]
                                [int_sequence[i+1]]
                                [int_sequence[i]] +
                    dangle_top [int_sequence[j-1]]
                                [int_sequence[i+1]]
                                [int_sequence[j]] +
                    misc.multi_helix_penalty +
                    2*misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i] = -1;
            pred_pairings[j] = -1;
            tmp = tmp - loss (i,i) - loss (j,j);
        }                                
        if (tmp < min)
        {
            min = tmp;
            best_row = 4;
        }
            
        tmp = VM->get_energy_WM (i+1,j) + misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[i] = -1;
            tmp = tmp - loss (i,i);
        }                    
        if (tmp < min)
        {
            min = tmp;
            best_row = 5;
        }
            
        tmp = VM->get_energy_WM (i,j-1) + misc.multi_free_base_penalty;
        // add the loss
        if (pred_pairings != NULL)
        {
            pred_pairings[j] = -1;
            tmp = tmp - loss (j,j);
        }            
        
        if (tmp < min)
        {
            min = tmp;
            best_row = 6;
        }
    
        for (int k=i; k < j; k++)
        {
            tmp = VM->get_energy_WM (i, k) + VM->get_energy_WM (k+1, j);
            if (tmp < min)
            {
                min = tmp;
                best_k = k;
                best_row = 7;
            }
        }
        switch (best_row)
        {
            case 1: insert_node (i, j, LOOP); break;
            case 2: insert_node (i+1, j, LOOP); break;
            case 3: insert_node (i, j-1, LOOP); break;
            case 4: insert_node (i+1, j-1, LOOP); break;
            case 5: 
                if (j-i-1 > TURN)
                    insert_node (i+1, j, M_WM);
                break;
            case 6:
                if (j-1-i > TURN)
                    insert_node (i, j-1, M_WM);
                break;            
            case 7: 
                if (best_k-i > TURN) 
                    insert_node (i, best_k, M_WM);
                if (j-best_k-1 > TURN)
                    insert_node (best_k+1, j, M_WM); 
                break;
        }        
    }
}


void s_min_folding::backtrack_restricted (seq_interval *cur_interval, str_features *fres)
// PRE:  All matrixes V, VM, WM and W have been filled
// POST: Discover the MFE path   
{
    char type;

    if(cur_interval->type == LOOP)
    {
        int i = cur_interval->i;
        int j = cur_interval->j;
        if (i >= j) 
            return;
        f[i].pair = j;
        f[j].pair = i;
        structure[i] = '(';
        structure[j] = ')';
      
        type = V->get_type (i,j);
        if (debug) 
            printf ("\t(%d,%d) LOOP - type %c\n", i,j,type);
        if (type == STACK)
        {
            f[i].type = STACK;
            f[j].type = STACK;
            if (i+1 < j-1)            
                insert_node (i+1, j-1, LOOP);
            else
            {
                printf ("NOT GOOD RESTR STACK, i=%d, j=%d\n", i, j);
                exit (0);
            }    
            
        }
        else if (type == HAIRP)
        {
            f[i].type = HAIRP;
            f[j].type = HAIRP;
        }
        else if (type == INTER) 
        {
            f[i].type = INTER;
            f[j].type = INTER;
            // detect the other closing pair
            int ip, jp, best_ip, best_jp, minq;
            PARAMTYPE tmp, min = INF;
            for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1); ip++) 
            {
                minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    
                for (jp = minq; jp < j; jp++)
                {
                    if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres))
                        continue;
                    tmp = VBI->get_energy_str (i,j,ip,jp);
                    if (tmp < min)
                    {
                        min = tmp;
                        best_ip = ip;
                        best_jp = jp;
                    }
                }
            }
            if (best_ip < best_jp)
                insert_node (best_ip, best_jp, LOOP);
            else
            {
                printf ("NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
                exit (0);
            }                    
        }
        else if (type == MULTI)
        {
            f[i].type = MULTI;
            f[j].type = MULTI;
            int k, best_k, best_row;
            PARAMTYPE tmp, min = INF;
            for (k = i+1; k <= j-1; k++)
              {
                tmp = VM->get_energy_WM (i+1,k) + VM->get_energy_WM (k+1, j-1);
                if (tmp < min)
                  {
                    min = tmp;
                    best_k = k;
                    best_row = 1;
                  }
                if (fres[i+1].pair <= -1)
                {  
                    tmp = VM->get_energy_WM (i+2,k) + VM->get_energy_WM (k+1, j-1) + 
                    dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] + 
                    misc.multi_free_base_penalty;
                    if (tmp < min)
                    {
                        min = tmp;
                        best_k = k;
                        best_row = 2;
                    }
                }
                if (fres[j-1].pair <= -1)
                {    
                    tmp = VM->get_energy_WM (i+1,k) + VM->get_energy_WM (k+1, j-2) + 
                    dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] + 
                    misc.multi_free_base_penalty;
                    if (tmp < min)
                    {
                        min = tmp;
                        best_k = k;
                        best_row = 3;
                    }
                }
                if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
                {    
                    tmp = VM->get_energy_WM (i+2,k) + VM->get_energy_WM (k+1, j-2) + 
                    dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] + 
                    dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] + 
                    2*misc.multi_free_base_penalty;
                    if (tmp < min)
                    {
                        min = tmp;
                        best_k = k;
                        best_row = 4;
                    }
                }    
              }
            switch (best_row)
              {
              case 1: insert_node (i+1, best_k, M_WM);
                insert_node (best_k+1, j-1, M_WM); break;
              case 2: insert_node (i+2, best_k, M_WM);
                insert_node (best_k+1, j-1, M_WM); break;
              case 3: insert_node (i+1, best_k, M_WM);
                insert_node (best_k+1, j-2, M_WM); break;
              case 4: insert_node (i+2, best_k, M_WM);
                insert_node (best_k+1, j-2, M_WM); break;
              }
        }
    }
    else if(cur_interval->type == FREE)
    {
        int j = cur_interval->j;

        if (j==0) return;
        //if (j <= TURN) return; 
        
        PARAMTYPE min = INF, tmp, acc, energy_ij;
        int best_row, i, best_i;
    
        if (debug)
            printf ("\t(0,%d) FREE\n", j);
            
        // this case is for j unpaired, so I have to check that.
        if (fres[j].pair <= -1)  
        {
        //printf ("j=%d\n", j);
            tmp = W[j-1];
            if (tmp < min)
            {
                min = tmp;          
                best_row = 0;
            }
        }      
        for (i=0; i<=j-1; i++)    // no TURN
        {

            // Don't need to make sure i and j don't have to pair with something else
            //  it's INF, done in fold_sequence_restricted
            acc = (i-1>0) ? W[i-1] : 0;
            energy_ij = V->get_energy(i,j);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
//                 if (fres[i].pair == j)
//                 {
//                     min = tmp;
//                     best_i = i;
//                     best_row = 1;
//                     continue;
//                 }
                if (tmp < min)
                {
                min = tmp;
                best_i = i;
                best_row = 1;
                }
            }
                        
            if (fres[i].pair <= -1)
            {  
                energy_ij = V->get_energy(i+1,j);            
                if (energy_ij < INF)
                {
                    tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
                    tmp += dangle_bot [int_sequence[j]]
                        [int_sequence[i+1]]                    
                        [int_sequence[i]];
    //                 if (fres[i+1].pair == j)
    //                 {
    //                     min = tmp;
    //                     best_i = i;
    //                     best_row = 2;
    //                     continue;
    //                 }                    
                    if (tmp < min)
                    {
                        min = tmp;
                        best_i = i;
                        best_row = 2;
                    }
                }          
            }
            if (fres[j].pair <= -1)
            {      
                energy_ij = V->get_energy(i,j-1);
                if (energy_ij < INF)
                {
                    tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j-1]) + acc;
                    tmp += dangle_top [int_sequence[j-1]]
                        [int_sequence[i]]
                        [int_sequence[j]];
    //                 if (fres[i].pair == j-1)
    //                 {
    //                     min = tmp;
    //                     best_i = i;
    //                     best_row = 3; 
    //                     continue;
    //                 }                    
                    if (tmp < min)
                    {
                        min = tmp;
                        best_i = i;
                        best_row = 3;
                    }
                }                  
            }
            if (fres[i].pair <= -1 && fres[j].pair <= -1)
            {      
                energy_ij = V->get_energy(i+1,j-1);
                if (energy_ij < INF)
                {
                    tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
                    tmp += dangle_bot [int_sequence[j-1]]
                        [int_sequence[i+1]]
                        [int_sequence[i]];
                    tmp += dangle_top [int_sequence[j-1]]
                        [int_sequence[i+1]]
                        [int_sequence[j]];
    //                 if (fres[i+1].pair == j-1)
    //                 {
    //                     min = tmp;
    //                     best_i = i;
    //                     best_row = 4;
    //                     continue;
    //                 }                    
                    if (tmp < min)
                    {
                        min = tmp;
                        best_i = i;
                        best_row = 4;
                    }
                }
            }      
        }
        switch (best_row)
        {
            case 0: insert_node (0, j-1, FREE); break;
            case 1: insert_node (best_i, j, LOOP); 
                if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
                    insert_node (0, best_i-1, FREE); 
                break;
            case 2: insert_node (best_i+1, j, LOOP); 
                if (best_i-1 > 0) 
                    insert_node (0, best_i-1, FREE); 
                break;
            case 3: insert_node (best_i, j-1, LOOP); 
                if (best_i-1 > 0) 
                    insert_node (0, best_i-1, FREE); 
                break;
            case 4: insert_node (best_i+1, j-1, LOOP); 
                if (best_i-1 > 0) 
                    insert_node (0, best_i-1, FREE); 
                break;
        }
    }
  else if(cur_interval->type == M_WM)
    {
      int i = cur_interval->i;
      int j = cur_interval->j;
      PARAMTYPE tmp, min = INF;
      int best_k, best_row;

      if (debug)
        printf ("\t (%d,%d) M_WM\n", i,j);

      tmp = V->get_energy(i,j) +
        AU_penalty (int_sequence[i], int_sequence[j]) +
        misc.multi_helix_penalty;
      if (tmp < min)
        {
          min = tmp;
          best_row = 1;
        }
      if (fres[i].pair <= -1)
      {
          tmp = V->get_energy(i+1,j) +
                AU_penalty (int_sequence[i+1], int_sequence[j]) +
                dangle_bot [int_sequence[j]]
                [int_sequence[i+1]]
                [int_sequence[i]] +
                misc.multi_helix_penalty +
                misc.multi_free_base_penalty;
          if (tmp < min)
          {
              min = tmp;
              best_row = 2;
          }
      }      
      if (fres[j].pair <= -1)
      {    
          tmp = V->get_energy(i,j-1) +
                AU_penalty (int_sequence[i], int_sequence[j-1]) +
                dangle_top [int_sequence[j-1]]
                            [int_sequence[i]]
                            [int_sequence[j]] +
                misc.multi_helix_penalty +
                misc.multi_free_base_penalty;
          if (tmp < min)
          {
              min = tmp;
              best_row = 3;
          }
      }  
      if (fres[i].pair <= -1 && fres[j].pair <= -1)
      {
          tmp = V->get_energy(i+1,j-1) +
                AU_penalty (int_sequence[i+1], int_sequence[j-1]) +
                dangle_bot [int_sequence[j-1]]
                            [int_sequence[i+1]]
                            [int_sequence[i]] +
                dangle_top [int_sequence[j-1]]
                            [int_sequence[i+1]]
                            [int_sequence[j]] +
                misc.multi_helix_penalty +
                2*misc.multi_free_base_penalty;
          if (tmp < min)
          {
              min = tmp;
              best_row = 4;
          }
      }   
      if (fres[i].pair <= -1)
      {
          tmp = VM->get_energy_WM (i+1,j) + misc.multi_free_base_penalty;
          if (tmp < min)
          {
              min = tmp;
              best_row = 5;
          }
      }      
      if (fres[j].pair <= -1)
      {          
          tmp = VM->get_energy_WM (i,j-1) + misc.multi_free_base_penalty;
          if (tmp < min)
          {
              min = tmp;
              best_row = 6;
          }
      }      

      for (int k=i; k < j; k++)
        {
            tmp = VM->get_energy_WM (i, k) + VM->get_energy_WM (k+1, j);
            if (tmp < min)
              {
                min = tmp;
                best_k = k;
                best_row = 7;
              }
        }
      switch (best_row)
        {
          case 1: insert_node (i, j, LOOP); break;
          case 2: insert_node (i+1, j, LOOP); break;
          case 3: insert_node (i, j-1, LOOP); break;
          case 4: insert_node (i+1, j-1, LOOP); break;
          case 5: 
            if (j-i-1 > 0)
              insert_node (i+1, j, M_WM);
            break;
          case 6:
            if (j-1-i > 0)
              insert_node (i, j-1, M_WM);
            break;            
          case 7: 
            if (best_k-i > 0) 
              insert_node (i, best_k, M_WM);
            if (j-best_k-1 > 0)
              insert_node (best_k+1, j, M_WM); 
            break;
          }        
    }
}





PARAMTYPE s_min_folding::compute_W_br2 (int j)
// Compute the second branch of W formula.
//       This branch has to consider the AU_penalties and the dangling energies.
{
    PARAMTYPE min = INF, tmp, energy_ij = INF, acc;
    int i;

    for (i=0; i<=j-TURN-1; i++)
    {  
        acc = (i-1>0) ? W[i-1]: 0;
        
        energy_ij = V->get_energy(i,j);
        if (energy_ij < INF)
        {
            tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
            if (tmp < min)
            {
                min = tmp;          
            }
        }
        
        energy_ij = V->get_energy(i+1,j);            
        if (energy_ij < INF)
        {
            tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
            tmp += dangle_bot [int_sequence[j]]
                              [int_sequence[i+1]]
                              [int_sequence[i]];
            // add the loss
            if (pred_pairings != NULL)
            {
                pred_pairings[i] = -1;
                tmp = tmp - loss (i,i);
            }            
                              
            if (tmp < min)
            {
                min = tmp;
            }
        }
        
        energy_ij = V->get_energy(i,j-1);
        if (energy_ij < INF)
        {
            tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j-1]) + acc;
            tmp += dangle_top [int_sequence [j-1]]
                              [int_sequence [i]]
                              [int_sequence [j]];
            // add the loss
            if (pred_pairings != NULL)
            {
                pred_pairings[j] = -1;
                tmp = tmp - loss (j,j);
            }            
                              
            if (tmp < min)
            {
                min = tmp;
            }
        }

        
        energy_ij = V->get_energy(i+1,j-1);
        if (energy_ij < INF)
        {
            tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
            tmp += dangle_bot [int_sequence[j-1]]
                              [int_sequence[i+1]]
                              [int_sequence[i]];
            tmp += dangle_top [int_sequence [j-1]]
                              [int_sequence [i+1]]
                              [int_sequence [j]];
            // add the loss
            if (pred_pairings != NULL)
            {
                pred_pairings[i] = -1;
                pred_pairings[j] = -1;
                tmp = tmp - loss (i,i) - loss (j,j);
            }            
                              
            if (tmp < min)
            {
                min = tmp;
            }
        }
    }
    return min;    
}





void s_min_folding::compute_W (int j)
// compute W(j)
{
    PARAMTYPE m1, m2;
    m1 = W[j-1];
    // add the loss
    if (pred_pairings != NULL)
    {
        pred_pairings[j] = -1;
        m1 = m1 - loss (j,j);
    }            
    
    m2 = compute_W_br2(j);
    if (m1 < m2)     // does this make sense?    || (m1 >= MAXENERGY && m2 >= MAXENERGY))
    {
        W[j] = m1;
    }
    else
    {
        W[j] = m2;        
    }
}


PARAMTYPE s_min_folding::compute_W_br2_restricted (int j, str_features *fres, int &must_choose_this_branch)
// Compute the second branch of W formula.
//       This branch has to consider the AU_penalties and the dangling energies.
// May 16, 2007: this function had some bugs. 
{
    PARAMTYPE min = INF, tmp, energy_ij = INF, acc;
    int i;
    int chosen = 0;
    int best_i = 0;

    // I have to check if j is restricted
    // added Jan 28, 2006

 
    // j does not HAVE to pair (or not with a base downstream)
    must_choose_this_branch = 0;
    for (i=0; i<=j-1; i++)    // TURN shouldn't be there
    {  
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

        // We don't need to make sure i and j don't have to pair with something else,
        //  because that would be INF - done in fold_sequence_restricted
        acc = (i-1>0) ? W[i-1]: 0;
            
        energy_ij = V->get_energy(i,j);
        if (energy_ij < INF)
        {
            tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
            // May 16, 2007: I don't think I need this
//             if (fres[i].pair == j)
//             {
//                 must_choose_this_branch = 1;
//                 printf ("i=%d, j=%d, en=%d, chose i j\n", i, j, tmp);
//                 return tmp;
//             }
            if (tmp < min)
            {
                min = tmp;
                chosen = 21;        best_i = i;
                if (fres[i].pair == j)  must_choose_this_branch = 1;
                else                    must_choose_this_branch = 0;
            }
        }

        // I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
        if (fres[i].pair <= -1 && i+1 < j)
        {
            energy_ij = V->get_energy(i+1,j);            
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
                tmp += dangle_bot [int_sequence[j]]
                                [int_sequence[i+1]]
                                [int_sequence[i]];
                // May 16, 2007: I don't think I need this                                
//                 if (fres[i+1].pair == j)
//                 {
//                     must_choose_this_branch = 1;
//                     printf ("i=%d, j=%d, en=%d, chose i+1 j\n", i, j, tmp);
//                     return tmp;
//                 }                                
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 22;  best_i = i;
                    if (fres[i+1].pair == j)  must_choose_this_branch = 1;
                    else                      must_choose_this_branch = 0;                    
                }
            }
        }
         
        // I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
        if (fres[j].pair <= -1 && i < j-1)
        {   
            energy_ij = V->get_energy(i,j-1);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j-1]) + acc;
                tmp += dangle_top [int_sequence [j-1]]
                                [int_sequence [i]]
                                [int_sequence [j]];
                // May 16, 2007: I don't think I need this
//                 if (fres[i].pair == j-1)
//                 {
//                     must_choose_this_branch = 1;
//                     printf ("i=%d, j=%d, en=%d, chose i j-1\n", i, j, tmp);
//                     return tmp;
//                 }                                                                
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 23;  best_i = i;
                    if (fres[i].pair == j-1)  must_choose_this_branch = 1;
                    else                      must_choose_this_branch = 0;                    
                }
            }
        }    

        if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
        {
            energy_ij = V->get_energy(i+1,j-1);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
                tmp += dangle_bot [int_sequence[j-1]]
                                [int_sequence[i+1]]
                                [int_sequence[i]];
                tmp += dangle_top [int_sequence [j-1]]
                                [int_sequence [i+1]]
                                [int_sequence [j]];
                // May 16, 2007: I don't think I need this
//                 if (fres[i+1].pair == j-1)
//                 {
//                     must_choose_this_branch = 1;
//                     printf ("i=%d, j=%d, en=%d, chose i+1 j-1\n", i, j, tmp);
//                     return tmp;
//                 }                                                                
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 24;  best_i = i;
                    if (fres[i+1].pair == j-1)  must_choose_this_branch = 1;
                    else                        must_choose_this_branch = 0;                                        
                }
            }
        }    
    }
    //printf ("Chosen=%d, best_i=%d\n", chosen, best_i);
    return min;    
}


void s_min_folding::compute_W_restricted (int j, str_features *fres)
// compute W(j)
{
    PARAMTYPE m1, m2;
    int must_choose_this_branch;
    //if (fres[j].pair <= -1)
        m1 = W[j-1];
    //else m1 = INF;    // not sure about this!! TO COME BACK HERE!!!!!!
    m2 = compute_W_br2_restricted (j, fres, must_choose_this_branch);
    if (must_choose_this_branch)
    {
        W[j] = m2;
        //printf ("j=%d, chose branch 2, W[%d]=%d\n", j, j, W[j]);
    }
    else
    {
        if (m1 < m2) // this is kind of stupid, dunno why it's here    || (m1 >= MAXENERGY && m2 >= MAXENERGY))
        {
            W[j] = m1;
            //printf ("j=%d, chose branch 1, W[%d]=%d\n", j, j, W[j]);
        }
        else
        {        
            W[j] = m2;
            //printf ("j=%d, chose branch 2, W[%d]=%d\n", j, j, W[j]);
        }
    }        
}



void s_min_folding::print_result ()
// PRE:  The matrix V has been calculated and the results written in f
// POST: Prints details of each elementary structure
{
    int i;
//  char type[5][10] = {"hairpin","stacked","internal", "multi", "bulge"};
    PARAMTYPE energy = INF, sum;

    printf ("Minimum energy: %d\n", W[nb_nucleotides-1]);    
    sum = 0;

    for (i=0; i< nb_nucleotides; i++)
    {
        if (f[i].pair > i)
        {
            if (f[i].type == HAIRP)
                energy = V->get_energy(i, f[i].pair);
            else if (f[i].type == STACK)
                energy = V->get_energy(i, f[i].pair) - V->get_energy(i+1, f[i+1].pair);
            /*
            else if (f[i].type == INTER)
            {
                V_node = V->get_node(i, f[i].pair);    
                energy = V_node->energy - V->get_energy (((internal_details *)V_node->details)->prime.i,
                                                         ((internal_details *)V_node->details)->prime.j);
            }

            else if (f[i].type == MULTI)
            {
                multi_details *m_node;
                V_node = V->get_node(i, f[i].pair);    
                m_node = (multi_details*)V_node->details;
                energy = V_node->energy;
                for (k=0; k< m_node->num_branches; k++)
                {
                    pair branch = m_node->branches[k];
                    energy -= V->get_energy (branch.i, branch.j);
                }
            }
            */
            printf ("Pair (%d,%d), type %c,\tenergy %6d\n", i, f[i].pair, f[i].type, energy);
            sum += energy;
        }
    }
    printf ("0....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8\n");
    printf ("%s\n", sequence);
    printf ("%s\n", structure);

}
