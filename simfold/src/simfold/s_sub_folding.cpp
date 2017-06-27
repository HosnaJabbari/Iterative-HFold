/***************************************************************************
                          s_sub_folding.h  -  description
                             -------------------
    begin                : Thu Apr 11 2002
    copyright            : (C) 2003 orginally started by Zhi Chuan Zhang, 
                           continued by Mirela Andronescu
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

// this class represents free energy minimization and backtracking for pairfold
// WITH suboptimal structures. Recurrences are from Wutchy complete suboptimal 
// folding paper. The dangling free energies for suboptimal structures are not like
// in the MFE folding, but the 5' and 3' dangling ends are added all in any case.


#include "s_sub_folding.h"
#include "constants.h"
#include "externs.h"
#include "common.h"
#include "simfold.h"
#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"
#include "s_multi_loop_sub.h"
#include "s_energy_matrix.h"


s_sub_folding::s_sub_folding(char *sequence, PARAMTYPE var)
// CONSTRUCTOR
{
	printf("in s_subfolding.cpp: calling allocate_space()\n");
    allocate_space (sequence, var);
}


s_sub_folding::s_sub_folding(char *sequence, char *restricted, PARAMTYPE var)
// CONSTRUCTOR, the restricted case
{
    this->restricted = restricted;
    allocate_space (sequence, var);
}


void s_sub_folding::allocate_space (char *sequence, PARAMTYPE var)
{
    /* Set energy variation, min_energy will be initilized when W table is done*/
	printf("in s_subfolding.cpp: beginning of allocate_space()\n");
    en_var = var;
    // M: added:
    num_partial_structures = 0;
    num_partial_structures_thrown_away = 0;
    num_complete_structures = 0;
    limit = MAXSUBSTR;
    nb_nucleotides = strlen(sequence);
    folding_list = NULL;
    tail_folding_list = NULL;    
    this->sequence = new char[nb_nucleotides+1];
    if (this->sequence == NULL) giveup ("Cannot allocate memory", "s_sub_folding");
    int i;
    for(i=0; i<nb_nucleotides; i++)
    {
        this->sequence[i] = sequence[i];
    }
    this->sequence[i] = '\0';

    // W(j) is the mfe of s0 to sj
    W = new PARAMTYPE [nb_nucleotides];
    if (W == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) W[i] = 0;

    // Integer representation of the sequence 
    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);

    H = new s_hairpin_loop (sequence, int_sequence, nb_nucleotides);
    if (H == NULL) giveup ("Cannot allocate memory", "energy");
    S = new s_stacked_pair (int_sequence, nb_nucleotides);
    if (S == NULL) giveup ("Cannot allocate memory", "energy");
    VBI = new s_internal_loop (int_sequence, nb_nucleotides);
    if (VBI == NULL) giveup ("Cannot allocate memory", "energy");
    VM_sub = new s_multi_loop_sub (int_sequence, nb_nucleotides);
    if (VM_sub == NULL) giveup ("Cannot allocate memory", "energy");
    V = new s_energy_matrix (int_sequence, nb_nucleotides);
    if (V == NULL) giveup ("Cannot allocate memory", "energy");

    S->set_energy_matrix (V);
    VBI->set_energy_matrix (V);
    VM_sub->set_energy_matrix (V);
    V->set_loops (H, S, VBI, NULL, VM_sub);
    result_list = NULL;
	printf("in s_subfolding.cpp: allocate_space() DONE\n");
}


s_sub_folding::~s_sub_folding()
// DESTRUCTOR
{          
    delete [] sequence;
    delete [] int_sequence;  

    delete [] W;
  
    delete V;
    delete VM_sub;
    delete VBI;
    delete S;
    delete H;
    
    // release the result list       
    struct_node *tmp;
    tmp = result_list;
    while(tmp != NULL)
    {     
        result_list = result_list->next;        
        release_struct(tmp);
        tmp = result_list;
    }     
}


double s_sub_folding::s_simfold (double &enthalpy)
// fold sequence
{
	printf("in s_sub_folding.cpp->s_simfold \n");
    double energy;
    check_sequence (sequence);
	printf("in s_sub_folding.cpp->s_simfold: check_sequence() DONE \n");
    energy = fold_sequence (enthalpy);
	printf("in s_sub_folding.cpp->s_simfold: fold_sequence() DONE  and energy is: %f\n", energy);
    return energy;
}


double s_sub_folding::s_simfold_restricted (double &enthalpy)
// fold sequence, the restricted case
{
    double energy;
    check_sequence (sequence);
    energy = fold_sequence_restricted (enthalpy);
    return energy;
} 



double s_sub_folding::fold_sequence (double &enthalpy)
// PRE: the p_sub_folding object has been constructed and the resources has been
//      allocated
// POST: the sequence has been folded and the result been put into result_list
{
	printf("in s_sub_folding.cpp -> fold_sequence() beginning\n");
    double energy;
    int i, j;
 
    /* 1). Fill table V and WM, FM1, FM*/
    for (j=0; j < nb_nucleotides; j++)
    {         
        for (i=0; i<j; i++)
        {      
            V->compute_energy_sub (i,j);
        }            
        //  Compute the FM Table
        // TODO: uncomment
        if (!ignore_multi)
        {
            VM_sub->compute_energy_FM1 (j); // it has to be after V(i,j), because FM1(i,j) needs V(i,j)
            VM_sub->compute_energy_FM (j);
        }
    }

    /* 2). Fill the table W */
    for (j=1; j < nb_nucleotides; j++)
    {
        compute_W(j);
    }

    
    // backtrack: 
    // 1). Initialize the struct_node list to point to the only struct_node now
    //     
    
    seq_interval* seq = new seq_interval;
    if(!seq){
        giveup("s_sub_folding", "no memory");
    }
    seq->energy =W[nb_nucleotides-1];
    seq->i=0;
    seq->j=nb_nucleotides-1;
    seq->next=NULL;
    seq->type='W';

    struct_node* st_node = new struct_node; // REMEMBER TO DELETE IT
    struct_node* tmp_node;
    if(!st_node){
        giveup("s_sub_folding", "no memory");
    }

    st_node->f = new minimum_fold [nb_nucleotides];    
    if (st_node->f == NULL) giveup ("Cannot allocate memory", "energy");
    // M: added: set all pairs to -1
    for (i=0; i < nb_nucleotides; i++)
    {
        st_node->f[i].pair = -1;
        st_node->f[i].type = NONE;
        st_node->f[i].filled = 'N';
    }
    
    st_node->structure = new char[nb_nucleotides+2];      // I may want to put a space inside too
    if (st_node->structure == NULL) giveup ("Cannot allocate memory", "s_sub_folding");
    for (i=0; i<nb_nucleotides; i++)
    {
        st_node->structure[i] = '.';
    }
    st_node->structure[nb_nucleotides] = '\0';
    
    st_node->intervals = seq;
    st_node->energy = W[nb_nucleotides-1];
    st_node->previous = NULL;
    st_node->next = NULL;

    cur_folding = NULL;
    cur_interval = NULL;
    // 2). change the code for backtrack

    folding_list = st_node; // No need to initialize in constructor    
    num_partial_structures++;
	printf("in s_sub_folding.cpp -> fold_sequence(): and num_partial_structures = %d \n",num_partial_structures);
    min_energy = W[nb_nucleotides-1];
   // if (debug) 
		printf ("in s_sub_folding.cpp -> fold_sequence(): Min energy: %d\n", min_energy);
    max_energy = min_energy + en_var;
	printf ("in s_sub_folding.cpp -> fold_sequence(): Max energy: %d\n", max_energy);
	
	if (folding_list == NULL){
		printf("in s_sub_folding.cpp -> fold_sequence(): folding_list is NULL when it should not be!!!\n");
	}

    // Take out one by one the partial structures and continue folding them.
    while(folding_list != NULL)  // means there are still partial structures on the stack of partial structures
    {
        if(folding_list->intervals != NULL)        
        {
            if (debug)
                printf ("========\nPop the next partial structure from folding_list\n=========\n");
            cur_folding = folding_list;
            tmp_node = folding_list;
            folding_list = folding_list->next;
            num_partial_structures--;
			//printf("in s_sub_folding.cpp-> fold_sequence(): now the num_partial_structures = %d \n",num_partial_structures);

            cur_interval = cur_folding->intervals;
            //cur_folding->intervals = cur_interval->next;    //(I added ->next to cur_folding->intervals)       

            backtrack();            
            release_struct(tmp_node);
            //delete cur_interval;
        }
		/*
		if (folding_list == NULL){
			printf("in s_sub_solding.cpp -> fold_sequence(), and folding_list is NULL when it should not be!!\n");
		}
		if (folding_list != NULL && folding_list->intervals != NULL){
			printf("in s_sub_solding.cpp -> fold_sequence(), and folding_list->intervals is NOT NULL when it should be!!\n");
		}
        */
		
        // Add result to result_list.
        if(folding_list != NULL && folding_list->intervals == NULL)
        {
            //if (debug)
            {
                printf ("------\nConstruction of one structure has finished ---- 2 \n------\n");
            }
            // The construction of this structure has been finished
            // SHould be put into the structure list for print out
            if (result_list == NULL) // THIS IS NOT CORRECT
            {
                result_list = folding_list;
                folding_list = folding_list->next;
                last_list = result_list;
                last_list->next = NULL;
            }
            else
            {
                last_list->next = folding_list;
                folding_list = folding_list->next;
                last_list = last_list->next;
                last_list->next = NULL;
            }
            num_complete_structures++;
            num_partial_structures--;
            if (num_complete_structures >= limit)  // we are done
                break;
        }
        
    }// outer while
    
    energy = W[nb_nucleotides-1];
	printf("in s_sub_folding.cpp ->fold_sequence(): energy=%f \n",energy);
    return energy;
}



double s_sub_folding::fold_sequence_restricted (double &enthalpy)
// PRE: the p_sub_folding object has been constructed and the resources has been
//      allocated
// POST: the sequence has been folded and the result been put into result_list
{
    double energy;
    int i, j;
 
    str_features *fres;
    if ((fres = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");   
    // detect the structure features
    detect_structure_features (restricted, fres);        
    
    /* 1). Fill table V and WM, FM1, FM*/
    for (j=0; j < nb_nucleotides; j++)
    {         
        for (i=0; i<j; i++)
        {      
            // V(i,j) = infinity if i restricted or j restricted and pair of i is not j
            if ((fres[i].pair > -1 && fres[i].pair !=j) || (fres[j].pair > -1 && fres[j].pair != i)) 
                continue;                                    
            V->compute_energy_sub_restricted (i, j, fres);
        }            
        //  Compute the FM Table
        VM_sub->compute_energy_FM1_restricted (j, fres); // it has to be after V(i,j), because FM1(i,j) needs V(i,j)
        VM_sub->compute_energy_FM_restricted (j, fres);  
    }

    /* 2). Fill the table W */
    for (j=1; j < nb_nucleotides; j++)
    {
        compute_W_restricted (j, fres);
    }

    
    // backtrack: 
    // 1). Initialize the struct_node list to point to the only struct_node now
    //     
    
    seq_interval* seq = new seq_interval;
    if(!seq){
        giveup("s_sub_folding", "no memory");
    }
    seq->energy =W[nb_nucleotides-1];
    seq->i=0;
    seq->j=nb_nucleotides-1;
    seq->next=NULL;
    seq->type='W';

    struct_node* st_node = new struct_node; // REMEMBER TO DELETE IT
    struct_node* tmp_node;
    if(!st_node){
        giveup("s_sub_folding", "no memory");
    }

    st_node->f = new minimum_fold [nb_nucleotides];    
    if (st_node->f == NULL) giveup ("Cannot allocate memory", "energy");
    // M: added: set all pairs to -1
    for (i=0; i < nb_nucleotides; i++)
    {
        st_node->f[i].pair = -1;
        st_node->f[i].type = NONE;
        st_node->f[i].filled = 'N';
    }
    
    st_node->structure = new char[nb_nucleotides+2];      // I may want to put a space inside too
    if (st_node->structure == NULL) giveup ("Cannot allocate memory", "s_sub_folding");
    for (i=0; i<nb_nucleotides; i++)
    {
        st_node->structure[i] = '.';
    }
    st_node->structure[nb_nucleotides] = '\0';
    
    st_node->intervals = seq;
    st_node->energy = W[nb_nucleotides-1];
    st_node->previous = NULL;
    st_node->next = NULL;

    cur_folding = NULL;
    cur_interval = NULL;
    // 2). change the code for backtrack

    folding_list = st_node; // No need to initialize in constructor    
    num_partial_structures++;
    min_energy = W[nb_nucleotides-1];
    if (debug) printf ("Min energy: %d\n", min_energy);
    max_energy = min_energy + en_var;

    // Take out one by one the partial structures and continue folding them.
    while(folding_list != NULL)  // means there are still partial structures on the stack of partial structures
    {
        if(folding_list->intervals != NULL)        
        {
            if (debug)
                printf ("========\nPop the next partial structure from folding_list\n=========\n");
            cur_folding = folding_list;
            tmp_node = folding_list;
            folding_list = folding_list->next;
            num_partial_structures--;

            cur_interval = cur_folding->intervals;
            //cur_folding->intervals = cur_interval->next;    //(I added ->next to cur_folding->intervals)       

            backtrack_restricted (fres);
            release_struct(tmp_node);
            //delete cur_interval;
        }

        
        // Add result to result_list.
        if(folding_list != NULL && folding_list->intervals == NULL)
        {
            if (debug)
            {
                printf ("------\nConstruction of one structure has finished ---- 2 \n------\n");
            }
            // The construction of this structure has been finished
            // SHould be put into the structure list for print out
            if (result_list == NULL) // THIS IS NOT CORRECT
            {
                result_list = folding_list;
                folding_list = folding_list->next;
                last_list = result_list;
                last_list->next = NULL;
            }
            else
            {
                last_list->next = folding_list;
                folding_list = folding_list->next;
                last_list = last_list->next;
                last_list->next = NULL;
            }
            num_complete_structures++;
            num_partial_structures--;
            if (num_complete_structures >= limit)  // we are done
                break;
        }        
    }// outer while
    
    energy = W[nb_nucleotides-1];
    delete [] fres;
    return energy;
}



void s_sub_folding::insert_node(struct_node* sn1)
// PRE: sn1 is a struct_node that contains the structure to be further folded
// POST: sn1 is inserted into folding_list s.t. the nodes are ordered after the energy.
// M: rewrote some parts
{
    num_partial_structures++;

    // folding_list is NULL
    if(folding_list == NULL)
    {
        folding_list = sn1;
        tail_folding_list = sn1;
    }

    // add at the head
    else if(sn1->energy < folding_list->energy)
    {
        sn1->next = folding_list;
        folding_list->previous = sn1;
        folding_list = sn1;
      // tail_folding_list is not modified     
    }

    // find the right place
    else
    {
        struct_node* tmp;
        tmp = folding_list;
        while (tmp->next != NULL)
        {
            if (sn1->energy < tmp->next->energy)
            {
                sn1->next = tmp->next;
                sn1->previous = tmp;
                tmp->next->previous = sn1;
                tmp->next = sn1;
                break;
            }
            tmp = tmp->next;
        }
        if (tmp->next == NULL)  // insert at the end
        {
            tmp->next = sn1;
            sn1->previous = tmp;
            tail_folding_list = sn1;
        }
    }

    // M: We shouldn't keep more nodes than limit on folding_list 
    if (num_complete_structures + num_partial_structures > limit)   
    // if limit exceeded, remove the last element; should be exceeded only by 1.
    {
        struct_node *tmp;
        tmp = tail_folding_list;
        tail_folding_list = tail_folding_list->previous;
        tail_folding_list->next = NULL;
        num_partial_structures--;
        release_struct(tmp);
        num_partial_structures_thrown_away++;    
    }

    // M: update max_energy
 
    if (num_complete_structures + num_partial_structures == limit && num_partial_structures > 0)
    {
        if (tail_folding_list->energy < max_energy)
            max_energy = tail_folding_list->energy;
    }    
  
    // print the stack of substructures
    if (debug)
    {
        printf ("The stack of substructures folding_list:\n%d partial structures\n%d complete structures:\nMax energy: %d\n", 
          num_partial_structures, num_complete_structures, max_energy);
        struct_node *tmp;
        seq_interval *tmpint;
        printf ("0....,....1....,....2....,....3....,....4....,....5\n");
        printf ("%s\n", sequence);
        for (tmp = folding_list; tmp != NULL; tmp = tmp->next)
        {
            printf ("%s\t%d\nInts: ", tmp->structure, tmp->energy);
            for (tmpint = tmp->intervals; tmpint != NULL; tmpint = tmpint->next)
            {
                printf ("(%d,%d,%c) ", tmpint->i, tmpint->j, tmpint->type);
            }
            printf ("\nPins: ");
            for (int i=0; i < nb_nucleotides; i++)
            {
                if (tmp->f[i].pair > i )
                    printf ("(%d,%d,%c)", i, tmp->f[i].pair, tmp->f[i].type);
            }
            printf ("\n");
        }
        printf ("\n");
      //  printf ("Tail:     %s\t%d\n", tail_folding_list->structure, tail_folding_list->energy);
    }
}



void s_sub_folding::backtrack()
// PRE:  All cells in Vs matrix have been computed
//       j is a column in the matrix Vs
// POST: Read info from Vs matrix, and write into f and structure
//       This function is called when we don't know the pair of j, so we search for the minimum
{
	//printf("in s_sub_folding.cpp ->backtrack()\n");
    int i,j;    
    if(cur_interval->type == LOOP)
    {
        i = cur_interval->i;
        j = cur_interval->j;

        // Each backtrack will add valid structures to the structure list
        backtrack_hairpin(i, j);
        // TODO: uncomment
        if (!ignore_internal)
            backtrack_VBI(i, j);
        backtrack_stack(i, j);
        // TODO: uncomment
        if (!ignore_multi)
            backtrack_multi(i, j);
    }
    else if(cur_interval->type == FREE)
    {
        i = cur_interval->i;
        j = cur_interval->j;
        backtrack_freebases(i, j);
    }
    else if(cur_interval->type == M_FM1 && !ignore_multi)
    {
        i = cur_interval->i;
        j = cur_interval->j;
        backtrack_MFM1(i, j);    
    }
    else if(cur_interval->type == M_FM && !ignore_multi)
    {
        i = cur_interval->i;
        j = cur_interval->j;
        backtrack_MFM(i, j);
    }
    else
    {
         printf("ERROR: INVALID TYPE FOR BACKTRACKING, type = %c \n", cur_interval->type);
    }
	//printf("in s_sub_folding.cpp ->backtrack() DONE!\n");
}


void s_sub_folding::backtrack_restricted (str_features *fres)
// PRE:  All cells in Vs matrix have been computed
//       j is a column in the matrix Vs
// POST: Read info from Vs matrix, and write into f and structure
//       This function is called when we don't know the pair of j, so we search for the minimum
{
    int i,j;    
    if(cur_interval->type == LOOP)
    {
        i = cur_interval->i;
        j = cur_interval->j;

        // Each backtrack will add valid structures to the structure list
        backtrack_hairpin_restricted (i, j, fres);
        backtrack_VBI_restricted (i, j, fres);
        backtrack_stack(i, j);
        backtrack_multi(i, j);    
    }
    else if(cur_interval->type == FREE)
    {
        i = cur_interval->i;
        j = cur_interval->j;
        backtrack_freebases_restricted (i, j, fres);
    }
    else if(cur_interval->type == M_FM1)
    {
        i = cur_interval->i;
        j = cur_interval->j;
        backtrack_MFM1_restricted (i, j, fres);
    }
    else if(cur_interval->type == M_FM)
    {
        i = cur_interval->i;
        j = cur_interval->j;
        backtrack_MFM_restricted (i, j, fres);
    }
    else
    {
         printf("ERROR: INVALID TYPE FOR BACKTRACKING, type = %c \n", cur_interval->type);
    }
}


void s_sub_folding::backtrack_hairpin(int i, int j)
// PRE: sequence[i] to sequence[j] forms a loop
// POST: the loop is backtracked, its decomposition is reflected in structure

{
    PARAMTYPE h_energy, V_en, increment;
    if(j<i || j<0)
    {
        return;
    }

    V_en = V->get_energy(i, j);
    h_energy = H->compute_energy (i, j); //OK
    increment = h_energy-V_en;
    if((h_energy < INF) && (V_en < INF) && (increment+cur_folding->energy) <= max_energy)
    {
        struct_node* sn = copy_struct();
        sn->f[i].pair = j;
        sn->f[i].filled = 'Y';
        sn->f[i].type = HAIRP;
        sn->f[j].pair = i;
        sn->f[j].filled = 'Y';
        sn->f[j].type = HAIRP;
        sn->structure[i] = '(';
        sn->structure[j] = ')';
        /*
        // all bases that are between i and j are free bases
        for (k=j-1; k>i; k--)
          {
        sn->f[k].pair = -1;
        sn->f[k].filled = 'Y';
        sn->f[k].type = NONE;
        sn->structure[k] = '.';
          }
        */
        // Put the cur_folding back into the struct_list
        // Need to delete the cur_interval, memory leak
        
        sn->energy += increment;
        sn->next=NULL;

        if (debug)
          printf ("Insert node in bt_hairpin, i=%d, j=%d\n", i, j);
        insert_node(sn);
    }
    else 
        return;
}

void s_sub_folding::backtrack_hairpin_restricted (int i, int j, str_features *fres)
// PRE: sequence[i] to sequence[j] forms a loop
// POST: the loop is backtracked, its decomposition is reflected in structure

{
    PARAMTYPE h_energy, V_en, increment;
    if(j<i || j<0)
    {
        return;
    }

    V_en = V->get_energy(i, j);
    h_energy = H->compute_energy_restricted (i, j, fres); //OK
    increment = h_energy-V_en;
    if((h_energy < INF) && (V_en < INF) && (increment+cur_folding->energy) <= max_energy)
    {
        struct_node* sn = copy_struct();
        sn->f[i].pair = j;
        sn->f[i].filled = 'Y';
        sn->f[i].type = HAIRP;
        sn->f[j].pair = i;
        sn->f[j].filled = 'Y';
        sn->f[j].type = HAIRP;
        sn->structure[i] = '(';
        sn->structure[j] = ')';      
        sn->energy += increment;
        sn->next=NULL;

        if (debug)
          printf ("Insert node in bt_hairpin, i=%d, j=%d\n", i, j);
        insert_node(sn);
    }
    else 
        return;
}

void s_sub_folding::backtrack_VBI(int i, int j)
// PRE: sequence[i] to sequence[j] forms a loop
// POST: the loop is backtracked, its decomposition is reflected in structure
{
    // 2). i and j forms a internal loop
    //     Calculate energy direclty: VBI->calculate_energy(i, j);
    int ip, jp, minq;
    PARAMTYPE i_energy;
    PARAMTYPE increment;

    for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1) ; ip++)  // j-2-TURN 
      /*** MAXLOOP is the max size of a loop that is penalized ***/
    {
        minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // ip+1+TURN);
        for (jp = minq; jp < j; jp++)
        {        
            // M: make sure this is not a stacked pair
            if (jp == j-1 && ip == i+1)
                continue;

            //printf ("Before VBI call: ip=%d, jp=%d\n", ip, jp);
            i_energy = VBI->get_energy_str(i, j, ip, jp); 
            //printf ("After VBI call: i_energy = %d\n", i_energy);
            
            // Update structure and f
            increment = i_energy - V->get_energy(i, j);
            if((i_energy <INF) && (increment + cur_folding->energy) <= max_energy)          
            {
                // Need to create new struct and interval, put it onto list
                struct_node* sn = copy_struct();
                sn->f[i].pair = j;
                sn->f[i].filled = 'Y';
                sn->f[i].type = INTER;
                sn->f[j].pair = i;
                sn->f[j].filled = 'Y';
                sn->f[j].type = INTER;
                sn->structure[i] = '(';
                sn->structure[j] = ')';

                sn->energy += increment;                
                sn->next = NULL;
                // Need to update cur_folding energy
                   
                seq_interval* si = new seq_interval;//();
                if(!si){
                    giveup("s_sub_folding", "no memeory");
                }                
                si->i = ip;
                si->j = jp;                
                si->type = LOOP;                

                si->next = sn->intervals;
                sn->intervals = si;           

                if (debug)
                    printf ("Insert node in bt_VBI, i=%d, j=%d, ip=%d, jp=%d, i_energy=%d, increment=%d, V(i,j)=%d, V(ip,jp)=%d, type=%c\n", 
                  i, j, ip, jp, i_energy, increment, V->get_energy(i,j), V->get_energy(ip,jp), V->get_type(ip,jp));
                insert_node(sn);
            }//if
        // add (m, n) onto intervals
        }//for           
    }//for
}


void s_sub_folding::backtrack_VBI_restricted (int i, int j, str_features *fres)
// PRE: sequence[i] to sequence[j] forms a loop
// POST: the loop is backtracked, its decomposition is reflected in structure
{
    // 2). i and j forms a internal loop
    //     Calculate energy direclty: VBI->calculate_energy(i, j);
    int ip, jp, minq;
    PARAMTYPE i_energy;
    PARAMTYPE increment;

    for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1) ; ip++)  // j-2-TURN 
      /*** MAXLOOP is the max size of a loop that is penalized ***/
    {
        minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // ip+1+TURN);
        for (jp = minq; jp < j; jp++)
        {        
            if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres))
                continue;        
            // M: make sure this is not a stacked pair
            if (jp == j-1 && ip == i+1)
                continue;

            //printf ("Before VBI call: ip=%d, jp=%d\n", ip, jp);
            i_energy = VBI->get_energy_str(i, j, ip, jp); 
            //printf ("After VBI call: i_energy = %d\n", i_energy);
            
            // Update structure and f
            increment = i_energy - V->get_energy(i, j);
            if((i_energy <INF) && (increment + cur_folding->energy) <= max_energy)          
            {
                // Need to create new struct and interval, put it onto list
                struct_node* sn = copy_struct();
                sn->f[i].pair = j;
                sn->f[i].filled = 'Y';
                sn->f[i].type = INTER;
                sn->f[j].pair = i;
                sn->f[j].filled = 'Y';
                sn->f[j].type = INTER;
                sn->structure[i] = '(';
                sn->structure[j] = ')';

                sn->energy += increment;                
                sn->next = NULL;
                // Need to update cur_folding energy
                   
                seq_interval* si = new seq_interval;//();
                if(!si){
                    giveup("s_sub_folding", "no memeory");
                }                
                si->i = ip;
                si->j = jp;                
                si->type = LOOP;                

                si->next = sn->intervals;
                sn->intervals = si;           

                if (debug)
                    printf ("Insert node in bt_VBI, i=%d, j=%d, ip=%d, jp=%d, i_energy=%d, increment=%d, V(i,j)=%d, V(ip,jp)=%d, type=%c\n", 
                  i, j, ip, jp, i_energy, increment, V->get_energy(i,j), V->get_energy(ip,jp), V->get_type(ip,jp));
                insert_node(sn);
            }//if
        // add (m, n) onto intervals
        }//for           
    }//for
}


void s_sub_folding::backtrack_stack(int i, int j)
// PRE: Sequence[i] to sequence[j] forms a loop
// POST: 
{
    PARAMTYPE s_energy, increment;
    s_energy = S->compute_energy(i, j); //
    increment = s_energy-V->get_energy(i, j);
    if((s_energy < INF)&&(increment+cur_folding->energy) <= max_energy)
    {    
        struct_node* sn1;
        sn1 = copy_struct();
        sn1->f[i].pair = j;
        sn1->f[i].filled = 'Y';
        sn1->f[i].type = STACK;
        sn1->f[j].pair = i;
        sn1->f[j].filled = 'Y';
        sn1->f[j].type = STACK;
        sn1->structure[i] = '(';
        sn1->structure[j] = ')';
        seq_interval* tmp_in = new seq_interval;    
        if(!tmp_in){
            giveup("s_sub_folding", "no memeory");
        }
        tmp_in->next = sn1->intervals;        
        sn1->intervals = tmp_in;                     
        
        sn1->intervals->i = i+1;
        sn1->intervals->j = j-1;
        sn1->intervals->type = LOOP;
        sn1->intervals->energy = V->get_energy(i+1, j-1);//
        sn1->energy += increment; // need to add the increment        
        sn1->next = NULL;
    
        if (debug)
            printf ("Insert node in bt_stack, i=%d, j=%d\n", i, j);
        insert_node(sn1);
    }
}


void s_sub_folding::backtrack_multi(int i, int j)
// PRE: sequence[i] to sequence[j] forms a loop
// POST: the loop is backtracked, its decomposition is reflected in structure
{

    // decompose (i, j) into FM(i, m-1) + FM1(m, n-1) + FM1(n, j)
    // On the top-level, insure that multi-loop is not a hairpin or internal loop
    int k;
    PARAMTYPE fm_en, fm1_en1,  m_en, v_en;
    PARAMTYPE increment;
    v_en = V->get_energy(i, j);
    for (k = i+1; k < j; k++)
    {
        fm_en = VM_sub->get_FM_energy(i+1, k);
        fm1_en1 = VM_sub->get_FM1_energy(k+1, j-1);
        // M: added dangling energies
        m_en = fm_en + fm1_en1+misc.multi_offset+misc.multi_helix_penalty+
            AU_penalty (int_sequence[i], int_sequence[j]) + 
            dangle_top[int_sequence[i]][int_sequence[j]][int_sequence[i+1]] + 
            dangle_bot[int_sequence[i]][int_sequence[j]][int_sequence[j-1]];  // M: changed sequence into int_sequence
        increment = m_en - v_en;
        if( //(fm_en < INF && VM_sub->check_decomp(i+1, m-1) != -1) && (fm1_en1 <INF) &&
            (cur_folding-> energy + increment <= max_energy))
        {
            struct_node* sn = copy_struct();
            // update the structure
            sn->f[i].pair = j;
            sn->f[i].filled = 'Y';
            sn->f[i].type = MULTI;
            sn->f[j].pair = i;
            sn->f[j].filled = 'Y';
            sn->f[j].type = MULTI;
            sn->structure[i] = '(';
            sn->structure[j] = ')';
            // Add two new intervals: FM + FM1
            seq_interval* si1 = new seq_interval;
            seq_interval* si2 = new seq_interval;
            //seq_interval* si3 = new seq_interval;
            if(!si1 || !si2){
                giveup("s_sub_folding", "no memeory");
            }
            si1->i = i+1;
            si1->j = k;
            si1->type = M_FM;
            si1->next = NULL;
            si1->energy = fm_en; 
            si2->i = k+1;
            si2->j = j-1;
            si2->type = M_FM1;
            si2->next = si1;
            si2->energy = fm1_en1;
            // Update the intervals of sn and insert sn tp the struc_list;
            si1->next = sn->intervals;
            sn->intervals = si2;
            sn->energy += increment;
            //if (debug)
                printf ("Insert node in bt_multi, i=%d, j=%d\n", i, j);
            insert_node(sn);
        }                
    }
}


void s_sub_folding::backtrack_freebases(int i, int j)
// PRE: sequence[i] to sequence[j] are free bases.
// POST: sequence[i] to sequence[j] are decomposed into
//      smaller structure
{
    int k;
    PARAMTYPE increment;
    PARAMTYPE tmp, tmp1, tmp2, acc, energy_ij;
    int tmp_i, tmp_j;
    
    //if(i==j)// i=j=0
    if (j-i <= 1)   // M: changed this: if i and j are at distance <= TURN, there is no point in continuing
    {
        // M: We need these
        struct_node* t_node=copy_struct();
        if (debug)
            printf ("Insert node in bt_freebases, i=%d, j=%d\n", i, j);
        insert_node(t_node);
        return;
    }

//    for( k=0; k<=j-TURN-1; k++)
    for( k=j-1; k>=0; k--)
    {
        // M: only one situation, we take V(k,j), and we always add the dangling ends, except boundaries
        tmp = tmp1 = tmp2 = INF;
        acc = (k-1>0) ? W[k-1]: 0;

        // the first situation is when dangle_top should not be added
        // i.e. j+1 is paired or j is the last base of the sequence
        energy_ij = V->get_energy(k,j);        
        if (energy_ij < INF)
        {
            tmp = AU_penalty (int_sequence[k],int_sequence[j])+acc;
            tmp += energy_ij;
            // always add the 5' dangling energy, except when k is the first base in sequence
            if (k > 0)
            {
                tmp += dangle_bot[int_sequence[j]]
                  [int_sequence[k]]
                  [int_sequence[k-1]];
            }
            // add the 3' dangling energy only if j is not the last base in seq
            if (j < nb_nucleotides-1)
            {
                tmp += dangle_top[int_sequence [j]]
                  [int_sequence [k]]
                  [int_sequence [j+1]];
            }
            tmp_i = k;
            tmp_j = j;
                      
            if((cur_folding->energy + tmp - W[j] ) <= max_energy)
            {
                struct_node* sn = copy_struct();  // the new partial structure
                increment = tmp - W[j]; 
                
                // M: if the (0.j) loop is too small, don't add it in the intervals list
                seq_interval* tmp_in1;
                //if (k-1 <= TURN) // Hosna, Sep 11, 2012, it only makes sense when i=0, shouldn't we have a general case here?
				if (k-1-i <= TURN)
                {
                    tmp_in1 = NULL;
                }
                else
                {
                    tmp_in1 = new seq_interval;
                    if(!tmp_in1){
                        giveup("s_sub_folding", "no memory");
                    }
                    tmp_in1->i = i; //0;    // M: should it be 0 or i??? // Hosna, Sep 11,2012 I think it should be i not 0
                    tmp_in1->j = k-1;
                    tmp_in1->next=NULL;              
                    tmp_in1->type=FREE;              
                    //tmp_in1->energy = W[i-1]; // not used
                }
                
                seq_interval* tmp_in2 = new seq_interval;
                if(!tmp_in2){
                    giveup("s_sub_folding", "no memory");
                }
                tmp_in2->i = tmp_i;
                tmp_in2->j = tmp_j;
                
                tmp_in2->next=tmp_in1;
                tmp_in2->type=LOOP;          
                
                if(tmp_in1 == NULL){
                    tmp_in2->next = sn->intervals;
                }
                else{
                    tmp_in1->next = sn->intervals;
                }
                sn->intervals = tmp_in2;          
                sn->energy += increment;          
                sn->previous = NULL;
                sn->next=NULL;
                sn->f[tmp_i].pair = tmp_j;
                sn->f[tmp_i].filled = 'Y';
                sn->f[tmp_i].type = LOOP;
                sn->f[tmp_j].pair = tmp_i;
                sn->f[tmp_j].filled = 'Y';
                sn->f[tmp_j].type = LOOP;
                sn->structure[tmp_i] = '(';
                sn->structure[tmp_j] = ')';
                if (debug)
                    printf ("Insert node in bt_freebases, i=%d, j=%d, k=%d\n", i, j, k);
                insert_node(sn);
            }//if (cur_folding->energy ...
        } // if energy_ij   
    }//for

    // TO COME BACK: do I need to check for restriction here?
    PARAMTYPE wj1 = W[j-1];
    if((wj1 <INF) && ((wj1 + (cur_folding->energy-W[j])) <= max_energy))
    {
        // The second way to decompose the sequence
        struct_node* sn;
        sn = copy_struct();
    
        seq_interval* tmp_in = new seq_interval();
        if(!tmp_in){
            giveup("s_sub_folding", "no memory");
        }
        sn->f[j].pair = -1;
        sn->f[j].filled = 'Y';
        sn->f[j].type = NONE;
        sn->structure[j] = '.';
        sn->energy = cur_folding->energy - W[j] + W[j-1]; 
    
        tmp_in->i = i; //0; // Hosna, Sep 11, 2012, again shouldn't it be i instead??
        tmp_in->j = j-1;
        tmp_in->energy = W[j-1];
        tmp_in->type = FREE;
        tmp_in->next = sn->intervals;
        sn->intervals = tmp_in;
        sn->next=NULL;
        if (debug)
            printf ("Insert node in bt_freebases last, i=%d, j=%d, en=%d\n", i, j, sn->energy);
        insert_node(sn);
    }
}


void s_sub_folding::backtrack_freebases_restricted (int i, int j, str_features *fres)
// PRE: sequence[i] to sequence[j] are free bases.
// POST: sequence[i] to sequence[j] are decomposed into
//      smaller structure
{
    int k;
    PARAMTYPE increment;
    PARAMTYPE tmp, tmp1, tmp2, acc, energy_ij;
    int tmp_i, tmp_j;
    
    //if(i==j)// i=j=0
    if (j-i <= 1)   // M: changed this: if i and j are at distance <= TURN, there is no point in continuing
    {
        // M: We need these
        struct_node* t_node=copy_struct();
        if (debug)
            printf ("Insert node in bt_freebases, i=%d, j=%d\n", i, j);
        insert_node(t_node);
        return;
    }

//    for( k=0; k<=j-TURN-1; k++)
    for( k=j-1; k>=0; k--)
    {
        // M: only one situation, we take V(k,j), and we always add the dangling ends, except boundaries
        tmp = tmp1 = tmp2 = INF;
        acc = (k-1>0) ? W[k-1]: 0;

        // the first situation is when dangle_top should not be added
        // i.e. j+1 is paired or j is the last base of the sequence
        energy_ij = V->get_energy(k,j);        
        if (energy_ij < INF)
        {
            tmp = AU_penalty (int_sequence[k],int_sequence[j])+acc;
            tmp += energy_ij;
            // always add the 5' dangling energy, except when k is the first base in sequence
            if (k > 0)
            {
                tmp += dangle_bot[int_sequence[j]]
                  [int_sequence[k]]
                  [int_sequence[k-1]];
            }
            // add the 3' dangling energy only if j is not the last base in seq
            if (j < nb_nucleotides-1)
            {
                tmp += dangle_top[int_sequence [j]]
                  [int_sequence [k]]
                  [int_sequence [j+1]];
            }
            tmp_i = k;
            tmp_j = j;
                      
            if((cur_folding->energy + tmp - W[j] ) <= max_energy)
            {
                struct_node* sn = copy_struct();  // the new partial structure
                increment = tmp - W[j]; 
                
                // M: if the (0.j) loop is too small, don't add it in the intervals list
                seq_interval* tmp_in1;
                if (k-1 <= TURN)       // TODO: should this TURN be replaced, because it's in a "restricted" function?
                {
                    tmp_in1 = NULL;
                }
                else
                {
                    tmp_in1 = new seq_interval;
                    if(!tmp_in1){
                        giveup("s_sub_folding", "no memory");
                    }
                    tmp_in1->i = 0;    // M: should it be 0 or i???
                    tmp_in1->j = k-1;
                    tmp_in1->next=NULL;              
                    tmp_in1->type=FREE;              
                    if (i > 0)
                        tmp_in1->energy = W[i-1]; // not used
                }
                
                seq_interval* tmp_in2 = new seq_interval;
                if(!tmp_in2){
                    giveup("s_sub_folding", "no memory");
                }
                tmp_in2->i = tmp_i;
                tmp_in2->j = tmp_j;
                
                tmp_in2->next=tmp_in1;
                tmp_in2->type=LOOP;          
                
                if(tmp_in1 == NULL){
                    tmp_in2->next = sn->intervals;
                }
                else{
                    tmp_in1->next = sn->intervals;
                }
                sn->intervals = tmp_in2;          
                sn->energy += increment;          
                sn->previous = NULL;
                sn->next=NULL;
                sn->f[tmp_i].pair = tmp_j;
                sn->f[tmp_i].filled = 'Y';
                sn->f[tmp_i].type = LOOP;
                sn->f[tmp_j].pair = tmp_i;
                sn->f[tmp_j].filled = 'Y';
                sn->f[tmp_j].type = LOOP;
                sn->structure[tmp_i] = '(';
                sn->structure[tmp_j] = ')';
                if (debug)
                    printf ("Insert node in bt_freebases, i=%d, j=%d, k=%d\n", i, j, k);
                insert_node(sn);
            }//if (cur_folding->energy ...
        } // if energy_ij    
    }//for

    // this case is for j unpaired, so check this is true
    if (fres[j].pair <= -1)  
    {
        PARAMTYPE wj1 = W[j-1];
        if((wj1 <INF) && ((wj1 + (cur_folding->energy-W[j])) <= max_energy))
        {
            // The second way to decompose the sequence
            struct_node* sn;
            sn = copy_struct();
        
            seq_interval* tmp_in = new seq_interval();
            if(!tmp_in){
                giveup("s_sub_folding", "no memory");
            }
            sn->f[j].pair = -1;
            sn->f[j].filled = 'Y';
            sn->f[j].type = NONE;
            sn->structure[j] = '.';
            sn->energy = cur_folding->energy - W[j] + W[j-1]; 
        
            tmp_in->i = 0;
            tmp_in->j = j-1;
            tmp_in->energy = W[j-1];
            tmp_in->type = FREE;
            tmp_in->next = sn->intervals;
            sn->intervals = tmp_in;
            sn->next=NULL;
            if (debug)
                printf ("Insert node in bt_freebases last, i=%d, j=%d, en=%d\n", i, j, sn->energy);
            insert_node(sn);
        }
    }    
}


void s_sub_folding::backtrack_MFM1(int i, int j)
// PRE:
// POST:
{
    int k;
    PARAMTYPE fm1_energy, increment;
    PARAMTYPE tmp;
    tmp=INF;
    fm1_energy = VM_sub->get_FM1_energy(i, j);
    /* Find the FM1 energy for each k */
    for(k=i+1; k<=j; k++)
    {
        tmp = V->get_energy(i,k) + AU_penalty(int_sequence[i], int_sequence[k]) + 
          misc.multi_helix_penalty + (j-k)*misc.multi_free_base_penalty;
        if (k < nb_nucleotides-1)
            tmp += dangle_top [int_sequence [k]][int_sequence [i]][int_sequence [k+1]];
        if(i>0)
            tmp += dangle_bot [int_sequence[k]][int_sequence[i]][int_sequence[i-1]];
        
        increment = tmp-fm1_energy;
        if(//(fm1_energy != INF ) && (tmp != INF) && 
           (increment + cur_folding->energy) <= max_energy)
        {    
            struct_node* sn = copy_struct();
            // Need to create new struct and interval, put it onto list
            // Need to update cur_folding energy
            seq_interval* si = new seq_interval;
            if(!si){
                giveup("s_sub_folding", "no memory");
            }
            
            si->i = i;
            si->j = k;
            si->type = LOOP;
            si->energy = 0; // no use here
            si->next = sn->intervals;
            sn->intervals = si;//cur_interval;            
            sn->energy += increment;        
            sn->next = NULL;
            if (debug)
                printf ("Insert node in bt_MFM1, i=%d, j=%d\n", i, j);
            insert_node(sn);
        }
    }           
}


void s_sub_folding::backtrack_MFM1_restricted (int i, int j, str_features *fres)
// PRE:
// POST:
{
    int k;
    PARAMTYPE fm1_energy, increment;
    PARAMTYPE tmp;
    tmp=INF;
    fm1_energy = VM_sub->get_FM1_energy(i, j);
    /* Find the FM1 energy for each k */
    for(k=i+1; k<=j; k++)
    {
        if (exists_restricted (k, j+1, fres))
            continue;         
        tmp = V->get_energy(i,k) + AU_penalty(int_sequence[i], int_sequence[k]) + 
          misc.multi_helix_penalty + (j-k)*misc.multi_free_base_penalty;
        if (k < nb_nucleotides-1)
            tmp += dangle_top [int_sequence [k]][int_sequence [i]][int_sequence [k+1]];
        if(i>0)
            tmp += dangle_bot [int_sequence[k]][int_sequence[i]][int_sequence[i-1]];
        
        increment = tmp-fm1_energy;
        if(//(fm1_energy != INF ) && (tmp != INF) && 
           (increment + cur_folding->energy) <= max_energy)
        {    
            struct_node* sn = copy_struct();
            // Need to create new struct and interval, put it onto list
            // Need to update cur_folding energy
            seq_interval* si = new seq_interval;
            if(!si){
                giveup("s_sub_folding", "no memory");
            }
            
            si->i = i;
            si->j = k;
            si->type = LOOP;
            si->energy = 0; // no use here
            si->next = sn->intervals;
            sn->intervals = si;//cur_interval;            
            sn->energy += increment;        
            sn->next = NULL;
            if (debug)
                printf ("Insert node in bt_MFM1, i=%d, j=%d\n", i, j);
            insert_node(sn);
        }
    }           
}



void s_sub_folding::backtrack_MFM(int i, int j)
// PRE: sequence[i] to sequence[j] is part of a multi loop
//      j>i
// POST: MFM from i to j is further decomposed into FM1 and FM
{    
    int m;
    PARAMTYPE fm_en, fm1_en, tmp_energy, vme, increment;

    vme = VM_sub->get_FM_energy(i, j);

    // 2). FM = FM + FM1
    //for( m = i; m < j-4; m++)
    for (m = i+1; m < j; m ++)   // M: changed
    {
        fm_en = VM_sub->get_FM_energy(i, m);
        fm1_en = VM_sub->get_FM1_energy(m+1, j);
        tmp_energy = fm_en+fm1_en; 
        increment = tmp_energy - vme;    
        
        if(//(fm_en < INF) && (fm1_en < INF) && 
           (cur_folding->energy+increment <= max_energy))
        {
            if(increment < 0)
            {
                printf("SO WEIRED, %d, %d, %d\n",fm_en, fm1_en, VM_sub->get_FM_energy(i, j));
            }
    

            // The second way to decompose the sequence
            struct_node* sn = copy_struct();//new struct_node(*cur_folding);
            seq_interval* tmp_in1 = new seq_interval;        
            if(!tmp_in1){
                giveup("s_sub_folding", "no memory");
            }
            tmp_in1->i = i;
            tmp_in1->j = m;
            tmp_in1->next=NULL;
            tmp_in1->type=M_FM;
            tmp_in1->energy = 0;

            seq_interval* tmp_in2 = new seq_interval;
            if(!tmp_in2){
                giveup("s_sub_folding", "no memory");
            }
            tmp_in2->i = m+1;
            tmp_in2->j = j;
            tmp_in2->next=tmp_in1;
            tmp_in2->type=M_FM1;
            tmp_in2->energy = 0;

            tmp_in1->next = sn->intervals;
            sn->intervals = tmp_in2;
            sn->energy += increment;
            sn->next=NULL;
            if (debug)
                printf ("Insert node in bt_MFM, i=%d, j=%d\n", i, j);
            insert_node(sn);
        }
    }

    
    // M: 1). added: backtrack over m

    for (m = i; m < j; m++)  // M: added
    {
    
        fm1_en = VM_sub->get_FM1_energy(m, j);
        increment =  fm1_en + (m-i)*misc.multi_free_base_penalty - vme;  // M: added (m-i)*...
        if((cur_folding->energy + increment) <= max_energy)
        {
            struct_node* sn = copy_struct();
            seq_interval* tmp_in1 = new seq_interval;
            if(!tmp_in1){
                giveup("s_sub_folding", "no memory");
            }
            tmp_in1->i = m;
            tmp_in1->j = j;
            tmp_in1->next=NULL;
            tmp_in1->type=M_FM1;
            tmp_in1->energy = 0;
            tmp_in1->next = sn->intervals;
            sn->intervals = tmp_in1;
            sn->energy += increment;
            sn->next=NULL;
            if (debug)
                printf ("Insert node in bt_MFM, i=%d, j=%d\n", i, j);
            insert_node(sn);
        }
    }
}


void s_sub_folding::backtrack_MFM_restricted (int i, int j, str_features *fres)
// PRE: sequence[i] to sequence[j] is part of a multi loop
//      j>i
// POST: MFM from i to j is further decomposed into FM1 and FM
{    
    int m;
    PARAMTYPE fm_en, fm1_en, tmp_energy, vme, increment;

    vme = VM_sub->get_FM_energy(i, j);

    // 2). FM = FM + FM1
    //for( m = i; m < j-4; m++)
    for (m = i+1; m < j; m ++)   // M: changed
    {
        fm_en = VM_sub->get_FM_energy(i, m);
        fm1_en = VM_sub->get_FM1_energy(m+1, j);
        tmp_energy = fm_en+fm1_en; 
        increment = tmp_energy - vme;    
        
        if(//(fm_en < INF) && (fm1_en < INF) && 
           (cur_folding->energy+increment <= max_energy))
        {
            if(increment < 0)
            {
                printf("SO WEIRED, %d, %d, %d\n",fm_en, fm1_en, VM_sub->get_FM_energy(i, j));
            }
    

            // The second way to decompose the sequence
            struct_node* sn = copy_struct();//new struct_node(*cur_folding);
            seq_interval* tmp_in1 = new seq_interval;        
            if(!tmp_in1){
                giveup("s_sub_folding", "no memory");
            }
            tmp_in1->i = i;
            tmp_in1->j = m;
            tmp_in1->next=NULL;
            tmp_in1->type=M_FM;
            tmp_in1->energy = 0;

            seq_interval* tmp_in2 = new seq_interval;
            if(!tmp_in2){
                giveup("s_sub_folding", "no memory");
            }
            tmp_in2->i = m+1;
            tmp_in2->j = j;
            tmp_in2->next=tmp_in1;
            tmp_in2->type=M_FM1;
            tmp_in2->energy = 0;

            tmp_in1->next = sn->intervals;
            sn->intervals = tmp_in2;
            sn->energy += increment;
            sn->next=NULL;
            if (debug)
                printf ("Insert node in bt_MFM, i=%d, j=%d\n", i, j);
            insert_node(sn);
        }
    }

    
    // M: 1). added: backtrack over m

    for (m = i; m < j; m++)  // M: added
    {
        if (exists_restricted (i-1, m, fres))
            continue;
        fm1_en = VM_sub->get_FM1_energy(m, j);
        increment =  fm1_en + (m-i)*misc.multi_free_base_penalty - vme;  // M: added (m-i)*...
        if((cur_folding->energy + increment) <= max_energy)
        {
            struct_node* sn = copy_struct();
            seq_interval* tmp_in1 = new seq_interval;
            if(!tmp_in1){
                giveup("s_sub_folding", "no memory");
            }
            tmp_in1->i = m;
            tmp_in1->j = j;
            tmp_in1->next=NULL;
            tmp_in1->type=M_FM1;
            tmp_in1->energy = 0;
            tmp_in1->next = sn->intervals;
            sn->intervals = tmp_in1;
            sn->energy += increment;
            sn->next=NULL;
            if (debug)
                printf ("Insert node in bt_MFM, i=%d, j=%d\n", i, j);
            insert_node(sn);
        }
    }
}


void s_sub_folding::set_limit(int limit)
// PRE: sn is an allocated mem clocation
// POST: the mem is released
{
    if (limit < 0)
        return;
    this->limit = limit;
	printf("in s_sub_folding.cpp: set_limit() DONE\n");
}


void s_sub_folding::release_struct(struct_node* sn)
// PRE: sn is an allocated mem clocation
// POST: the mem is released
// M:added a few more lines 
{
    seq_interval*  tmp_seq;

    if(sn->f != NULL)
        delete [] sn->f;
    // delete intervals 
    
    tmp_seq = sn->intervals;
    while(tmp_seq != NULL)
    {        
        sn->intervals = sn->intervals->next;
        delete tmp_seq;
        tmp_seq = sn->intervals;
    }
    if(sn->structure != NULL)
        delete [] sn->structure; 

    delete sn;    
}


void s_sub_folding::copy_list (seq_interval* from, seq_interval *& to)
  // PRE: from is a linked list
  // POST: copy all linked list from "from" to "to"
  // M: added
{
    seq_interval *to_tail, *from_tmp, *to_tmp;
    to = NULL;
    to_tail = NULL;

    for (from_tmp = from; from_tmp != NULL; from_tmp = from_tmp->next)
    {
        to_tmp = new seq_interval;
        if(to_tmp == NULL){
            giveup("copy_list", "no memory");
        }
        from_tmp->copy(to_tmp);
        to_tmp->next = NULL;
        if (to == NULL)
        {
            to = to_tmp;
            to_tail = to_tmp;
        }
        else
        {
            to_tail->next = to_tmp;
            to_tail = to_tmp;
        }
    }
}


struct_node* s_sub_folding::copy_struct()
// PRE: Cur_folding is initialized
// POST: Cur_folding is copied and returned.
// M: added \0 at the end of structure -> it was showing weird characters
// M: this function is correct, but not too elegant: to be rewritten
{
    struct_node* sn = new struct_node;
    if(!sn){
        giveup("s_sub_folding", "no memory");
    }
    sn->intervals = NULL;
    sn->previous = NULL;
    sn->next = NULL;
    sn->energy = cur_folding->energy;
    sn->structure = new char[nb_nucleotides+1];   // M: added +1
    if(!sn->structure){
        giveup("s_sub_folding", "no memory");
    }
    int m;
    for(m = 0; m< nb_nucleotides; m++)
    {
        sn->structure[m] = cur_folding->structure[m];
    }
    sn->structure[m] = '\0';   // M: added

    sn->f = new minimum_fold[nb_nucleotides];
    if(!sn->f){
        giveup("s_sub_folding", "no memory");
    }

    for(m = 0; m< nb_nucleotides; m++)
    {
        sn->f[m].filled = cur_folding->f[m].filled;
        sn->f[m].pair = cur_folding->f[m].pair;
        sn->f[m].type = cur_folding->f[m].type;   
    }

    copy_list (cur_folding->intervals->next, sn->intervals);  // M: added

    return sn;
    //End of Copy 
}




PARAMTYPE s_sub_folding::compute_W_br2 (int j)
// PRE:  the nucleotide index j and the word it is in (bj) - are given
//      j >= 1
//       min_h is the minimum i that was found
// POST: The second branch of compute_Ws formula.
//       This branch has to consider the AU_penalties and the dangling energies.
{
    //next_back points to the next smaller W[j] 
    // type seems not used here 
    PARAMTYPE min, tmp, energy_ij, acc;
    int i;

    min = INF; energy_ij = INF;

    for (i=0; i < j; i++)
    {  
        acc = (i>0) ? W[i-1]:0;
        
        energy_ij = V->get_energy(i,j);
        if (energy_ij < INF)
        {
            tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
            if (i > 0)
            {
                tmp += dangle_bot[int_sequence[j]]
                  [int_sequence[i]]
                  [int_sequence[i-1]];
            }
            // add the 3' dangling energy only if j is not the last base in seq
            if (j < nb_nucleotides-1)
            {
                tmp += dangle_top[int_sequence [j]]
                  [int_sequence [i]]
                  [int_sequence [j+1]];
            }
            if (tmp < min)
            {
                min = tmp;
            }
        }        
    }
    return min;    
}



void s_sub_folding::compute_W (int j)
// PRE:  V, FM and FM1 tables have been initialized 
//       j >= 1
// POST: compute Ws
{
    PARAMTYPE m1, m2;
    m1 = W[j-1];
    m2 = compute_W_br2(j);
    if (m1 < m2)
    {
        W[j] = m1;
    }
    else
    {
        W[j] = m2;
    }
	//sprintf("W(%d) = %d \n",j,W[j]);
}

void s_sub_folding::compute_W_restricted (int j, str_features *fres)
// compute W(j)
{
    PARAMTYPE m1, m2;
    if (fres[j].pair <= -1)
        m1 = W[j-1];
    else m1 = INF;    // not sure about this!! 
    m2 = compute_W_br2 (j);
    if (m1 < m2)
    {
        W[j] = m1;
    }
    else
    {
        W[j] = m2;        
    }
}



int s_sub_folding::return_structures (char structures[][MAXSLEN], double energies[])
//int s_sub_folding::return_structures (char **structures, double energies[])
// PRE:  The sequence has been folded and the result is in the result_list.
//       The memory for structures and energies is already allocated.
// POST: Return the resulting structures and the corresponding energies.
//       Also return the number of structures actually found.
// (added by Mirela on Nov 14 2003)
{
	//printf("in s_sub_folding.cpp -> retun_structures and nb_nucleotides= %d\n", nb_nucleotides);
    int i, j;
    struct_node* tmp;
    tmp = result_list;
    
    i = 1;
	/*
	if (tmp == NULL){
		printf("result_list is NULL, and it should not be!!!\n");
	}
	 */
    while(tmp != NULL)
    {
        energies[i-1] = tmp->energy/100.00;
        tmp->structure[nb_nucleotides] = '\0';
        for (j=0; j < nb_nucleotides; j++)
            structures[i-1][j] = tmp->structure[j];
        structures[i-1][j] = '\0';    
        //strcpy (structures[i-1],tmp->structure);

        i++;
        tmp=tmp->next;
    }
	//printf("in s_sub_folding.cpp -> retun_structures() DONE and the number of structures found is %d \n",i-1);
    return i-1;
}


void s_sub_folding::print_result (int flag)
// PRE:  The sequence has been folded and the result is in the result_list
// POST: Prints results
{ 
    // SET THIS FLAG TO PRINT OUT THE INTERVALS TRAVERSED
    struct_node* tmp;
    tmp = result_list;
    int i;

    printf("Min energy is:        %5.3f\n\n", W[nb_nucleotides-1]/100.00);
    //seq_interval* seq_in;
    i = 1;
    while(tmp != NULL)
    {
        printf ("The No.%d Structure (%.2f Kcal/mol):\n",i, tmp->energy/100.00);    
        printf ("0....,....1....,....2....,....3....,....4....,....5....\n");
        printf ("%s\n", this->sequence);
        printf ("%s\n\n", tmp->structure);       
        i++;           
        tmp=tmp->next;
      }      
}
