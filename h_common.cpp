#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "externs.h"
#include "h_externs.h"
#include "common.h"
#include "h_struct.h"
#include "h_common.h"
#include "params.h"

// Hosna feb 12, 2008
#include "W_final.h"
#include "hfold.h"
// Hosna, May 3rd, 2012
#include "hfold_pkonly.h"
// Hosna, November 16, 2015

/*
 * This function is just the same as detect_original_pairs
 * but it also calculates the arcs for each base and saves them in arc_table
 *
 * The algorithm for finding the arcs is as follows:
 * When checking each base in the input structure,
 *
 * case 1) if (structure[i] == '.' or ' ' or '_')
 * 		1-a) if stack is not empty, then put arc_table[i] = top element on the stack
 * 		1-b) otherwise arc_table[i] = -1
 *
 * case 2) if (structure[i] == '(')
 * 		2-a) if stack is not empty, then put arc_table[i] = top element on the stack and push i in the stack
 * 		2-b) otherwise arc_table[i] = -1 and push i in the stack
 *
 * case 3) if (strcuture[i] == ')'), pop the stack and match it with j
 * 		3-a) if stack is not empty, then put arc_table[i] = top element on the stack
 * 		3-b) otherwise put arc_table[i] = -1
 *
 */

void detect_original_pairs_arcs(char *structure, int *p_table, int *arc_table)
// PRE:  structure contains the desired structure
// POST: pairs will contain the index of each base pair
//               or -1 if it does not pair
//
{
		//printf("structure: %s\n", structure);
        int i, j, struct_len;
        stack_ds st;
        h_init (&st);
        remove_space (structure);
		//printf("The given structure is: \n %s\n", structure);
        struct_len = strlen (structure);
	// Hosna March 8, 2012
	// since index i starts at 0 and stack top is also set to 0 to show stack is empty, if i=0 is paired then we have incorrect arc values!
	// So I am introducing STACK_EMPTY = -1 to h_common.h and change h_init and h_pop accordingly

        for (i=0; i < struct_len; i++)
          {
			  // Hosna March 8, 2012
			  // changing nested ifs to switch for optimality
			  switch (structure[i])
				{
					case '.':
					{
					  p_table[i] = RESTRICTED_UNPAIR;
					  if (st.top > STACK_EMPTY){//0){
		//              	if (debug)
		//	              	printf("base %d: stack has %d elements in it and its top element is %d \n",i, st.top, st.elem[st.top]);
						arc_table[i] = st.elem[st.top];
					  }else{
						arc_table[i] = NOT_COVERED;
					  }
					}
						break;
					case ' ':
					case '_':
					{
					  p_table[i] = FREE_TO_PAIR;
					  if (st.top > STACK_EMPTY){//0){
		//              	if (debug)
		//	              	printf("base %d: stack has %d elements in it and its top element is %d \n",i, st.top, st.elem[st.top]);
						arc_table[i] = st.elem[st.top];
					  }else{
						arc_table[i] = NOT_COVERED;
					  }
					}
						break;
					case '(':
						{
							if (st.top > STACK_EMPTY){//0){
			//              	if (debug)
			//	              	printf("base %d: stack has %d elements in it and its top element is %d \n",i, st.top, st.elem[st.top]);
								arc_table[i] = st.elem[st.top];
							  }else{
								arc_table[i] = NOT_COVERED;
							  }
							  h_push (&st, i);
						}
						break;
					case ')':
				//else if (structure[i] == ')')
					  {
						j = h_pop (&st);
						p_table[i] = j;
						p_table[j] = i;
						  if (st.top > STACK_EMPTY){//0){
		//                	if (debug)
		//	                	printf("base %d: stack has %d elements in it and its top element is %d \n",i, st.top, st.elem[st.top]);
							arc_table[i] = st.elem[st.top];
							}else{
								arc_table[i] = NOT_COVERED;
							}
					  }
						break;
			  }
          }

	/* for (i=137; i<struct_len; i++){
			printf("p_table[%d] = %d AND arc_table[%d]=%d \n\n",i,p_table[i],i,arc_table[i]);
		}
		printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"); */
        if (st.top != STACK_EMPTY)//0)
        {
            fprintf (stderr, "The given structure is not valid: %d more left parentheses than right parentheses\n", st.top);
            exit (1);
        }

}

/**
 * This function is just like the above function except the case that
 * it can handle the pseudoknotted structures of at most density 2
 *
 * the density 2 structures can be presented with ( and [ in dot parenthesis format
 * so we need at most 2 stacks to keep track of the base pairings
 * we call the stacks st and st_brack for ( and [ respectively.
 *
 *
 */

void detect_original_PKed_pairs(char *structure, int *p_table)
// PRE:  structure contains the desired structure
// POST: pairs will contain the index of each base pair
//               or -1 if it does not pair
{
	int i, j, struct_len;
	stack_ds st; //stach used for (
    stack_ds st_brack; // stack used for [
	h_init (&st);
	h_init (&st_brack);
	remove_space (structure);
	struct_len = strlen (structure);
	for (i=0; i < struct_len; i++) {
		// Hosna March 8, 2012
		// changing nested ifs to switch for optimality
		switch (structure[i]) {
			case '.':
				p_table[i] = -1;
				break;

			case ' ':
			case '_':
				p_table[i] = -2;
				break;

			case '(':
				h_push (&st, i);
				break;

			case '[':
				h_push (&st_brack, i);
				break;

			case ')':
				j = h_pop (&st);
				p_table[i] = j;
				p_table[j] = i;
				break;

			case ']':
				j = h_pop (&st_brack);
				p_table[i] = j;
				p_table[j] = i;
				break;
		}
	}

	if (st.top != STACK_EMPTY || st_brack.top != STACK_EMPTY) { //0 || st_brack.top != 0)
    	fprintf (stderr, "The given structure is not valid: %d more left parenthesis than right parentheses\n", st.top);
        exit (1);
    }
}

/*
 * The algorithm for finding the weakly closed regions for a given pseudoknot free structure is as follows
 *
 * initialization:
 * 				open = -1;
 * 				weakly_closed array initialized to 0
 *
 *
 * case 1) if the region is of the form [i,i] and i is not paired,
 * then it is considered as a weakly closed region and we put weakly_closed[ii] = 1
 * otherwise we put open = i
 *
 * case 2) if we are considering region [i,j] where i!=j and we know region [i, j-1] is weakly closed
 *		2-a) if j is not paired, then [i,j] is also weakly closed and we put weakly_closed[ij] = 1
 * 		2-b) if j is paired and we have i <= bp(j) < j then [i,j] is weakly_closed and we put weakly_closed[ij] = 1
 * 		THIS CASE CAN NEVER HAPPERN ==> REMOVE FROM CODE
 * 		2-c) otherwise [i,j] is NOT weakly closed
 * 			if (open == -1) then open = j
 *
 * case 3) if we are considering region [i,j] where i!=j and we know region [i,j-1] is NOT weakly closed
 * 		3-a) if j is paired and bp(j) == open then [i,j] is weakly closed and we put weakly_closed[ij] = 1 and open = -1
 * 		3-b) otherwise the region is still NOT weakly closed
 *
 */

void detect_weakly_closed(h_str_features *fres, int *weakly_closed, int nb_nucleotides, int *index){
	int i,j;
	int open = -1;

	for (i = 0; i < nb_nucleotides; i++){
		open = -1;
		for (j = i; j < nb_nucleotides; j++){
			int ij = index[i]+j-i; // index[i]+j-i gives the index ij
//			if ( i == 7 && j == 27){
//				printf("weakly_closed[%d,%d] = %d and open = %d \n", i,j-1,weakly_closed[ij-1], open);
//			}
			if (i == j ){
				if (fres[j].pair < 0){ // j is not paired
					weakly_closed[ij]= 1;
//					if (debug ){//_weakly_closed){
//						printf("%d==%d and %d is not paired => weakly closed\n",i,j,j);
//					}
				}else{
					open = j;
//					if (debug){//_weakly_closed){
//						printf("%d==%d and %d is paired => not weakly closed, and open = %d\n",i,j,j,open);
//					}
				}
			}else if (weakly_closed[ij-1]){
				if (fres[j].pair < 0){ // j is not paired
					weakly_closed[ij] = 1;
//					if (debug ){//_weakly_closed){
//						printf("[%d,%d] is weakly closed and %d is not paired => weakly closed\n",i,j-1,j);
//					}
				}else if (open == -1){
					open = j;
//					if (debug ){//_weakly_closed){
//						printf("[%d,%d] is weakly closed but %d is paired => not weakly closed, and open = %d\n",i,j-1,j,open);
//					}
				}
			}else if(fres[j].pair >= 0 && open == fres[j].pair){
				weakly_closed[ij] = 1;
				open = -1;
//				if (debug && i == 27 && j == 68){//_weakly_closed){
//					printf("[%d,%d] is not weakly closed but %d is paired with %d => weakly closed\n",i,j-1,j,fres[j].pair);
//				}
			}
			//if (debug){//_weakly_closed){
//				if (i == 27 && j == 68){//_weakly_closed){
//					printf("weakly_closed[%d,%d] = %d\n",i,j,weakly_closed[ij]);
//				}
//				if (weakly_closed[ij] == 1)
//					printf("Region [%d,%d] is Weakly closed\n",i,j);
			//}


		}
	}
}
/*
 * Hosna: Feb 12, 2007
 * algorithm for finding if region [i,j] is an empty region
 * We define an empty region as follows:
 * region [i,j] is an empty region if for all k i<k<j, k is unpaired
 *
 * in our program we only need to know empty regions based on the original structure
 *
 * for every base i and j, such that i <= j
 * 1) if (i == j && i is not paired in G)
 * then region [i,i] is an empty region and
 * we put not_paired_all[i,j] = 1
 *
 * 2) else if j is not paired in G and we know that region[i,j-1] is an empty region
 * then region [i,j] is an empty region and we put
 * not_paired_all[i,j] = 1
 *
 * 3) otherwise region[i,j] is not an empty region and we put
 * not_paired_all[i,j] = 0
 *
 */

void detect_not_paired_all(h_str_features *fres, int *not_paired_all, int nb_nucleotides, int *index){
	int i, j;
	for(i = 0; i < nb_nucleotides; i++){
		for(j = i; j < nb_nucleotides; j++){
			int ij = index[i]+j-i;
			if (i == j){
				if (fres[i].pair < 0){
					not_paired_all[ij]=1;
//					if (debug){//_empty){
//						printf("region [%d,%d] is empty \n",i,j);
//					}
				}else{
					not_paired_all[ij] = 0;
//					if (debug){//_empty){
//						printf("region [%d,%d] is NOT empty \n",i,j);
//					}
				}
			}else if (not_paired_all[ij-1] == 1 && fres[j].pair < 0){
				not_paired_all[ij] = 1;
//				if (debug){//_empty){
//					printf("region [%d,%d] is empty \n",i,j);
//				}
			}
			else{
				not_paired_all[ij] = 0;
//				if (debug){//_empty){
//					printf("region [%d,%d] is NOT empty \n",i,j);
//				}
			}
		}
	}
}



/* Hosna:
 * The algorithm for finding b(i,l) and B(l,j) is as follows:
 *
 *  for every base l, if l is in an arc and l is not paired (i.e. f[l].arc != -1)
 *  and every i (0 < i < nb_nucleotides)
 *
 *  1) if (i <= arc(l)) then we are finding b(i,l):
 * 		1-a) temp = arc(arc(l))
 * 		1-b) if (temp == -1 || temp <= i) then we put
 * 			b(i,l) = b'(i,l) (i.e. border_bs[l][i] = arc(l))
 * 		1-c) otherwise
 * 				1-c-i) while (arc(temp) ! = -1 && arc(temp) > i)
 * 							temp = arc(temp)
 * 				1-c-ii) b(i,l) = temp and we put border_bs[l][i] = temp
 *
 *  2) if (i >= pair(arc(l))) then we are finding B(l,j):
 * 		1-a) temp = arc(arc(l))
 * 		1-b) if (temp == -1 || pair(temp) >= i) then we put
 * 			B(l,i) = B'(l,i) (i.e. border_bs[l][i] = pair(arc(l)))
 * 		1-c) otherwise
 * 				1-c-i) while (arc(temp) ! = -1 && pair(arc(temp)) < i)
 * 							temp = arc(temp)
 * 				1-c-ii) B(l,i) = pair(temp) and we put border_bs[l][i] = pair(temp)
 *
 *  3) if (i > arc(l) && i < pair(arc(l))) then border_bs[l][i] = -1
 */

void detect_border_bs(h_str_features *fres, int** border_bs, int nb_nucleotides){

	int l,i;
	for (l = 0; l < nb_nucleotides; l++){
		for (i = 0; i < nb_nucleotides; i++){
			int cover_l = fres[l].arc, pair_l=fres[l].pair; // Hosna March 8, 2012, using local varibales for optimality
			if (cover_l == -1 || pair_l >= 0){
				border_bs[l][i] = -2;
//				if (debug_border_bs){
//					printf("fres[%d].arc = %d, fres[%d].pair = %d \n",l, fres[l].arc,l,fres[l].pair);
//					printf("border_bs[%d][%d] = %d \n",l,i,border_bs[l][i]);
//				}

			}else{
				if (i <= cover_l){
					int temp = fres[cover_l].arc;
					if (temp == -1 || temp < i){ // Hosna: Jan 31, 2007: temp <= i changed to < to include i itself too
						border_bs[l][i] = cover_l;
					}else{
						while(fres[temp].arc != -1 && fres[temp].arc >= i){ // Hosna: Jan 31, 2007: fres[temp].arc > i changed to >= to include i itself too
							temp = fres[temp].arc;
						}
						border_bs[l][i] = temp;
					}
				}
				if (i >= fres[cover_l].pair){
					int temp = fres[cover_l].arc;
					if (temp == -1 || fres[temp].pair > i){ //Hosna: Jan 31, 2007: fres[temp].pair >= i changed to > to include i itself too
						border_bs[l][i] = fres[cover_l].pair;
					}else{
						while(fres[temp].arc != -1 && fres[fres[temp].arc].pair <= i){ //Hosna: Jan 31, 2007: fres[fres[temp].arc].pair < i changed to <= to include i itself too
							temp = fres[temp].arc;
						}
						border_bs[l][i] = fres[temp].pair;
					}
				}
				if (i > cover_l && i < fres[cover_l].pair){
					border_bs[l][i] = -1;
				}
//				if (debug_border_bs){
//					printf("fres[%d].arc = %d, fres[%d].pair = %d \n",l, fres[l].arc,l,fres[l].pair);
//					printf("border_bs[%d][%d] = %d \n",l,i,border_bs[l][i]);
//				}
			}
		}
	}
}

/* Hosna:
 * The algorithm for finding b'(i,l) and B'(l,j) is as follows:
 *
 *  for every base l, if l is in an arc and l is not paired (i.e. f[l].arc != -1)
 *  and every i (0 < i < nb_nucleotides)
 *
 * 	1) if (i < arc(l)) then b'(i,l) = arc(l)
 *  and we put the corresponding value in border_bps[l][i]
 *
 *  2) if (i > pair(arc(l))) then B'(l,i) = pair(arc(l))
 *  and we put the corresponding value in border_bps[l][i]
 *
 *  3) if (i >= arc(l) && i<= pair(arc(l))) then border[l][i] = -1
 *
 */

void detect_border_bps(h_str_features *fres, int** border_bps, int nb_nucleotides){

	int l, i;
	for (l = 0; l < nb_nucleotides; l++){
		for (i = 0; i < nb_nucleotides ; i++){
			int cover_l = fres[l].arc, pair_l=fres[l].pair; // Hosna March 8, 2012, using local varibales for optimality
			if (cover_l  == -1 || pair_l >= 0){
				border_bps[l][i] = -2;//INF;//-2;
//				if (debug_border_bps){
//					printf("fres[%d].arc = %d, fres[%d].pair = %d \n",l, fres[l].arc,l,fres[l].pair);
//					printf("border_bps[%d][%d] = %d \n",l,i,border_bps[l][i]);
//				}
				//break;
			}else{
				if (i <= cover_l ){		//Hosna: Jan 31, 2007: < changed to <= to include i itself too
					border_bps[l][i] = cover_l ;
				}
				if (i >= fres[cover_l ].pair){ //Hosna: Jan 31, 2007: > changed to >= to include i itself too
					border_bps[l][i] = fres[fres[l].arc].pair;
				}
				if ( i > cover_l  && i < fres[cover_l].pair){ //Hosna: Jan 31, 2007: >= and <= changed to > and < to include i itself too
					border_bps[l][i] = -1;
				}
//				if (debug_border_bps){
//					printf("fres[%d].arc = %d, fres[%d].pair = %d \n",l, fres[l].arc,l,fres[l].pair);
//					printf("border_bps[%d][%d] = %d \n",l,i,border_bps[l][i]);
//				}
			}
		}
	}
}

void h_init (stack_ds *st)
// PRE:  None
// POST: Initialize the stack st
{
	st->top = STACK_EMPTY;//0;
}

void h_push (stack_ds *st, int el)
// PRE:  st is a valid stack
// POST: Push an element to the stack
{
	st->top = st->top +1;
    st->elem[st->top] = el;
}

int h_pop (stack_ds *st)
// PRE:  st is a valid stack, that is not empty
// POST: pop an element from the stack and return it
{
    if (st->top <= STACK_EMPTY)//0)
    {
        fprintf (stderr, "The given structure is not valid: more right parentheses than left parentheses\n");
        exit (1);
    }
    int result = st->elem[st->top];
    st->top = st->top -1 ;
    return result;
}

h_str_features *convert_str_features_to_h_str_features(str_features *f){
	h_str_features *fres;
	fres->arc = -1;
	fres->num_branches = f->num_branches;
	fres->pair = f->pair;
	fres->type = f->type;
	int i;
	for (i = 0; i < fres->num_branches; i++){
		fres->bri[i] = f->bri[i];
	}
	return fres;

}

str_features *convert_h_str_features_to_str_features(h_str_features *f){
	str_features *fres;
	fres->num_branches = f->num_branches;
	fres->pair = f->pair;
	fres->type = f->type;
	int i;
	for (i = 0; i < fres->num_branches; i++){
		fres->bri[i] = f->bri[i];
	}
	return fres;
}

void detect_h_structure_features (char *structure, h_str_features *f)
// PRE:  None
// POST: The variable f is filled with structure features, i.e. what type of elementary structure
//       this base is closing (such as stacked pair, hairpin loop etc.)
{
    int num_branches, i, j;
    int p_table[MAXSLEN];
    int arc_table[MAXSLEN];
    int bri[MAX_BRANCHES];
    int nb_nucleotides;
    nb_nucleotides = strlen(structure);
    detect_original_pairs_arcs (structure, p_table, arc_table);
    for (i=0; i < nb_nucleotides; i++)
    {
		// Hosna, March 8, 2012
		// use local variables instead of getting array values all the time
		int i_pair = p_table[i];
		f[i].pair = i_pair;//p_table[i];
        f[i].arc = arc_table[i];


        if (i_pair>i)//p_table[i] > i)
        {
            //f[i].pair = p_table[i]; //Hosna March 8, 2012, this statement seemed redundant! so removed
            f[i_pair].pair = i;//f[p_table[i]].pair = i;
            // check if it is stacked pair
			int i_pair_plus1 = p_table[i+1];
            if (i_pair_plus1 == i_pair-1 && i_pair_plus1 > i+1)//(p_table[i+1] == p_table[i]-1 && p_table[i+1] > i+1)
            {
                f[i].type = STACK;
                f[i_pair].type = STACK;//f[p_table[i]].type = STACK;
				/* if (i>=136){
						printf("END OF THE LOOP \n");
						printf("p_table[%d]=%d \n",i,p_table[i]);
						printf("f[%d].pair = %d, f[%d].type = %c, f[%d].arc = %d \n",i,f[i].pair,i,f[i].type,i,f[i].arc);

					} */
                continue;
            }
            // check if it is hairpin, internal loop or multi-loop
            num_branches = 0;
            for (j=i+1; j<i_pair; j++)//(j=i+1; j < p_table[i]; j++)
            {
                if (p_table[j] > j)
                {
                    bri[num_branches] = j;
                    num_branches++;
                    j = p_table[j];
                }
            }
            if (num_branches == 0)  // hairpin
            {
                f[i].type = HAIRP;
                f[i_pair].type = HAIRP; //f[p_table[i]].type = HAIRP;
            }
            else if (num_branches == 1) // internal loop
            {
                f[i].type = INTER;
                f[i_pair].type = INTER; //f[p_table[i]].type = INTER;
                f[i].num_branches = 1;
                f[i].bri[0] = bri[0];
            }
            else    // multi loop
            {
                f[i].type = MULTI;
                f[i_pair].type = MULTI; //f[p_table[i]].type = MULTI;
                f[i].num_branches = num_branches;
                for (j=0; j < num_branches; j++)
                    f[i].bri[j] = bri[j];
            }
        }

		/* if (i>=136){
							printf("END OF THE LOOP \n");
							printf("p_table[%d]=%d \n",i,p_table[i]);
							printf("f[%d].pair = %d, f[%d].type = %c, f[%d].arc = %d \n",i,f[i].pair,i,f[i].type,i,f[i].arc);

						} */
    }
    if (debug){
    	printf("h_str_features was successful! \n");
    }
}

/*
 * Hosna: January 10, 2008
 * The following two functions are modified versions of
 * the functions found in simfold/src/common/common.cpp
 * the modifications are to make them work for density-2 structures
 *
 */


double compute_h_sensitivity (char *ref_structure, char *pred_structure)
// returns 0 if undefined (denominator is 0)
{
    int ptable_ref[MAXSLEN];
    int ptable_pred[MAXSLEN];
    int distance;
    int len, i;
    double sens;
    int num_correct_bp;
    int num_true_bp;

    len = strlen(ref_structure);
    detect_original_PKed_pairs (ref_structure, ptable_ref);
    detect_original_PKed_pairs (pred_structure, ptable_pred);
    num_correct_bp = 0;
    num_true_bp = 0;
    for (i=0; i < len; i++)
    {
        if (ptable_ref[i] > -1)    // paired base
        {
            num_true_bp++;
            if (ptable_pred[i] == ptable_ref[i])
                num_correct_bp++;
        }
    }
    if (num_true_bp == 0)
        return -1.0;
    sens = num_correct_bp*1.0/num_true_bp;
    return sens;
}


double compute_h_ppv (char *ref_structure, char *pred_structure)
// returns 0 if undefined (denominator is 0)
{
    int ptable_ref[MAXSLEN];
    int ptable_pred[MAXSLEN];
    int distance;
    int len, i;
    double ppv;
    int num_correct_bp;
    int num_pred_bp;

    len = strlen(ref_structure);
    detect_original_PKed_pairs (ref_structure, ptable_ref);
    detect_original_PKed_pairs (pred_structure, ptable_pred);
    num_correct_bp = 0;
    num_pred_bp = 0;
    for (i=0; i < len; i++)
    {
        if (ptable_ref[i] > -1 && ptable_pred[i] == ptable_ref[i])    // paired base
            num_correct_bp++;
        if (ptable_pred[i] > -1)    // paired base
            num_pred_bp++;
    }
    if (num_pred_bp == 0)
        return -1.0;
    ppv = num_correct_bp*1.0/num_pred_bp;
    return ppv;
}

/*
 * Hosna Jan 10, 2008
 * The following function is used to tune the parameters
 * using Andronescu's GC program
 * This function is supposed to reset the pseudoknotted parameters
 */

void h_fill_data_structures_with_new_parameters (char *filename){
    FILE *file;
    char buffer[100];
    double param;
    int line = 0;

    //printf ("FILENAME: %s\n", filename);
	if ((file = fopen (filename, "r")) == NULL)
	{
	    giveup ("Cannot open file", filename);
	}

	// PS_penalty: exterior pseudoloop initiation penalty (originally 9.6 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
    param *= 100;
    PS_penalty = (int) param;

	//PSM_penalty: penalty for introducing pseudoknot inside a multiloop (originally 15 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    PSM_penalty = (int) param;

	//PSP_penalty: penalty for introducing pseudoknot inside a pseudoloop (originally 15 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    PSP_penalty = (int) param;

	//PB_penalty: band penalty (originally 0.2 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    PB_penalty = (int) param;

	//PUP_penalty: penalty for an un-paired base in a pseudoloop or a band (originally 0.1 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    PUP_penalty = (int) param;

	//PPS_penalty: penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band(originally 0.1 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    PPS_penalty = (int) param;


	//a_penalty: penalty for introducing a multiloop (originally 3.4 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    a_penalty = (int) param;

	//b_penalty: penalty for base pair in a multiloop (originally 0.4 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    b_penalty = (int) param;


	//c_penalty: penalty for un-paired base in a multi-loop (originally 0)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
    param *= 100;
    c_penalty = (int) param;



	// e_stP = 0.83 * e_s
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
    e_stP_penalty = param;


	// e_intP = 0.83 * e_int
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
    e_intP_penalty = param;


	//ap_penalty: penalty for introducing a multiloop that spans a band (originally 3.4 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    ap_penalty = (int) param;

	//bp_penalty: base pair penalty for a multiloop that spans a band (originally 0.4 Kcal/mol)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    bp_penalty = (int) param;


	//cp_penalty: penalty for unpaired base in a multiloop that spans a band (originally 0)
    fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%lf", &param);
	param *= 100;
    cp_penalty = (int) param;


	fclose (file);
	//printf ("****** we must have 14 lines by now: LINES = %d\n", line);


}


void h_fill_data_structures_with_new_parameters (double *param)
{

	if (param == NULL){
		giveup ("Incorrect parameter length", "h_fill_data_structure_with_new_parameters");
	}

	// PS_penalty: exterior pseudoloop initiation penalty (originally 9.6 Kcal/mol)
    param[0] *= 100;
    PS_penalty = (int) param[0];

	//PSM_penalty: penalty for introducing pseudoknot inside a multiloop (originally 15 Kcal/mol)
	param[1] *= 100;
    PSM_penalty = (int) param[1];

	//PSP_penalty: penalty for introducing pseudoknot inside a pseudoloop (originally 15 Kcal/mol)
	param[2] *= 100;
    PSP_penalty = (int) param[2];

	//PB_penalty: band penalty (originally 0.2 Kcal/mol)
	param[3] *= 100;
    PB_penalty = (int) param[3];

	//PUP_penalty: penalty for an un-paired base in a pseudoloop or a band (originally 0.1 Kcal/mol)
	param[4] *= 100;
    PUP_penalty = (int) param[4];

	//PPS_penalty: penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band(originally 0.1 Kcal/mol)
    param[5] *= 100;
    PPS_penalty = (int) param[5];


	//a_penalty: penalty for introducing a multiloop (originally 3.4 Kcal/mol)
	param[6] *= 100;
    a_penalty = (int) param[6];

	//b_penalty: penalty for base pair in a multiloop (originally 0.4 Kcal/mol)
    param[7] *= 100;
    b_penalty = (int) param[7];


	//c_penalty: penalty for un-paired base in a multi-loop (originally 0)
    param[8] *= 100;
    c_penalty = (int) param[8];



	// e_stP = 0.83 * e_s
    e_stP_penalty = param[9];


	// e_intP = 0.83 * e_int
    e_intP_penalty = param[10];


	//ap_penalty: penalty for introducing a multiloop that spans a band (originally 3.4 Kcal/mol)
    param[11] *= 100;
    ap_penalty = (int) param[11];

	//bp_penalty: base pair penalty for a multiloop that spans a band (originally 0.4 Kcal/mol)
    param[12] *= 100;
    bp_penalty = (int) param[12];


	//cp_penalty: penalty for unpaired base in a multiloop that spans a band (originally 0)
    param[13] *= 100;
    cp_penalty = (int) param[13];

}

int h_create_string_params(){
	int index = create_string_params();

	sprintf (string_params[index++], "PS_penalty");
	sprintf (string_params[index++], "PSM_penalty");
	sprintf (string_params[index++], "PSP_penalty");
	sprintf (string_params[index++], "PB_penalty");
	sprintf (string_params[index++], "PUP_penalty");
	sprintf (string_params[index++], "PPS_penalty");


	sprintf (string_params[index++], "a_penalty");
	sprintf (string_params[index++], "b_penalty");
	sprintf (string_params[index++], "c_penalty");


	sprintf (string_params[index++], "e_stP_penalty");
	sprintf (string_params[index++], "e_intP_penalty");


	sprintf (string_params[index++], "ap_penalty");
	sprintf (string_params[index++], "bp_penalty");
	sprintf (string_params[index++], "cp_penalty");
	return index;
}

double hfold(char *sequence, char *restricted, char *structure){
	W_final *min_fold = new W_final (sequence, restricted);
	if (min_fold == NULL) giveup ("Cannot allocate memory", "HFold");
	double energy = min_fold->hfold();
    min_fold->return_structure (structure);
    delete min_fold;
    return energy;
}

// Hosna, May 3rd, 2012
// added function for the pkonly version
double hfold_pkonly(char *sequence, char *restricted, char *structure){
	W_final *min_fold = new W_final (sequence, restricted);
	if (min_fold == NULL) giveup ("Cannot allocate memory", "HFoldPKonly");
	double energy = min_fold->hfold_pkonly();
    min_fold->return_structure (structure);
    delete min_fold;
    return energy;
}
