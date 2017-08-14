#include "pseudo_loop.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "constants.h"
#include "h_struct.h"
#include "externs.h"
#include "h_externs.h"
#include "h_common.h"
#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"
#include "VM_final.h"
#include "V_final.h"
#include "s_specific_functions.h"

// Ian Wark July 19 2017
// constant that defines what fres[i].pair will be compared against (>=) for impossible cases
// set to -1 because >= 0 means there is already a base pair there,
// and -1 means restricted struture says there is no base pair there.
#define FRES_RESTRICTED_MIN -1

pseudo_loop::pseudo_loop(char *seq, char* restricted, V_final *V, s_hairpin_loop *H, s_stacked_pair *S, s_internal_loop *VBI, VM_final *VM)
{
	this->sequence = seq;
	this->restricted = restricted;
	this->V = V;
	this->H = H;
	this->S = S;
	this->VBI = VBI;
	this->VM = VM;
    allocate_space();
    if (debug){
    	printf("an object of pseudo_loop was successfully created! \n");
    }
}

void pseudo_loop::allocate_space()
{
    int i;
    nb_nucleotides = strlen(sequence);
    needs_computation = 0; // Hosna, March 14, 2012 I need to remove this variable!! to make everything a lot faster

    index = new int [nb_nucleotides];
    int total_length = (nb_nucleotides *(nb_nucleotides+1))/2;
    index[0] = 0;
    for (int i=1; i < nb_nucleotides; i++)
        index[i] = index[i-1]+nb_nucleotides-i+1;

    WI = new int [total_length];
    if (WI == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < total_length; i++) WI[i] = 0; // if i == j -> p_up

    weakly_closed = new int[total_length];
    if (weakly_closed == NULL) giveup ("Cannot allocate memory", "weakly_closed");
    for (i=0; i < total_length; i++) weakly_closed[i] = 0;

    not_paired_all = new int[total_length];
    if (not_paired_all == NULL) giveup ("Cannot allocate memory", "not_paired_all");
    for (i=0; i < total_length; i++) not_paired_all[i] = 0;

    VP = new int[total_length];
    if (VP == NULL) giveup ("Cannot allocate memory", "VP");
    for (i=0; i < total_length; i++) VP[i] = INF;

    WMB = new int[total_length];
    if (WMB == NULL) giveup ("Cannot allocate memory", "WMB");
    for (i=0; i < total_length; i++) WMB[i] = INF;

    WMBP = new int[total_length];
	if (WMBP == NULL) giveup("Cannot allocate memory","WMBP");
	for (i=0; i < total_length; i++) WMBP[i] = INF;

    WIP = new int[total_length];
    if (WIP == NULL) giveup ("Cannot allocate memory", "WIP");
    for (i=0; i < total_length; i++) WIP[i] = INF;


    VPP = new int[total_length];
    if (VPP == NULL) giveup ("Cannot allocate memory", "VPP");
    for (i=0; i < total_length; i++) VPP[i] = INF;

    BE = new int[total_length];
    if (BE == NULL) giveup ("Cannot allocate memory", "BE");
    for (i=0; i < total_length; i++) BE[i] = 0; //check

    border_bs = new int*[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) border_bs[i] = new int[nb_nucleotides];

    border_bps = new int*[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) border_bps[i] = new int[nb_nucleotides];


    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);

}

pseudo_loop::~pseudo_loop()
{
    delete [] WI;
    delete [] WIP;
    delete [] VP;
    delete [] VPP;
    delete [] WMB;
    delete [] WMBP;

    delete [] BE;
    delete [] weakly_closed;
    delete [] not_paired_all;

    for(int i = 0; i < nb_nucleotides; i++)
    {
        delete [] border_bs[i];
        delete []  border_bps[i];
    }
    delete [] border_bs;
    delete [] border_bps;


    delete [] index;
    delete [] int_sequence;
}

void pseudo_loop::set_features(h_str_features *f){
	fres = f;
}

int pseudo_loop::is_weakly_closed(int i, int j){
	// base case: if i > j then the region is weakly closed
	if (i>j){
		return 1;
	}
	int ij = index[i]+j-i;
	if (weakly_closed[ij] == 1)
		return 1;
	return 0;
}


int pseudo_loop::is_empty_region(int i, int j){
	//base case: if i> j then the region is empty
	if (i>j){
		return 1;
	}
	int ij = index[i]+j-i;
	if (not_paired_all[ij] == 1){
		return 1;
	}
	return 0;
}

void pseudo_loop::initialize(){

	int i, j;

    //Hosna: before going further, we should fill up the weakly closed array
    detect_weakly_closed(fres, weakly_closed, nb_nucleotides, index);
    detect_not_paired_all(fres, not_paired_all, nb_nucleotides, index);
    detect_border_bs(fres,border_bs, nb_nucleotides);
    detect_border_bps(fres,border_bps, nb_nucleotides);

    //testing to see if all the structures were found correctly
//    if (debug){
//	    for (i=0; i < nb_nucleotides; i++){
	        //if (fres[i].pair != -1 && fres[i].pair != -2)
//	        if(debug)
//	        printf ("%d pairs %d, is in arc %d, and has type %c\n", i, fres[i].pair, fres[i].arc, fres[i].type);

//	    }
//	    for (i = 0; i < nb_nucleotides; i++){
//	    	int j;
//	    	for (j = i; j < nb_nucleotides; j++){
//	    		int ij = index[i]+j -i;
//	    		if (weakly_closed[ij] == 1){
//		    		printf("region [%d,%d] is weakly closed. \n", i, j);
//	    		}
//	    		else{
//	    			printf("region [%d,%d] is NOT weakly closed. \n", i, j);
//	    		}
//				if (not_paired_all[ij] == 1){::compute_WI
//					printf("region [%d,%d] is EMPTY \n",i,j);
//				}
				//checking WI, VP and WMB to see if they have been initialized correctly
//				if (debug)
//				if(WI[ij] != 0){
//					printf("WI[%d,%d] NOT initialized correctly!\n",i,j);
//				}
//				if(VP[ij] != 0){
//					printf("VP[%d,%d] NOT initialized correctly!\n",i,j);
//				}
//				if(WMB[ij] != 0){
//					printf("WMB[%d,%d] NOT initialized correctly!\n",i,j);
//				}
//	    	}
//	    }

//    }
//    printf("WMB was initialized successfully! \n");

}

void pseudo_loop::compute_energies(int i, int j)
{

	// Hosna, April 18th, 2007
	// based on discussion with Anne, we changed WMB to case 2 and WMBP(containing the rest of the recurrences)

	//	if(debug){
	//		printf("calculating VP(%d,%d) \n",i,j);
	//	}
	compute_VP(i,j,fres); // Hosna, March 14, 2012, changed the positionof computing VP from after BE to befor WMBP


//	if(debug){
//		printf("calculating WMBP(%d,%d) \n",i,j);
//	}::compute_WI
	compute_WMBP(i,j,fres);
//	if(debug){
//		printf("calculating WMB(%d,%d) \n",i,j);
//	}
    compute_WMB(i,j,fres);
//	if(debug){
//		printf("calculating WI(%d,%d) \n",i,j);
//	}
    compute_WI(i,j,fres);
//    if(debug){
//		printf("calculating WIP(%d,%d) \n",i,j);
//	}
    compute_WIP(i,j,fres);
//    if(debug){
//		printf("calculating VPP(%d,%d) \n",i,j);
//	}
    compute_VPP(i,j,fres);
//	if(debug){
//		printf("calculating BE(%d,%d) \n",i,j);
//	}

	compute_BE(fres[j].pair,j,fres[i].pair,i,fres);

//    if (debug){
//    	printf("WI(%d,%d) = %d \n", i,j, get_WI(i,j));
//    	printf("WIP(%d,%d) = %d \n", i,j, get_WIP(i,j));
//    	printf("VPP(%d,%d) = %d \n", i,j, get_VPP(i,j));
//    	printf("BE(%d,%d,%d,%d) = %d \n", i,fres[i].pair,j,fres[j].pair, get_BE(i,fres[i].pair,j,fres[j].pair));
//    	printf("VP(%d,%d) = %d \n", i,j, get_VP(i,j));
//    	printf("WMBP(%d,%d) = %d \n", i,j, get_WMBP(i,j));
//    	printf("WMB(%d,%d) = %d \n", i,j, get_WMB(i,j));
//    }
}

void pseudo_loop::compute_WI(int i, int j , h_str_features *fres){
	int min = INF, m1 = INF, m2= INF, m3= INF;
	int ij = index[i]+j-i;
	if (WI[ij] != 0){ //calculated before
//		if (debug)
//		{
//			printf("WI(%d,%d) was calculated before ==> WI(%d,%d) = %d \n",i,j,i,j,get_WI(i,j));
//		}
		return;
	}

	//base cases
	// if [i,j] is not weakly closed then WI[i,j] = INF
	if (is_weakly_closed(i,j) == 0){
		WI[ij] = INF;
//		if (debug)
//		{
//			printf("[%d,%d] is not weakly closed ==> WI(%d,%d) = %d \n",i,j,i,j,get_WI(i,j));
//		}
		return;
	}

	// branch 4, one base
	if (i == j){
		WI[ij] = PUP_penalty;
//		if (debug){
//			printf("i ==j => WI[%d,%d]= %d \n",i,j,WI[ij]);
//		}
		return;
	}
	// Hosna: Feb 12, 2007
	// changed this part to see if it works better

	// Hosna: Feb 16, 2007:
	// we don't need to check to see if i and j are inside an arc
	// because they are not in an arc in G but they will be in an arc in G'
	if (fres[i].arc != fres[j].arc){
		WI[ij] = INF;
//		if (debug){
//			printf("i and j not in the same arc => WI[%d,%d]= %d \n",i,j,WI[ij]);
//		}
		return;
	}

// Hosna: July 2nd, 2007
// in branch 1 of WI, we can have a case like
// ((..))((...))
// such that both i and j are paired but we can chop them

	// branch 1:
//	if (fres[i].pair < 0 && fres[j].pair < 0)
//	{
	int t;
	for (t = i; t< j; t++){
		int wi_1 = get_WI(i,t);
		int wi_2 = get_WI(t+1,j);
		int energy = wi_1 + wi_2;
		m1 = (m1 > energy)? energy : m1;
//		if (debug_WI){
//			printf("WI branch 1: WI[%d,%d] = %d and WI[%d,%d] = %d => energy = %d and m1 = %d \n",i,t,wi_1,(t+1),j,wi_2,energy, m1);
//		}
	}
//	if (debug){
//		printf("WI(%d,%d) branch 1: m1 = %d\n",i,j,m1);
//	}
//	}
	// branch 2:

	// Ian Wark July 19 2017
	// fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if ((fres[i].pair == j && fres[j].pair == i)
	||(fres[i].pair < FRES_RESTRICTED_MIN && fres[j].pair < FRES_RESTRICTED_MIN)){
		// Hosna, April 16th, 2007
		// changed back to see if it will work fine
		// Hosna: April 19th, 2007
		// I think we should call the restricted version

		int v_ener = (i>j)? INF: V->get_energy(i,j);
		m2 = v_ener + PPS_penalty;
//		if (debug){
//			printf("WI(%d,%d) branch 2: m2 = %d\n",i,j,m2);
//		}
	}
	// branch 3:
	// Hosna: April 20, 2007
	// removed the penalty of PPS

	// Hosna: July 5th, 2007
	// Anne said we should put PPS back
	// change PSM to PSP
//	if (debug && i == 6 && j == 15){
//		printf("WI(6,15) is calling WMB \n");
//	}

	m3 = get_WMB(i,j) + PSP_penalty + PPS_penalty;

//	if (debug){
//		printf("WI(%d,%d) branch 3: m3 = %d\n",i,j,m3);
//	}

	min = MIN(m1,MIN(m2,m3));
	WI[ij] = min;
//	if (debug ){
//		printf("WI[%d,%d]: m1 = %d, m2 = %d and m3 = %d ==> min = %d \n",i,j,m1,m2,m3,WI[ij]);
//	}
}


void pseudo_loop::compute_WI_pkonly(int i, int j , h_str_features *fres){
	int min = INF, m1 = INF, m2= INF, m3= INF;
	int ij = index[i]+j-i;
	if (WI[ij] != 0){ //calculated before
		return;
	}

	//base cases
	// if [i,j] is not weakly closed then WI[i,j] = INF
	if (is_weakly_closed(i,j) == 0){
		WI[ij] = INF;
		return;
	}

	// branch 4, one base
	if (i == j){
		WI[ij] = PUP_penalty;
		return;
	}
	if (fres[i].arc != fres[j].arc){
		WI[ij] = INF;
		return;
	}


	// branch 1:
	int t;
	for (t = i; t< j; t++){
		int wi_1 = get_WI(i,t);
		int wi_2 = get_WI(t+1,j);
		int energy = wi_1 + wi_2;
		m1 = (m1 > energy)? energy : m1;
	}
	// branch 2:
	if (fres[i].pair == j && fres[j].pair == i){

		int v_ener = (i>j)? INF: V->get_energy(i,j);
		m2 = v_ener + PPS_penalty;

	}
	// branch 3:
	m3 = get_WMB(i,j) + PSP_penalty + PPS_penalty;


	min = MIN(m1,MIN(m2,m3));
	WI[ij] = min;
}


void pseudo_loop::compute_VP(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (VP[ij] != INF){//has been calculated before
//		if (debug)
//		{
//			printf("VP(%d,%d) was calculated before ==> VP(%d,%d) = %d \n",i,j,i,j,VP[ij]);
//		}
		return;
	}
	// base cases:
	// a) i == j => VP[ij] = INF
	// b) [i,j] is a weakly_closed region => VP[ij] = INF
	// c) i or j is paired in original structure => VP[ij] = INF

	// Ian Wark July 19 2017
	// fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says it should not be modified
	if (i == j || weakly_closed[ij] == 1 || fres[i].pair >= FRES_RESTRICTED_MIN || fres[j].pair >= FRES_RESTRICTED_MIN || can_pair(int_sequence[i],int_sequence[j]) != 1)	{
		VP[ij] = INF;
//		if (debug){
//			printf("VP[%d,%d] = %d and can_pair(%d,%d) = %d\n", i,j, VP[ij],int_sequence[i],int_sequence[j],can_pair(int_sequence[i],int_sequence[j]));
//		}
		return;
	}
	else{
		int m1 = INF, m2 = INF, m3 = INF, m4= INF, m5 = INF, m6 = INF, m7 = INF, m8 = INF; //different branches
		//branchs:
		// 1) inArc(i) and NOT_inArc(j)
		// WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// need to check the borders as they may be negative
		if(fres[i].arc > -1 && fres[j].arc == -1 && get_Bp(i,j) >= 0 && get_Bp(i,j)< nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides){
			int Bp_i = get_Bp(i,j);
			int B_i = get_B(i,j);
			int WI_ipus1_BPminus = get_WI(i+1,Bp_i - 1) ;
			int WI_Bplus_jminus = get_WI(B_i + 1,j-1);
			m1 =   WI_ipus1_BPminus + WI_Bplus_jminus;
//			if(debug){
//				printf("VP[%d,%d] branch 1: WI(%d+1)(BP(%d)-1) = %d and WI(B(%d)+1)(%d-1) = %d => m1 = %d \n",i,j,i,i,WI_ipus1_BPminus,i,j,WI_Bplus_jminus, m1);
//			}
		}

		// 2) NOT_inArc(i) and inArc(j)
		// WI(i+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// checking the borders as they may be negative
		if (fres[i].arc == -1 && fres[j].arc > -1 && get_b(i,j)>= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j) >= 0 && get_bp(i,j) < nb_nucleotides){
			int b_i = get_b(i,j);
			int bp_i = get_bp(i,j);
			int WI_i_plus_b_minus = get_WI(i+1,b_i - 1);
			int WI_bp_plus_j_minus = get_WI(bp_i + 1,j-1);
			m2 = WI_i_plus_b_minus + WI_bp_plus_j_minus;
//			if(debug){
//				printf("VP[%d,%d] branch 2: WI(%d+1)(b(%d)-1) = %d and WI(bp(%d)+1)(%d-1) = %d => m2 = %d \n",i,j,i,i,WI_i_plus_b_minus,i,j,WI_bp_plus_j_minus, m2);
//			}
		}

		// 3) inArc(i) and inArc(j)
		// WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// checking the borders as they may be negative
		if(fres[i].arc > -1 && fres[j].arc > -1 && get_Bp(i,j) >= 0 && get_Bp(i,j) < nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides && get_b(i,j) >= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j)>= 0 && get_bp(i,j) < nb_nucleotides){
			int Bp_i = get_Bp(i,j);
			int B_i = get_B(i,j);
			int b_i = get_b(i,j);
			int bp_i = get_bp(i,j);
			int WI_i_plus_Bp_minus = get_WI(i+1,Bp_i - 1);
			int WI_B_plus_b_minus = get_WI(B_i + 1,b_i - 1);
			int WI_bp_plus_j_minus = get_WI(bp_i +1,j - 1);
			m3 = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
//			if(debug){
//				printf("VP[%d,%d] branch 3: WI(%d+1)(B'(%d)-1) = %d, WI(B(%d)+1)(b(%d)-1) = %d and WI(b'(%d)+1)(%d-1) = %d => m3 = %d \n",i,j,i,i, WI_i_plus_Bp_minus,i,i,WI_B_plus_b_minus,i,j,WI_bp_plus_j_minus, m3);
//			}
		}

        // Ian Wark July 19 2017
        // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs


		// 4) NOT_paired(i+1) and NOT_paired(j-1) and they can pair together
		// e_stP(i,i+1,j-1,j) + VP(i+1)(j-1)
		if(fres[i+1].pair < FRES_RESTRICTED_MIN && fres[j-1].pair < FRES_RESTRICTED_MIN && can_pair(int_sequence[i+1],int_sequence[j-1])){
			m4 = get_e_stP(i,j)+ get_VP(i+1,j-1);
//			if (debug){
//				printf("VP[%d,%d] branch 4: S->get_energy(%d,%d) = %d and VP[%d,%d] = %d  so m4 = %d\n", i,j,i,j, S->get_energy(i,j,int_sequence), i+1, j-1, get_VP(i+1,j-1), m4);
//			}
		}

		// 5) NOT_paired(r) and NOT_paired(rp)
		//  VP(i,j) = e_intP(i,ip,jp,j) + VP(ip,jp)
		int ip, jp;
		int max_borders;
		// Hosna, April 6th, 2007
		// whenever we use get_borders we have to check for the correct values
		int min_borders = 0; // what if both are negative
		if (get_Bp(i,j)> 0 && get_Bp(i,j) < nb_nucleotides && get_b(i,j) >0 && get_b(i,j) < nb_nucleotides){
			min_borders = MIN(get_Bp(i,j),get_b(i,j));
		}else if (get_b(i,j) > 0 && get_b(i,j) < nb_nucleotides && (get_Bp(i,j) < 0 || get_Bp(i,j) > nb_nucleotides)){
			min_borders = get_b(i,j);
		}else if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && (get_b(i,j) < 0 || get_b(i,j) > nb_nucleotides)){
			min_borders = get_Bp(i,j);
		}
//		printf("B'(%d,%d) = %d, b(%d,%d) = %d, min_borders = %d\n",i,j,get_Bp(i,j),i,j,get_b(i,j), min_borders);
		for (ip = i+1; ip < min_borders; ip++){
			// Hosna: April 20, 2007
			// i and ip and j and jp should be in the same arc
			// also it should be the case that [i+1,ip-1] && [jp+1,j-1] are empty regions

			// Ian Wark July 19 2017
            // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
            // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
			if (fres[ip].pair < FRES_RESTRICTED_MIN && (fres[i].arc == fres[ip].arc) && is_empty_region(i+1,ip-1) == 1){
				// Hosna, April 6th, 2007
				// whenever we use get_borders we have to check for the correct values
				max_borders= 0;
				if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides){
					max_borders = MAX(get_bp(i,j),get_B(i,j));
				}else if (get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides && (get_bp(i,j) < 0 || get_bp(i,j) > nb_nucleotides)){
					max_borders = get_B(i,j);
				}else if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && (get_B(i,j) < 0 || get_B(i,j) > nb_nucleotides)){
					max_borders = get_bp(i,j);
				}
//				printf("b'(%d,%d) = %d, B(%d,%d) = %d, max_borders = %d\n",i,j,get_bp(i,j),i,j,get_B(i,j), max_borders);
				for (jp = max_borders+1; jp < j ; jp++){
                    // Ian Wark July 19 2017
                    // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
                    // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
					if (fres[jp].pair < FRES_RESTRICTED_MIN && can_pair(int_sequence[ip],int_sequence[jp]) && is_empty_region(jp+1,j-1) == 1){
						// Hosna: April 20, 2007
						// i and ip and j and jp should be in the same arc
						if (fres[j].arc == fres[jp].arc ){
							int temp = get_e_intP(i,ip,jp,j) + get_VP(ip,jp);
//							if (debug){
//								printf("VP(%d,%d) branch 5: e_intP(%d,%d,%d,%d) = %d, VP(%d,%d) = %d, temp = %d \n",i,j,i,ip,jp,j,get_e_intP(i,ip,jp,j),ip,jp,get_VP(ip,jp),temp);
//							}
//							printf("m5 = %d \n", m5);
							if (m5 > temp){
								m5 = temp;
							}
						}
					}
				}
			}
		}

		// 6) VP(i,j) = WIP(i+1,r-1) + VPP(r,j-1)
		int r;
		// Hosna April 9th, 2007
		// checking the borders as they may be negative numbers
		int min_Bp_j = j;
		if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_Bp(i,j) < min_Bp_j){
			min_Bp_j = get_Bp(i,j);
		}
		for (r = i+1; r < min_Bp_j ; r++){
            // Ian Wark July 19 2017
            // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
            // otherwise it will create pairs in spots where the restricted structure says there should be no pairs

			if (fres[r].pair < FRES_RESTRICTED_MIN){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				int tmp = get_WIP(i+1,r-1) + get_VPP(r,j-1) + ap_penalty + 2*bp_penalty;
//				if (debug)
//				{
//					printf("VP(%d,%d) branch 6: WIP(%d,%d) = %d, VPP(%d,%d) = %d ==> tmp = %d and m6 = %d \n",i,j,i+1,r-1,get_WIP(i+1,r-1),r,j-1,get_VPP(r,j-1),tmp,m6);
//				}
				if (tmp < m6){
					m6 = tmp;
				}
			}
		}


		// 7) VP(i,j) = VPP(i+1,r) + WIP(r+1,j-1)
		// Hosna April 9th, 2007
		// checking the borders as they may be negative numbers
		int max_i_bp = i;
		if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
			max_i_bp = get_bp(i,j);
		}
		for (r = max_i_bp+1; r < j ; r++){
            // Ian Wark July 19 2017
            // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
            // otherwise it will create pairs in spots where the restricted structure says there should be no pairs

			if (fres[r].pair < FRES_RESTRICTED_MIN){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				int tmp = get_VPP(i+1,r) + get_WIP(r+1,j-1)+ ap_penalty + 2* bp_penalty;
//				if (debug)
//				{
//					printf("VP(%d,%d) branch 7: VPP(%d,%d) = %d, WIP(%d,%d) = %d ==> tmp = %d and m7 = %d \n",i,j,i+1,r,get_VPP(i+1,r),r+1,j-1,get_WIP(r+1,j-1),tmp,m7);
//				}
				if (tmp < m7){
					m7 = tmp;
				}
			}
		}


// Hosna: June 29, 2007
// As case 6 and 7 handle what we wanted from case 8, then I think we don't need this any more

//		// Hosna: April 20, 2007
//		// based on the discussion with Anne, we decided we need another case for VP
//		// which is like case 5 with the difference that it can have some structure
//		// between [i+1,ip-1] and [jp+1,j-1]
//
//		// 8) VP(i,j) = WI(i+1,ip-1) + VP(ip,jp) + WI(jp+1,j-1) + P_ps + ap_penalty
//		// Hosna, April 6th, 2007
//		// whenever we use get_borders we have to check for the correct values
//		if (get_Bp(i,j)> 0 && get_Bp(i,j) < nb_nucleotides && get_b(i,j) >0 && get_b(i,j) < nb_nucleotides){
//			min_borders = MIN(get_Bp(i,j),get_b(i,j));
//		}else if (get_b(i,j) > 0 && get_b(i,j) < nb_nucleotides && (get_Bp(i,j) < 0 || get_Bp(i,j) > nb_nucleotides)){
//			min_borders = get_b(i,j);
//		}else if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && (get_b(i,j) < 0 || get_b(i,j) > nb_nucleotides)){
//			min_borders = get_Bp(i,j);
//		}
//		for (ip = i+1; ip < min_borders; ip++){
//			// Hosna: April 20, 2007
//			// i and ip and j and jp should be in the same arc
//			if (fres[ip].pair < 0 && (fres[i].arc == fres[ip].arc)){
//				// Hosna, April 6th, 2007
//				// whenever we use get_borders we have to check for the correct values
//				int max_borders= 0;
//				if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides){
//					max_borders = MAX(get_bp(i,j),get_B(i,j));
//				}else if (get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides && (get_bp(i,j) < 0 || get_bp(i,j) > nb_nucleotides)){
//					max_borders = get_B(i,j);
//				}else if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && (get_B(i,j) < 0 || get_B(i,j) > nb_nucleotides)){
//					max_borders = get_bp(i,j);
//				}
//				for (jp = max_borders+1; jp < j ; jp++){
//					if (fres[jp].pair < 0 && can_pair(int_sequence[ip],int_sequence[jp])){
//						// Hosna: April 20, 2007
//						// i and ip and j and jp should be in the same arc
//						if (fres[j].arc == fres[jp].arc){
//							int temp = get_WI(i+1,ip-1) + get_VP(ip,jp) + get_WI(jp+1,j-1) + PPS_penalty + ap_penalty;
//							if (m8 > temp){
//								m8 = temp;
//							}
//						}
//					}
//				}
//			}
//		}

		//finding the min energy
		int min = MIN(MIN(m1,m8),MIN(m2,m3));
		min = (min > MIN(m4,m5))? MIN(m4,m5) : min;
		min = (min > MIN(m6,m7))? MIN(m6,m7) : min;
//		if (debug ){
//			printf("VP[%d,%d]: m1 = %d, m2 = %d, m3 = %d, m4 = %d, m5 = %d, m6 = %d, m7 = %d and min = %d \n",i,j,m1,m2,m3,m4,m5,m6,m7,min);
//		}
		VP[ij] = min;

	}
}

void pseudo_loop::compute_WMBP(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (WMBP[ij] != INF){
		return;
	}
	//base case
	if (i == j){
		WMBP[ij] = INF;
		return;
	}
	// Hosna: July 6th, 2007
	// added impossible cases

	// Ian Wark July 19 2017
	// fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if ((fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair > j)
	||  (fres[j].pair >= FRES_RESTRICTED_MIN && fres[j].pair < i)
	||  (fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair < i )
	||  (fres[j].pair >= FRES_RESTRICTED_MIN && j < fres[j].pair)){
		WMB[ij] = INF;
		return;
	}
	else{
		int m1 = INF, m3 = INF, m4 = INF, m5 = INF;
		// if not paired(j) and paired(i) then
		// WMBP(i,j) = 2*Pb + min_{i<l<bp(i)}(BE(i,bp(i),b'(i,l),bp(b'(i,l)))+WI(b'+1,l-1)+VP(l,j))
		if(fres[j].pair < 0 && fres[i].pair >= 0){
			int tmp = INF, l, l_min=-1;
			// Hosna: June 29, 2007
			// if j is inside i's arc then the l should be
			// less than j not bp(i)
			// check with Anne
//			for (l = i+1; l < MIN(fres[i].pair,j); l++){
			// Hosna: July 5th, 2007:
			// if we have bp(i)> j then we should not have come to the WMBP
			for (l = i+1; l < j; l++){
				// Hosna, March 14, 2007
				// fixing the for loop

				// Hosna, April 9th, 2007
				// checking the borders as they may be negative
//				if(fres[l].pair < 0 && get_bp(i,l) >= 0 && get_bp(i,l) < nb_nucleotides && l+TURN <= j){
				// Hosna: July 5th, 2007:
				// removed bp(l)<0 as VP should handle that
				if(get_bp(i,l) >= 0 && get_bp(i,l) < nb_nucleotides && l+TURN <= j){
					int bp_i_l = get_bp(i,l);
					int BE_energy = get_BE(i,fres[i].pair,bp_i_l,fres[bp_i_l].pair);
					int WI_energy = get_WI(bp_i_l +1,l-1);
					int VP_energy = get_VP(l,j);
					int sum = BE_energy + WI_energy + VP_energy;
					if (tmp > sum){
						tmp = sum;
						l_min = l;
					}
//					if (debug  && i == 6 && bp_i_l == 11 ){
//						printf("***************\n");
//						printf("WMBP(%d,%d) inside branch 1: %d is paired with %d  and b'(%d,%d) = %d  and bp(b')= %d \n",i,j, i,fres[i].pair, i,l,bp_i_l, fres[bp_i_l].pair);
//						printf("l = %d, BE(%d,%d,%d,%d) = %d, WI = %d, VP = %d, --> sum = %d \n",l,i,fres[i].pair,bp_i_l,fres[bp_i_l].pair,BE_energy,WI_energy,VP_energy, sum);
//						printf("***************\n");
//					}
				}
			}
			m1 = 2*PB_penalty + tmp;
//			if (debug ){
//				printf("WMBP(%d,%d) branch 1:  l = %d  ==> m1 = %d \n",i,j,l_min, m1);
//			}

		}

		// 3)
		if (fres[j].pair < 0){
			int l, temp = INF, l_min=-1;
			for (l = i+1; l<j ; l++)	{
				// Hosna, April 6th, 2007
				// whenever we use get_borders we have to check for the correct values

				if (fres[l].arc > -1 && get_B(l,j) >= 0 && get_B(l,j) < nb_nucleotides && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
					// Hosna: April 19th, 2007
					// the chosen l should be less than border_b(i,j)
					if (get_b(i,j) >= 0 && get_b(i,j) < nb_nucleotides && l < get_b(i,j)){

						// Hosna: June 29 2007
						// after going over the program with Cristina, we noticed that
						// l should be < B'(i,j)
	//					if (l < get_Bp(i,j) && l+TURN <= j){

						// Hosna: July 5th, 2007:
						// as long as we have i <= arc(l)< j we are fine
						if (i <= fres[l].arc && fres[l].arc < j && l+TURN <=j){
							int sum = get_BE(fres[get_B(l,j)].pair,get_B(l,j),fres[get_Bp(l,j)].pair,get_Bp(l,j))+ get_WMBP(i,l-1)+ get_VP(l,j);
							if (temp > sum){
								temp = sum;
								l_min = l;
							}
	//						if (debug && fres[get_B(l,j)].pair == 6 && fres[get_Bp(l,j)].pair == 11){
	//							printf("***************\n");
	//							printf("WMBP(%d,%d) inside branch 3: l = %d, BE = %d, WMBP = %d, VP = %d ~~> sum = %d \n",i,j,l,get_BE(fres[get_B(l,j)].pair,get_B(l,j),fres[get_Bp(l,j)].pair,get_Bp(l,j)),get_WMBP(i,l-1),get_VP(l,j), sum);
	//							printf("***************\n");
	//						}
						}
					}
				}
				// Hosna: April 5th
				// after going over the WMB recurrence with Anne, we think we should add another p_b penalty
				// to the 3rd case ==> 2*P_b
				m3 = 2*PB_penalty + temp;
	//			if (debug){
	//				printf("WMBP(%d,%d) branch 3:  l = %d So m3 = %d\n",i,j,l_min, m3);
	//			}
			}
		}

		// 4) WMB(i,j) = VP(i,j) + P_b
		int temp = get_VP(i,j) + PB_penalty;
		if (temp < m4){
			m4 = temp;
		}
//		if (debug){
//			printf("WMBP(%d,%d) branch 4:  VP = %d So m4 = %d\n",i,j,get_VP(i,j), m4);
//		}
		// 5) WMB(i,j) = min_{i<l<j}{WMB(i,l)+WI(l+1,j)} if bp(j)<j
		// Hosna: Feb 5, 2007
		if(fres[j].pair < j){
			int l,l_min =-1;
			for(l = i+1; l<j; l++){
				// Hosna: March 14th, 2007
				// I think l cannot be paired

				// Hosna: April 18th, 2007
				// l and j should be in the same arc
				if (fres[l].pair < 0 && fres[l].arc > -1 && fres[j].arc > -1 && fres[j].arc == fres[l].arc){
					int temp = get_WMBP(i,l) + get_WI(l+1,j);
					if (temp < m5){
						m5 = temp;
						l_min = l;
					}

//					if (debug_WMB){
//						printf("***************\n");
//						printf("WMB(%d,%d) inside branch 5: l = %d, WMB = %d, WI = %d \n",i,j,l,get_WMB(i,l),get_WI(l+1,j));
//						printf("***************\n");
//					}
				}
			}
//			if (debug ){
//				printf("WMB(%d,%d) branch 5:  l = %d So m5 = %d\n",i,j,l_min, m5);
//			}
		}

		// get the min for WMB
		WMBP[ij] = MIN(MIN(m1,m3),MIN(m4,m5));
//		if (debug && i == 1 && j == 87){
//			printf("m1 = %d, m3 = %d, m4 = %d and m5 = %d ==> WMBP[%d,%d] = %d\n",m1,m3,m4,m5,i,j,WMBP[ij]);
//		}
	}


}

void pseudo_loop::compute_WMB(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (WMB[ij] != INF){
		return;
	}
	//base case
	if (i == j){
		WMB[ij] = INF;
		return;
	}
	// Hosna: July 6th, 2007
	// added impossible cases

    // Ian Wark July 19 2017
	// fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if ((fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair > j)
	 || (fres[j].pair >= FRES_RESTRICTED_MIN && fres[j].pair < i)
	 || (fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair < i )
	 || (fres[j].pair >= FRES_RESTRICTED_MIN && j < fres[j].pair)){
		WMB[ij] = INF;
		return;
	}
	else{
		int m2 = INF, mWMBP = INF;
		// 2)
		if (fres[j].pair >= 0 && j > fres[j].pair){
			int l, l_min=-1;
			int bp_j = fres[j].pair;
			int temp = INF;
//			if (debug_WMB){
//				printf("\n INSIDE WMB BRANCH 2 \n where bp_j = %d and j = %d \n\n", bp_j,j);
//			}
			for (l = (bp_j +1); (l < j); l++){
				// Hosna: April 24, 2007
				// correct case 2 such that a multi-pseudoknotted
				// loop would not be treated as case 2

				// Hosna: July 5th, 2007
				// this restriction was removed as it is not needed here
//				if (l > fres[i].pair){
	//				printf("l = %d and bp_l = %d \n",l,fres[l].pair);
					// Hosna April 9th,
					// checking the borders as they may be negative numbers

					// Hosna: July 5th, 2007:
					// we don't need to check that l is unpaired here
//					if (fres[l].pair < 0 && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
					if (get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
						int sum = get_BE(bp_j,j,fres[get_Bp(l,j)].pair,get_Bp(l,j)) + get_WMBP(i,l) + get_WI(l+1,get_Bp(l,j)-1);
						if (l == 600 && i == 522 && j == 615) {
							int t = 0;
						}
						if (temp > sum){
							temp = sum;
							l_min = l;
						}
//						if (debug && bp_j == 6 && fres[get_Bp(l,j)].pair == 11){
//							printf("***************\n");
//							printf("WMB(%d,%d) inside branch 2: l = %d, BE(%d,%d) = %d, WMBP = %d, WI = %d ==> sum = %d while temp = %d\n",i,j,l,bp_j,j,get_BE(bp_j,j,fres[get_Bp(l,j)].pair,get_Bp(l,j)),get_WMBP(i,l),get_WI(l+1,get_Bp(l,j)-1), sum, temp);
//							printf("***************\n");
//						}
					}
//				}

			}
			m2 = PB_penalty + temp;
//			if (debug){
//				printf("WMB(%d,%d) branch 2:  l = %d So m2 = %d\n",i,j,l_min, m2);
//			}
		}
		// check the WMBP value
		mWMBP =  get_WMBP(i,j);

		// get the min for WMB
		WMB[ij] = MIN(m2,mWMBP);
//		if (debug && i == 1 && j == 87){
//			printf("m2 = %d, mWMBP = %d ==> WMB[%d,%d] = %d\n",m2,mWMBP,i,j,WMB[ij]);
//		}
	}
}

void pseudo_loop::compute_WIP(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (WIP[ij] < INF/2){ // was calculated before
		return;
	}
//	if (debug){
//		if (i >= 27 && i <= 28 && j <= 68 && j >= 67){
//			printf("\n ************************* \n");
//			printf("Computing WIP(%d,%d) when arc(%d) = %d and arc(%d) = %d and weakly_closed(%d,%d) = %d \n",i,j,i,fres[i].arc,j,fres[j].arc,i,j,weakly_closed[ij]);
//		}
//	}
	if (fres[i].arc != fres[j].arc || i == j || weakly_closed[ij]== 0){
		WIP[ij] = INF;
		return;
	}
	int m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;

    // Ian Wark July 19 2017
	// fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs

	// branch 1:
	if (fres[i].pair < FRES_RESTRICTED_MIN){
		m1 = get_WIP(i+1,j) + cp_penalty;
//		if (debug && (i == 27 || i == 28) && (j == 68 || j == 67)){
//			printf("\n ************************* \n");
//			printf("Computing WIP(%d,%d) when  WIP(%d,%d) = %d and m1 = %d\n",i,j,i+1,j,get_WIP(i+1,j),m1);
//		}
	}
	// branch 2:
	if (fres[j].pair < FRES_RESTRICTED_MIN){
		m2 = get_WIP(i,j-1) + cp_penalty;
//		if (debug && (i == 27 || i == 28) && (j == 68 || j == 67)){
//			printf("\n ************************* \n");
//			printf("Computing WIP(%d,%d) when  WIP(%d,%d) = %d and m2 = %d\n",i,j,i,j-1,get_WIP(i,j-1),m2);
//		}
	}
	//branch 3:
	int t;
	for (t = i; t <j; t++){
		int tmp = get_WIP(i,t) + get_WIP(t+1,j);
		if (tmp < m3){
			m3 = tmp;
		}
	}

	// branch 4:
	if (fres[i].pair == j
	|| (fres[i].pair < FRES_RESTRICTED_MIN && fres[j].pair < FRES_RESTRICTED_MIN && can_pair(int_sequence[i],int_sequence[j]))){

		m4 = V->get_energy(i,j)	+ bp_penalty;
//		if (debug && i ==15 && j == 20 ){
//			printf("*************************************\n");
//			printf("WIP(%d,%d) branch 4: V(%d,%d) = %d => m4 = %d \n",i,j,i,j,V->get_energy(i,j),m4);
//			printf("*************************************\n");
//		}

	}

	// branch 5:
	m5 = get_WMB(i,j) + PSM_penalty + bp_penalty;
//	if (debug && i == 6 && j == 15){
//		printf("WIP(6,15) is calling WMB \n");
//	}

	WIP[ij] = MIN(MIN(m1,MIN(m2,m3)),MIN(m4,m5));

//	if (debug){
//		printf("WIP(%d,%d): m1 = %d, m2 = %d, m3 = %d, m4 = %d, m5 = %d ==> min = %d \n",i,j,m1,m2,m3,m4,m5,WIP[ij]);
//	}
}

void pseudo_loop::compute_WIP_pkonly(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (WIP[ij] < INF/2){ // was calculated before
		return;
	}
	if (fres[i].arc != fres[j].arc || i == j || weakly_closed[ij]== 0){
		WIP[ij] = INF;
		return;
	}
	int m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;

	// branch 1:
	if (fres[i].pair < 0){
		m1 = get_WIP(i+1,j) + cp_penalty;
	}
	// branch 2:
	if (fres[j].pair < 0){
		m2 = get_WIP(i,j-1) + cp_penalty;
	}
	//branch 3:
	int t;
	for (t = i; t <j; t++){
		int tmp = get_WIP(i,t) + get_WIP(t+1,j);
		if (tmp < m3){
			m3 = tmp;
		}
	}

	// branch 4:
	if (fres[i].pair == j && fres[j].pair == i){

		m4 = V->get_energy(i,j)	+ bp_penalty;

	}

	// branch 5:
	m5 = get_WMB(i,j) + PSM_penalty + bp_penalty;

	WIP[ij] = MIN(MIN(m1,MIN(m2,m3)),MIN(m4,m5));

}





void pseudo_loop::compute_VPP(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (VPP[ij] != INF){ // computed before
//		if (debug){
//			printf("VPP(%d,%d) was calculated before ==> VPP(%d,%d)= %d \n",i,j,i,j,get_VPP(i,j));
//		}
		return;
	}
	if (i == j  || this->is_weakly_closed(i,j)){
		VPP[ij] = INF;
//		if (debug){
//			printf("VPP(%d,%d): i == j ==> VPP = %d \n",i,j,VPP[ij]);
//		}
		return;
	}
	int m1 = INF, m2 = INF;

	// Hosna: July 4th, 2007
	// After discussion with Anne, we figured out that we need to add
	// two more cases to VPP so that it can handle cases that in only one side
	// we have some structure and the other side is empty
	int m3 = INF, m4 = INF;

	//branch 1:
	int r=-1 ;
	// Hosna April 9th, 2007
	// checking the borders as they may be negative numbers
	int max_i_bp = i;
	if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
		max_i_bp = get_bp(i,j);
	}
	for (r = max_i_bp+1; r < j; r++ ){
        // Ian Wark July 19 2017
        // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (fres[r].pair < FRES_RESTRICTED_MIN){
			int tmp = get_VP(i,r) + get_WIP(r+1,j);
//			if (debug){
//				printf("VPP(%d,%d) branch 1: VP(%d,%d) = %d, WIP(%d,%d)= %d ==> tmp = %d  and m1 = %d\n",i,j,i,r,get_VP(i,r),r+1,j,get_WIP(r+1,j),tmp, m1);
//			}
			if (tmp < m1){
				m1 = tmp;
			}
		}
	}

	//branch 2:
	// Hosna April 9th, 2007
	// checking the borders as they may be negative numbers
	int min_Bp_j = j;
	if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_bp(i,j) < min_Bp_j){
		min_Bp_j = get_Bp(i,j);
	}
	for (r = i+1; r < min_Bp_j; r++){
        // Ian Wark July 19 2017
        // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (fres[r].pair < FRES_RESTRICTED_MIN){
			int tmp = get_WIP(i,r-1) + get_VP(r,j);
//			if (debug){
//				printf("VPP(%d,%d) branch 2: WIP(%d,%d) = %d, VP(%d,%d)= %d ==> tmp = %d  and m2 = %d\n",i,j,i,r-1,get_WIP(i,r-1),r,j,get_VP(r,j),tmp, m2);
//			}
			if (tmp < m2){
				m2 = tmp;
			}
		}
	}

	// Branch 3:
//	max_i_bp = i;
//	if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
//		max_i_bp = get_bp(i,j);
//	}
	for (r = max_i_bp+1; r < j; r++ ){
        // Ian Wark July 19 2017
        // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (fres[r].pair < FRES_RESTRICTED_MIN && this->is_empty_region(r+1,j)){
			int tmp = get_VP(i,r) + (cp_penalty *(j-r)); // check the (j-r) part
//			if (debug){
//				printf("VPP(%d,%d) branch 3: VP(%d,%d) = %d, %d *(%d-%d)= %d ==> tmp = %d  and m3 = %d\n",i,j,i,r,get_VP(i,r),cp_penalty,j,r,cp_penalty *(j-r),tmp, m3);
//			}
			if (tmp < m3){
				m3 = tmp;
			}
		}
	}

	// Branch 4:

//	min_Bp_j = j;
//	if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_bp(i,j) < min_Bp_j){
//		min_Bp_j = get_Bp(i,j);
//	}
	for (r = i+1; r < min_Bp_j; r++){
        // Ian Wark July 19 2017
        // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (fres[r].pair < FRES_RESTRICTED_MIN && this->is_empty_region(i,r-1)){
			int tmp = (cp_penalty * (r-i)) + get_VP(r,j);
//			if (debug){
//				printf("VPP(%d,%d) branch 4: %d *(%d-%d) = %d, VP(%d,%d)= %d ==> tmp = %d  and m4 = %d\n",i,j,cp_penalty,r,i,cp_penalty * (r-i),r,j,get_VP(r,j),tmp, m4);
//			}
			if (tmp < m4){
				m4 = tmp;
			}
		}
	}
	int min_branches = m1;
	if (m2 < min_branches){
		min_branches = m2;
	}
	if (m3 < min_branches){
		min_branches = m3;
	}
	if (m4 < min_branches){
		min_branches = m4;
	}
	VPP[ij] = min_branches; //MIN(MIN(m1,m2),MIN(m3,m4));
//	if (debug){
//		printf("VPP(%d,%d): m1 = %d, m2 = %d, m3 = %d and m4 = %d ==> min = %d \n", i,j,m1,m2,m3,m4,VPP[ij]);
//	}
}

void pseudo_loop::compute_BE(int i, int j, int ip, int jp, h_str_features * fres){

//	if (debug && i == 6 && ip == 11){
//		printf("coming to BE to calculate BE(6,69,11,24) \n");
//	}

    // Ian Wark July 19 2017
    // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
    // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if (!( i >= 0 && i <= ip && ip < jp && jp <= j && j < nb_nucleotides && fres[i].pair >= FRES_RESTRICTED_MIN && fres[j].pair >= FRES_RESTRICTED_MIN && fres[ip].pair >= FRES_RESTRICTED_MIN && fres[jp].pair >= FRES_RESTRICTED_MIN && fres[i].pair == j && fres[j].pair == i && fres[ip].pair == jp && fres[jp].pair == ip)){ //impossible cases
//		if (debug && i == 1 && ip == 11){
//			printf("BE(%d,%d,%d,%d): Impossible case!! \n",i,j,ip,jp);
//		}
		return;
	}
	int iip = index[i]+ip-i;
	if (BE[iip] != 0){ // computed before
//		if (debug && i == 6 && ip == 11){
//			printf("BE(%d,%d,%d,%d) was calculated before ==> BE=%d\n",i,j,ip,jp,BE[iip]);
//		}
		return;
	}
	// base case: i.j and ip.jp must be in G
	if (fres[i].pair != j || fres[ip].pair != jp){
//		if (debug ){
//			printf("BE(%d,%d,%d,%d) = INF \n",i,j,ip,jp);
//		}
		BE[iip] = INF;
		return;
	}

	// base case:
	if(i == ip && j == jp && i<j){
//		if (debug ){
//			printf("BE(%d,%d,%d,%d) = 0 \n",i,j,ip,jp);
//		}
		BE[iip] = 0;
		return;
	}

	int m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
	// 1) bp(i+1) == j-1
	if (fres[i+1].pair == j-1){
		m1 = get_e_stP(i,j) + get_BE(i+1,j-1,ip,jp);
//		if(debug ){
//			printf("BE(%d,%d,%d,%d) Case 1: e_stP(%d,%d) = %d and BE(%d,%d,%d,%d) = %d ==> m1 = %d \n",i,j,ip,jp,i,j,get_e_stP(i,j),i+1,j-1,ip,jp,get_BE(i+1,j-1,ip,jp),m1);
//		}
	}

	// cases 2-5 are all need an l s.t. i<l<=ip and jp<=bp(l)<j
	int l;
	for (l = i+1; l<= ip ; l++){
        // Ian Wark July 19 2017
        // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs

		// Hosna: March 14th, 2007
		if (fres[l].pair >= FRES_RESTRICTED_MIN && jp <= fres[l].pair && fres[l].pair < j){
			// Hosna, March 15, 2007
			// since not_paired_all[i,l] includes i and l themselves
			// and in BE energy calculation we are looking for the oepn region (i,l)
			// we have to look at not_paired_all[i+1,l-1]
			int lp = fres[l].pair;
			int il = index[i]+l-i;
			int lpj = index[lp]+j-lp;
			// 2)
			// Hosna June 29, 2007
			// when we pass a stacked pair instead of an internal loop to e_int, it returns underflow,
			// so I am checking explicitely that we won't have stems instead of internal loop
			if (is_empty_region(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1 ){//&& !(ip == (i+1) && jp==(j-1)) && !(l == (i+1) && lp == (j-1))){
				int temp = get_e_intP(i,l,lp,j)+ get_BE(l,lp,ip,jp);
//				if (debug){
//					printf("BE(%d,%d,%d,%d) branch 2: e_intP(%d,%d,%d,%d) = %d, BE(%d,%d,%d,%d)= %d ==> temp = %d  and m2 = %d\n",i,j,ip,jp,i,l,lp,j,get_e_intP(i,l,lp,j),l,lp,ip,jp,get_BE(l,lp,ip,jp));
//				}
				if (m2 > temp){
					m2 = temp;
				}
			}

			// 3)
			if (is_weakly_closed(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				int temp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp) + get_WIP(lp+1,j-1)+ ap_penalty + 2* bp_penalty;
//				if (debug){
//					printf("BE(%d,%d,%d,%d) branch 3: WIP(%d,%d) = %d, BE(%d,%d,%d,%d)= %d, WIP(%d,%d)= %d ==> temp = %d  and m3 = %d\n",i,j,ip,jp,i+1,l-1,get_WIP(i+1,l-1),l,lp,ip,jp,get_BE(l,lp,ip,jp),lp+1,j-1,get_WIP(lp+1,j-1),temp,m3);
//				}
				if (m3 > temp){
					m3 = temp;
				}
			}

			// 4)
			if (is_weakly_closed(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				int temp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp) + cp_penalty * (j-lp+1) + ap_penalty + 2*bp_penalty;
//				if (debug){
//					printf("BE(%d,%d,%d,%d) branch 4: WIP(%d,%d) = %d, BE(%d,%d,%d,%d)= %d ==> temp = %d  and m4 = %d\n",i,j,ip,jp,i+1,l-1,get_WIP(i+1,l-1),l,lp,ip,jp,get_BE(l,lp,ip,jp),temp,m4);
//				}
				if (m4 > temp){
					m4 = temp;
				}
			}

			// 5)
			if (is_empty_region(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				int temp = ap_penalty + 2*bp_penalty + (cp_penalty * (l-i+1)) + get_BE(l,lp,ip,jp) + get_WIP(lp+1,j-1);
//				if (debug && l == 9 && ip == 11){
//					printf("BE(%d,%d,%d,%d) branch 5: BE(%d,%d,%d,%d)= %d and WIP(%d,%d) = %d ==> temp = %d \n", i,j,ip,jp,l,lp,ip,jp,get_BE(l,lp,ip,jp),lp+1,j-1,get_WIP(lp+1,j-1),temp);
//				}
				if (m5 > temp){
					m5 = temp;
				}
			}
		}
	}

	// finding the min and putting it in BE[iip]
	BE[iip] = MIN(m1,MIN(MIN(m2,m3),MIN(m4,m5)));
//	if (debug && i == 6 && ip == 11){
//		printf("BE[%d,%d,%d,%d]: m1 = %d, m2 = %d, m3 = %d, m4 = %d, m5 = %d ==> min = %d \n",i,j,ip,jp,m1,m2,m3,m4,m5,BE[iip]);
//	}
}

int pseudo_loop::get_WI(int i, int j){
	if (i>j){
		return 0;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WI[ij] == 0){
		//printf("get_WI(%d,%d), and we need to compute WI (i.e it's 0)!\n",i,j);
		compute_WI(i,j,fres);
	}
	 */
	//printf("get_WI(%d,%d), after computation its value = %d!\n",i,j, WI[ij]);
	return WI[ij];


}

// Hosna, May 1st, 2012
// I don't think we need specific getter function for pkonly case
/*
int pseudo_loop::get_WI_pkonly(int i, int j){
	if (i>j){
		return 0;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	//if (needs_computation == 1 && WI[ij] == 0){
		//printf("get_WI(%d,%d), and we need to compute WI (i.e it's 0)!\n",i,j);
	//	compute_WI_pkonly(i,j,fres);
	//}

	//printf("get_WI(%d,%d), after computation its value = %d!\n",i,j, WI[ij]);
	return WI[ij];


}
*/

int pseudo_loop::get_VP(int i, int j){
	// Hosna, March 16, 2012
	// two bases should be at least 3 bases apart
	if (j-i < TURN || i >= j || fres[i].pair >= 0 || fres[j].pair >= 0 || this->is_weakly_closed(i,j) == 1){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && VP[ij] == INF){
		//printf("get_VP(%d,%d), and we need to compute VP (i.e. it's INF)!\n",i,j);
		compute_VP(i,j,fres);
	}
	 */
	//printf("get_VP(%d,%d), after computation its value = %d!\n",i,j, VP[ij]);
	return VP[ij];

}
int pseudo_loop::get_WMB(int i, int j){
	// Hosna: July 6th, 2007
	// added impossible cases
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN ||(fres[i].pair >= 0 && fres[i].pair > j) || (fres[j].pair >= 0 && fres[j].pair < i) || (fres[i].pair >= 0 && fres[i].pair < i ) || (fres[j].pair >= 0 && j < fres[j].pair)){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WMB[ij] == INF){
		//printf("get_WMB(%d,%d), and we need to compute WMB (i.e. it's INF)!\n",i,j);
		compute_WMB(i,j,fres);
	}
	 */
	//printf("get_WMB(%d,%d), after computation its value = %d!\n",i,j, WMB[ij]);
	return WMB[ij];
}

// Hosna: April 18th, 2007
// changed WMB to case 2 and WMBP
int pseudo_loop::get_WMBP(int i, int j){
	// Hosna: July 6th, 2007
	// added impossible cases
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i< TURN || (fres[i].pair >= 0 && fres[i].pair > j) || (fres[j].pair >= 0 && fres[j].pair < i) || (fres[i].pair >= 0 && fres[i].pair < i ) || (fres[j].pair >= 0 && j < fres[j].pair)){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WMBP[ij] == INF){
		//printf("get_WMBP(%d,%d), and we need to compute WMBP (i.e. it's INF)!\n",i,j);
		compute_WMBP(i,j,fres);
	}
	 */
	//printf("get_WMBP(%d,%d), after computation its value = %d!\n",i,j, WMBP[ij]);
	return WMBP[ij];
}

int pseudo_loop::get_BE(int i, int j, int ip, int jp){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i>= TURN && i >= 0 && i <= ip && ip < jp && jp <= j && j < nb_nucleotides && fres[i].pair >=0 && fres[j].pair >= 0 && fres[ip].pair >= 0 && fres[jp].pair >= 0 && fres[i].pair == j && fres[j].pair == i && fres[ip].pair == jp && fres[jp].pair == ip){
		if(i == ip && j == jp && i<j){
			return 0;
		}
		int iip = index[i]+ip-i;
		// Hosna, May 1st , 2012
		// these parts are not needed any more
		/*
		if (needs_computation == 1 && BE[iip] == 0){
//			if (debug ){
//				printf("computing BE(%d,%d,%d,%d) \n",i,j,ip,jp);
//			}
			//printf("get_BE(%d,%d), and we need to compute BE (i.e. it's 0)!\n",i,ip);
			compute_BE(i,j,ip,jp,fres);
		}
		 */
//		if (debug ){
//			printf("BE[%d,%d]=%d \n",i,ip,BE[iip]);
//		}
		//printf("get_BE(%d,%d), after computation its value = %d!\n",i,ip, BE[iip]);
		return BE[iip];
	}else{
//		if (debug && i == 6 && ip == 11){
//			printf("returning INF from BE(6,69,11,24) \n");
//		}
		return INF;
	}
}

int pseudo_loop::get_WIP(int i, int j){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN || i >= j || this->is_weakly_closed(i,j) != 1){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WIP[ij]== INF){
		//printf("get_WIP(%d,%d), and we need to compute WIP (i.e. it's INF)!\n",i,j);
		compute_WIP(i,j,fres);
	}
	 */
	//printf("get_WIP(%d,%d), after computation its value = %d!\n",i,j, WIP[ij]);
	return WIP[ij];
}

// Hosna, May 1st, 2012
// I don't think we need specific getter function for pkonly case
/*
int pseudo_loop::get_WIP_pkonly(int i, int j){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN || i >= j || this->is_weakly_closed(i,j) != 1){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	//if (needs_computation == 1 && WIP[ij]== INF){
		//printf("get_WIP(%d,%d), and we need to compute WIP (i.e. it's INF)!\n",i,j);
	//	compute_WIP_pkonly(i,j,fres);
	//}

	//printf("get_WIP(%d,%d), after computation its value = %d!\n",i,j, WIP[ij]);
	return WIP[ij];
}
*/

int pseudo_loop::get_VPP(int i, int j){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN || i >= j || this->is_weakly_closed(i,j) == 1){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && VPP[ij] == INF)
	{
		//printf("get_VPP(%d,%d), and we need to compute VPP (i.e. it's INF)!\n",i,j);
		compute_VPP(i,j,fres);
	}
	 */
	//printf("get_VPP(%d,%d), after computation its value = %d!\n",i,j, VPP[ij]);
	return VPP[ij];

}

// PRE: i< j
int pseudo_loop::get_b(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j){
		return INF;
	}
	int border = MIN(border_bs[j][i],INF);
//	if(debug){
//		printf("b(%d,%d) = %d \n",i,j,border);
//	}
	return border;
}

// PRE: i<j
int pseudo_loop::get_bp(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j ){
		return -1;
	}
	int border = MAX(border_bps[j][i],-1);
//	if(debug){
//		printf("bp(%d,%d) = %d \n",i,j,border);
//	}
	return border;
}
//PRE: i<j
int pseudo_loop::get_B(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j ){
		return -1;
	}
	int border = MAX(border_bs[i][j],-1);
//	if(debug){
//		printf("B(%d,%d) = %d \n",i,j,border);
//	}
	return border;
}
//PRE: i<j
int pseudo_loop::get_Bp(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j ){
		return INF;
	}
	int border = MIN(border_bps[i][j],INF);
//	if(debug){
//		printf("Bp(%d,%d) = %d \n",i,j,border);
//	}
	return border;
}

int pseudo_loop::get_e_stP(int i, int j){
	if (i+1 == j-1){ // TODO: do I need something like that or stack is taking care of this?
		return INF;
	}
	int ss = S->get_energy(i,j,int_sequence);
//	if (debug){
//		printf("stack energy got from simfold is %d \n", ss);
//	}
	return (int)round(e_stP_penalty * (double)ss);
}

int pseudo_loop::get_e_intP(int i, int ip, int jp, int j){
	// Hosna Feb 12th, 2007:
	// this function is only being called in branch 5 of VP
	// and branch 2 of BE
	// in both cases regions [i,ip] and [jp,j] are closed regions
	int e_int = VBI->get_energy(i,j,ip,jp,int_sequence);
//	int uPup = ((ip-i-1)+(j-jp-1)) * PUP_penalty;
//	return MIN(VBI_energy,uPup);

	// Hosna April 3rd, 2007
	// based on the discussion with Anne, we decided to have
	// e_intP = 0.83 * e_int
//	printf("test: e_int(5,30,7,29) = %d \n",VBI->get_energy(5,30,7,29,int_sequence));
//	printf("e_int(%d,%d,%d,%d) = %d \n",i,j,ip,jp,e_int);
	int energy = (int)round(e_intP_penalty * (double)e_int);
//	printf("e_intP(%d,%d,%d,%d) = %d \n",i,ip,jp,j,energy);
	return energy;
}

int pseudo_loop::get_energy(int i, int j){
	return get_WMB(i,j);
}

void pseudo_loop::back_track(char *structure, minimum_fold *f, seq_interval *cur_interval)
{
	this->structure = structure;
	this->f = f;
	needs_computation = 0;
	// Hosna March 8, 2012
	// changing the nested if structure to switch for optimality
	switch (cur_interval->type)
	{
			case P_WMB:
			{
				int i = cur_interval->i;
				int j = cur_interval->j;
				if (debug) {
					printf ("\t(%d,%d) P_WMB\n", i,j);
				}
				if (i >= j)
					return;
				int tmp = INF, best_l = -1, best_row = -1, min = INF;

				// case 1
				if (fres[j].pair >= 0 && j > fres[j].pair){
					int l, acc = INF;
					int bp_j = fres[j].pair;
					for (l = bp_j +1; l < j; l++){
						// Hosna: April 24, 2007
						// correct case 2 such that a multi-pseudoknotted
						// loop would not be treated as case 2

						// Hosna: July 5th, 2007:
						// We don't need this restriction
		//				if (l > fres[i].pair){
							// Hosna April 9th,
							// checking the borders as they may be negative numbers

							// Hosna: July 5th, 2007:
							// We don't need to check for l not being paired
		//					if (fres[l].pair < 0 && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
							if (get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
								int sum = get_BE(bp_j,j,fres[get_Bp(l,j)].pair,get_Bp(l,j)) + get_WMBP(i,l) + get_WI(l+1,get_Bp(l,j)-1);
								if (acc > sum){
									acc = sum;
									best_l = l;
								}
							}
		//				}
					}
					tmp = PB_penalty + acc;
					if (tmp < min){
						min = tmp;
						best_row = 1;
		//				if (debug){
		//					printf("P_WMB: case 2: min = %d, best_l = %d and best_row = %d \n", min, best_l, best_row);
		//				}
					}

				}
				// case WMBP
				tmp = get_WMBP(i,j);
				if (tmp < min){
					min = tmp;
					best_row = 2;
		//			if (debug && i == 1 && j == 88){
		//				printf("P_WMB: case WMBP: min = %d, best_row = %d \n", min, best_row);
		//			}
				}

				switch (best_row)
				{
					case 1:
						if (best_l > -1){
							if (debug){
								printf("WMB(%d,%d)(1): Pushing WMBP(%d,%d), WI(%d,%d) and BE(%d,%d) \n",i,j,i,best_l,best_l+1,get_Bp(best_l,j)-1,fres[j].pair,fres[get_Bp(best_l,j)].pair);
							}
							insert_node(i,best_l,P_WMBP);
							insert_node(best_l +1,get_Bp(best_l,j)-1,P_WI);
							insert_node(fres[j].pair,fres[get_Bp(best_l,j)].pair, P_BE);
						}
						break;
					case 2:
						if (debug){
							printf("WMB(%d,%d)(2): Pushing WMBP(%d,%d) \n",i,j,i,j);
						}
						insert_node(i,j,P_WMBP);
						break;
				}
			}
				break;
			case P_WMBP:
			{
				int i = cur_interval->i;
				int j = cur_interval->j;
				if (debug) {
					printf ("\t(%d,%d) P_WMBP\n", i,j);
				}
				if (i >= j)
					return;
				int tmp = INF, best_l = -1, best_row = -1, min = INF;

				// case 1
				if(fres[j].pair < 0 && fres[i].pair >= 0){
					int l, l1 = -1;
					// Hosna: June 29, 2007
					// if j is inside i's arc then the l should be
					// less than j not bp(i)
					// check with Anne
		//			for (l = i+1; l < MIN(fres[i].pair,j); l++){
					// Hosna: July 5th, 2007:
					// if we have bp(i)> j then we should not have come to the WMBP
					for (l = i+1; l < j; l++){
						// Hosna, April 9th, 2007
						// checking the borders as they may be negative
		//				if(fres[l].pair < 0 && get_bp(i,l) >= 0 && get_bp(i,l) < nb_nucleotides){
						// Hosna: July 5th, 2007:
						// removed bp(l)<0 as VP should handle that
						if(get_bp(i,l) >= 0 && get_bp(i,l) < nb_nucleotides && l+TURN <= j){
							// Hosna: April 19th, 2007
							// the chosen l should be less than border_b(i,j)
		//					if (get_b(i,j) >= 0 && l < get_b(i,j)){
								int bp_i_l = get_bp(i,l);
								int BE_energy = get_BE(i,fres[i].pair,bp_i_l,fres[bp_i_l].pair);
								int WI_energy = get_WI(bp_i_l +1,l-1);
								int VP_energy = get_VP(l,j);
								int sum = BE_energy + WI_energy + VP_energy;
								if (tmp > sum){
									tmp = sum;
									l1 = l;
								}
		//					}
						}
					}
					tmp = 2*PB_penalty + tmp;
					if (tmp < min){
						min = tmp;
						best_row = 1;
						best_l = l1;
		//				if (debug){
		//					printf("P_WMBP: case 1: min = %d, best_l = %d and best_row = %d \n", min, best_l, best_row);
		//				}
					}
				}
		//		if (debug){
		//			printf("P_WMBP: after case 1: best_l = %d \n", best_l);
		//		}
				// case 3
				if (fres[j].pair < 0){
					int l, acc = INF;
					int l3 = -1;
					for (l = i+1; l<j ; l++)	{
						// Hosna, April 6th, 2007
						// whenever we use get_borders we have to check for the correct values
						if (fres[l].arc > -1 && get_B(l,j) >= 0 && get_B(l,j) < nb_nucleotides && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
							// Hosna: April 19th, 2007
							// the chosen l should be less than border_b(i,j)
							if (get_b(i,j) >= 0 && get_b(i,j)<nb_nucleotides && l < get_b(i,j)){
								// Hosna: June 29 2007
								// after going over the program with Cristina, we noticed that
								// l should be < B'(i,j)
			//					if (l < get_Bp(i,j)){
								// Hosna: July 5th, 2007:
								// as long as we have i <= arc(l)< j we are fine
								if (i <= fres[l].arc && fres[l].arc < j && l+TURN <=j){
									int sum = get_BE(fres[get_B(l,j)].pair,get_B(l,j),fres[get_Bp(l,j)].pair,get_Bp(l,j))+ get_WMBP(i,l-1)+ get_VP(l,j);
									if (acc > sum){
										acc = sum;
										l3 = l;
									}
								}
							}
						}
					}
					// Hosna: April 5th
					// after going over the WMB recurrence with Anne, we think we should add another p_b penalty
					// to the 3rd case ==> 2*P_b
					tmp = 2 *PB_penalty + acc;
					if (tmp < min){
						min = tmp;
						best_row = 3;
						best_l = l3;
		//				if (debug){
		//					printf("P_WMBP: case 3: min = %d, best_l = %d and best_row = %d \n", min, best_l, best_row);
		//				}
					}
				}
		//		if (debug){
		//			printf("P_WMBP: after case 3 best_l = %d \n", best_l);
		//		}
				// case 4
				tmp = get_VP(i,j) + PB_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 4;
		//			if (debug){
		//				printf("P_WMBP: case 4: min = %d and best_row = %d \n", min, best_row);
		//			}
				}

				// case 5
				if(fres[j].pair < j){
					int l, acc = INF;
					for(l = i+1; l<j; l++){
						// Hosna: April 18th, 2007
						// l and j should be in the same arc
						if (fres[l].pair < 0 && fres[l].arc > -1 && fres[j].arc > -1 && fres[l].arc == fres[j].arc){
						tmp = get_WMBP(i,l) + get_WI(l+1,j);
							if (tmp < min){
								min = tmp;
								best_l = l;
								best_row = 5;
		//						if (debug){
		//							printf("P_WMBP: case 5: min = %d, best_l = %d and best_row = %d \n", min, best_l, best_row);
		//						}
							}
						}
					}
				}

		//		if (debug){
		//			printf("P_WMBP: best_row = %d, best_l = %d and min = %d \n",best_row, best_l,min);
		//		}


				switch (best_row)
				{
					case 1:
						if (best_l > -1){
							if (debug){
								printf("WMBP(%d,%d)(1): Pushing BE(%d,%d), WI(%d,%d) and VP(%d,%d) \n",i,j,i,get_bp(i,best_l),get_bp(i,best_l)+1,best_l-1,best_l,j);
							}
							insert_node(i,get_bp(i,best_l),P_BE);
							insert_node(get_bp(i,best_l)+1,best_l-1,P_WI);
							insert_node(best_l,j,P_VP);
						}
						break;
					case 3:
						if (best_l > -1){
							if (debug){
								printf("WMBP(%d,%d)(2): Pushing WMBP(%d,%d), VP(%d,%d) and BE(%d,%d) \n",i,j,i,best_l-1,best_l,j,fres[get_B(best_l,j)].pair,fres[get_Bp(best_l,j)].pair);
							}
							insert_node(i,best_l -1,P_WMBP);
							insert_node(best_l,j,P_VP);
							insert_node(fres[get_B(best_l,j)].pair,fres[get_Bp(best_l,j)].pair,P_BE);
						}
						break;
					case 4:
						if (debug){
							printf("WMBP(%d,%d)(3): Pushing VP(%d,%d) \n",i,j,i,j);
						}
						insert_node(i,j,P_VP);
						break;
					case 5:
						if (best_l > -1){
							if (debug){
								printf("WMBP(%d,%d)(4): Pushing WMBP(%d,%d) and WI(%d,%d)\n",i,j,i,best_l,best_l+1,j);
							}
							insert_node(i,best_l,P_WMBP);
							insert_node(best_l +1,j,P_WI);
						}
						break;
				}

			}
				break;
			case P_VP:
			{
				int i = cur_interval->i;
				int j = cur_interval->j;
				if (debug) {
					printf ("\t(%d,%d) P_VP\n", i,j);
				}
				if (i>=j){
					return;
				}
				f[i].pair = j;
				f[j].pair = i;
				structure[i] = '[';
				structure[j] = ']';
				//printf("----> original VP: adding (%d,%d) <-------\n",i,j);
				f[i].type = P_VP;
				f[j].type = P_VP;

				int min = INF, tmp = INF, best_ip = INF, best_jp = INF, best_row = -1, best_r = INF;

				//case 1
				// Hosna April 9th, 2007
				// need to check the borders as they may be negative
				if(fres[i].arc > -1 && fres[j].arc == -1 && get_Bp(i,j) >= 0 && get_Bp(i,j)< nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides){
					int Bp_i = get_Bp(i,j);
					int B_i = get_B(i,j);
					int WI_ipus1_BPminus = get_WI(i+1,Bp_i - 1) ;
					int WI_Bplus_jminus = get_WI(B_i + 1,j-1);
					tmp =   WI_ipus1_BPminus + WI_Bplus_jminus;
					if (tmp < min){
						min = tmp;
						best_row = 1;
						if(debug){
		//					printf("P_VP: case 1: min = %d and best_row = %d \n",min,best_row);
						}
					}
				}
				//case 2
				// Hosna April 9th, 2007
				// checking the borders as they may be negative
				if (fres[i].arc == -1 && fres[j].arc > -1 && get_b(i,j)>= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j) >= 0 && get_bp(i,j) < nb_nucleotides){
					int b_i = get_b(i,j);
					int bp_i = get_bp(i,j);
					int WI_i_plus_b_minus = get_WI(i+1,b_i - 1);
					int WI_bp_plus_j_minus = get_WI(bp_i + 1,j-1);
					tmp = WI_i_plus_b_minus + WI_bp_plus_j_minus;
					if (tmp < min){
						min = tmp;
						best_row = 2;
						if(debug){
		//					printf("P_VP: case 2: min = %d and best_row = %d \n",min,best_row);
						}
					}
				}
				//case 3
				// Hosna April 9th, 2007
				// checking the borders as they may be negative
				if(fres[i].arc > -1 && fres[j].arc > -1 && get_Bp(i,j) >= 0 && get_Bp(i,j) < nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides && get_b(i,j) >= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j)>= 0 && get_bp(i,j) < nb_nucleotides){
					int Bp_i = get_Bp(i,j);
					int B_i = get_B(i,j);
					int b_i = get_b(i,j);
					int bp_i = get_bp(i,j);
					int WI_i_plus_Bp_minus = get_WI(i+1,Bp_i - 1);
					int WI_B_plus_b_minus = get_WI(B_i + 1,b_i - 1);
					int WI_bp_plus_j_minus = get_WI(bp_i +1,j - 1);
					tmp = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
					if (tmp < min){
						min = tmp;
						best_row = 3;
						if(debug){
		//					printf("P_VP: case 3: min = %d and best_row = %d \n",min,best_row);
						}
					}
				}
				//case 4
				if(fres[i+1].pair < 0 && fres[j-1].pair < 0 && can_pair(int_sequence[i+1],int_sequence[j-1])){
					tmp = get_e_stP(i,j)+ get_VP(i+1,j-1);
					if (tmp < min){
						min = tmp;
						best_row = 4;
						if(debug){
		//					printf("P_VP: case 4: min = %d and best_row = %d \n",min,best_row);
						}
					}
				}
				//case 5
				int ip, jp;
				// Hosna, April 9th, 2007
				// whenever we use get_borders we have to check for the correct values
				int min_borders = 0; // what if both are negative
				if (get_Bp(i,j)> 0 && get_Bp(i,j) < nb_nucleotides && get_b(i,j) >0 && get_b(i,j) < nb_nucleotides){
					min_borders = MIN(get_Bp(i,j),get_b(i,j));
				}else if (get_b(i,j) > 0 && get_b(i,j) < nb_nucleotides && (get_Bp(i,j) < 0 || get_Bp(i,j) > nb_nucleotides)){
					min_borders = get_b(i,j);
				}else if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && (get_b(i,j) < 0 || get_b(i,j) > nb_nucleotides)){
					min_borders = get_Bp(i,j);
				}
				for (ip = i+1; ip < min_borders; ip++){
					// Hosna: April 20, 2007
					// i and ip and j and jp should be in the same arc
					// it should also be the case that [i+1,ip-1] && [jp+1,j-1] are empty regions
					if (fres[ip].pair < 0 && fres[i].arc == fres[ip].arc && is_empty_region(i+1,ip-1)==1){
						// Hosna, April 9th, 2007
						// whenever we use get_borders we have to check for the correct values
						int max_borders= 0;
						if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides){
							max_borders = MAX(get_bp(i,j),get_B(i,j));
						}else if (get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides && (get_bp(i,j) < 0 || get_bp(i,j) > nb_nucleotides)){
							max_borders = get_B(i,j);
						}else if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && (get_B(i,j) < 0 || get_B(i,j) > nb_nucleotides)){
							max_borders = get_bp(i,j);
						}
						for (jp = max_borders+1 ; jp < j ; jp++){
							if (fres[jp].pair < 0 && can_pair(int_sequence[ip],int_sequence[jp]) && is_empty_region(jp+1,j-1)==1){
								// Hosna: April 20, 2007
								// i and ip and j and jp should be in the same arc
								if (fres[j].arc == fres[jp].arc){
									tmp = get_e_intP(i,ip,jp,j) + get_VP(ip,jp);
									if (tmp < min){
										min = tmp;
										best_row = 5;
										best_ip = ip;
										best_jp = jp;
										if(debug){
		//									printf("P_VP: case 5: min = %d, best_row = %d, best_ip = %d and best_jp = %d \n",min,best_row,best_ip,best_jp);
										}
									}
								}
							}
						}
					}
				}
				//case 6
				int r;
				// Hosna April 9th, 2007
				// checking the borders as they may be negative numbers
				int min_Bp_j = j;
				if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_Bp(i,j) < min_Bp_j){
					min_Bp_j = get_Bp(i,j);
				}
				for (r = i+1; r < min_Bp_j ; r++){
					if (fres[r].pair < 0){
						// Hosna: July 5th, 2007
						// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
						// in a multiloop that spans a band
						tmp = get_WIP(i+1,r-1) + get_VPP(r,j-1) + ap_penalty + 2* bp_penalty;
						if (tmp < min){
							min = tmp;
							best_row = 6;
							best_r = r;
							if(debug){
		//						printf("P_VP: case 6: min = %d, best_row = %d and best_r = %d \n",min,best_row,best_r);
							}
						}
					}
				}
				//case 7
				// Hosna April 9th, 2007
				// checking the borders as they may be negative numbers
				int max_i_bp = i;
				if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
					max_i_bp = get_bp(i,j);
				}
				for (r = max_i_bp+1; r < j ; r++){
					if (fres[r].pair < 0){
						// Hosna: July 5th, 2007
						// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
						// in a multiloop that spans a band
						tmp = get_VPP(i+1,r) + get_WIP(r+1,j-1) + ap_penalty + 2* bp_penalty;
						if (tmp < min){
							min = tmp;
							best_row = 7;
							best_r = r;
							if(debug){
		//						printf("P_VP: case 7: min = %d, best_row = %d and best_r = %d \n",min,best_row,best_r);
							}
						}
					}
				}

				switch (best_row)
				{
					case 1:
						if(debug){
							printf("P_VP(%d,%d)(1): Pushing WI(%d,%d) and WI(%d,%d)\n",i,j,i+1,get_Bp(i,j)-1,get_B(i,j)+1,j-1);
						}
						if (i+1 <= get_Bp(i,j)-1){
							if(debug){
								printf("P_VP(%d,%d)(1): Pushing WI(%d,%d) \n",i,j,i+1,get_Bp(i,j)-1);
							}
							insert_node(i+1,get_Bp(i,j)-1,P_WI);
						}
						if (get_B(i,j)+1 <= j-1){
							if(debug){
								printf("P_VP(%d,%d)(1): Pushing WI(%d,%d) \n",i,j,get_B(i,j)+1,j-1);
							}
							insert_node(get_B(i,j)+1,j-1,P_WI);
						}
						break;
					case 2:
						if (debug){
							printf("P_VP(%d,%d)(2): Pushing WI(%d,%d) and WI(%d,%d)\n",i,j,i+1,get_b(i,j)-1,get_bp(i,j)+1,j-1);
						}
						if (i+1 <= get_b(i,j)-1){
							if(debug){
								printf("P_VP(%d,%d)(2): Pushing WI(%d,%d) \n",i,j,i+1,get_b(i,j)-1);
							}
							insert_node(i+1,get_b(i,j)-1,P_WI);
						}
						if (get_bp(i,j)+1 <= j-1){
							if(debug){
								printf("P_VP(%d,%d)(2): Pushing WI(%d,%d) \n",i,j,get_bp(i,j)+1,j-1);
							}
							insert_node(get_bp(i,j)+1,j-1,P_WI);
						}
						break;
					case 3:
						if(debug){
							printf("P_VP(%d,%d)(3): Pushing WI(%d,%d), WI(%d,%d) and WI(%d,%d) \n",i,j,i+1,get_Bp(i,j)-1,get_B(i,j)+1,get_b(i,j)-1,get_bp(i,j)+1,j-1);
						}
						if (i+1 <= get_Bp(i,j)-1){
							if(debug){
								printf("P_VP(%d,%d)(3): Pushing WI(%d,%d) \n",i,j,i+1,get_Bp(i,j)-1);
							}
							insert_node(i+1,get_Bp(i,j)-1,P_WI);
						}
						if (get_B(i,j)+1 <= get_b(i,j)-1){
							if(debug){
								printf("P_VP(%d,%d)(3): Pushing WI(%d,%d) \n",i,j,get_B(i,j)+1,get_b(i,j)-1);
							}
							insert_node(get_B(i,j)+1,get_b(i,j)-1,P_WI);
						}
						if (get_bp(i,j)+1 <= j-1){
							if(debug){
								printf("P_VP(%d,%d)(3): Pushing WI(%d,%d) \n",i,j,get_bp(i,j)+1,j-1);
							}
							insert_node(get_bp(i,j)+1,j-1,P_WI);
						}
						break;
					case 4:
						if (i+1 <= j-1){
							if(debug){
								printf("P_VP(%d,%d)(4): Pushing VP(%d,%d) \n",i,j,i+1,j-1);
							}
							insert_node(i+1,j-1,P_VP);
						}
						break;
					case 5:
						/*if(debug){
							printf("P_VP(%d,%d)(5): Pushing VP(%d,%d) and WI(%d,%d) \n",i,j,best_ip,best_jp,best_ip+1,best_jp-1);
						}*/
						if (best_ip <= best_jp){
							if(debug){
								printf("P_VP(%d,%d)(5): Pushing VP(%d,%d) \n",i,j,best_ip,best_jp);
							}
							insert_node(best_ip,best_jp,P_VP);
						}
						// Hosna, May 28, 2012
						// adding WI recurrence is totally wrong!!
						// as any structures in region [ip+1,jp-1] would be considered in VP recurrence
						/*
						// Hosna: April 20, 2007
						// we should consider any structure that can form in region [ip+1,jp-1]
						if (best_ip+1 < best_jp-1){
							if (debug){
								printf("P_VP(%d,%d)(5): Pushing WI(%d,%d) \n",i,j,best_ip+1,best_jp-1);
							}
							insert_node(best_ip+1,best_jp-1,P_WI);
						}
						 */
						break;
					case 6:
						if (debug){
							printf("P_VP(%d,%d)(6): Pushing WIP(%d,%d) and  VPP(%d,%d) \n",i,j,i+1,best_r-1,best_r,j-1);
						}
						if (i+1 <= best_r-1){
							if(debug){
								printf("P_VP(%d,%d)(6): Pushing WIP(%d,%d) \n",i,j,i+1,best_r-1);
							}
							insert_node(i+1,best_r-1,P_WIP);
						}
						if (best_r <= j-1){
							if(debug){
								printf("P_VP(%d,%d)(6): Pushing VPP(%d,%d) \n",i,j,best_r,j-1);
							}
							insert_node(best_r,j-1,P_VPP);
						}
						break;
					case 7:
						if (debug){
							printf("P_VP(%d,%d)(7): Pushing VPP(%d,%d) and WIP(%d,%d) \n",i,j,i+1,best_r,best_r+1,j-1);
						}
						if (i+1 <= best_r){
							if(debug){
								printf("P_VP(%d,%d)(7): Pushing VPP(%d,%d) \n",i,j,i+1,best_r);
							}
							insert_node(i+1,best_r,P_VPP);
						}
						if (best_r+1 <= j-1){
							if(debug){
								printf("P_VP(%d,%d)(7): Pushing WIP(%d,%d) \n",i,j,best_r+1,j-1);
							}
							insert_node(best_r+1,j-1,P_WIP);
						}
						break;
				}
			}
				break;
			case P_VPP:
			{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (debug) {
				printf ("\t(%d,%d) P_VPP\n", i,j);
			}
			if (i >= j){
				return;
			}
			int min = INF, tmp = INF, best_r = INF, best_row = -1;


			//case 1
			int r ;
			// Hosna April 9th, 2007
			// checking the borders as they may be negative numbers
			int max_i_bp = i;
			if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
				max_i_bp = get_bp(i,j);
			}
			for (r = max_i_bp+1; r < j; r++ ){
				if (fres[r].pair < 0){
					tmp = get_VP(i,r) + get_WIP(r+1,j);
					if (tmp < min){
						min = tmp;
						best_row = 1;
						best_r = r;
						if(debug){
	//						printf("P_VPP: case 1: min = %d, best_row = %d and best_r = %d \n",min,best_row,best_r);
						}
					}
				}
			}
			//case 2:
			// Hosna April 9th, 2007
			// checking the borders as they may be negative numbers
			int min_Bp_j = j;
			if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_bp(i,j) < min_Bp_j){
				min_Bp_j = get_Bp(i,j);
			}
			for (r = i+1; r < min_Bp_j; r++){
				if (fres[r].pair < 0){
					tmp = get_WIP(i,r-1) + get_VP(r,j);
					if (tmp < min){
						min = tmp;
						best_row = 2;
						best_r = r;
						if(debug){
	//						printf("P_VPP: case 2: min = %d, best_row = %d and best_r = %d \n",min,best_row,best_r);
						}
					}
				}
			}

			// Hosna: July 4th, 2007
			// After discussion with Anne, we figured out that we need to add
			// two more cases to VPP so that it can handle cases that in only one side
			// we have some structure and the other side is empty

			// Branch 3:
	//		max_i_bp = i;
	//		if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
	//			max_i_bp = get_bp(i,j);
	//		}
			for (r = max_i_bp+1; r < j; r++ ){
				if (fres[r].pair < 0 && this->is_empty_region(r+1,j)){
					tmp = get_VP(i,r) + cp_penalty *(j-r); // check the (j-r) part
	//				if (debug){
	//					printf("VPP(%d,%d) branch 3: VP(%d,%d) = %d, %d *(%d-%d)= %d ==> tmp = %d  and m3 = %d\n",i,j,i,r,get_VP(i,r),cp_penalty,j,r,cp_penalty *(j-r),tmp, m3);
	//				}
					if (tmp < min){
						min = tmp;
						best_row = 3;
						best_r = r;
					}
				}
			}

			// Branch 4:

	//		min_Bp_j = j;
	//		if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_bp(i,j) < min_Bp_j){
	//			min_Bp_j = get_Bp(i,j);
	//		}
			for (r = i+1; r < min_Bp_j; r++){
				if (fres[r].pair < 0 && this->is_empty_region(i,r-1)){
					tmp = cp_penalty * (r-i) + get_VP(r,j);
	//				if (debug){
	//					printf("VPP(%d,%d) branch 4: %d *(%d-%d) = %d, VP(%d,%d)= %d ==> tmp = %d  and m4 = %d\n",i,j,cp_penalty,r,i,cp_penalty * (r-i),r,j,get_VP(r,j),tmp, m4);
	//				}
					if (tmp < min){
						min = tmp;
						best_row = 4;
						best_r = r;
					}
				}
			}
			switch(best_row)
			{
				case 1:
					if (best_r != INF){
						if (debug){
							printf("P_VPP(%d,%d)(1): Pushing VP(%d,%d) and WIP(%d,%d) \n",i,j,i,best_r,best_r+1,j);
						}
						if (i <= best_r){
							if(debug){
								printf("P_VPP(%d,%d)(1): Pushing VP(%d,%d) \n",i,j,i,best_r);
							}
							insert_node(i,best_r,P_VP);
						}
						if (best_r+1 <= j){
							if(debug){
								printf("P_VPP(%d,%d)(1): Pushing WIP(%d,%d) \n",i,j,best_r+1,j);
							}
							insert_node(best_r +1,j,P_WIP);
						}
					}
					break;
				case 2:
					if (best_r != INF){
						if (debug){
							printf("P_VPP(%d,%d)(2): Pushing WIP(%d,%d) and VP(%d,%d) \n",i,j,i,best_r-1,best_r,j);
						}
						if (i <= best_r-1){
							if(debug){
								printf("P_VPP(%d,%d)(2): Pushing WIP(%d,%d) \n",i,j,i,best_r-1);
							}
							insert_node(i,best_r-1,P_WIP);
						}
						if (best_r <= j){
							if(debug){
								printf("P_VPP(%d,%d)(2): Pushing VP(%d,%d) \n",i,j,best_r,j);
							}
							insert_node(best_r,j,P_VP);
						}
					}
					break;
				case 3:
					if (best_r != INF){
						if (i <= best_r){
							if(debug){
								printf("P_VPP(%d,%d)(3): Pushing VP(%d,%d) \n",i,j,i,best_r);
							}
							insert_node(i,best_r,P_VP);
						}
					}
					break;
				case 4:
					if (best_r != INF){
						if (best_r <= j){
							if(debug){
								printf("P_VPP(%d,%d)(4): Pushing VP(%d,%d) \n",i,j,best_r,j);
							}
							insert_node(best_r,j,P_VP);
						}
					}
					break;
			}
		}
			break;
		case P_WI:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (debug) {
				printf ("\t(%d,%d) P_WI\n", i,j);
			}
			if (i >= j){
				return;
			}
			int min = INF, tmp = INF, best_row = -1, best_t= -1;

	// Hosna: July 2nd, 2007
	// in branch 1 of WI, we can have a case like
	// ((..))((...))
	// such that both i and j are paired but we can chop them
			//case 1
	//    	if (fres[i].pair < 0 && fres[j].pair < 0)
	//		{
				int t;
				for (t = i; t< j; t++){
					int wi_1 = get_WI(i,t);
					int wi_2 = get_WI(t+1,j);
					tmp = wi_1 + wi_2;
					if(tmp < min){
						min = tmp;
						best_row = 1;
						best_t = t;
						if(debug){
	//						printf("P_WI: case 1: min = %d, best_row = %d and best_t = %d \n",min,best_row,best_t);
						}
					}
				}
	//		}
			//case 2
			if ((fres[i].pair == j && fres[j].pair == i)||(fres[i].pair < 0 && fres[j].pair < 0)){
				// Hosna, April 16th, 2007
				// changed back to see if it will work fine
				// Hosna: April 19th, 2007
				// I think we should call the restricted version
	//			tmp = V->get_energy_restricted(i,j,fres) + PPS_penalty;

				tmp = V->get_energy(i,j) + PPS_penalty;
				if(tmp < min){
					min = tmp;
					best_row = 2;
					if(debug){
	//					printf("P_WI: case 2: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 3
			// Hosna: April 20, 2007
			// removed the penalty of PPS

			// Hosna: July 5th, 2007
			// Anne said we should put PPS back
			// change PSM to PSP
			tmp = get_WMB(i,j) + PSP_penalty + PPS_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 3;
				if(debug){
	//				printf("P_WI: case 3: min = %d and best_row = %d \n",min,best_row);
				}
			}
			switch (best_row)
			{
				case 1:
					if (best_t != -1){
						if (debug){
							printf("P_WI(%d,%d)(1): Pushing WI(%d,%d) and WI(%d,%d)\n",i,j,i,best_t,best_t+1,j);
						}
						if (i <= best_t){
							if(debug){
								printf("P_WI(%d,%d)(1): Pushing WI(%d,%d) \n",i,j,i,best_t);
							}
							insert_node(i,best_t,P_WI);
						}
						if (best_t+1 <= j){
							if(debug){
								printf("P_WI(%d,%d)(1): Pushing WI(%d,%d) \n",i,j,best_t+1,j);
							}
							insert_node(best_t+1,j,P_WI);
						}
					}
					break;
				case 2:
					if (i < j){
						if(debug){
							printf("P_WI(%d,%d)(2): Pushing LOOP(%d,%d) \n",i,j,i,j);
						}
						insert_node(i,j,LOOP);
					}
					break;
				case 3:
					if (i < j){
						if(debug){
							printf("P_WI(%d,%d)(3): Pushing WMB(%d,%d) \n",i,j,i,j);
						}
						insert_node(i,j,P_WMB);
					}
					break;
			}
		}
			break;
		case P_BE:
		{
			int i = cur_interval->i;
			int j = fres[i].pair;
			int ip = cur_interval->j;
			int jp = fres[ip].pair;
			if (debug) {
				printf ("\t(%d,%d) P_BE\n", i,ip);
			}
			if (i > ip || i > j || ip > jp || jp > j){
				return;
			}

			f[i].pair = j;
			f[j].pair = i;
			structure[i] = '(';
			structure[j] = ')';
			f[i].type = P_BE;
			f[j].type = P_BE;
			f[ip].pair = jp;
			f[jp].pair = ip;
			structure[ip] = '(';
			structure[jp] = ')';
			f[ip].type = P_BE;
			f[jp].type = P_BE;

	//    	if (debug){
	//    		printf("P_BE: paired (%d,%d) and (%d,%d) in structure \n",i,j,ip,jp);
	//    	}

			int min = INF, tmp = INF, best_row = -1, best_l = INF;
			//case 1
			if (fres[i+1].pair == j-1){
				tmp = get_e_stP(i,j) + get_BE(i+1,j-1,ip,jp);
				if(tmp < min){
					min = tmp;
					best_row = 1;
					/* if(debug){
											printf("P_BE: case 1: min = %d and best_row = %d \n",min,best_row);
										} */
				}
			}
			int l;
			for (l = i+1; l<= ip ; l++){
				if (fres[l].pair >= 0 && jp <= fres[l].pair && fres[l].pair < j){
				int lp = fres[l].pair;
				int il = index[i]+l-i;
				int lpj = index[lp]+j-lp;

	//			if (debug){
	//				printf("BE: checking cases 2-5 for i = %d, j= %d, ip = %d, jp = %d, l = %d and lp = %d \n",i,j,ip,jp,l,lp);
	//
	//			}

				// case 2
	//			if (is_empty_region(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1){
				// Hosna June 29, 2007
				// when we pass a stacked pair instead of an internal loop to e_int, it returns underflow,
				// so I am checking explicitely that we won't have stems instead of internal loop
				if (is_empty_region(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1 ){
					tmp = get_e_intP(i,l,lp,j)+ get_BE(l,lp,ip,jp);
					if (min > tmp){
						min = tmp;
						best_row = 2;
						best_l = l;
						if(debug){
	//						printf("P_BE: case 2: min = %d, best_row = %d and best_l = %d \n",min,best_row,best_l);
						}
					}
				}

				// case 3
				if (is_weakly_closed(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
					tmp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp) + get_WIP(lp+1,j-1);
					if (min > tmp){
						min = tmp;
						best_row = 3;
						best_l = l;
						if(debug){
	//						printf("P_BE: case 3: min = %d, best_row = %d and best_l = %d \n",min,best_row,best_l);
						}
					}
				}

				// case 4
				if (is_weakly_closed(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1){
					// Hosna: July 5th, 2007
					// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
					// in a multiloop that spans a band
					tmp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp) + c_penalty * (j-lp+1) + ap_penalty + 2* bp_penalty;
					if (min > tmp){
						min = tmp;
						best_row = 4;
						best_l = l;
						if(debug){
	//						printf("P_BE: case 4: min = %d, best_row = %d and best_l = %d \n",min,best_row,best_l);
						}
					}
				}

				// case 5
				if (is_empty_region(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
					// Hosna: July 5th, 2007
					// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
					// in a multiloop that spans a band
					tmp = ap_penalty + 2* bp_penalty+ c_penalty * (l-i+1) + get_BE(l,lp,ip,jp) + get_WIP(lp+1,j-1);
					if (min > tmp){
						min = tmp;
						best_row = 5;
						best_l = l;
						if(debug){
	//						printf("P_BE: case 5: min = %d, best_row = %d and best_l = %d \n",min,best_row,best_l);
						}
					}
				}
				}
			}
			switch(best_row)
			{
				case 1:
					if (i+1 <= ip){
						if (debug){
							printf("P_BE(%d,%d)(1): Pushing BE(%d,%d) \n",i,ip,i+1,ip);
						}
						insert_node(i+1,ip,P_BE);
					}
					break;
				case 2:
					if (best_l != INF && best_l <= ip){
						if (debug){
							printf("P_BE(%d,%d)(2): Pushing BE(%d,%d) \n",i,ip,best_l,ip);
						}
						insert_node(best_l,ip,P_BE);
					}
					break;
				case 3:
					if (best_l != INF){
						if (debug){
							printf("P_BE(%d,%d)(3): Pushing WIP(%d,%d), BE(%d,%d) and WIP(%d,%d) \n",i,ip,i+1,best_l-1,best_l,ip,fres[best_l].pair+1,j-1);
						}
						if (i+1 <= best_l-1){
							if (debug){
								printf("P_BE(%d,%d)(3): Pushing WIP(%d,%d) \n",i,ip,i+1,best_l-1);
							}
							insert_node(i+1,best_l-1,P_WIP);
						}
						if (best_l <= ip){
							if (debug){
								printf("P_BE(%d,%d)(3): Pushing BE(%d,%d) \n",i,ip,best_l,ip);
							}
							insert_node(best_l,ip,P_BE);
						}
						if (fres[best_l].pair +1 <= j-1){
							if (debug){
								printf("P_BE(%d,%d)(3): Pushing WIP(%d,%d) \n",i,ip,fres[best_l].pair+1,j-1);
							}
							insert_node(fres[best_l].pair +1,j-1,P_WIP);
						}
					}
					break;
				case 4:
					if (best_l != INF){
						if (debug){
							printf("P_BE(%d,%d)(4): Pushing WIP(%d,%d) and BE(%d,%d)\n",i,ip,i+1,best_l-1,best_l,ip);
						}
						if (i+1 <= best_l-1){
							if (debug){
								printf("P_BE(%d,%d)(4): Pushing WIP(%d,%d) \n",i,ip,i+1,best_l-1);
							}
							insert_node(i+1,best_l-1,P_WIP);
						}
						if (best_l <= ip){
							if (debug){
								printf("P_BE(%d,%d)(4): Pushing BE(%d,%d) \n",i,ip,best_l,ip);
							}
							insert_node(best_l,ip,P_BE);
						}
					}
					break;
				case 5:
					if (best_l != INF){
						if (debug){
							printf("P_BE(%d,%d)(5): Pushing BE(%d,%d) and WIP(%d,%d) \n",i,ip,best_l,ip,fres[best_l].pair+1,j-1);
						}
						if (best_l <= ip){
							if (debug){
								printf("P_BE(%d,%d)(5): Pushing BE(%d,%d) \n",i,ip,best_l,ip);
							}
							insert_node(best_l,ip,P_BE);
						}
						if (fres[best_l].pair +1 <= j-1){
							if (debug){
								printf("P_BE(%d,%d)(5): Pushing WIP(%d,%d) \n",i,ip,fres[best_l].pair+1,j-1);
							}
							insert_node(fres[best_l].pair +1,j-1,P_WIP);
						}
					}
					break;
			}

		}
			break;
		case P_WIP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (debug) {
				printf ("\t(%d,%d) P_WIP\n", i,j);
			}
			if (i == j){
				return;
			}
			int min = INF, tmp = INF, best_row = -1, best_t = INF;
			//case 1
			if (fres[i].pair < 0){
				tmp = get_WIP(i+1,j) + cp_penalty;
				if(tmp < min){
					min = tmp;
					best_row = 1;
					if(debug){
	//					printf("P_WIP: case 1: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 2
			if (fres[j].pair < 0){
				tmp = get_WIP(i,j-1) + cp_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 2;
					if(debug){
	//					printf("P_WIP: case 2: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 3
			int t;
			for (t = i; t <j; t++){
				tmp = get_WIP(i,t) + get_WIP(t+1,j);
	//			if (debug){
	//				printf("WIP, case 3: breaking WIP(%d,%d) to WIP(%d,%d) = %d and WIP(%d,%d) = %d ==> tmp = %d, min = %d \n",i,j,i,t,get_WIP(i,t),t+1,j,get_WIP(t+1,j),tmp,min);
	//			}
				if (tmp < min){
					min = tmp;
					best_row = 3;
					best_t = t;
					if(debug){
	//					printf("P_WIP: case 3: min = %d, best_row = %d and best_t = %d \n",min,best_row,best_t);
					}
				}
			}
			//case 4
			if (fres[i].pair == j || (fres[i].pair < 0 && fres[j].pair < 0 && can_pair(int_sequence[i],int_sequence[j]))){
				tmp = V->get_energy(i,j)+ bp_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 4;
					if(debug){
	//					printf("P_WIP: case 4: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 5
			tmp = get_WMB(i,j) + PSM_penalty + bp_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 5;
				if(debug){
	//				printf("P_WIP: case 5: min = %d and best_row = %d \n",min,best_row);
				}
			}
			switch(best_row)
			{
				case 1:
					if (i+1 <= j){
						if(debug){
							printf("P_WIP(%d,%d)(1): Pushing WIP(%d,%d) \n",i,j,i+1,j);
						}
						insert_node(i+1,j,P_WIP);
					}
					break;
				case 2:
					if (i <= j-1){
						if(debug){
							printf("P_WIP(%d,%d)(2): Pushing WIP(%d,%d) \n",i,j,i,j-1);
						}
						insert_node(i,j-1,P_WIP);
					}
					break;
				case 3:
					if (best_t != INF){
						if (debug){
							printf("P_WIP(%d,%d)(3): Pushing WIP(%d,%d) and WIP(%d,%d) \n",i,j,i,best_t,best_t+1,j);
						}
						if (i <= best_t){
							if(debug){
								printf("P_WIP(%d,%d)(3): Pushing WIP(%d,%d) \n",i,j,i,best_t);
							}
							insert_node(i,best_t,P_WIP);
						}
						if (best_t+1 <= j){
							if(debug){
								printf("P_WIP(%d,%d)(3): Pushing WIP(%d,%d) \n",i,j,best_t+1,j);
							}
							insert_node(best_t +1,j,P_WIP);
						}
					}
					break;
				case 4:
					if (i < j){
						if(debug){
							printf("P_WIP(%d,%d)(4): Pushing LOOP(%d,%d) \n",i,j,i,j);
						}
						insert_node(i,j,LOOP);
					}
					break;
				case 5:
					if (i < j){
						if(debug){
							printf("P_WIP(%d,%d)(5): Pushing WMB(%d,%d) \n",i,j,i,j);
						}
						insert_node(i,j,P_WMB);
					}
					break;
			}
		}
			break;

	//default:
		//printf("Should not happen!!!");
	}
}

void pseudo_loop::insert_node(int i, int j, char type)
{
//	if (debug)
//	{
//		printf("\n*************************\n");
//		printf("WMB insert_node: i = %d j = %d and type = %c \n",i,j,type);
//		printf("*************************\n");
//
//	}
	seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
//    if (debug)
//    {
//    	printf("###################\n");
//    	printf("WMB node inserted: stack_interval.i = %d, stack_interval.j = %d and atck_interval.type = %c \n",stack_interval->i,stack_interval->j,stack_interval->type);
//    	printf("###################\n\n");
//    }

}

void pseudo_loop::set_stack_interval(seq_interval *stack_interval){
	this->stack_interval = stack_interval;
}

void pseudo_loop::back_track_pkonly(char *structure, minimum_fold *f, seq_interval *cur_interval)
{
	this->structure = structure;
	this->f = f;
	needs_computation = 0;
	// Hosna March 8, 2012
	// changing the nested if structure to switch for optimality
	switch (cur_interval->type)
	{
		case P_WMB:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (debug) {
				printf ("\t(%d,%d) P_WMB\n", i,j);
			}
			if (i >= j)
				return;
			int tmp = INF, best_l = -1, best_row = -1, min = INF;

			// case 1
			if (fres[j].pair >= 0 && j > fres[j].pair){
				int l, acc = INF;
				int bp_j = fres[j].pair;
				for (l = bp_j +1; l < j; l++){
					// Hosna: April 24, 2007
					// correct case 2 such that a multi-pseudoknotted
					// loop would not be treated as case 2

					// Hosna: July 5th, 2007:
					// We don't need this restriction
					//				if (l > fres[i].pair){
					// Hosna April 9th,
					// checking the borders as they may be negative numbers

					// Hosna: July 5th, 2007:
					// We don't need to check for l not being paired
					//					if (fres[l].pair < 0 && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
					if (get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
						int sum = get_BE(bp_j,j,fres[get_Bp(l,j)].pair,get_Bp(l,j)) + get_WMBP(i,l) + get_WI(l+1,get_Bp(l,j)-1);
						if (acc > sum){
							acc = sum;
							best_l = l;
						}
					}
					//				}
				}
				tmp = PB_penalty + acc;
				if (tmp < min){
					min = tmp;
					best_row = 1;
					//				if (debug){
					//					printf("P_WMB: case 2: min = %d, best_l = %d and best_row = %d \n", min, best_l, best_row);
					//				}
				}

			}
			// case WMBP
			tmp = get_WMBP(i,j);
			if (tmp < min){
				min = tmp;
				best_row = 2;
				//			if (debug && i == 1 && j == 88){
				//				printf("P_WMB: case WMBP: min = %d, best_row = %d \n", min, best_row);
				//			}
			}

			switch (best_row)
			{
				case 1:
					if (best_l > -1){
						if (debug){
							printf("WMB(%d,%d)(1): Pushing WMBP(%d,%d), WI(%d,%d) and BE(%d,%d) \n",i,j,i,best_l,best_l+1,get_Bp(best_l,j)-1,fres[j].pair,fres[get_Bp(best_l,j)].pair);
						}
						insert_node(i,best_l,P_WMBP);
						insert_node(best_l +1,get_Bp(best_l,j)-1,P_WI);
						insert_node(fres[j].pair,fres[get_Bp(best_l,j)].pair, P_BE);
					}
					break;
				case 2:
					if (debug){
						printf("WMB(%d,%d)(2): Pushing WMBP(%d,%d) \n",i,j,i,j);
					}
					insert_node(i,j,P_WMBP);
					break;
			}
		}
			break;
		case P_WMBP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (debug) {
				printf ("\t(%d,%d) P_WMBP\n", i,j);
			}
			if (i >= j)
				return;
			int tmp = INF, best_l = -1, best_row = -1, min = INF;

			// case 1
			if(fres[j].pair < 0 && fres[i].pair >= 0){
				int l, l1 = -1;
				// Hosna: June 29, 2007
				// if j is inside i's arc then the l should be
				// less than j not bp(i)
				// check with Anne
				//			for (l = i+1; l < MIN(fres[i].pair,j); l++){
				// Hosna: July 5th, 2007:
				// if we have bp(i)> j then we should not have come to the WMBP
				for (l = i+1; l < j; l++){
					// Hosna, April 9th, 2007
					// checking the borders as they may be negative
					//				if(fres[l].pair < 0 && get_bp(i,l) >= 0 && get_bp(i,l) < nb_nucleotides){
					// Hosna: July 5th, 2007:
					// removed bp(l)<0 as VP should handle that
					if(get_bp(i,l) >= 0 && get_bp(i,l) < nb_nucleotides && l+TURN <= j){
						// Hosna: April 19th, 2007
						// the chosen l should be less than border_b(i,j)
						//					if (get_b(i,j) >= 0 && l < get_b(i,j)){
						int bp_i_l = get_bp(i,l);
						int BE_energy = get_BE(i,fres[i].pair,bp_i_l,fres[bp_i_l].pair);
						int WI_energy = get_WI(bp_i_l +1,l-1);
						int VP_energy = get_VP(l,j);
						int sum = BE_energy + WI_energy + VP_energy;
						if (tmp > sum){
							tmp = sum;
							l1 = l;
						}
						//					}
					}
				}
				tmp = 2*PB_penalty + tmp;
				if (tmp < min){
					min = tmp;
					best_row = 1;
					best_l = l1;
					//				if (debug){
					//					printf("P_WMBP: case 1: min = %d, best_l = %d and best_row = %d \n", min, best_l, best_row);
					//				}
				}
			}
			//		if (debug){
			//			printf("P_WMBP: after case 1: best_l = %d \n", best_l);
			//		}
			// case 3
			if (fres[j].pair < 0){
				int l, acc = INF;
				int l3 = -1;
				for (l = i+1; l<j ; l++)	{
					// Hosna, April 6th, 2007
					// whenever we use get_borders we have to check for the correct values
					if (fres[l].arc > -1 && get_B(l,j) >= 0 && get_B(l,j) < nb_nucleotides && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
						// Hosna: April 19th, 2007
						// the chosen l should be less than border_b(i,j)
						if (get_b(i,j) >= 0 && get_b(i,j)<nb_nucleotides && l < get_b(i,j)){
							// Hosna: June 29 2007
							// after going over the program with Cristina, we noticed that
							// l should be < B'(i,j)
							//					if (l < get_Bp(i,j)){
							// Hosna: July 5th, 2007:
							// as long as we have i <= arc(l)< j we are fine
							if (i <= fres[l].arc && fres[l].arc < j && l+TURN <=j){
								int sum = get_BE(fres[get_B(l,j)].pair,get_B(l,j),fres[get_Bp(l,j)].pair,get_Bp(l,j))+ get_WMBP(i,l-1)+ get_VP(l,j);
								if (acc > sum){
									acc = sum;
									l3 = l;
								}
							}
						}
					}
				}
				// Hosna: April 5th
				// after going over the WMB recurrence with Anne, we think we should add another p_b penalty
				// to the 3rd case ==> 2*P_b
				tmp = 2 *PB_penalty + acc;
				if (tmp < min){
					min = tmp;
					best_row = 3;
					best_l = l3;
					//				if (debug){
					//					printf("P_WMBP: case 3: min = %d, best_l = %d and best_row = %d \n", min, best_l, best_row);
					//				}
				}
			}
			//		if (debug){
			//			printf("P_WMBP: after case 3 best_l = %d \n", best_l);
			//		}
			// case 4
			tmp = get_VP(i,j) + PB_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 4;
				//			if (debug){
				//				printf("P_WMBP: case 4: min = %d and best_row = %d \n", min, best_row);
				//			}
			}

			// case 5
			if(fres[j].pair < j){
				int l, acc = INF;
				for(l = i+1; l<j; l++){
					// Hosna: April 18th, 2007
					// l and j should be in the same arc
					if (fres[l].pair < 0 && fres[l].arc > -1 && fres[j].arc > -1 && fres[l].arc == fres[j].arc){
						tmp = get_WMBP(i,l) + get_WI(l+1,j);
						if (tmp < min){
							min = tmp;
							best_l = l;
							best_row = 5;
							//						if (debug){
							//							printf("P_WMBP: case 5: min = %d, best_l = %d and best_row = %d \n", min, best_l, best_row);
							//						}
						}
					}
				}
			}

			//		if (debug){
			//			printf("P_WMBP: best_row = %d, best_l = %d and min = %d \n",best_row, best_l,min);
			//		}


			switch (best_row)
			{
				case 1:
					if (best_l > -1){
						if (debug){
							printf("WMBP(%d,%d)(1): Pushing BE(%d,%d), WI(%d,%d) and VP(%d,%d) \n",i,j,i,get_bp(i,best_l),get_bp(i,best_l)+1,best_l-1,best_l,j);
						}
						insert_node(i,get_bp(i,best_l),P_BE);
						insert_node(get_bp(i,best_l)+1,best_l-1,P_WI);
						insert_node(best_l,j,P_VP);
					}
					break;
				case 3:
					if (best_l > -1){
						if (debug){
							printf("WMBP(%d,%d)(2): Pushing WMBP(%d,%d), VP(%d,%d) and BE(%d,%d) \n",i,j,i,best_l-1,best_l,j,fres[get_B(best_l,j)].pair,fres[get_Bp(best_l,j)].pair);
						}
						insert_node(i,best_l -1,P_WMBP);
						insert_node(best_l,j,P_VP);
						insert_node(fres[get_B(best_l,j)].pair,fres[get_Bp(best_l,j)].pair,P_BE);
					}
					break;
				case 4:
					if (debug){
						printf("WMBP(%d,%d)(3): Pushing VP(%d,%d) \n",i,j,i,j);
					}
					insert_node(i,j,P_VP);
					break;
				case 5:
					if (best_l > -1){
						if (debug){
							printf("WMBP(%d,%d)(4): Pushing WMBP(%d,%d) and WI(%d,%d)\n",i,j,i,best_l,best_l+1,j);
						}
						insert_node(i,best_l,P_WMBP);
						insert_node(best_l +1,j,P_WI);
					}
					break;
			}

		}
			break;
		case P_VP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (debug) {
				printf ("\t(%d,%d) P_VP\n", i,j);
			}
			if (i>=j){
				return;
			}
			f[i].pair = j;
			f[j].pair = i;
			structure[i] = '[';
			structure[j] = ']';
			//printf("----> pkonly VP: adding (%d,%d) <-------\n",i,j);
			f[i].type = P_VP;
			f[j].type = P_VP;

			int min = INF, tmp = INF, best_ip = INF, best_jp = INF, best_row = -1, best_r = INF;

			//case 1
			// Hosna April 9th, 2007
			// need to check the borders as they may be negative
			if(fres[i].arc > -1 && fres[j].arc == -1 && get_Bp(i,j) >= 0 && get_Bp(i,j)< nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides){
				int Bp_i = get_Bp(i,j);
				int B_i = get_B(i,j);
				int WI_ipus1_BPminus = get_WI(i+1,Bp_i - 1) ;
				int WI_Bplus_jminus = get_WI(B_i + 1,j-1);
				tmp =   WI_ipus1_BPminus + WI_Bplus_jminus;
				if (tmp < min){
					min = tmp;
					best_row = 1;
					if(debug){
						//					printf("P_VP: case 1: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 2
			// Hosna April 9th, 2007
			// checking the borders as they may be negative
			if (fres[i].arc == -1 && fres[j].arc > -1 && get_b(i,j)>= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j) >= 0 && get_bp(i,j) < nb_nucleotides){
				int b_i = get_b(i,j);
				int bp_i = get_bp(i,j);
				int WI_i_plus_b_minus = get_WI(i+1,b_i - 1);
				int WI_bp_plus_j_minus = get_WI(bp_i + 1,j-1);
				tmp = WI_i_plus_b_minus + WI_bp_plus_j_minus;
				if (tmp < min){
					min = tmp;
					best_row = 2;
					if(debug){
						//					printf("P_VP: case 2: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 3
			// Hosna April 9th, 2007
			// checking the borders as they may be negative
			if(fres[i].arc > -1 && fres[j].arc > -1 && get_Bp(i,j) >= 0 && get_Bp(i,j) < nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides && get_b(i,j) >= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j)>= 0 && get_bp(i,j) < nb_nucleotides){
				int Bp_i = get_Bp(i,j);
				int B_i = get_B(i,j);
				int b_i = get_b(i,j);
				int bp_i = get_bp(i,j);
				int WI_i_plus_Bp_minus = get_WI(i+1,Bp_i - 1);
				int WI_B_plus_b_minus = get_WI(B_i + 1,b_i - 1);
				int WI_bp_plus_j_minus = get_WI(bp_i +1,j - 1);
				tmp = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
				if (tmp < min){
					min = tmp;
					best_row = 3;
					if(debug){
						//					printf("P_VP: case 3: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 4
			if(fres[i+1].pair < 0 && fres[j-1].pair < 0 && can_pair(int_sequence[i+1],int_sequence[j-1])){
				tmp = get_e_stP(i,j)+ get_VP(i+1,j-1);
				if (tmp < min){
					min = tmp;
					best_row = 4;
					if(debug){
						//					printf("P_VP: case 4: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 5
			int ip, jp;
			// Hosna, April 9th, 2007
			// whenever we use get_borders we have to check for the correct values
			int min_borders = 0; // what if both are negative
			if (get_Bp(i,j)> 0 && get_Bp(i,j) < nb_nucleotides && get_b(i,j) >0 && get_b(i,j) < nb_nucleotides){
				min_borders = MIN(get_Bp(i,j),get_b(i,j));
			}else if (get_b(i,j) > 0 && get_b(i,j) < nb_nucleotides && (get_Bp(i,j) < 0 || get_Bp(i,j) > nb_nucleotides)){
				min_borders = get_b(i,j);
			}else if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && (get_b(i,j) < 0 || get_b(i,j) > nb_nucleotides)){
				min_borders = get_Bp(i,j);
			}
			for (ip = i+1; ip < min_borders; ip++){
				// Hosna: April 20, 2007
				// i and ip and j and jp should be in the same arc
				// it should also be the case that [i+1,ip-1] && [jp+1,j-1] are empty regions
				if (fres[ip].pair < 0 && fres[i].arc == fres[ip].arc && is_empty_region(i+1,ip-1)==1){
					// Hosna, April 9th, 2007
					// whenever we use get_borders we have to check for the correct values
					int max_borders= 0;
					if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides){
						max_borders = MAX(get_bp(i,j),get_B(i,j));
					}else if (get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides && (get_bp(i,j) < 0 || get_bp(i,j) > nb_nucleotides)){
						max_borders = get_B(i,j);
					}else if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && (get_B(i,j) < 0 || get_B(i,j) > nb_nucleotides)){
						max_borders = get_bp(i,j);
					}
					for (jp = max_borders+1 ; jp < j ; jp++){
						if (fres[jp].pair < 0 && can_pair(int_sequence[ip],int_sequence[jp]) && is_empty_region(jp+1,j-1)==1){
							// Hosna: April 20, 2007
							// i and ip and j and jp should be in the same arc
							if (fres[j].arc == fres[jp].arc){
								tmp = get_e_intP(i,ip,jp,j) + get_VP(ip,jp);
								if (tmp < min){
									min = tmp;
									best_row = 5;
									best_ip = ip;
									best_jp = jp;
									if(debug){
										//									printf("P_VP: case 5: min = %d, best_row = %d, best_ip = %d and best_jp = %d \n",min,best_row,best_ip,best_jp);
									}
								}
							}
						}
					}
				}
			}
			//case 6
			int r;
			// Hosna April 9th, 2007
			// checking the borders as they may be negative numbers
			int min_Bp_j = j;
			if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_Bp(i,j) < min_Bp_j){
				min_Bp_j = get_Bp(i,j);
			}
			for (r = i+1; r < min_Bp_j ; r++){
				if (fres[r].pair < 0){
					// Hosna: July 5th, 2007
					// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
					// in a multiloop that spans a band
					tmp = get_WIP(i+1,r-1) + get_VPP(r,j-1) + ap_penalty + 2* bp_penalty;
					if (tmp < min){
						min = tmp;
						best_row = 6;
						best_r = r;
						if(debug){
							//						printf("P_VP: case 6: min = %d, best_row = %d and best_r = %d \n",min,best_row,best_r);
						}
					}
				}
			}
			//case 7
			// Hosna April 9th, 2007
			// checking the borders as they may be negative numbers
			int max_i_bp = i;
			if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
				max_i_bp = get_bp(i,j);
			}
			for (r = max_i_bp+1; r < j ; r++){
				if (fres[r].pair < 0){
					// Hosna: July 5th, 2007
					// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
					// in a multiloop that spans a band
					tmp = get_VPP(i+1,r) + get_WIP(r+1,j-1) + ap_penalty + 2* bp_penalty;
					if (tmp < min){
						min = tmp;
						best_row = 7;
						best_r = r;
						if(debug){
							//						printf("P_VP: case 7: min = %d, best_row = %d and best_r = %d \n",min,best_row,best_r);
						}
					}
				}
			}

			switch (best_row)
			{
				case 1:
					if(debug){
						printf("P_VP(%d,%d)(1): Pushing WI(%d,%d) and WI(%d,%d)\n",i,j,i+1,get_Bp(i,j)-1,get_B(i,j)+1,j-1);
					}
					if (i+1 <= get_Bp(i,j)-1){
						if(debug){
							printf("P_VP(%d,%d)(1): Pushing WI(%d,%d) \n",i,j,i+1,get_Bp(i,j)-1);
						}
						insert_node(i+1,get_Bp(i,j)-1,P_WI);
					}
					if (get_B(i,j)+1 <= j-1){
						if(debug){
							printf("P_VP(%d,%d)(1): Pushing WI(%d,%d) \n",i,j,get_B(i,j)+1,j-1);
						}
						insert_node(get_B(i,j)+1,j-1,P_WI);
					}
					break;
				case 2:
					if (debug){
						printf("P_VP(%d,%d)(2): Pushing WI(%d,%d) and WI(%d,%d)\n",i,j,i+1,get_b(i,j)-1,get_bp(i,j)+1,j-1);
					}
					if (i+1 <= get_b(i,j)-1){
						if(debug){
							printf("P_VP(%d,%d)(2): Pushing WI(%d,%d) \n",i,j,i+1,get_b(i,j)-1);
						}
						insert_node(i+1,get_b(i,j)-1,P_WI);
					}
					if (get_bp(i,j)+1 <= j-1){
						if(debug){
							printf("P_VP(%d,%d)(2): Pushing WI(%d,%d) \n",i,j,get_bp(i,j)+1,j-1);
						}
						insert_node(get_bp(i,j)+1,j-1,P_WI);
					}
					break;
				case 3:
					if(debug){
						printf("P_VP(%d,%d)(3): Pushing WI(%d,%d), WI(%d,%d) and WI(%d,%d) \n",i,j,i+1,get_Bp(i,j)-1,get_B(i,j)+1,get_b(i,j)-1,get_bp(i,j)+1,j-1);
					}
					if (i+1 <= get_Bp(i,j)-1){
						if(debug){
							printf("P_VP(%d,%d)(3): Pushing WI(%d,%d) \n",i,j,i+1,get_Bp(i,j)-1);
						}
						insert_node(i+1,get_Bp(i,j)-1,P_WI);
					}
					if (get_B(i,j)+1 <= get_b(i,j)-1){
						if(debug){
							printf("P_VP(%d,%d)(3): Pushing WI(%d,%d) \n",i,j,get_B(i,j)+1,get_b(i,j)-1);
						}
						insert_node(get_B(i,j)+1,get_b(i,j)-1,P_WI);
					}
					if (get_bp(i,j)+1 <= j-1){
						if(debug){
							printf("P_VP(%d,%d)(3): Pushing WI(%d,%d) \n",i,j,get_bp(i,j)+1,j-1);
						}
						insert_node(get_bp(i,j)+1,j-1,P_WI);
					}
					break;
				case 4:
					if (i+1 <= j-1){
						if(debug){
							printf("P_VP(%d,%d)(4): Pushing VP(%d,%d) \n",i,j,i+1,j-1);
						}
						insert_node(i+1,j-1,P_VP);
					}
					break;
				case 5:
					/*if(debug){
						printf("P_VP(%d,%d)(5): Pushing VP(%d,%d) and WI(%d,%d) \n",i,j,best_ip,best_jp,best_ip+1,best_jp-1);
					}*/
					if (best_ip <= best_jp){
						if(debug){
							printf("P_VP(%d,%d)(5): Pushing VP(%d,%d) \n",i,j,best_ip,best_jp);
						}
						insert_node(best_ip,best_jp,P_VP);
					}
					// Hosna, May 28, 2012
					// adding WI recurrence is totally wrong!!
					// as any structures in region [ip+1,jp-1] would be considered in VP recurrence
					/*
					// Hosna: April 20, 2007
					// we should consider any structure that can form in region [ip+1,jp-1]
					if (best_ip+1 < best_jp-1){
						if (debug){
							printf("P_VP(%d,%d)(5): Pushing WI(%d,%d) \n",i,j,best_ip+1,best_jp-1);
						}
						insert_node(best_ip+1,best_jp-1,P_WI);
					}
					 */
					break;
				case 6:
					if (debug){
						printf("P_VP(%d,%d)(6): Pushing WIP(%d,%d) and  VPP(%d,%d) \n",i,j,i+1,best_r-1,best_r,j-1);
					}
					if (i+1 <= best_r-1){
						if(debug){
							printf("P_VP(%d,%d)(6): Pushing WIP(%d,%d) \n",i,j,i+1,best_r-1);
						}
						insert_node(i+1,best_r-1,P_WIP);
					}
					if (best_r <= j-1){
						if(debug){
							printf("P_VP(%d,%d)(6): Pushing VPP(%d,%d) \n",i,j,best_r,j-1);
						}
						insert_node(best_r,j-1,P_VPP);
					}
					break;
				case 7:
					if (debug){
						printf("P_VP(%d,%d)(7): Pushing VPP(%d,%d) and WIP(%d,%d) \n",i,j,i+1,best_r,best_r+1,j-1);
					}
					if (i+1 <= best_r){
						if(debug){
							printf("P_VP(%d,%d)(7): Pushing VPP(%d,%d) \n",i,j,i+1,best_r);
						}
						insert_node(i+1,best_r,P_VPP);
					}
					if (best_r+1 <= j-1){
						if(debug){
							printf("P_VP(%d,%d)(7): Pushing WIP(%d,%d) \n",i,j,best_r+1,j-1);
						}
						insert_node(best_r+1,j-1,P_WIP);
					}
					break;
			}
		}
			break;
		case P_VPP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (debug) {
				printf ("\t(%d,%d) P_VPP\n", i,j);
			}
			if (i >= j){
				return;
			}
			int min = INF, tmp = INF, best_r = INF, best_row = -1;


			//case 1
			int r ;
			// Hosna April 9th, 2007
			// checking the borders as they may be negative numbers
			int max_i_bp = i;
			if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
				max_i_bp = get_bp(i,j);
			}
			for (r = max_i_bp+1; r < j; r++ ){
				if (fres[r].pair < 0){
					tmp = get_VP(i,r) + get_WIP(r+1,j);
					if (tmp < min){
						min = tmp;
						best_row = 1;
						best_r = r;
						if(debug){
							//						printf("P_VPP: case 1: min = %d, best_row = %d and best_r = %d \n",min,best_row,best_r);
						}
					}
				}
			}
			//case 2:
			// Hosna April 9th, 2007
			// checking the borders as they may be negative numbers
			int min_Bp_j = j;
			if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_bp(i,j) < min_Bp_j){
				min_Bp_j = get_Bp(i,j);
			}
			for (r = i+1; r < min_Bp_j; r++){
				if (fres[r].pair < 0){
					tmp = get_WIP(i,r-1) + get_VP(r,j);
					if (tmp < min){
						min = tmp;
						best_row = 2;
						best_r = r;
						if(debug){
							//						printf("P_VPP: case 2: min = %d, best_row = %d and best_r = %d \n",min,best_row,best_r);
						}
					}
				}
			}

			// Hosna: July 4th, 2007
			// After discussion with Anne, we figured out that we need to add
			// two more cases to VPP so that it can handle cases that in only one side
			// we have some structure and the other side is empty

			// Branch 3:
			//		max_i_bp = i;
			//		if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
			//			max_i_bp = get_bp(i,j);
			//		}
			for (r = max_i_bp+1; r < j; r++ ){
				if (fres[r].pair < 0 && this->is_empty_region(r+1,j)){
					tmp = get_VP(i,r) + cp_penalty *(j-r); // check the (j-r) part
					//				if (debug){
					//					printf("VPP(%d,%d) branch 3: VP(%d,%d) = %d, %d *(%d-%d)= %d ==> tmp = %d  and m3 = %d\n",i,j,i,r,get_VP(i,r),cp_penalty,j,r,cp_penalty *(j-r),tmp, m3);
					//				}
					if (tmp < min){
						min = tmp;
						best_row = 3;
						best_r = r;
					}
				}
			}

			// Branch 4:

			//		min_Bp_j = j;
			//		if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_bp(i,j) < min_Bp_j){
			//			min_Bp_j = get_Bp(i,j);
			//		}
			for (r = i+1; r < min_Bp_j; r++){
				if (fres[r].pair < 0 && this->is_empty_region(i,r-1)){
					tmp = cp_penalty * (r-i) + get_VP(r,j);
					//				if (debug){
					//					printf("VPP(%d,%d) branch 4: %d *(%d-%d) = %d, VP(%d,%d)= %d ==> tmp = %d  and m4 = %d\n",i,j,cp_penalty,r,i,cp_penalty * (r-i),r,j,get_VP(r,j),tmp, m4);
					//				}
					if (tmp < min){
						min = tmp;
						best_row = 4;
						best_r = r;
					}
				}
			}
			switch(best_row)
			{
				case 1:
					if (best_r != INF){
						if (debug){
							printf("P_VPP(%d,%d)(1): Pushing VP(%d,%d) and WIP(%d,%d) \n",i,j,i,best_r,best_r+1,j);
						}
						if (i <= best_r){
							if(debug){
								printf("P_VPP(%d,%d)(1): Pushing VP(%d,%d) \n",i,j,i,best_r);
							}
							insert_node(i,best_r,P_VP);
						}
						if (best_r+1 <= j){
							if(debug){
								printf("P_VPP(%d,%d)(1): Pushing WIP(%d,%d) \n",i,j,best_r+1,j);
							}
							insert_node(best_r +1,j,P_WIP);
						}
					}
					break;
				case 2:
					if (best_r != INF){
						if (debug){
							printf("P_VPP(%d,%d)(2): Pushing WIP(%d,%d) and VP(%d,%d) \n",i,j,i,best_r-1,best_r,j);
						}
						if (i <= best_r-1){
							if(debug){
								printf("P_VPP(%d,%d)(2): Pushing WIP(%d,%d) \n",i,j,i,best_r-1);
							}
							insert_node(i,best_r-1,P_WIP);
						}
						if (best_r <= j){
							if(debug){
								printf("P_VPP(%d,%d)(2): Pushing VP(%d,%d) \n",i,j,best_r,j);
							}
							insert_node(best_r,j,P_VP);
						}
					}
					break;
				case 3:
					if (best_r != INF){
						if (i <= best_r){
							if(debug){
								printf("P_VPP(%d,%d)(3): Pushing VP(%d,%d) \n",i,j,i,best_r);
							}
							insert_node(i,best_r,P_VP);
						}
					}
					break;
				case 4:
					if (best_r != INF){
						if (best_r <= j){
							if(debug){
								printf("P_VPP(%d,%d)(4): Pushing VP(%d,%d) \n",i,j,best_r,j);
							}
							insert_node(best_r,j,P_VP);
						}
					}
					break;
			}
		}
			break;
		case P_WI:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (debug) {
				printf ("\t(%d,%d) P_WI\n", i,j);
			}
			if (i >= j){
				return;
			}
			int min = INF, tmp = INF, best_row = -1, best_t= -1;

			// Hosna: July 2nd, 2007
			// in branch 1 of WI, we can have a case like
			// ((..))((...))
			// such that both i and j are paired but we can chop them
			//case 1
			//    	if (fres[i].pair < 0 && fres[j].pair < 0)
			//		{
			int t;
			for (t = i; t< j; t++){
				int wi_1 = get_WI(i,t);
				int wi_2 = get_WI(t+1,j);
				tmp = wi_1 + wi_2;
				if(tmp < min){
					min = tmp;
					best_row = 1;
					best_t = t;
					if(debug){
						//						printf("P_WI: case 1: min = %d, best_row = %d and best_t = %d \n",min,best_row,best_t);
					}
				}
			}
			//		}
			//case 2
			// Hosna, May 1st, 2012
			// changed for pkonly case
			if (fres[i].pair == j && fres[j].pair == i){
				// Hosna, April 16th, 2007
				// changed back to see if it will work fine
				// Hosna: April 19th, 2007
				// I think we should call the restricted version
				//			tmp = V->get_energy_restricted(i,j,fres) + PPS_penalty;

				tmp = V->get_energy(i,j) + PPS_penalty;
				if(tmp < min){
					min = tmp;
					best_row = 2;
					if(debug){
						//					printf("P_WI: case 2: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 3
			// Hosna: April 20, 2007
			// removed the penalty of PPS

			// Hosna: July 5th, 2007
			// Anne said we should put PPS back
			// change PSM to PSP
			tmp = get_WMB(i,j) + PSP_penalty + PPS_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 3;
				if(debug){
					//				printf("P_WI: case 3: min = %d and best_row = %d \n",min,best_row);
				}
			}
			switch (best_row)
			{
				case 1:
					if (best_t != -1){
						if (debug){
							printf("P_WI(%d,%d)(1): Pushing WI(%d,%d) and WI(%d,%d)\n",i,j,i,best_t,best_t+1,j);
						}
						if (i <= best_t){
							if(debug){
								printf("P_WI(%d,%d)(1): Pushing WI(%d,%d) \n",i,j,i,best_t);
							}
							insert_node(i,best_t,P_WI);
						}
						if (best_t+1 <= j){
							if(debug){
								printf("P_WI(%d,%d)(1): Pushing WI(%d,%d) \n",i,j,best_t+1,j);
							}
							insert_node(best_t+1,j,P_WI);
						}
					}
					break;
				case 2:
					if (i < j){
						if(debug){
							printf("P_WI(%d,%d)(2): Pushing LOOP(%d,%d) \n",i,j,i,j);
						}
						insert_node(i,j,LOOP);
					}
					break;
				case 3:
					if (i < j){
						if(debug){
							printf("P_WI(%d,%d)(3): Pushing WMB(%d,%d) \n",i,j,i,j);
						}
						insert_node(i,j,P_WMB);
					}
					break;
			}
		}
			break;
		case P_BE:
		{
			int i = cur_interval->i;
			int j = fres[i].pair;
			int ip = cur_interval->j;
			int jp = fres[ip].pair;
			if (debug) {
				printf ("\t(%d,%d) P_BE\n", i,ip);
			}
			if (i > ip || i > j || ip > jp || jp > j){
				return;
			}

			f[i].pair = j;
			f[j].pair = i;
			structure[i] = '(';
			structure[j] = ')';
			f[i].type = P_BE;
			f[j].type = P_BE;
			f[ip].pair = jp;
			f[jp].pair = ip;
			structure[ip] = '(';
			structure[jp] = ')';
			f[ip].type = P_BE;
			f[jp].type = P_BE;

			//    	if (debug){
			//    		printf("P_BE: paired (%d,%d) and (%d,%d) in structure \n",i,j,ip,jp);
			//    	}

			int min = INF, tmp = INF, best_row = -1, best_l = INF;
			//case 1
			if (fres[i+1].pair == j-1){
				tmp = get_e_stP(i,j) + get_BE(i+1,j-1,ip,jp);
				if(tmp < min){
					min = tmp;
					best_row = 1;
					/* if(debug){
					 printf("P_BE: case 1: min = %d and best_row = %d \n",min,best_row);
					 } */
				}
			}
			int l;
			for (l = i+1; l<= ip ; l++){
				if (fres[l].pair >= 0 && jp <= fres[l].pair && fres[l].pair < j){
					int lp = fres[l].pair;
					int il = index[i]+l-i;
					int lpj = index[lp]+j-lp;

					//			if (debug){
					//				printf("BE: checking cases 2-5 for i = %d, j= %d, ip = %d, jp = %d, l = %d and lp = %d \n",i,j,ip,jp,l,lp);
					//
					//			}

					// case 2
					//			if (is_empty_region(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1){
					// Hosna June 29, 2007
					// when we pass a stacked pair instead of an internal loop to e_int, it returns underflow,
					// so I am checking explicitely that we won't have stems instead of internal loop
					if (is_empty_region(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1 ){
						tmp = get_e_intP(i,l,lp,j)+ get_BE(l,lp,ip,jp);
						if (min > tmp){
							min = tmp;
							best_row = 2;
							best_l = l;
							if(debug){
								//						printf("P_BE: case 2: min = %d, best_row = %d and best_l = %d \n",min,best_row,best_l);
							}
						}
					}

					// case 3
					if (is_weakly_closed(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
						tmp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp) + get_WIP(lp+1,j-1);
						if (min > tmp){
							min = tmp;
							best_row = 3;
							best_l = l;
							if(debug){
								//						printf("P_BE: case 3: min = %d, best_row = %d and best_l = %d \n",min,best_row,best_l);
							}
						}
					}

					// case 4
					if (is_weakly_closed(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1){
						// Hosna: July 5th, 2007
						// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
						// in a multiloop that spans a band
						tmp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp) + c_penalty * (j-lp+1) + ap_penalty + 2* bp_penalty;
						if (min > tmp){
							min = tmp;
							best_row = 4;
							best_l = l;
							if(debug){
								//						printf("P_BE: case 4: min = %d, best_row = %d and best_l = %d \n",min,best_row,best_l);
							}
						}
					}

					// case 5
					if (is_empty_region(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
						// Hosna: July 5th, 2007
						// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
						// in a multiloop that spans a band
						tmp = ap_penalty + 2* bp_penalty+ c_penalty * (l-i+1) + get_BE(l,lp,ip,jp) + get_WIP(lp+1,j-1);
						if (min > tmp){
							min = tmp;
							best_row = 5;
							best_l = l;
							if(debug){
								//						printf("P_BE: case 5: min = %d, best_row = %d and best_l = %d \n",min,best_row,best_l);
							}
						}
					}
				}
			}
			switch(best_row)
			{
				case 1:
					if (i+1 <= ip){
						if (debug){
							printf("P_BE(%d,%d)(1): Pushing BE(%d,%d) \n",i,ip,i+1,ip);
						}
						insert_node(i+1,ip,P_BE);
					}
					break;
				case 2:
					if (best_l != INF && best_l <= ip){
						if (debug){
							printf("P_BE(%d,%d)(2): Pushing BE(%d,%d) \n",i,ip,best_l,ip);
						}
						insert_node(best_l,ip,P_BE);
					}
					break;
				case 3:
					if (best_l != INF){
						if (debug){
							printf("P_BE(%d,%d)(3): Pushing WIP(%d,%d), BE(%d,%d) and WIP(%d,%d) \n",i,ip,i+1,best_l-1,best_l,ip,fres[best_l].pair+1,j-1);
						}
						if (i+1 <= best_l-1){
							if (debug){
								printf("P_BE(%d,%d)(3): Pushing WIP(%d,%d) \n",i,ip,i+1,best_l-1);
							}
							insert_node(i+1,best_l-1,P_WIP);
						}
						if (best_l <= ip){
							if (debug){
								printf("P_BE(%d,%d)(3): Pushing BE(%d,%d) \n",i,ip,best_l,ip);
							}
							insert_node(best_l,ip,P_BE);
						}
						if (fres[best_l].pair +1 <= j-1){
							if (debug){
								printf("P_BE(%d,%d)(3): Pushing WIP(%d,%d) \n",i,ip,fres[best_l].pair+1,j-1);
							}
							insert_node(fres[best_l].pair +1,j-1,P_WIP);
						}
					}
					break;
				case 4:
					if (best_l != INF){
						if (debug){
							printf("P_BE(%d,%d)(4): Pushing WIP(%d,%d) and BE(%d,%d)\n",i,ip,i+1,best_l-1,best_l,ip);
						}
						if (i+1 <= best_l-1){
							if (debug){
								printf("P_BE(%d,%d)(4): Pushing WIP(%d,%d) \n",i,ip,i+1,best_l-1);
							}
							insert_node(i+1,best_l-1,P_WIP);
						}
						if (best_l <= ip){
							if (debug){
								printf("P_BE(%d,%d)(4): Pushing BE(%d,%d) \n",i,ip,best_l,ip);
							}
							insert_node(best_l,ip,P_BE);
						}
					}
					break;
				case 5:
					if (best_l != INF){
						if (debug){
							printf("P_BE(%d,%d)(5): Pushing BE(%d,%d) and WIP(%d,%d) \n",i,ip,best_l,ip,fres[best_l].pair+1,j-1);
						}
						if (best_l <= ip){
							if (debug){
								printf("P_BE(%d,%d)(5): Pushing BE(%d,%d) \n",i,ip,best_l,ip);
							}
							insert_node(best_l,ip,P_BE);
						}
						if (fres[best_l].pair +1 <= j-1){
							if (debug){
								printf("P_BE(%d,%d)(5): Pushing WIP(%d,%d) \n",i,ip,fres[best_l].pair+1,j-1);
							}
							insert_node(fres[best_l].pair +1,j-1,P_WIP);
						}
					}
					break;
			}

		}
			break;
		case P_WIP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (debug) {
				printf ("\t(%d,%d) P_WIP\n", i,j);
			}
			if (i == j){
				return;
			}
			int min = INF, tmp = INF, best_row = -1, best_t = INF;
			//case 1
			if (fres[i].pair < 0){
				tmp = get_WIP(i+1,j) + cp_penalty;
				if(tmp < min){
					min = tmp;
					best_row = 1;
					if(debug){
						//					printf("P_WIP: case 1: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 2
			if (fres[j].pair < 0){
				tmp = get_WIP(i,j-1) + cp_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 2;
					if(debug){
						//					printf("P_WIP: case 2: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 3
			int t;
			for (t = i; t <j; t++){
				tmp = get_WIP(i,t) + get_WIP(t+1,j);
				//			if (debug){
				//				printf("WIP, case 3: breaking WIP(%d,%d) to WIP(%d,%d) = %d and WIP(%d,%d) = %d ==> tmp = %d, min = %d \n",i,j,i,t,get_WIP(i,t),t+1,j,get_WIP(t+1,j),tmp,min);
				//			}
				if (tmp < min){
					min = tmp;
					best_row = 3;
					best_t = t;
					if(debug){
						//					printf("P_WIP: case 3: min = %d, best_row = %d and best_t = %d \n",min,best_row,best_t);
					}
				}
			}
			//case 4
			// Hosna, May 1st, 2012
			// changed for pkonly case
			if (fres[i].pair == j && fres[j].pair==i){
				tmp = V->get_energy(i,j)+ bp_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 4;
					if(debug){
						//					printf("P_WIP: case 4: min = %d and best_row = %d \n",min,best_row);
					}
				}
			}
			//case 5
			tmp = get_WMB(i,j) + PSM_penalty + bp_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 5;
				if(debug){
					//				printf("P_WIP: case 5: min = %d and best_row = %d \n",min,best_row);
				}
			}
			switch(best_row)
			{
				case 1:
					if (i+1 <= j){
						if(debug){
							printf("P_WIP(%d,%d)(1): Pushing WIP(%d,%d) \n",i,j,i+1,j);
						}
						insert_node(i+1,j,P_WIP);
					}
					break;
				case 2:
					if (i <= j-1){
						if(debug){
							printf("P_WIP(%d,%d)(2): Pushing WIP(%d,%d) \n",i,j,i,j-1);
						}
						insert_node(i,j-1,P_WIP);
					}
					break;
				case 3:
					if (best_t != INF){
						if (debug){
							printf("P_WIP(%d,%d)(3): Pushing WIP(%d,%d) and WIP(%d,%d) \n",i,j,i,best_t,best_t+1,j);
						}
						if (i <= best_t){
							if(debug){
								printf("P_WIP(%d,%d)(3): Pushing WIP(%d,%d) \n",i,j,i,best_t);
							}
							insert_node(i,best_t,P_WIP);
						}
						if (best_t+1 <= j){
							if(debug){
								printf("P_WIP(%d,%d)(3): Pushing WIP(%d,%d) \n",i,j,best_t+1,j);
							}
							insert_node(best_t +1,j,P_WIP);
						}
					}
					break;
				case 4:
					if (i < j){
						if(debug){
							printf("P_WIP(%d,%d)(4): Pushing LOOP(%d,%d) \n",i,j,i,j);
						}
						insert_node(i,j,LOOP);
					}
					break;
				case 5:
					if (i < j){
						if(debug){
							printf("P_WIP(%d,%d)(5): Pushing WMB(%d,%d) \n",i,j,i,j);
						}
						insert_node(i,j,P_WMB);
					}
					break;
			}
		}
			break;

			//default:
			//printf("Should not happen!!!");
	}
}
