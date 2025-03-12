
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pseudo_loop.h"
#include "V_final.h"
#include "W_final.h"
#include "constants.h"
#include "h_struct.h"
#include "externs.h"
#include "h_externs.h"
#include "h_common.h"




// Hosna June 20th, 2007
// calls the constructor for s_min_folding
// to create all the matrixes required for simfold
// and then calls allocate_space in here to allocate
// space for WMB and V_final
W_final::W_final(char *seq, char *res):s_min_folding(seq,res)
{
	this->nb_nucleotides = strlen(seq);
	/*
        this->int_sequence = new int[this->nb_nucleotides];
	if (int_sequence == NULL) giveup ("Cannot allocate memory", "W_final");
	int i;
        for (i=0; i < this->nb_nucleotides; i++) int_sequence[i] = nuc_to_int(seq[i]);
        */
	space_allocation();
}


W_final::~W_final()
{
	delete vm;
	delete v;
	delete WMB;
	//delete [] int_sequence;

}

// Hosna June 20th, 2007
// allocates space for WMB object and V_final
void W_final::space_allocation(){

	// Hosna June 20th, 2007
	vm = new VM_final(this->int_sequence,this->nb_nucleotides);
	if (vm == NULL) giveup ("Cannot allocate memory", "W_final");
	if (debug){
		printf("nb_nucleotides = %d \n",this->nb_nucleotides);
	}



	// Hosna June 20th, 2007
	// I don't think we need the following line
	//vm->set_energy_matrix(s_min_folding::V);

	// Hosna June 20th, 2007
	v = new V_final(nb_nucleotides);
	if (v == NULL) giveup ("Cannot allocate memory", "W_final");
	//s_min_folding::V, s_min_folding::H, s_min_folding::S, s_min_folding::VBI, vm);
	v->setloops(this->V,vm);

	// Hosna: June 20th 2007
    WMB = new pseudo_loop (sequence,restricted,v,this->H,this->S,this->VBI,vm);
    if (WMB == NULL) giveup ("Cannot allocate memory", "W_final");

    // Hosna: June 20th 2007
    vm->set_V_matrix(v);
    vm->set_WMB_matrix(WMB);


}


double W_final::hfold_pkonly(){
	double energy;
    int i, j;

	h_str_features *h_fres;
    if ((h_fres = new h_str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "h_str_features");
    // detect the structure features
    detect_h_structure_features (restricted, h_fres);
    WMB->set_features(h_fres);
    WMB->initialize();

    str_features *fres;
    if ((fres = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");
    // detect the structure features
    detect_structure_features (restricted, fres);



    // Hosna: June 28, 2007
    // set the features for checking
    v->set_features(fres);

    // Hosna: July 2nd, 2007
    // set the VM matrix for VM_final
    vm->set_VM_matrix(VM);


	// TODO:
	// I think I shoud fill simfold tables here, before filling the HFold tables
	// Hosna, March 8, 2012

	// 1) fill all th ematrices simfold needs (i.e. pk free ones)
	// This is done similar to s_min_folding::fold_sequence_restricted()
	for (j=0; j < nb_nucleotides; j++)
    {
        for (i=0; i<j; i++)
        {
            // V(i,j) = infinity if i restricted or j restricted and pair of i is not j
            if ((fres[i].pair > -1 && fres[i].pair !=j) || (fres[j].pair > -1 && fres[j].pair != i))
                continue;
            if (fres[i].pair == -1 || fres[j].pair == -1)   // i or j MUST be unpaired
                continue;
            V->compute_energy_restricted_pkonly (i, j, fres); // in s_energy_matrix in simfold package


        }
        // if I put this before V calculation, WM(i,j) cannot be calculated, because it returns infinity
        VM->compute_energy_WM_restricted_pkonly (j, fres); // added April 18, 2012

    }


	for (j=0; j < nb_nucleotides; j++)
    {
		// Hosna, March 19, 2012
        for (i =j; i >= 0; i--)//for (i=0; i<=j; i++)
        {
			WMB->compute_energies(i,j);

        	vm->WM_compute_energy(i,j);
        }

	}


	// end of addition at March 8, 2012, Hosna

	for (j= 1; j < nb_nucleotides; j++)
    {
    	this->compute_W_restricted_pkonly(j,fres);
    }

	energy = this->W[nb_nucleotides-1]/100.0;
//    printf("energy = %f \n", energy);


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
        backtrack_restricted_pkonly (cur_interval, fres); // added April 30, 2012 Hosna
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }


    if (debug)
    {
        print_result ();
    }
    delete [] h_fres;
    delete [] fres;
    return energy;

}

double W_final::hfold(){

	double energy;
    int i, j;

	h_str_features *h_fres;
    if ((h_fres = new h_str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "h_str_features");
    // detect the structure features
    detect_h_structure_features (restricted, h_fres);
    WMB->set_features(h_fres);
    WMB->initialize();

    str_features *fres;
    if ((fres = new str_features[nb_nucleotides]) == NULL) giveup ("Cannot allocate memory", "str_features");
    // detect the structure features
    detect_structure_features (restricted, fres);



    // Hosna: June 28, 2007
    // set the features for checking
    v->set_features(fres);

    // Hosna: July 2nd, 2007
    // set the VM matrix for VM_final
    vm->set_VM_matrix(VM);


	// TODO:
	// I think I shoud fill simfold tables here, before filling the HFold tables
	// Hosna, March 8, 2012

	// 1) fill all th ematrices simfold needs (i.e. pk free ones)
	// This is done similar to s_min_folding::fold_sequence_restricted()
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

		// test V values
		/*
		 for (i=0; i<j; i++)
		 {
		 if (fres[i].pair ==j && fres[j].pair ==i){
		 printf("---->> V(%d,%d) = %d \n",i,j, V->get_energy(i,j));

		 }
		 }
		 */
    }



	for (j=0; j < nb_nucleotides; j++)
    {
        for (i =j; i >= 0; i--)//for (i=0; i<=j; i++)
        {
			WMB->compute_energies(i,j);

			vm->WM_compute_energy(i,j);
			//        	if (debug){
			//        		printf("WM_final(%d,%d) = %d \n",i,j,vm->get_energy_WM(i,j));
			//        	}
        }

	}


	// end of addition at March 8, 2012, Hosna

	for (j= 1; j < nb_nucleotides; j++)
    {
    	this->compute_W_restricted(j,fres);
    }
    energy = this->W[nb_nucleotides-1]/100.0;
	//    printf("energy = %f \n", energy);





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
    delete [] h_fres;
    delete [] fres;
    return energy;

}





void W_final::return_structure(char *structure){
	strcpy (structure, this->structure);
	//s_min_folding::return_structure(structure);
}

void W_final::compute_W_restricted (int j, str_features *fres)
// compute W(j)
{
    int m1, m2, m3;
    int must_choose_this_branch;
    m1 = W[j-1];
    m2 = compute_W_br2_restricted (j, fres, must_choose_this_branch);
    m3 = compute_W_br3_restricted (j, fres);
    if (WMB->is_weakly_closed(0,j) < 0){
    	W[j] = INF;
    	return;
    }

    if (must_choose_this_branch)
    {
        W[j] = MIN(m2,m3);
    }
    else
    {
        W[j] = MIN(m1,MIN(m2,m3));
    }
}


void W_final::compute_W_restricted_pkonly (int j, str_features *fres)
// compute W(j)
{
    int m1, m2, m3;
    int must_choose_this_branch;
    m1 = W[j-1];
    m2 = compute_W_br2_restricted_pkonly (j, fres, must_choose_this_branch);
    m3 = compute_W_br3_restricted (j, fres);
    if (WMB->is_weakly_closed(0,j) < 0){
    	W[j] = INF;
    	return;
    }

    if (must_choose_this_branch)
    {
        W[j] = MIN(m2,m3);
    }
    else
    {
        W[j] = MIN(m1,MIN(m2,m3));
    }

}

int W_final::compute_W_br2_restricted (int j, str_features *fres, int &must_choose_this_branch)
{
	int min = INF, tmp, energy_ij = INF, acc;
    int i;
    int chosen = 0;
    int best_i = 0;

	must_choose_this_branch = 0;
    for (i=0; i<=j-1; i++)    // TURN shouldn't be there
    {
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

        // We don't need to make sure i and j don't have to pair with something else,
        //  because that would be INF - done in fold_sequence_restricted
        acc = (i-1>0) ? W[i-1]: 0;

        energy_ij = v->get_energy(i,j);

        if (energy_ij < INF)
        {
            tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
            if (tmp < min)
            {
                min = tmp;
                chosen = 21;        best_i = i;
                if (fres[i].pair == j){
					must_choose_this_branch = 1;
				}
                else                    must_choose_this_branch = 0;
            }
        }

        // I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
        if (fres[i].pair <= -1 && i+1 < j)
        {
            energy_ij = v->get_energy(i+1,j);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
                PARAMTYPE dan = dangle_bot [int_sequence[j]]
											[int_sequence[i+1]]
											[int_sequence[i]];
				//Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i+1].pair == j && fres[j].pair==i+1){
					dan = MIN(0,dan);
				}
				tmp += dan;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 22;  best_i = i;
                    if (fres[i+1].pair == j){
						must_choose_this_branch = 1;
					}
                    else                      must_choose_this_branch = 0;
                }
            }
        }

        // I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
        if (fres[j].pair <= -1 && i < j-1)
        {
            energy_ij = v->get_energy(i,j-1);
            if (energy_ij < INF)
            {
				PARAMTYPE AU_pen=AU_penalty (int_sequence[i],int_sequence[j-1]);
                tmp = energy_ij + AU_pen+ acc;
				PARAMTYPE dan = dangle_top  [int_sequence [j-1]]
											[int_sequence [i]]
											[int_sequence [j]];
				//Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i].pair == j-1 && fres[j-1].pair==i){
					dan = MIN(0,dan);
				 }
                tmp += dan;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 23;  best_i = i;
                    if (fres[i].pair == j-1){
						must_choose_this_branch = 1;
					}
                    else                      must_choose_this_branch = 0;
                }
            }
        }

        if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
        {
            energy_ij = v->get_energy(i+1,j-1);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
				PARAMTYPE dan_bot = dangle_bot [int_sequence[j-1]]
												[int_sequence[i+1]]
												[int_sequence[i]];

				PARAMTYPE dan_top = dangle_top [int_sequence [j-1]]
									[int_sequence [i+1]]
									[int_sequence [j]];
                //Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i+1].pair == j-1 && fres[j-1].pair==i+1){
					dan_bot = MIN(0,dan_bot);
					dan_top = MIN(0,dan_top);
				}
				tmp += dan_bot;
                tmp += dan_top;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 24;  best_i = i;
                    if (fres[i+1].pair == j-1){
						must_choose_this_branch = 1;
					}
                    else                        must_choose_this_branch = 0;
                }
            }
        }
    }
    //printf ("Chosen=%d, best_i=%d\n", chosen, best_i);
    return min;
}

int W_final::compute_W_br2_restricted_pkonly (int j, str_features *fres, int &must_choose_this_branch)
{
	int min = INF, tmp, energy_ij = INF, acc;
    int i;
    int chosen = 0;
    int best_i = 0;

	must_choose_this_branch = 0;
    for (i=0; i<=j-1; i++)    // TURN shouldn't be there
    {
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

        // We don't need to make sure i and j don't have to pair with something else,
        //  because that would be INF - done in fold_sequence_restricted
        acc = (i-1>0) ? W[i-1]: 0;

        energy_ij = (fres[j].pair == i && fres[i].pair ==j)? v->get_energy(i,j) : INF; //v->get_energy(i,j);

        if (energy_ij < INF)
        {
            tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
            if (tmp < min)
            {
                min = tmp;
                chosen = 21;        best_i = i;
                if (fres[i].pair == j){
					must_choose_this_branch = 1;
				}
                else                    must_choose_this_branch = 0;
            }
        }

        // I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
		// in the pk_only version we don't add any pseudoknot free base pairs
        if (fres[i].pair <= -1 && i+1 < j && fres[i+1].pair ==j && fres[j].pair==i+1)
        {
            energy_ij = v->get_energy(i+1,j);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
                PARAMTYPE dan = dangle_bot [int_sequence[j]]
				[int_sequence[i+1]]
				[int_sequence[i]];
				//Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i+1].pair == j && fres[j].pair==i+1){
					dan = MIN(0,dan);
				}
				tmp += dan;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 22;  best_i = i;
                    if (fres[i+1].pair == j){
						must_choose_this_branch = 1;
					}
                    else                      must_choose_this_branch = 0;
                }
            }
        }

        // I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
		// in the pkonly version we don't add any pseudoknot free base pairs
        if (fres[j].pair <= -1 && i < j-1 && fres[i].pair ==j-1 && fres[j-1].pair ==i)
        {
            energy_ij = v->get_energy(i,j-1);
            if (energy_ij < INF)
            {
				PARAMTYPE AU_pen=AU_penalty (int_sequence[i],int_sequence[j-1]);
                tmp = energy_ij + AU_pen+ acc;
				PARAMTYPE dan = dangle_top  [int_sequence [j-1]]
				[int_sequence [i]]
				[int_sequence [j]];
				//Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i].pair == j-1 && fres[j-1].pair==i){
					dan = MIN(0,dan);
				}
                tmp += dan;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 23;  best_i = i;
                    if (fres[i].pair == j-1){
						must_choose_this_branch = 1;
					}
                    else                      must_choose_this_branch = 0;
                }
            }
        }

		// in the pkonly version we don't add any pseudoknot free base pairs
        if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1 && fres[i+1].pair==j-1 && fres[j-1].pair==i+1)
        {
            energy_ij = v->get_energy(i+1,j-1);
            if (energy_ij < INF)
            {
                tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
				PARAMTYPE dan_bot = dangle_bot [int_sequence[j-1]]
				[int_sequence[i+1]]
				[int_sequence[i]];

				PARAMTYPE dan_top = dangle_top [int_sequence [j-1]]
				[int_sequence [i+1]]
				[int_sequence [j]];
                //Hosna, March 27, 2012
				// dangle is INF if the bases are non-canonical and \leq 0 otherwise
				// to accommodate non-canonical base pairing in input structure I add the MIN
				if (fres[i+1].pair == j-1 && fres[j-1].pair==i+1){
					dan_bot = MIN(0,dan_bot);
					dan_top = MIN(0,dan_top);
				}
				tmp += dan_bot;
                tmp += dan_top;
                if (tmp < min)
                {
                    min = tmp;
                    chosen = 24;  best_i = i;
                    if (fres[i+1].pair == j-1){
						must_choose_this_branch = 1;
					}
                    else                        must_choose_this_branch = 0;
                }
            }
        }
    }
    //printf ("Chosen=%d, best_i=%d\n", chosen, best_i);
    return min;
}

int W_final::compute_W_br3_restricted(int j, str_features *fres){
	// Hosna June 30, 2007
	// The following would not take care of when
	// we have some unpaired bases before the start of the WMB
	//return WMB->get_energy(0,j) + PS_penalty;
	int min = INF, tmp, energy_ij = INF, acc;
    int i;
    int chosen = 0;
    int best_i = 0;

//	must_choose_this_branch = 0;
    for (i=0; i<=j-1; i++)    // TURN shouldn't be there
    {
        // don't allow pairing with restricted i's
        // added Jan 28, 2006

		// Hosna: July 9, 2007
		// We only chop W to W + WMB when the bases before WMB are free
		if (i == 0 || (WMB->is_weakly_closed(0,i-1) && WMB->is_weakly_closed(i,j))){

	        // We don't need to make sure i and j don't have to pair with something else,
	        //  because that would be INF - done in fold_sequence_restricted
	        acc = (i-1>0) ? W[i-1]: 0;

	        energy_ij = WMB->get_energy(i,j);
	        if (energy_ij < INF)
	        {
	            tmp = energy_ij + PS_penalty + acc;
	            if (tmp < min)
	            {
	                min = tmp;
	                chosen = 31;
	                best_i = i;
	//                if (fres[i].pair == j)  must_choose_this_branch = 1;
	//                else                    must_choose_this_branch = 0;
	            }
	        }

	        // I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
	        if (fres[i].pair <= -1 && i+1 < j)
	        {
	            energy_ij = WMB->get_energy(i+1,j);
	            if (energy_ij < INF)
	            {
	                tmp = energy_ij + PS_penalty + acc;
	//                tmp += dangle_bot [int_sequence[j]]
	//                                [int_sequence[i+1]]
	//                                [int_sequence[i]];
	                if (tmp < min)
	                {
	                    min = tmp;
	                    chosen = 32;
	                    best_i = i;
	//                    if (fres[i+1].pair == j)  must_choose_this_branch = 1;
	//                    else                      must_choose_this_branch = 0;
	                }
	            }
	        }

	        // I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
	        if (fres[j].pair <= -1 && i < j-1)
	        {
	            energy_ij = WMB->get_energy(i,j-1);
	            if (energy_ij < INF)
	            {
	                tmp = energy_ij + PS_penalty + acc;
	//                tmp += dangle_top [int_sequence [j-1]]
	//                                [int_sequence [i]]
	//                                [int_sequence [j]];
	                if (tmp < min)
	                {
	                    min = tmp;
	                    chosen = 33;
	                    best_i = i;
	//                    if (fres[i].pair == j-1)  must_choose_this_branch = 1;
	//                    else                      must_choose_this_branch = 0;
	                }
	            }
	        }

	        if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
	        {
	            energy_ij = WMB->get_energy(i+1,j-1);
	            if (energy_ij < INF)
	            {
	                tmp = energy_ij + PS_penalty + acc;
	//                tmp += dangle_bot [int_sequence[j-1]]
	//                                [int_sequence[i+1]]
	//                                [int_sequence[i]];
	//                tmp += dangle_top [int_sequence [j-1]]
	//                                [int_sequence [i+1]]
	//                                [int_sequence [j]];
	                if (tmp < min)
	                {
	                    min = tmp;
	                    chosen = 34;
	                    best_i = i;
	//                    if (fres[i+1].pair == j-1)  must_choose_this_branch = 1;
	//                    else                        must_choose_this_branch = 0;
	                }
	            }
	        }
		}
		if (min == -8376) {
			int temp = 0;;
		}
    }
    //printf ("Chosen=%d, best_i=%d\n", chosen, best_i);
    return min;
}

void W_final::backtrack_restricted(seq_interval *cur_interval, str_features *fres){
    char type;



	//Hosna, March 8, 2012
	// changing nested if to switch for optimality
	switch (cur_interval->type){
		case LOOP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (i >= j)
				return;
			f[i].pair = j;
			f[j].pair = i;

			// Hosna Jun. 28 2007
			// if the pairing is part of original structure, put '(' and ')' in the structure
			// otherwise make it '[' and ']'
			if (fres[i].pair == j){
				structure[i] = '(';
				structure[j] = ')';
			}else{
				structure[i] = '[';
				structure[j] = ']';
			}

			type = v->get_type (i,j);
			if (debug)
				printf ("\t(%d,%d) LOOP - type %c\n", i,j,type);
			// Hosna, March 8, 2012
			// changing nested ifs to switch for optimality
			switch (type){
				case STACK:
			//if (type == STACK)
				{
					f[i].type = STACK;
					f[j].type = STACK;
					if (i+1 < j-1)
						this->insert_node(i+1,j-1, LOOP);
						//insert_node (i+1, j-1, LOOP);
					else
					{
						fprintf (stderr, "NOT GOOD RESTR STACK, i=%d, j=%d\n", i, j);
						exit (0);
					}

				}
					break;
				case HAIRP:
			//else if (type == HAIRP)
				{
					f[i].type = HAIRP;
					f[j].type = HAIRP;
				}
					break;
				case INTER:
			//else if (type == INTER)
				{
					f[i].type = INTER;
					f[j].type = INTER;
					// detect the other closing pair
					int ip, jp, best_ip, best_jp, minq;
					int tmp, min = INF;
					for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1); ip++)
					{
						// Hosna, August 28, 2012
						// TODO: cannot understand why we have th efollowing calculations, as it makes the following case be missed!
						// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
						// (____(_____((__(_______)__))_______________________________)____)
						// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
						// in this example int(5,59,11,27) is falsely missed and is equal to INF
						// So I am changing it to the be jp=ip+1; jp<j; jp++ instead
						//minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
						minq = ip+1;
						for (jp = minq; jp < j; jp++)
						{
							if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres))
								continue;
							//tmp = VBI->get_energy_str (i,j,ip,jp);
							// Hosna, March 26, 2012
							// modified to accommodate non-canonical base pairing in restricted structure
							tmp = VBI->get_energy_str_restricted(i,j,ip,jp,fres);
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
						fprintf (stderr, "NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
						exit (0);
					}
				}
					break;
				case MULTI:
			//else if (type == MULTI)
				{
					f[i].type = MULTI;
					f[j].type = MULTI;
					int k, best_k = -1, best_row = -1;
					int tmp= INF, min = INF;
					for (k = i+1; k <= j-1; k++)
					  {
						tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-1);
						if (tmp < min)
						  {
							min = tmp;
							best_k = k;
							best_row = 1;
						  }
						  // TODO:
						  // Hosna, May 1st, 2012
						  // do I need to check for non-canonical base pairings here as well so the dangle values not be INF??
						if (fres[i+1].pair <= -1)
						{
							tmp = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-1) +
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
							tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-2) +
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
							tmp = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-2) +
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

						// Hosna: June 28, 2007
						// the last branch of VM, which is WMB_(i+1),(j-1)
						// Hosna: July 5th, 2007
						// add a b_penalty to this case to match the previous cases
						tmp = WMB->get_energy(i+1,j-1)+ a_penalty + PSM_penalty+ b_penalty;
						if (tmp < min)
						{
							min = tmp;
							best_row = 5;
						}

					  }
					switch (best_row)
					  {
					  case 1:
		//              	printf("M_WM(%d,%d) branch 1: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-1);
						insert_node (i+1, best_k, M_WM);
						insert_node (best_k+1, j-1, M_WM);
						break;
					  case 2:
		//              	printf("M_WM(%d,%d) branch 2: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-1);
						insert_node (i+2, best_k, M_WM);
						insert_node (best_k+1, j-1, M_WM);
						break;
					  case 3:
		//              	printf("M_WM(%d,%d) branch 3: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-2);
						insert_node (i+1, best_k, M_WM);
						insert_node (best_k+1, j-2, M_WM);
						break;
					  case 4:
		//              	printf("M_WM(%d,%d) branch 4: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-2);
						insert_node (i+2, best_k, M_WM);
						insert_node (best_k+1, j-2, M_WM);
						break;
					  // Hosna: June 28, 2007
					  // the last branch of VM, which is WMB_(i+1),(j-1)
					  case 5:
		//              	printf("M_WM(%d,%d) branch 3: pushing P_WMB(%d,%d)\n", i,j,i+1,j-1);
						insert_node(i+1,j-1, P_WMB);
						break;
					  }
				}
					break;
			}
		}
			break;
		case FREE:
		{
			int j = cur_interval->j;

			if (j==0) return;

			int min = INF, tmp, best_row, i, best_i, acc, energy_ij;

			if (debug)
				printf ("\t(0,%d) FREE\n", j);

			// this case is for j unpaired, so I have to check that.
			if (fres[j].pair <= -1)
			{
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
				energy_ij = v->get_energy(i,j);

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

				if (fres[i].pair <= -1)
				{
					energy_ij = v->get_energy(i+1,j);
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
						tmp += dangle_bot [int_sequence[j]]
							[int_sequence[i+1]]
							[int_sequence[i]];
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
					energy_ij = v->get_energy(i,j-1);
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j-1]) + acc;
						tmp += dangle_top [int_sequence[j-1]]
							[int_sequence[i]]
							[int_sequence[j]];
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
					energy_ij = v->get_energy(i+1,j-1);
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
						tmp += dangle_bot [int_sequence[j-1]]
							[int_sequence[i+1]]
							[int_sequence[i]];
						tmp += dangle_top [int_sequence[j-1]]
							[int_sequence[i+1]]
							[int_sequence[j]];
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 4;
						}
					}
				}
			}

			// Hosna: June 28, 2007
			// the last branch of W, which is WMB_i,j
	//        energy_ij = WMB->get_energy(0,j);
	//        if (energy_ij < INF){
	//          	tmp = energy_ij + PS_penalty;
	//           	if (tmp < min){
	//           		min = tmp;
	//           		best_row = 5;
	//           	}
	//        }
		// Hosna June 30, 2007
		// The following would not take care of when
		// we have some unpaired bases before the start of the WMB
		for (i=0; i<=j-1; i++)
		{
			// Hosna: July 9, 2007
			// We only chop W to W + WMB when the bases before WMB are free
			if (i == 0 || (WMB->is_weakly_closed(0,i-1) && WMB->is_weakly_closed(i,j))){

				acc = (i-1>0) ? W[i-1]: 0;

				energy_ij = WMB->get_energy(i,j);

				if (energy_ij < INF)
				{
					tmp = energy_ij + PS_penalty + acc;

					if (tmp < min)
					{
						min = tmp;
						best_row = 5;
						best_i = i;
					}
				}

				// I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
				if (fres[i].pair <= -1 && i+1 < j)
				{
					energy_ij = WMB->get_energy(i+1,j);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 6;
							best_i = i;
						}
					}
				}

				// I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
				if (fres[j].pair <= -1 && i < j-1)
				{
					energy_ij = WMB->get_energy(i,j-1);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 7;
							best_i = i;
						}
					}
				}

				if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
				{
					energy_ij = WMB->get_energy(i+1,j-1);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 8;
							best_i = i;
						}
					}
				}
			}
		}
			switch (best_row)
			{
				case 0:
					//printf("W(%d) case 0: inserting Free (0,%d)\n",j,j-1);
					insert_node (0, j-1, FREE); break;
				case 1:
					//printf("W(%d) case 1: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, LOOP);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 2:
					//printf("W(%d) case 2: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, LOOP);
					if (best_i >= 0)// Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 3:
					//printf("W(%d) case 3: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, LOOP);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 4:
					//printf("W(%d) case 4: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, LOOP);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				// Hosna: June 28, 2007
				// the last branch of W, which is WMB_i,j
				case 5:
					//printf("W(%d) case 5: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, P_WMB);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 6:
					//printf("W(%d) case 6: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 7:
					//printf("W(%d) case 7: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, P_WMB);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 8:
					//printf("W(%d) case 8: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
			}
		}
			break;
		case M_WM:
//  else if(cur_interval->type == M_WM)
		{
			  int i = cur_interval->i;
			  int j = cur_interval->j;
			  int tmp, min = INF;
			  int best_k, best_row;

			  if (debug)
				printf ("\t (%d,%d) M_WM\n", i,j);

			  tmp = v->get_energy(i,j) +
				AU_penalty (int_sequence[i], int_sequence[j]) +
				misc.multi_helix_penalty;
			  if (tmp < min)
				{
				  min = tmp;
				  best_row = 1;
				}
			  if (fres[i].pair <= -1)
			  {
				  tmp = v->get_energy(i+1,j) +
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
				  tmp = v->get_energy(i,j-1) +
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
				  tmp = v->get_energy(i+1,j-1) +
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
				  tmp = vm->get_energy_WM (i+1,j) + misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 5;
				  }
			  }
			  if (fres[j].pair <= -1)
			  {
				  tmp = vm->get_energy_WM (i,j-1) + misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 6;
				  }
			  }

			  for (int k=i; k < j; k++)
				{
					tmp = vm->get_energy_WM (i, k) + vm->get_energy_WM (k+1, j);
					if (tmp < min)
					  {
						min = tmp;
						best_k = k;
						best_row = 7;
					  }
				}
			  // Hosna: June 28, 2007
			  // the last branch of WW, which is WMB_i,j
			  tmp = WMB->get_energy(i,j)+PSM_penalty;
			  if (tmp < min){
				min = tmp;
				best_row = 8;
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
				  // Hosna: June 28, 2007
				  // the last branch of W, which is WMB_i,j
				  case 8:
					insert_node(i,j,P_WMB);
					break;
				  }
			}
			break;
    // Hosna: Feb 19th 2007
		case P_WMB:
   // else if(cur_interval->type == P_WMB)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WMBP:
    // Hosna: April 18th, 2007
    // changed WMB to case 2 and WMBP
   // else if(cur_interval->type == P_WMBP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VP:
    //else if(cur_interval->type == P_VP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VPP:
    //else if(cur_interval->type == P_VPP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WI:
    //else if(cur_interval->type == P_WI)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_BE:
    //else if(cur_interval->type == P_BE)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WIP:
    //else if(cur_interval->type == P_WIP)
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		default:
			printf("Should not be here!\n");
	}


}

void W_final::backtrack_restricted_pkonly (seq_interval *cur_interval, str_features *fres){
    char type;
	//Hosna, March 8, 2012
	// changing nested if to switch for optimality

	switch (cur_interval->type){
			// TODO:
			// April 3, 2012
			// for the pk only case, I don't think I need to change any part of the LOOP case
			// April 18, 2012
			// I think the closing base pairs of the loops are set previously, but the pkonly condition needs to be checked
			// for the inner base pairs.
		case LOOP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (i >= j)
				return;
			f[i].pair = j;
			f[j].pair = i;

			// Hosna Jun. 28 2007
			// if the pairing is part of original structure, put '(' and ')' in the structure
			// otherwise make it '[' and ']'
			if (fres[i].pair == j){
				structure[i] = '(';
				structure[j] = ')';
			}else{
				fprintf(stderr, "Base pairing between %d and %d is pseudoknot free and should not happen here!! \n",i,j);
				exit(0);
				//structure[i] = '[';
				//structure[j] = ']';
			}

			type = v->get_type (i,j);
			if (debug)
				printf ("\t(%d,%d) LOOP - type %c\n", i,j,type);
			// Hosna, March 8, 2012
			// changing nested ifs to switch for optimality
			switch (type){
				case STACK:
					//if (type == STACK)
				{
					f[i].type = STACK;
					f[j].type = STACK;
					if (i+1 < j-1 && fres[i+1].pair == j-1) //check for pkonly
						this->insert_node(i+1,j-1, LOOP);
					//insert_node (i+1, j-1, LOOP);
					else
					{
						fprintf (stderr, "NOT GOOD RESTR STACK (i+1 and j-1 are not paired!), i=%d, j=%d\n", i, j);
						exit (0);
					}

				}
					break;
				case HAIRP:
					//else if (type == HAIRP)
				{
					f[i].type = HAIRP;
					f[j].type = HAIRP;
				}
					break;
				case INTER:
					//else if (type == INTER)
				{
					f[i].type = INTER;
					f[j].type = INTER;
					// detect the other closing pair
					int ip, jp, best_ip, best_jp, minq;
					int tmp, min = INF;
					// Hosna, August 31, 2012
					// The following restriction misses the long restricted loops, so I am chaning it
					//for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1) ; ip++)  // the -TURN shouldn't be there
					for (ip = i+1; ip <= j-2 ; ip++)  // the -TURN shouldn't be there
					{
						// Hosna, August 28, 2012
						// TODO: cannot understand why we have th efollowing calculations, as it makes the following case be missed!
						// GACAUAUAAUCGCGUGGAUAUGGCACGCAAGUUUCUACCGGGCACCGUAAAUGUCCGAUUAUGUC
						// (____(_____((__(_______)__))_______________________________)____)
						// 0....5....1....5....2....5....3....5....4....5....5....5....6...4
						// in this example int(5,59,11,27) is falsely missed and is equal to INF
						// So I am changing it to the be jp=ip+1; jp<j; jp++ instead
						//minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
						minq = ip+1;
						for (jp = minq; jp < j; jp++)
						{
							if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres))
								continue;
							//tmp = VBI->get_energy_str (i,j,ip,jp);
							// Hosna, March 26, 2012
							// modified to accommodate non-canonical base pairing in restricted structure
							// April 18, 2012
							// added a condition for pkonly
							tmp = (fres[ip].pair == jp && fres[jp].pair == ip)? VBI->get_energy_str_restricted(i,j,ip,jp,fres):INF;
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
						fprintf (stderr, "NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
						exit (0);
					}
				}
					break;
				case MULTI:
					//else if (type == MULTI)
				{
					f[i].type = MULTI;
					f[j].type = MULTI;
					int k, best_k = -1, best_row = -1;
					int tmp= INF, min = INF;
					for (k = i+1; k <= j-1; k++)
					{
						tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-1);
						if (tmp < min)
						{
							min = tmp;
							best_k = k;
							best_row = 1;
						}
						if (fres[i+1].pair <= -1)
						{
							tmp = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-1) +
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
							tmp = vm->get_energy_WM (i+1,k) + vm->get_energy_WM (k+1, j-2) +
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
							tmp = vm->get_energy_WM (i+2,k) + vm->get_energy_WM (k+1, j-2) +
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

						// Hosna: June 28, 2007
						// the last branch of VM, which is WMB_(i+1),(j-1)
						// Hosna: July 5th, 2007
						// add a b_penalty to this case to match the previous cases
						tmp = WMB->get_energy(i+1,j-1)+ a_penalty + PSM_penalty+ b_penalty;
						if (tmp < min)
						{
							min = tmp;
							best_row = 5;
						}

					}
					switch (best_row)
					{
						case 1:
							//              	printf("M_WM(%d,%d) branch 1: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-1);
							insert_node (i+1, best_k, M_WM);
							insert_node (best_k+1, j-1, M_WM);
							break;
						case 2:
							//              	printf("M_WM(%d,%d) branch 2: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-1);
							insert_node (i+2, best_k, M_WM);
							insert_node (best_k+1, j-1, M_WM);
							break;
						case 3:
							//              	printf("M_WM(%d,%d) branch 3: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-2);
							insert_node (i+1, best_k, M_WM);
							insert_node (best_k+1, j-2, M_WM);
							break;
						case 4:
							//              	printf("M_WM(%d,%d) branch 4: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-2);
							insert_node (i+2, best_k, M_WM);
							insert_node (best_k+1, j-2, M_WM);
							break;
							// Hosna: June 28, 2007
							// the last branch of VM, which is WMB_(i+1),(j-1)
						case 5:
							//              	printf("M_WM(%d,%d) branch 3: pushing P_WMB(%d,%d)\n", i,j,i+1,j-1);
							insert_node(i+1,j-1, P_WMB);
							break;
					}
				}
					break;
			}
		}
			break;
		case FREE:
		{
			int j = cur_interval->j;

			if (j==0) return;

			int min = INF, tmp, best_row, i, best_i, acc, energy_ij;

			if (debug)
				printf ("\t(0,%d) FREE\n", j);

			// this case is for j unpaired, so I have to check that.
			if (fres[j].pair <= -1)
			{
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
				// pkonly condition:
				energy_ij = (fres[j].pair == i && fres[i].pair== j)? v->get_energy(i,j): INF;
				if (debug){
						printf("V(%d,%d)=%d \n",i,j,energy_ij);
				}
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

				if (fres[i].pair <= -1)
				{
					energy_ij = (fres[j].pair == i+1 && fres[i+1].pair== j) ? v->get_energy(i+1,j) : INF;
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
						tmp += dangle_bot [int_sequence[j]]
						[int_sequence[i+1]]
						[int_sequence[i]];
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
					energy_ij = (fres[j-1].pair == i && fres[i].pair== j-1)? v->get_energy(i,j-1) : INF;
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j-1]) + acc;
						tmp += dangle_top [int_sequence[j-1]]
						[int_sequence[i]]
						[int_sequence[j]];
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
					energy_ij = (fres[j-1].pair == i+1 && fres[i+1].pair == j-1) ? v->get_energy(i+1,j-1) : INF;
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
						tmp += dangle_bot [int_sequence[j-1]]
						[int_sequence[i+1]]
						[int_sequence[i]];
						tmp += dangle_top [int_sequence[j-1]]
						[int_sequence[i+1]]
						[int_sequence[j]];
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 4;
						}
					}
				}
			}

			// Hosna: June 28, 2007
			// the last branch of W, which is WMB_i,j
			//        energy_ij = WMB->get_energy(0,j);
			//        if (energy_ij < INF){
			//          	tmp = energy_ij + PS_penalty;
			//           	if (tmp < min){
			//           		min = tmp;
			//           		best_row = 5;
			//           	}
			//        }
			// Hosna June 30, 2007
			// The following would not take care of when
			// we have some unpaired bases before the start of the WMB
			for (i=0; i<=j-1; i++)
			{
				// Hosna: July 9, 2007
				// We only chop W to W + WMB when the bases before WMB are free
				if (i == 0 || (WMB->is_weakly_closed(0,i-1) && WMB->is_weakly_closed(i,j))){

					acc = (i-1>0) ? W[i-1]: 0;

					energy_ij = WMB->get_energy(i,j);

					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;

						if (tmp < min)
						{
							min = tmp;
							best_row = 5;
							best_i = i;
						}
					}

					// I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
					if (fres[i].pair <= -1 && i+1 < j)
					{
						energy_ij = WMB->get_energy(i+1,j);
						if (energy_ij < INF)
						{
							tmp = energy_ij + PS_penalty + acc;
							if (tmp < min)
							{
								min = tmp;
								best_row = 6;
								best_i = i;
							}
						}
					}

					// I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
					if (fres[j].pair <= -1 && i < j-1)
					{
						energy_ij = WMB->get_energy(i,j-1);
						if (energy_ij < INF)
						{
							tmp = energy_ij + PS_penalty + acc;
							if (tmp < min)
							{
								min = tmp;
								best_row = 7;
								best_i = i;
							}
						}
					}

					if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
					{
						energy_ij = WMB->get_energy(i+1,j-1);
						if (energy_ij < INF)
						{
							tmp = energy_ij + PS_penalty + acc;
							if (tmp < min)
							{
								min = tmp;
								best_row = 8;
								best_i = i;
							}
						}
					}
				}
			}
			switch (best_row)
			{
				case 0:
					//printf("W(%d) case 0: inserting Free (0,%d)\n",j,j-1);
					insert_node (0, j-1, FREE); break;
				case 1:
					//printf("W(%d) case 1: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, LOOP);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 2:
					//printf("W(%d) case 2: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, LOOP);
					if (best_i >= 0)// Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 3:
					//printf("W(%d) case 3: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, LOOP);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 4:
					//printf("W(%d) case 4: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, LOOP);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
					// Hosna: June 28, 2007
					// the last branch of W, which is WMB_i,j
				case 5:
					//printf("W(%d) case 5: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, P_WMB);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 6:
					//printf("W(%d) case 6: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 7:
					//printf("W(%d) case 7: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, P_WMB);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 8:
					//printf("W(%d) case 8: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
			}
		}
			break;
		case M_WM:
			//  else if(cur_interval->type == M_WM)
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			int tmp, min = INF;
			int best_k, best_row;

			if (debug)
				printf ("\t (%d,%d) M_WM\n", i,j);

			// the if clause added April 3, 2012
			if (fres[j].pair == i && fres[i].pair == j){
				tmp = v->get_energy(i,j) +
				AU_penalty (int_sequence[i], int_sequence[j]) +
				misc.multi_helix_penalty;
				if (tmp < min)
				{
					min = tmp;
					best_row = 1;
				}
			}
			if (fres[i].pair <= -1 && fres[j].pair == i+1 && fres[i+1].pair == j)
			{
				tmp = v->get_energy(i+1,j) +
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
			if (fres[j].pair <= -1 && fres[j-1].pair == i && fres[i].pair == j-1)
			{
				tmp = v->get_energy(i,j-1) +
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
			if (fres[i].pair <= -1 && fres[j].pair <= -1 && fres[j-1].pair == i+1 && fres[i+1].pair == j-1)
			{
				tmp = v->get_energy(i+1,j-1) +
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

			// TODO: April 3, 2012
			// do I need to change WM to pk_only as well?
			if (fres[i].pair <= -1)
			{
				tmp = vm->get_energy_WM (i+1,j) + misc.multi_free_base_penalty;
				if (tmp < min)
				{
					min = tmp;
					best_row = 5;
				}
			}
			if (fres[j].pair <= -1)
			{
				tmp = vm->get_energy_WM (i,j-1) + misc.multi_free_base_penalty;
				if (tmp < min)
				{
					min = tmp;
					best_row = 6;
				}
			}

			for (int k=i; k < j; k++)
			{
				tmp = vm->get_energy_WM (i, k) + vm->get_energy_WM (k+1, j);
				if (tmp < min)
				{
					min = tmp;
					best_k = k;
					best_row = 7;
				}
			}
			// Hosna: June 28, 2007
			// the last branch of W, which is WMB_i,j
			tmp = WMB->get_energy(i,j)+PSM_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 8;
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
					// Hosna: June 28, 2007
					// the last branch of W, which is WMB_i,j
				case 8:
					insert_node(i,j,P_WMB);
					break;
			}
		}
			break;
			// Hosna: Feb 19th 2007
		case P_WMB:
			// else if(cur_interval->type == P_WMB)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WMBP:
			// Hosna: April 18th, 2007
			// changed WMB to case 2 and WMBP
			// else if(cur_interval->type == P_WMBP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VP:
			//else if(cur_interval->type == P_VP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VPP:
			//else if(cur_interval->type == P_VPP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WI:
			//else if(cur_interval->type == P_WI)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_BE:
			//else if(cur_interval->type == P_BE)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WIP:
			//else if(cur_interval->type == P_WIP)
		{
			WMB->set_stack_interval(stack_interval);
			// Hosna, May 1st, 2012
			// changed to back_track_pkonly instead of back_track to accommodate pkonly case
			WMB->back_track_pkonly(structure,f,cur_interval);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		default:
			printf("Should not be here!\n");
	}


}

void W_final::print_result ()
// PRE:  The matrix V has been calculated and the results written in f
// POST: Prints details of each elementary structure
{
    int i;
    int energy = INF, sum;

    printf ("Minimum energy: %d\n", W[nb_nucleotides-1]);
    sum = 0;

    for (i=0; i< nb_nucleotides; i++)
    {
        if (f[i].pair > i)
        {
			//Hosna March 8, 2012
			// changing nested ifs to switch for optimality
			switch (f[i].type){
				case HAIRP:
				//if (f[i].type == HAIRP)
					energy = V->get_energy(i, f[i].pair);
					break;
				case STACK:
				//else if (f[i].type == STACK)
					energy = V->get_energy(i, f[i].pair) - V->get_energy(i+1, f[i+1].pair);
					break;
				case P_VP:
				// Hosna: June 28th, 2007
				//else if (f[i].type == P_VP){
					energy = WMB->get_VP(i,f[i].pair);
					break;
				case P_VPP:
				//}else if(f[i].type == P_VPP){
					energy = WMB->get_VPP(i,f[i].pair);
				//}
					break;
			}
            printf ("Pair (%d,%d), type %c,\tenergy %6d\n", i, f[i].pair, f[i].type, energy);
            sum += energy;
        }
    }
    printf ("0....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8\n");
    printf ("%s\n", sequence);
    printf ("%s\n", structure);

}
