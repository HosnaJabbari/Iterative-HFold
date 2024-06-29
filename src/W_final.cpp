#include "W_final.hh"
#include "h_struct.hh"
#include "h_externs.hh"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>


// Hosna June 20th, 2007
// calls the constructor for s_min_folding
// to create all the matrixes required for simfold
// and then calls allocate_space in here to allocate
// space for WMB and V_final
W_final::W_final(std::string seq,std::string res,bool pk_free, bool pk_only, int dangle) : params_(scale_parameters())
{
	seq_ = seq;
	this->res = res;
	this->n = seq.length();
	make_pair_matrix();
	params_->model_details.dangles = dangle;
    S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);
	this->pk_free = pk_free;
	this->pk_only = pk_only;
	W.resize(n+1,0);
	space_allocation();
}


W_final::~W_final()
{
	delete WMB;
	delete V;
	delete [] f;
	free(params_);
	free(S_);
	free(S1_);
}

// Hosna June 20th, 2007
// allocates space for WMB object and V_final
void W_final::space_allocation(){

	// From simfold
	f = new minimum_fold [n+1];

    V = new s_energy_matrix (seq_, n,S_,S1_,params_);
	structure = std::string (n+1,'.');

	// Hosna: June 20th 2007
    WMB = new pseudo_loop (seq_,res,V,S_,S1_,params_);

}

energy_t W_final::hfold(sparse_tree &tree){

		for (int i = n; i >=1; --i)
		{	
			for (int j =i; j<=n; ++j)//for (i=0; i<=j; i++)
			{
				const bool evaluate = tree.weakly_closed(i,j);
				const pair_type ptype_closing = pair[S_[i]][S_[j]];
				const bool restricted = tree.tree[i].pair == -1 || tree.tree[j].pair == -1;
				const bool paired = (tree.tree[i].pair == j && tree.tree[j].pair == i);

				const bool pkonly = (!pk_only || paired);

				if(ptype_closing> 0 && evaluate && !restricted && pkonly)
				V->compute_energy_restricted (i,j,tree);

				if(!pk_free) WMB->compute_energies(i,j,tree);


				V->compute_WMv_WMp(i,j,WMB->get_WMB(i,j),tree.tree);
				V->compute_energy_WM_restricted(i,j,tree,WMB->WMB);
			}

		}

	for (cand_pos_t j= TURN+1; j <= n; j++){
		energy_t m1 = INF;
		energy_t m2 = INF;
		energy_t m3 = INF;
		if(tree.tree[j].pair < 0) m1 = W[j-1];
		
		
		for (cand_pos_t k=1; k<=j-TURN-1; ++k){
		 	// m2 = compute_W_br2_restricted (j, fres, must_choose_this_branch);
			energy_t acc = (k>1) ? W[k-1]: 0;
			m2 = std::min(m2,acc + E_ext_Stem(V->get_energy(k,j),V->get_energy(k+1,j),V->get_energy(k,j-1),V->get_energy(k+1,j-1),S_,params_,k,j,n,tree.tree));
			if (k == 1 || (tree.weakly_closed(1,k-1) && tree.weakly_closed(k,j))) m3 = std::min(m3,acc + WMB->get_WMB(k,j) + PS_penalty);
			}
		W[j] = std::min({m1,m2,m3});
	}

    energy_t energy = W[n];

    // backtrack
    // first add (1,n) on the stack
    stack_interval = new seq_interval;
    stack_interval->i = 1;
    stack_interval->j = n;
    stack_interval->energy = W[n];
    stack_interval->type = FREE;
    stack_interval->next = NULL;

    seq_interval *cur_interval = stack_interval;

    while ( cur_interval != NULL)
    {
        stack_interval = stack_interval->next;
        backtrack_restricted (cur_interval,tree);
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }
	this->structure = structure.substr(1,n);
    return energy;

}

/**
 * @brief Gives the W(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * Until the changes to fres, I am adding +1 to the ptype closing and Si and Sj's to make them match - Mateo 2024
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t W_final::E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, paramT* params, const cand_pos_t i,const cand_pos_t j, cand_pos_t n, std::vector<Node> &tree){

	energy_t e = INF,en = INF;
  	pair_type tt  = pair[S[i]][S[j]];
	
    if ((tree[i].pair <-1 && tree[j].pair <-1) || (tree[i].pair == j && tree[j].pair == i)) {
				en = vij; // i j

				if (en != INF) {
					if (params->model_details.dangles == 2){
						base_type si1 = i>1 ? S[i-1] : -1;
                		base_type sj1 = j<n ? S[j+1] : -1;
                        en += vrna_E_ext_stem(tt, si1, sj1, params);
					}
                    else{
                        en += vrna_E_ext_stem(tt, -1, -1, params);
					}

                    e = MIN2(e, en);
					
				}

	}

	if(params->model_details.dangles  == 1){
        tt  = pair[S[i+1]][S[j]];
        if (((tree[i+1].pair <-1 && tree[j].pair <-1) || (tree[i+1].pair == j)) && tree[i].pair<0) {
            en = (j-i-1>TURN) ? vi1j : INF; //i+1 j

            if (en != INF) {

                base_type si1 = S[i];
                en += vrna_E_ext_stem(tt, si1, -1, params);
            }

            e = MIN2(e,en);

        }
        tt  = pair[S[i]][S[j-1]];
        if (((tree[i].pair <-1 && tree[j-1].pair <-1) || (tree[i].pair == j-1)) && tree[j].pair<0) {
            en = (j-1-i>TURN) ? vij1 : INF; // i j-1
            if (en != INF) {

                base_type sj1 = S[j];

                en += vrna_E_ext_stem(tt, -1, sj1, params);
            }
            e = MIN2(e,en);

        }
        tt  = pair[S[i+1]][S[j-1]];
        if (((tree[i+1].pair <-1 && tree[j-1].pair <-1) || (tree[i+1].pair == j-1)) && tree[i].pair < 0 && tree[j].pair<0) {
            en = (j-1-i-1>TURN) ? vi1j1 : INF; // i+1 j-1

            if (en != INF) {

                base_type si1 = S[i];
                base_type sj1 = S[j];

                en += vrna_E_ext_stem(tt, si1, sj1, params);
            }
            e = MIN2(e,en);
        }
	}
	return e;
}

void W_final::backtrack_restricted(seq_interval *cur_interval, sparse_tree &tree){
    char type;


	// printf("type is %c and i is %d and j is %d\n",cur_interval->type,cur_interval->i,cur_interval->j);
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
			// otherwise make it '[' and ']' -- changed to () if pseudoknot-free and [] if pseudoknotted -Mateo
			structure[i] = '(';
			structure[j] = ')';		

			type = V->get_type (i,j);
			
			// Hosna, March 8, 2012
			// changing nested ifs to switch for optimality
			switch (type){
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
					cand_pos_t best_ip=j, best_jp=i;
					energy_t min = INF;
					cand_pos_t max_ip = std::min(j-TURN-2,i+MAXLOOP+1);
					for (cand_pos_t k = i+1; k <= max_ip; ++k)
					{
						if (tree.up[k-1]>=(k-i-1)){
							cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
							for (cand_pos_t l = j-1; l >= min_l; --l)
							{
								
								if(tree.up[j-1]>=(j-l-1)){
							
									energy_t tmp = V->compute_int(i,j,k,l,params_);
									if (tmp < min)
									{
										min = tmp;
										best_ip = k;
										best_jp = l;
									}
								}
							}
						}
					}

					if (best_ip < best_jp)
						insert_node (best_ip, best_jp, LOOP);
					else
					{
						fprintf(stderr,"NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
						exit (0);
					}
				}
					break;
				case MULTI:
			//else if (type == MULTI)
				{
					f[i].type = MULTI;
					f[j].type = MULTI;
					int best_k = -1, best_row = -1;
					int tmp= INF, min = INF;
					for (cand_pos_t k = i+1; k <= j-1; k++){
						tmp = V->get_energy_WM (i+1,k-1) + std::min(V->get_energy_WMv(k, j-1),V->get_energy_WMp(k, j-1)) + E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_);
							
						if (tmp < min)
						  {
							min = tmp;
							best_k = k;
							best_row = 1;
						  }
						  // TODO:
						  // Hosna, May 1st, 2012
						  // do I need to check for non-canonical base pairings here as well so the dangle values not be INF??
						if (tree.tree[i+1].pair <= -1)
						{
							tmp = V->get_energy_WM (i+2,k-1) + std::min(V->get_energy_WMv(k, j-1),V->get_energy_WMp(k, j-1)) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_);
							
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 2;
							}
						}
						if (tree.tree[j-1].pair <= -1)
						{
							tmp = V->get_energy_WM (i+1,k-1) + std::min(V->get_energy_WMv(k, j-2),V->get_energy_WMp(k, j-2)) + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_);
							
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 3;
							}
						}
						if (tree.tree[i+1].pair <= -1 && tree.tree[j-1].pair <= -1)
						{
							tmp = V->get_energy_WM (i+2,k-1) + std::min(V->get_energy_WMv(k, j-2),V->get_energy_WMp(k, j-2));
							
							tmp = tmp + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_);
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 4;
							}
						}

						tmp = static_cast<energy_t>((k-i-1)*params_->MLbase + V->get_energy_WMp(k,j-1))+ E_MLstem(pair[S_[j]][S_[i]],-1,-1,params_);
						if (tmp < min)
						  {
							min = tmp;
							best_k = k;
							best_row = 5;
						  }
						  // TODO:
						  // Hosna, May 1st, 2012
						  // do I need to check for non-canonical base pairings here as well so the dangle values not be INF??
						if (tree.tree[i+1].pair <= -1)
						{
							if((k-(i+1)-1) >=0) tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + V->get_energy_WMp(k,j-1) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],-1,params_);
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 6;
							}
						}
						if (tree.tree[j-1].pair <= -1)
						{
							tmp = static_cast<energy_t>((k-i-1)*params_->MLbase) + V->get_energy_WMp(k,j-2) + E_MLstem(pair[S_[j]][S_[i]],-1,S_[i+1],params_);
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 7;
							}
						}
						if (tree.tree[i+1].pair <= -1 && tree.tree[j-1].pair <= -1)
						{
							if((k-(i+1)-1) >=0) tmp = static_cast<energy_t>((k-(i+1)-1)*params_->MLbase) + V->get_energy_WMp(k,j-2) + E_MLstem(pair[S_[j]][S_[i]],S_[j-1],S_[i+1],params_);
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 8;
							}
						}						
					  }
					switch (best_row)
					  {
					  case 1:
						insert_node (i+1, best_k-1, M_WM);
						insert_node (best_k, j-1, M_WM);
						break;
					  case 2:
						insert_node (i+2, best_k-1, M_WM);
						insert_node (best_k, j-1, M_WM);
						break;
					  case 3:
		             	// printf("M_WM(%d,%d) branch 3: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k-1,best_k,j-2);
						insert_node (i+1, best_k-1, M_WM);
						insert_node (best_k, j-2, M_WM);
						break;
					  case 4:
		             	// printf("M_WM(%d,%d) branch 4: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-2);
						insert_node (i+2, best_k-1, M_WM);
						insert_node (best_k, j-2, M_WM);
						break;
					  case 5:
						insert_node (best_k, j-1, M_WM);
						break;
					  case 6:
						insert_node (best_k, j-1, M_WM);
						break;
					  case 7:
						insert_node (best_k, j-2, M_WM);
						break;
					  case 8:
						insert_node (best_k, j-2, M_WM);
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

			if (j==1) return;

			int min = INF, tmp, best_row, i, best_i, acc, energy_ij;

			// this case is for j unpaired, so I have to check that.
			if (tree.tree[j].pair <= -1)
			{
				tmp = W[j-1];
				if (tmp < min)
				{
					min = tmp;
					best_row = 0;
				}
			}
			for (i=1; i<=j-1; i++)    // no TURN
			{

				// Don't need to make sure i and j don't have to pair with something else
				//  it's INF, done in fold_sequence_restricted
				acc = (i>1) ? W[i-1] : 0;
				energy_ij = V->get_energy(i,j);

				if (energy_ij < INF)
				{	
					if(params_->model_details.dangles == 2){
						base_type si1 = i>1 ? S_[i-1] : -1;
						base_type sj1 = j<n ? S_[j+1] : -1;
						tmp = energy_ij + E_ExtLoop(pair[S_[i]][S_[j]],si1,sj1,params_) + acc;
					} else 
						tmp = energy_ij + E_ExtLoop(pair[S_[i]][S_[j]],-1,-1,params_) + acc; 
					if (tmp < min)
					{
					min = tmp;
					best_i = i;
					best_row = 1;
					}
					
				}
				if(params_->model_details.dangles ==1){
					if (tree.tree[i].pair <= -1)
					{
						energy_ij = V->get_energy(i+1,j);
						if (energy_ij < INF)
						{
							tmp = energy_ij + E_ExtLoop(pair[S_[i+1]][S_[j]],S_[i],-1,params_) + acc;
							
							if (tmp < min)
							{
								min = tmp;
								best_i = i;
								best_row = 2;
							}
							
						}
					}
					if (tree.tree[j].pair <= -1)
					{
						energy_ij = V->get_energy(i,j-1);
						if (energy_ij < INF)
						{
							tmp = energy_ij + E_ExtLoop(pair[S_[i]][S_[j-1]],-1,S_[j],params_) + acc;
							if (tmp < min)
							{
								min = tmp;
								best_i = i;
								best_row = 3;
							}
						}
					}
					if (tree.tree[i].pair <= -1 && tree.tree[j].pair <= -1)
					{
						energy_ij = V->get_energy(i+1,j-1);
						if (energy_ij < INF)
						{
							tmp = energy_ij + E_ExtLoop(pair[S_[i+1]][S_[j-1]],S_[i],S_[j],params_) + acc;
							if (tmp < min)
							{
								min = tmp;
								best_i = i;
								best_row = 4;
							}
						}
					}
				}
			}
		// Hosna June 30, 2007
		// The following would not take care of when
		// we have some unpaired bases before the start of the WMB
		for (i=1; i<=j-1; i++)
		{
			// Hosna: July 9, 2007
			// We only chop W to W + WMB when the bases before WMB are free
			if (i == 1 || (tree.weakly_closed(1,i-1) && tree.weakly_closed(i,j))){

				acc = (i-1>0) ? W[i-1]: 0;

				energy_ij = WMB->get_WMB(i,j);

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

				if (tree.tree[i].pair <= -1 && i+1 < j)
				{
					energy_ij = WMB->get_WMB(i+1,j);
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

				if (tree.tree[j].pair <= -1 && i < j-1)
				{
					energy_ij = WMB->get_WMB(i,j-1);
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

				if (tree.tree[i].pair <= -1 && tree.tree[j].pair <= -1 && i+1 < j-1)
				{
					energy_ij = WMB->get_WMB(i+1,j-1);
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
					insert_node (1, j-1, FREE); break;
				case 1:
					//printf("W(%d) case 1: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, LOOP);
					if (best_i-1 > 1)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (1, best_i-1, FREE);
					break;
				case 2:
					//printf("W(%d) case 2: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, LOOP);
					if (best_i >= 1)// Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (1, best_i, FREE);
					break;
				case 3:
					//printf("W(%d) case 3: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, LOOP);
					if (best_i-1 > 1)
						insert_node (1, best_i-1, FREE);
					break;
				case 4:
					//printf("W(%d) case 4: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, LOOP);
					if (best_i >= 1) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (1, best_i, FREE);
					break;
				// Hosna: June 28, 2007
				// the last branch of W, which is WMB_i,j
				case 5:
					//printf("W(%d) case 5: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, P_WMB);
					if (best_i-1 > 1)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (1, best_i-1, FREE);
					break;
				case 6:
					//printf("W(%d) case 6: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, P_WMB);
					if (best_i >= 1) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (1, best_i, FREE);
					break;
				case 7:
					//printf("W(%d) case 7: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, P_WMB);
					if (best_i-1 > 1)
						insert_node (1, best_i-1, FREE);
					break;
				case 8:
					//printf("W(%d) case 8: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, P_WMB);
					if (best_i >= 1) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (1, best_i, FREE);
					break;
			}
		}
			break;
		case M_WM:
		{
			  cand_pos_t i = cur_interval->i;
			  cand_pos_t j = cur_interval->j;
			  energy_t min = INF;
			  cand_pos_t best_k = j, best_row;

			  if(tree.tree[j].pair<0){
				if(V->get_energy_WM(i,j-1)< min){
					min = V->get_energy_WM(i,j-1);
					best_row = 5;
				}
			}

			  for (cand_pos_t k=i; k <= j-TURN-1; k++)
				{	energy_t m1 = INF,m2 = INF;
					bool can_pair = tree.up[k-1] >= (k-(i));
					if(can_pair) m1 = static_cast<energy_t>((k-i)*params_->MLbase) + V->get_energy_WMv (k, j);
					if (m1 < min){
						min = m1;
						best_k = k;
						best_row = 1;
					}
					if(can_pair) m2 = static_cast<energy_t>((k-i)*params_->MLbase) + V->get_energy_WMp(k, j);
					if (m2 < min){
						min = m2;
						best_k = k;
						best_row = 2;
					}
					energy_t m3 = V->get_energy_WM (i, k-1) + V->get_energy_WMv (k, j);
					if (m3 < min){
						min = m3;
						best_k = k;
						best_row = 3;
					}
					energy_t m4 = V->get_energy_WM (i, k-1) + V->get_energy_WMp(k, j);
					if (m4 < min){
						min = m4;
						best_k = k;
						best_row = 4;
					}
				}
			  switch (best_row)
				{
				  case 1: insert_node (best_k, j, M_WMv); break;
				  case 2: insert_node (best_k, j, M_WMp); break;
				  case 3:
					  insert_node (i, best_k-1, M_WM);
					  insert_node (best_k, j, M_WMv);
					break;
				  case 4:
					  insert_node (i, best_k-1, M_WM);
					  insert_node (best_k+1, j, M_WMp);
					break;
				  case 5:
				  	  insert_node (i,j-1,M_WM); break;
				  }
			}
			break;
		case M_WMv:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = cur_interval->j;
			energy_t min = INF;
			cand_pos_t best_row;
			cand_pos_t si = S_[i];
			cand_pos_t sj = S_[j];
			cand_pos_t si1 = (i>1) ? S_[i-1] : -1;
			cand_pos_t sj1 = (j<n) ? S_[j+1] : -1;
			pair_type tt = pair[S_[i]][S_[j]];
			min = V->get_energy(i,j) + ((params_->model_details.dangles == 2) ? E_MLstem(tt,si1,sj1,params_) : E_MLstem(tt,-1,-1,params_));
			best_row = 1;
			if(params_->model_details.dangles == 1){
				if(tree.tree[i].pair<0){
					tt = pair[S_[i+1]][S_[j]];
					energy_t tmp = V->get_energy(i+1,j) + E_MLstem(tt,si,-1,params_);
					if(tmp<min){
						min = tmp;
						best_row = 2;
					}
				}
				if(tree.tree[j].pair<0){
					tt = pair[S_[i]][S_[j-1]];
					energy_t tmp = V->get_energy(i,j-1) + E_MLstem(tt,-1,sj,params_);
					if(tmp<min){
						min = tmp;
						best_row = 3;
					}
				}
				if(tree.tree[i].pair<0 && tree.tree[j].pair<0){
					tt = pair[S_[i+1]][S_[j-1]];
					energy_t tmp = V->get_energy(i+1,j-1) + E_MLstem(tt,si,sj,params_);
					if(tmp<min){
						min = tmp;
						best_row = 4;
					}
				}
			}
			if(tree.tree[j].pair<0){
				if(V->get_energy_WMv(i,j-1)< min){
					min = V->get_energy_WMv(i,j-1);
					best_row = 5;
				}
			}
			switch (best_row){
				case 1: insert_node (i, j, LOOP); break;
				case 2: insert_node (i+1, j, LOOP); break;
				case 3: insert_node (i, j-1, LOOP); break;
				case 4: insert_node (i+1, j-1, LOOP); break;
				case 5: insert_node (i, j-1, M_WMv); break;
			}
		}
		break;
		case M_WMp:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			int min = INF;
			int best_row;

			min = WMB->get_WMB(i,j) + PSM_penalty + b_penalty;
			best_row = 1;
			if(tree.tree[j].pair<0){
				if(V->get_energy_WMp(i,j-1)< min){
					min = V->get_energy_WMp(i,j-1);
					best_row = 2;
				}
			}
			switch (best_row){
				case 1: insert_node (i, j, P_WMB); break;
				case 2: insert_node (i, j-1, M_WMp); break;
			}
		}
		break;
		case P_WMB:
		case P_WMBP:
		case P_WMBW:
		case P_VP:
		case P_VPR:
		case P_VPL:
		case P_WI:
		case P_BE:
		case P_WIP:
		{
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(structure,f,cur_interval,tree);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		default:
			printf("Should not be here!\n");
	}


}

void W_final::insert_node (int i, int j, char type)
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


//Mateo 13 Sept 2023
//return number of bases in between the two inclusive index
int distance(int left, int right){
    return (right-left-1);
}

//Mateo 13 Sept 2023
//given a initial hotspot which is a hairpin loop, keep trying to add a arc to form a larger stack
void expand_hotspot(s_energy_matrix *V, Hotspot &hotspot, int n){
    //printf("\nexpanding hotspot: i: %d j: %d\n",hotspot->get_left_inner_index(),hotspot->get_right_inner_index());
    //calculation for the hairpin that is already in there
    V->compute_hotspot_energy(hotspot.get_left_outer_index(),hotspot.get_right_outer_index(),0);


    //try to expand by adding a arc right beside the current out most arc
    while(hotspot.get_left_outer_index()-1 >= 1 && hotspot.get_right_outer_index()+1 <= n){
		base_type sim1 = V->S_[hotspot.get_left_outer_index()-1];
		base_type sjp1 = V->S_[hotspot.get_right_outer_index()+1];
		pair_type ptype_closing = pair[sim1][sjp1];
        if(ptype_closing>0){
            hotspot.move_left_outer_index();
            hotspot.move_right_outer_index();
            hotspot.increment_size();
            V->compute_hotspot_energy(hotspot.get_left_outer_index(),hotspot.get_right_outer_index(),1);
        }else{
            break;
        }
    }
	base_type i = hotspot.get_left_outer_index();
	base_type j = hotspot.get_right_outer_index();
	pair_type tt = pair[V->S_[i]][V->S_[j]];
	base_type si1 = i>1 ? V->S_[i-1] : -1;
	base_type sj1 = j<=n ? V->S_[j+1] : -1;
	energy_t dangle_penalty = vrna_E_ext_stem(tt, si1, sj1, V->params_);


    double energy = V->get_energy(hotspot.get_left_outer_index(),hotspot.get_right_outer_index());

    // printf("here and %d\n",energy);
    //printf("energy: %lf, AU_total: %d, dangle_top_total: %d, dangle_bot_total: %d\n",energy,non_gc_penalty,dangle_top_penalty,dangle_bot_penalty);

    energy = (energy + dangle_penalty) / 100;

    hotspot.set_energy(energy);
    //printf("done: %d %d %d %d\n",hotspot->get_left_outer_index(),hotspot->get_left_inner_index(),hotspot->get_right_inner_index(),hotspot->get_right_outer_index());
    return;
}

//Mateo 13 Sept 2023
//look for every possible hairpin loop, and try to add a arc to form a larger stack with at least min_stack_size bases
void get_hotspots(std::string seq,std::vector<Hotspot> &hotspot_list,int max_hotspot, vrna_param_s *params){
    
	int n = seq.length();
	s_energy_matrix *V;
	make_pair_matrix();
	short *S_ = encode_sequence(seq.c_str(),0);
	short *S1_ = encode_sequence(seq.c_str(),1);
	V = new s_energy_matrix (seq,n,S_,S1_,params);
    int min_bp_distance = 3;
    int min_stack_size = 3; //the hotspot must be a stack of size >= 3
    // Hotspot current_hotspot;
    //start at min_stack_size-1 and go outward to try to add more arcs to form bigger stack because we cannot expand more than min_stack_size from there anyway
    for(int i = min_stack_size; i <= n; i++){
        for(int j = i; j <= n; j++){
			int ptype_closing = pair[V->S_[i]][V->S_[j]];
            if(ptype_closing>0 && distance(i,j) >= min_bp_distance){
                // current_hotspot = new Hotspot(i,j,nb_nucleotides);
				
                Hotspot current_hotspot(i,j,n);

                expand_hotspot(V,current_hotspot,n);


                if(current_hotspot.get_size() < min_stack_size || current_hotspot.is_invalid_energy()){

                }else{
                    
                    current_hotspot.set_structure();
                    hotspot_list.push_back(current_hotspot);

                }
            }
        }
    }

    //make sure we only keep top 20 hotspot with lowest energy
    std::sort(hotspot_list.begin(), hotspot_list.end(),compare_hotspot_ptr);
    while(hotspot_list.size() > max_hotspot){
        hotspot_list.pop_back();
    }

    //if no hotspot found, add all _ as restricted
    if(hotspot_list.size() == 0){
        Hotspot hotspot(1,n,n+1);
        hotspot.set_default_structure();
        hotspot_list.push_back(hotspot);
    }
	delete V;
	free(S_);
	free(S1_);

    return;
}

bool compare_hotspot_ptr(Hotspot &a, Hotspot &b) { 
    return (a.get_energy() < b.get_energy()); 
}