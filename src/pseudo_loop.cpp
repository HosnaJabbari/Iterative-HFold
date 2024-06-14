#include "pseudo_loop.hh"
#include "h_externs.hh"
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>

pseudo_loop::pseudo_loop(std::string seq, std::string res, s_energy_matrix *V, short *S, short *S1, vrna_param_t *params)
{
	this->seq = seq;
	this->res = res;
	this->V = V;
	S_ = S;
	S1_ = S1;
	params_ = params;
	make_pair_matrix();
    allocate_space();
}

void pseudo_loop::allocate_space()
{
    n = seq.length();

    index.resize(n+1);
    cand_pos_t total_length = ((n+1) *(n+2))/2;
    index[1] = 0;
    for (cand_pos_t i=2; i <= n; i++)
        index[i] = index[i-1]+(n+1)-i+1;

    WI.resize(total_length,0);

    VP.resize(total_length,INF);

	VPL.resize(total_length,INF);

	VPR.resize(total_length,INF);

    WMB.resize(total_length,INF);

	WMBW.resize(total_length,INF);

    WMBP.resize(total_length,INF);

    WIP.resize(total_length,INF);

    BE.resize(total_length,0);

}

pseudo_loop::~pseudo_loop()
{
}

void pseudo_loop::compute_energies(cand_pos_t i, cand_pos_t j, sparse_tree &tree)
{
	cand_pos_t ij = index[i]+j-i;
	const pair_type ptype_closing = pair[S_[i]][S_[j]];
	bool weakly_closed_ij = tree.weakly_closed(i,j);
	// base cases:
	// a) i == j => VP[ij] = INF
	// b) [i,j] is a weakly_closed region => VP[ij] = INF
	// c) i or j is paired in original structure => VP[ij] = INF
	if ((i == j || j-i<4 || weakly_closed_ij))	{
		VP[ij] = INF;
		VPL[ij] = INF;
		VPR[ij] = INF;
	}
	else{
		if(ptype_closing>0 && tree.tree[i].pair < -1 && tree.tree[j].pair < -1) compute_VP(i,j,tree);
		
		if(tree.tree[j].pair < -1) compute_VPL(i,j,tree);

		if(tree.tree[j].pair < j) compute_VPR(i,j,tree);
	}

	if (!((j-i-1) <= TURN || (tree.tree[i].pair >= -1 && tree.tree[i].pair > j) || (tree.tree[j].pair >= -1 && tree.tree[j].pair < i) || (tree.tree[i].pair >= -1 && tree.tree[i].pair < i ) || (tree.tree[j].pair >= -1 && j < tree.tree[j].pair))){
		compute_WMBW(i,j,tree);
		
		compute_WMBP(i,j,tree);

		compute_WMB(i,j,tree);
	}

	if(!weakly_closed_ij){
		WI[ij] = INF;
		WIP[ij] = INF;
	}
	else{
		compute_WI(i,j,tree);
		compute_WIP(i,j,tree);
	}

	cand_pos_t ip = tree.tree[i].pair; // i's pair ip should be right side so ip = )
	cand_pos_t jp = tree.tree[j].pair; // j's pair jp should be left side so jp = (

	compute_BE(i,ip,jp,j,tree);

}
// Added +1 to fres/tree indices as they are 1 ahead at the moment
void pseudo_loop::compute_WI(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	energy_t min = INF, m1 = INF, m2= INF, m3= INF, m4= INF, m5 = INF;
	cand_pos_t ij = index[i]+j-i;
	// branch 4, one base
	if (i == j){
		WI[ij] = PUP_penalty;
		return;
	}
	
	for (cand_pos_t k = i+1; k < j-TURN-1; ++k){
		energy_t wi_1 = get_WI(i,k-1);
		energy_t v_energy = wi_1 + V->get_energy(k,j);
		energy_t wmb_energy = wi_1 + get_WMB(k,j);
		m1 = std::min(m1,v_energy);
		m2 = std::min(m2,wmb_energy);
	}
	m1 += PPS_penalty;
	m2 += PSP_penalty + PPS_penalty;
	if (tree.tree[j].pair < 0) m3 = get_WI(i,j-1) + PUP_penalty; 
	m4 = V->get_energy(i,j) + PPS_penalty;
	m5 = get_WMB(i,j) + PSP_penalty + PPS_penalty;

	WI[ij] = std::min({m1,m2,m3,m4,m5});
}

void pseudo_loop::compute_WIP(cand_pos_t  i, cand_pos_t  j, sparse_tree &tree){
	cand_pos_t ij = index[i]+j-i;

	energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF, m6 = INF, m7 = INF;

	// branch 1:
	for (cand_pos_t k = i+1; k < j-TURN-1; ++k){
		bool can_pair = tree.up[k-1] >= (k-i);
		energy_t wi_1 = get_WIP(i,k-1);
		energy_t v_energy = V->get_energy(k,j);
		energy_t wmb_energy = get_WMB(k,j);
		m1 = std::min(m1,wi_1 + v_energy);
		m2 = std::min(m2,wi_1 + wmb_energy);
		if(can_pair) m3 = std::min(m3, static_cast<energy_t>((k-i)*cp_penalty) + v_energy);
		if(can_pair) m4 = std::min(m4, static_cast<energy_t>((k-i)*cp_penalty) + wmb_energy);
	}
	m1 += bp_penalty;
	m2 += PSM_penalty + bp_penalty;
	m3 += bp_penalty;
	m4 += PSM_penalty + bp_penalty;
	// branch 2:
	if (tree.tree[j].pair < 0){
		m5 = get_WIP(i,j-1) + cp_penalty;
	}
	m6 = V->get_energy(i,j) + bp_penalty;
	m7 = get_WMB(i,j) + PSM_penalty + bp_penalty;

	WIP[ij] = std::min({m1,m2,m3,m4,m5,m6,m7});

}

void pseudo_loop::compute_VPL(cand_pos_t i, cand_pos_t j, sparse_tree &tree){

	cand_pos_t ij = index[i]+j-i;
	energy_t m1 = INF;

	cand_pos_t min_Bp_j = std::min((cand_pos_tu) tree.b(i,j), (cand_pos_tu) tree.Bp(i,j));
	for(cand_pos_t k = i+1; k<min_Bp_j; ++k){
		bool can_pair = tree.up[k-1] >= (k-i);
		if(can_pair) m1 = std::min(m1, static_cast<energy_t>((k-i)*cp_penalty) + get_VP(k,j));
	}


	VPL[ij] = m1;

}

void pseudo_loop::compute_VPR(cand_pos_t i, cand_pos_t j, sparse_tree &tree){

	cand_pos_t ij = index[i]+j-i;
	energy_t m1 = INF, m2 = INF;

	cand_pos_t max_i_bp = std::max(tree.B(i,j),tree.bp(i,j));

	for(cand_pos_t k = max_i_bp+1; k<j; ++k){
		energy_t VP_energy = get_VP(i,k);
		bool can_pair = tree.up[j-1] >= (j-k);
		m1 = std::min(m1, VP_energy + get_WIP(k+1,j));
		if(can_pair) m2 = std::min(m2,VP_energy + static_cast<energy_t>((j-k)*cp_penalty));
	}

	VPR[ij] = std::min(m1,m2);
}


void pseudo_loop::compute_VP(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	cand_pos_t ij = index[i]+j-i;

	const pair_type ptype_closing = pair[S_[i]][S_[j]];	
	
	energy_t m1 = INF, m2 = INF, m3 = INF, m4= INF, m5 = INF, m6 = INF, m7 = INF, m8 = INF, m9 = INF; //different branches
	
	// Borders -- added one to i and j to make it fit current bounds but also subtracted 1 from answer as the tree bounds are shifted as well
	cand_pos_t Bp_ij = tree.Bp(i,j);
	cand_pos_t B_ij = tree.B(i,j);
	cand_pos_t b_ij = tree.b(i,j);
	cand_pos_t bp_ij = tree.bp(i,j);
	
	//branchs:
	// 1) inArc(i) and NOT_inArc(j)
	// WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(j-1)

	// Hosna April 9th, 2007
	// need to check the borders as they may be negative
	if((tree.tree[i].parent->index) > 0 && (tree.tree[j].parent->index) < (tree.tree[i].parent->index) && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0){
		energy_t WI_ipus1_BPminus = get_WI(i+1,Bp_ij - 1) ;
		energy_t WI_Bplus_jminus = get_WI(B_ij + 1,j-1);
		m1 =   WI_ipus1_BPminus + WI_Bplus_jminus;
	}

	// 2) NOT_inArc(i) and inArc(j)
	// WI(i+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)

	// Hosna April 9th, 2007
	// checking the borders as they may be negative
	if ((tree.tree[i].parent->index) < (tree.tree[j].parent->index) && (tree.tree[j].parent->index) > 0 && b_ij>= 0 && bp_ij >= 0 && Bp_ij < 0){
		energy_t WI_i_plus_b_minus = get_WI(i+1,b_ij - 1);
		energy_t WI_bp_plus_j_minus = get_WI(bp_ij + 1,j-1);
		m2 = WI_i_plus_b_minus + WI_bp_plus_j_minus;
	}

	// 3) inArc(i) and inArc(j)
	// WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)

	// Hosna April 9th, 2007
	// checking the borders as they may be negative
	if((tree.tree[i].parent->index) > 0 && (tree.tree[j].parent->index) > 0 && Bp_ij >= 0 && B_ij >= 0  && b_ij >= 0 && bp_ij>= 0){
		energy_t WI_i_plus_Bp_minus = get_WI(i+1,Bp_ij - 1);
		energy_t WI_B_plus_b_minus = get_WI(B_ij + 1,b_ij - 1);
		energy_t WI_bp_plus_j_minus = get_WI(bp_ij +1,j - 1);
		m3 = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
	}

	// 4) NOT_paired(i+1) and NOT_paired(j-1) and they can pair together
	// e_stP(i,i+1,j-1,j) + VP(i+1)(j-1)
	pair_type ptype_closingip1jm1 = pair[S_[i+1]][S_[j-1]];
	if((tree.tree[i+1].pair) < -1 && (tree.tree[j-1].pair) < -1 && ptype_closingip1jm1>0){
		m4 = get_e_stP(i,j)+ get_VP(i+1,j-1);
	}

	// 5) NOT_paired(r) and NOT_paired(rp)
	//  VP(i,j) = e_intP(i,ip,jp,j) + VP(ip,jp)
	// Hosna, April 6th, 2007
	// whenever we use get_borders we have to check for the correct values
	cand_pos_t min_borders = std::min((cand_pos_tu) Bp_ij, (cand_pos_tu) b_ij);
	cand_pos_t edge_i = std::min(i+MAXLOOP+1,j-TURN-1);
	min_borders = std::min({min_borders,edge_i});
//		printf("B'(%d,%d) = %d, b(%d,%d) = %d, min_borders = %d\n",i,j,get_Bp(i,j),i,j,get_b(i,j), min_borders);
	for (cand_pos_t k = i+1; k < min_borders; ++k){
		// Hosna: April 20, 2007
		// i and ip and j and jp should be in the same arc
		// also it should be the case that [i+1,ip-1] && [jp+1,j-1] are empty regions

		if (tree.tree[k].pair < -1 && (tree.up[(k)-1] >= ((k)-(i)-1))){
			// Hosna, April 6th, 2007
			// whenever we use get_borders we have to check for the correct values
			cand_pos_t max_borders = std::max(bp_ij,B_ij)+1;
			cand_pos_t edge_j = k+j-i-MAXLOOP-2;
			max_borders = std::max({max_borders,edge_j});
			for (cand_pos_t l = j-1; l > max_borders ; --l){

				pair_type ptype_closingkj = pair[S_[k]][S_[l]];
				if (tree.tree[l].pair < -1 && ptype_closingkj>0 && (tree.up[(j)-1] >= ((j)-(l)-1))){
					// Hosna: April 20, 2007
					// i and ip and j and jp should be in the same arc -- If it's unpaired between them, they have to be
					energy_t tmp = get_e_intP(i,k,l,j) + get_VP(k,l);
					m5 = std::min(m5,tmp);
					
				}
			}
		}
	}

		cand_pos_t min_Bp_j = std::min((cand_pos_tu) tree.b(i,j), (cand_pos_tu) tree.Bp(i,j));
		cand_pos_t max_i_bp = std::max(tree.B(i,j),tree.bp(i,j));

		for(cand_pos_t k = i+1; k<min_Bp_j; ++k){
			m6 = get_WIP(i+1,k-1) + get_VP(k,j-1);
		}
		
		m6 += ap_penalty + 2*bp_penalty;

		for(cand_pos_t k = max_i_bp+1; k<j; ++k){
			m7 = std::min(m7,get_VP(i+1,k) + get_WIP(k+1,j-1));
		}

		m7 += ap_penalty + 2*bp_penalty;

		for(cand_pos_t k = i+1; k<min_Bp_j; ++k){
			m8 = std::min(m8,get_WIP(i+1,k-1) + get_VPR(k,j-1));
		}

		m8 += ap_penalty + 2*bp_penalty;

		for(cand_pos_t k = max_i_bp+1; k<j; ++k){
			m9 = std::min(m9,get_VPL(i+1,k) + get_WIP(k+1,j-1));
		}

		m9 += ap_penalty + 2*bp_penalty;





	//finding the min energy
	energy_t vp_h = std::min({m1,m2,m3});
	energy_t vp_iloop = std::min({m4,m5});
	energy_t vp_split = std::min({m6,m7,m8,m9});
	energy_t min = std::min({vp_h,vp_iloop,vp_split});
	

	VP[ij] = min;

	
}

void pseudo_loop::compute_WMBW(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	cand_pos_t ij = index[i]+j-i;

	energy_t m1 = INF;

	if(tree.tree[j].pair < j){
		for(cand_pos_t l = i+1; l<j; l++){
			if (tree.tree[l].pair < 0 && tree.tree[l].parent->index > -1 && tree.tree[j].parent->index > -1 && tree.tree[j].parent->index == tree.tree[l].parent->index){
				energy_t tmp = get_WMBP(i,l) + get_WI(l+1,j);
				m1 = std::min(m1,tmp);
			}
		}
	}
	WMBW[ij] = m1;
}

void pseudo_loop::compute_WMBP(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	cand_pos_t ij = index[i]+j-i;

	energy_t m1 = INF, m2 = INF, m4 = INF;	

	// 1)
	if (tree.tree[j].pair < 0){
		energy_t tmp = INF;
		cand_pos_t b_ij = tree.b(i,j);
		for (cand_pos_t l = i+1; l<j ; l++)	{
			// Hosna, April 6th, 2007
			// whenever we use get_borders we have to check for the correct values
			cand_pos_t bp_il = tree.bp(i,l);
			cand_pos_t Bp_lj = tree.Bp(l,j);
			// Hosna: April 19th, 2007
			// the chosen l should be less than border_b(i,j) -- should be greater than border_b(i,l)
			if(b_ij > 0 && l < b_ij){
				if (bp_il >= 0 && l>bp_il && Bp_lj > 0 && l<Bp_lj){ // bp(i,l) < l < Bp(l,j)
					cand_pos_t B_lj = tree.B(l,j);

					// Hosna: July 5th, 2007:
					// as long as we have i <= arc(l)< j we are fine
					if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l+TURN <=j){
						energy_t sum = get_BE(tree.tree[B_lj].pair,B_lj,tree.tree[Bp_lj].pair,Bp_lj,tree)+ get_WMBP(i,l-1)+ get_VP(l,j);
						tmp = std::min(tmp,sum);
					}
				}
			}
			
			m1 = 2*PB_penalty + tmp;
		}
	}
	// 2) WMB(i,j) = min_{i<l<j}{WMB(i,l)+WI(l+1,j)} if bp(j)<j
	// Hosna: Feb 5, 2007
	if (tree.tree[j].pair < 0){
		energy_t tmp = INF;
		cand_pos_t b_ij = tree.b(i,j);
		for (cand_pos_t l = i+1; l<j ; l++)	{
			// Hosna, April 6th, 2007
			// whenever we use get_borders we have to check for the correct values
			cand_pos_t bp_il = tree.bp(i,l);
			cand_pos_t Bp_lj = tree.Bp(l,j);
			// Hosna: April 19th, 2007
			// the chosen l should be less than border_b(i,j) -- should be greater than border_b(i,l)
			if(b_ij>0 && l<b_ij){
				if (bp_il >= 0 && l>bp_il && Bp_lj > 0 && l<Bp_lj){ // bp(i,l) < l < Bp(l,j)
					cand_pos_t B_lj = tree.B(l,j);

					// Hosna: July 5th, 2007:
					// as long as we have i <= arc(l)< j we are fine
					if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l+TURN <=j){
						energy_t sum = get_BE(tree.tree[B_lj].pair,B_lj,tree.tree[Bp_lj].pair,Bp_lj,tree)+ get_WMBW(i,l-1)+ get_VP(l,j);
						tmp = std::min(tmp,sum);
					}
				}
			}
			m2 = 2*PB_penalty + tmp;
		}
	}
	// 3) WMB(i,j) = VP(i,j) + P_b
	energy_t m3 = get_VP(i,j) + PB_penalty;

	// if not paired(j) and paired(i) then
	// WMBP(i,j) = 2*Pb + min_{i<l<bp(i)}(BE(i,bp(i),b'(i,l),bp(b'(i,l)))+WI(b'+1,l-1)+VP(l,j))
	if(tree.tree[j].pair < 0 && tree.tree[i].pair >= 0){
		energy_t tmp = INF;
		// Hosna: June 29, 2007
		// if j is inside i's arc then the l should be
		// less than j not bp(i)
		// check with Anne
		for (cand_pos_t l = i+1; l < j; l++){


			// Hosna, April 9th, 2007
			// checking the borders as they may be negative
			cand_pos_t bp_il = tree.bp(i,l);
			if(bp_il >= 0 && bp_il < n && l+TURN <= j){
				energy_t BE_energy = get_BE(i,tree.tree[i].pair,bp_il,tree.tree[bp_il].pair,tree);
				energy_t WI_energy = get_WI(bp_il +1,l-1);
				energy_t VP_energy = get_VP(l,j);
				energy_t sum = BE_energy + WI_energy + VP_energy;
				tmp = std::min(tmp,sum);
			}
		}
		m4 = 2*PB_penalty + tmp;
	}

		// get the min for WMB
		WMBP[ij] = std::min({m1,m2,m3,m4});
	


}

void pseudo_loop::compute_WMB(cand_pos_t  i, cand_pos_t  j, sparse_tree &tree){
	cand_pos_t ij = index[i]+j-i;
	//base case
	if (i == j){
		WMB[ij] = INF;
		return;
	}
	// Hosna: July 6th, 2007
	// added impossible cases
	energy_t m2 = INF, mWMBP = INF;
	// 2)
	if (tree.tree[j].pair >= 0 && j > tree.tree[j].pair && tree.tree[j].pair > i){
		cand_pos_t bp_j = tree.tree[j].pair;
		energy_t temp = INF;
		for (cand_pos_t l = (bp_j +1); (l < j); l++){
			// Hosna: April 24, 2007
			// correct case 2 such that a multi-pseudoknotted
			// loop would not be treated as case 2
			cand_pos_t Bp_lj = tree.Bp(l,j);

			if (Bp_lj >= 0 && Bp_lj<n){

				energy_t sum = get_BE(bp_j,j,tree.tree[Bp_lj].pair,Bp_lj,tree) + get_WMBP(i,l) + get_WI(l+1,Bp_lj-1);
				temp = std::min(temp,sum);
			}

		}
		m2 = PB_penalty + temp;
	}
	// check the WMBP value
	mWMBP =  get_WMBP(i,j);

	// get the min for WMB
	WMB[ij] = std::min(m2,mWMBP);
	
}

void pseudo_loop::compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, sparse_tree &tree){


    // Ian Wark July 19 2017
    // otherwise it will create pairs in spots where the restricted structure says there should be no pairs

	if (!( i >= 1 && i <= ip && ip < jp && jp <= j && j <= n && tree.tree[i].pair > 0 && tree.tree[j].pair > 0 && tree.tree[ip].pair > 0 && tree.tree[jp].pair > 0 && tree.tree[i].pair == j && tree.tree[j].pair == i && tree.tree[ip].pair == jp && tree.tree[jp].pair == ip)){ //impossible cases
		return;
	}
	cand_pos_t iip = index[i]+ip-i;
	// base case: i.j and ip.jp must be in G
	if (tree.tree[i].pair != j || tree.tree[ip].pair != jp){
		BE[iip] = INF;
		return;
	}

	// base case:
	if(i == ip && j == jp && i<j){
		BE[iip] = 0;
		return;
	}

	energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
	// 1) bp(i+1) == j-1
	if (tree.tree[i+1].pair == j-1){
		m1 = get_e_stP(i,j) + get_BE(i+1,j-1,ip,jp,tree);

	}

	// cases 2-5 are all need an l s.t. i<l<=ip and jp<=bp(l)<j
	for (cand_pos_t l = i+1; l<= ip ; l++){

		// Hosna: March 14th, 2007
		if (tree.tree[l].pair >= -1 && jp <= tree.tree[l].pair && tree.tree[l].pair < j){
			// Hosna, March 15, 2007
			// since not_paired_all[i,l] includes i and l themselves
			// and in BE energy calculation we are looking for the oepn region (i,l)
			// we have to look at not_paired_all[i+1,l-1]
			cand_pos_t lp = tree.tree[l].pair;
			cand_pos_t il = index[i]+l-i;
			cand_pos_t lpj = index[lp]+j-lp;
			// 2)
			// Hosna June 29, 2007
			// when we pass a stacked pair instead of an internal loop to e_int, it returns underflow,
			// so I am checking explicitely that we won't have stems instead of internal loop
			bool empty_region_il = (tree.up[(l)-1] >= l-i-1); //empty between i+1 and lp-1
			bool empty_region_lpj = (tree.up[(j)-1] >= j-lp-1); // empty between l+1 and ip-1
			bool weakly_closed_il = tree.weakly_closed(i+1,l-1); // weakly closed between i+1 and lp-1
			bool weakly_closed_lpj = tree.weakly_closed(lp+1,j-1); // weakly closed between l+1 and ip-1


			if (empty_region_il && empty_region_lpj){//&& !(ip == (i+1) && jp==(j-1)) && !(l == (i+1) && lp == (j-1))){
				energy_t tmp = get_e_intP(i,l,lp,j)+ get_BE(l,lp,ip,jp,tree);
				m2 = std::min(m2,tmp);
			}

			// 3)
			if (weakly_closed_il && weakly_closed_lpj){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				energy_t tmp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp,tree) + get_WIP(lp+1,j-1)+ ap_penalty + 2* bp_penalty;
				m3 = std::min(m3,tmp);
			}

			// 4)
			if (weakly_closed_il && empty_region_lpj){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				energy_t tmp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp,tree) + cp_penalty * (j-lp+1) + ap_penalty + 2*bp_penalty;
				m4 = std::min(m4,tmp);
			}

			// 5)
			if (empty_region_il && weakly_closed_lpj){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				energy_t tmp = ap_penalty + 2*bp_penalty + (cp_penalty * (l-i+1)) + get_BE(l,lp,ip,jp,tree) + get_WIP(lp+1,j-1);
				m5 = std::min(m5,tmp);
			}
		}
	}

	// finding the min and putting it in BE[iip]
	BE[iip] = std::min({m1,m2,m3,m4,m5});
}

energy_t pseudo_loop::get_WI(cand_pos_t i, cand_pos_t j){
	if (i>j) return 0;
	cand_pos_t ij = index[i]+j-i;
	return WI[ij];
}

energy_t pseudo_loop::get_WIP(cand_pos_t i, cand_pos_t j){
	if(i>=j) return INF;
	cand_pos_t ij = index[i]+j-i;
	return WIP[ij];
}

energy_t pseudo_loop::get_VP(cand_pos_t i, cand_pos_t j){
	if(i>=j) return INF;
	cand_pos_t ij = index[i]+j-i;
	return VP[ij];
}
energy_t pseudo_loop::get_VPL(cand_pos_t i, cand_pos_t j){
	if(i>=j) return INF;
	cand_pos_t ij = index[i]+j-i;
	return VPL[ij];
}
energy_t pseudo_loop::get_VPR(cand_pos_t i, cand_pos_t j){
	if(i>=j) return INF;
	cand_pos_t ij = index[i]+j-i;
	return VPR[ij];
}
energy_t pseudo_loop::get_WMB(cand_pos_t i, cand_pos_t j){
	if(i>=j) return INF;
	cand_pos_t ij = index[i]+j-i;
	return WMB[ij];
}

energy_t pseudo_loop::get_WMBW(cand_pos_t i, cand_pos_t j){
	if(i>=j) return INF;
	cand_pos_t ij = index[i]+j-i;
	return WMBW[ij];
}

energy_t pseudo_loop::get_WMBP(cand_pos_t i, cand_pos_t j){
	if(i>=j) return INF;
	cand_pos_t ij = index[i]+j-i;
	return WMBP[ij];
}

energy_t pseudo_loop::get_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, sparse_tree &tree){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i>= TURN && i >= 1 && i <= ip && ip < jp && jp <= j && j <=n && tree.tree[i].pair >=0 && tree.tree[j].pair >= 0 && tree.tree[ip].pair >= 0 && tree.tree[jp].pair >= 0 && tree.tree[i].pair == j && tree.tree[j].pair == i && tree.tree[ip].pair == jp && tree.tree[jp].pair == ip){
		if(i == ip && j == jp && i<j){
			return 0;
		}
		cand_pos_t iip = index[i]+ip-i;

		return BE[iip];
	}else{
		return INF;
	}
}

energy_t pseudo_loop::compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, const paramT *params){

	const pair_type ptype_closing = pair[S_[i]][S_[j]];
    return E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<paramT *>(params));
}

energy_t pseudo_loop::get_e_stP(cand_pos_t i, cand_pos_t j){
	if (i+1 == j-1){ // TODO: do I need something like that or stack is taking care of this?
		return INF;
	}
	energy_t ss = compute_int(i,j,i+1,j-1,params_);
	return lrint(e_stP_penalty * ss);
}

energy_t pseudo_loop::get_e_intP(cand_pos_t i, cand_pos_t ip, cand_pos_t jp, cand_pos_t j){
	// Hosna Feb 12th, 2007:
	// this function is only being called in branch 5 of VP
	// and branch 2 of BE
	// in both cases regions [i,ip] and [jp,j] are closed regions
	energy_t e_int = compute_int(i,j,ip,jp,params_);
	energy_t energy = lrint(e_intP_penalty * e_int);
	return energy;
}

void pseudo_loop::back_track(std::string structure, minimum_fold *f, seq_interval *cur_interval, sparse_tree &tree)
{
	this->structure = structure;
	this->f = f;

	// changing the nested if structure to switch for optimality
	switch (cur_interval->type)
	{
			case P_WMB:
			{
				cand_pos_t i = cur_interval->i;
				cand_pos_t j = cur_interval->j;
				if (i >= j)
					return;
				cand_pos_t best_l = -1, best_row = -1;
				energy_t tmp = INF, min = INF;

				// case 1
				if (tree.tree[j].pair >= 0 && j > tree.tree[j].pair && tree.tree[j].pair > i){
					energy_t acc = INF;
					cand_pos_t bp_j = tree.tree[j].pair;
					for (cand_pos_t l = bp_j +1; l < j; l++){
						// Hosna: April 24, 2007
						// correct case 2 such that a multi-pseudoknotted
						// loop would not be treated as case 2
						cand_pos_t Bp_lj = tree.Bp(l,j);

						if (Bp_lj >= 0 && Bp_lj<n){
							energy_t sum = get_BE(bp_j,j,tree.tree[Bp_lj].pair,Bp_lj,tree) + get_WMBP(i,l) + get_WI(l+1,Bp_lj-1);
							if (acc > sum){
								acc = sum;
								best_l = l;
							}
						}
					}
					tmp = PB_penalty + acc;
					if (tmp < min){
						min = tmp;
						best_row = 1;
					}

				}
				// case WMBP
				tmp = get_WMBP(i,j);
				if (tmp < min){
					min = tmp;
					best_row = 2;
				}

				switch (best_row)
				{
					case 1:
						if (best_l > -1){
							insert_node(i,best_l,P_WMBP);
							insert_node(best_l +1,tree.Bp(best_l,j)-1,P_WI);
							insert_node(tree.tree[j].pair,tree.tree[tree.Bp(best_l,j)].pair, P_BE);
						}
						break;
					case 2:
						insert_node(i,j,P_WMBP);
						break;
				}
			}
				break;

			case P_WMBW:
			{
				cand_pos_t i = cur_interval->i;
				cand_pos_t j = cur_interval->j;
				if (i >= j)
					return;
				cand_pos_t best_l = -1;
				energy_t min = INF;

				if(tree.tree[j].pair < j){
					for(cand_pos_t l = i+1; l<j; l++){
						if (tree.tree[l].pair < 0 && tree.tree[l].parent->index > -1 && tree.tree[j].parent->index > -1 && tree.tree[j].parent->index == tree.tree[l].parent->index){
							energy_t tmp = get_WMBP(i,l) + get_WI(l+1,j);
							if(tmp<min){
								min = tmp;
								best_l = l;
							}

						}
					}
				}

				if (best_l > -1){
					insert_node(i,best_l,P_WMBP);
					insert_node(best_l +1,j,P_WI);
				}
			}
				break;
			case P_WMBP:
			{
				cand_pos_t i = cur_interval->i;
				cand_pos_t j = cur_interval->j;
				if (i >= j)
					return;
				cand_pos_t best_l = -1, best_row = -1; 
				energy_t min = INF;

				// case 1
				cand_pos_t b_ij = tree.b(i,j);
				if (tree.tree[j].pair < 0){
					energy_t acc = INF;
					cand_pos_t l3 = -1;
					cand_pos_t b_ij = tree.b(i,j);
					for (cand_pos_t l = i+1; l<j ; l++)	{
						cand_pos_t bp_il = tree.bp(i,l);
						cand_pos_t Bp_lj = tree.Bp(l,j);
						if(b_ij > 0 && l < b_ij){
							if (bp_il >= 0 && l>bp_il && Bp_lj > 0 && l<Bp_lj){ // bp(i,l) < l < Bp(l,j)
		
								cand_pos_t B_lj = tree.B(l,j);
								if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l+TURN <=j){
									energy_t sum = get_BE(tree.tree[B_lj].pair,B_lj,tree.tree[Bp_lj].pair,Bp_lj,tree)+ get_WMBP(i,l-1)+ get_VP(l,j);
									if (acc > sum){
										acc = sum;
										l3 = l;
									}
								}
							}
						}
					}
					energy_t tmp = 2 *PB_penalty + acc;
					if (tmp < min){
						min = tmp;
						best_row = 1;
						best_l = l3;
					}
				}

				// case 2
				if (tree.tree[j].pair < 0){
					energy_t acc = INF;
					cand_pos_t l3 = -1;
					cand_pos_t b_ij = tree.b(i,j);
					for (cand_pos_t l = i+1; l<j ; l++)	{
						cand_pos_t bp_il = tree.bp(i,l);
						cand_pos_t Bp_lj = tree.Bp(l,j);
						if(b_ij>0 && l<b_ij){
							if (bp_il >= 0 && l>bp_il && Bp_lj > 0 && l<Bp_lj){ // bp(i,l) < l < Bp(l,j)
		
								cand_pos_t B_lj = tree.B(l,j);
								if (i <= tree.tree[l].parent->index && tree.tree[l].parent->index < j && l+TURN <=j){
									energy_t sum = get_BE(tree.tree[B_lj].pair,B_lj,tree.tree[Bp_lj].pair,Bp_lj,tree)+ get_WMBW(i,l-1)+ get_VP(l,j);
									if (acc > sum){
										acc = sum;
										l3 = l;
									}
								}
							}
						}
					}
					energy_t tmp = 2 *PB_penalty + acc;
					if (tmp < min){
						min = tmp;
						best_row = 2;
						best_l = l3;
					}
				}

				// case 3
				energy_t temp = get_VP(i,j) + PB_penalty;
				if (temp < min){
					min = temp;
					best_row = 3;
				}

				// case 4
				if(tree.tree[j].pair < 0 && tree.tree[i].pair >= 0){
					cand_pos_t l1 = -1;
					energy_t acc = INF;
					for (cand_pos_t l = i+1; l < j; l++){
						// Hosna, April 9th, 2007
						// checking the borders as they may be negative
						// Hosna: July 5th, 2007:
						cand_pos_t bp_il = tree.bp(i,l);
						// removed bp(l)<0 as VP should handle that
						if( bp_il >= 0 &&  bp_il < n && l+TURN <= j){
							// Hosna: April 19th, 2007
							// the chosen l should be less than border_b(i,j)
							energy_t BE_energy = get_BE(i,tree.tree[i].pair,bp_il,tree.tree[bp_il].pair,tree);
							energy_t WI_energy = get_WI(bp_il +1,l-1);
							energy_t VP_energy = get_VP(l,j);
							energy_t sum = BE_energy + WI_energy + VP_energy;
							if (acc > sum){
								acc = sum;
								l1 = l;
							}
						}
					}
					energy_t tmp = 2*PB_penalty + acc;
					if (tmp < min){
						min = tmp;
						best_row = 4;
						best_l = l1;
					}
				}

				switch (best_row)
				{
					case 1:
						if (best_l > -1){
							insert_node(i,best_l -1,P_WMBP);
							insert_node(best_l,j,P_VP);
							insert_node(tree.tree[tree.B(best_l,j)].pair,tree.tree[tree.Bp(best_l,j)].pair,P_BE);
						}
						break;
					case 2:
						if (best_l > -1){
							insert_node(tree.tree[tree.B(best_l,j)].pair,tree.tree[tree.Bp(best_l,j)].pair,P_BE);
							insert_node(i,best_l-1,P_WMBW);
							insert_node(best_l,j,P_VP);
						}
						break;
					case 3:
						insert_node(i,j,P_VP);
						break;
					case 4:
						if (best_l > -1){
							insert_node(i,tree.bp(i,best_l),P_BE);
							insert_node(tree.bp(i,best_l)+1,best_l-1,P_WI);
							insert_node(best_l,j,P_VP);
						}
						break;
				}

			}
				break;
			case P_VP:
			{
				cand_pos_t i = cur_interval->i;
				cand_pos_t j = cur_interval->j;
				if (i>=j){
					return;
				}
				f[i].pair = j;
				f[j].pair = i;
				this->structure[i] = '[';
				this->structure[j] = ']';
				//printf("----> original VP: adding (%d,%d) <-------\n",i,j);
				f[i].type = P_VP;
				f[j].type = P_VP;

				int min = INF, tmp = INF, best_ip = INF, best_jp = INF, best_row = -1, best_r = INF;
				cand_pos_t Bp_ij = tree.Bp(i,j);
				cand_pos_t B_ij = tree.B(i,j);
				cand_pos_t b_ij = tree.b(i,j);
				cand_pos_t bp_ij = tree.bp(i,j);
				//case 1
				// Hosna April 9th, 2007
				// need to check the borders as they may be negative
				if((tree.tree[i].parent->index) > 0 && (tree.tree[j].parent->index) < (tree.tree[i].parent->index) && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0){
					energy_t WI_ipus1_BPminus = get_WI(i+1,Bp_ij - 1) ;
					energy_t WI_Bplus_jminus = get_WI(B_ij + 1,j-1);
					tmp =   WI_ipus1_BPminus + WI_Bplus_jminus;
					if (tmp < min){
						min = tmp;
						best_row = 1;
					}
				}
				//case 2
				// Hosna April 9th, 2007
				// checking the borders as they may be negative
				if ((tree.tree[i].parent->index) < (tree.tree[j].parent->index) && (tree.tree[j].parent->index) > 0 && b_ij>= 0 && bp_ij >= 0 && Bp_ij < 0){
					energy_t WI_i_plus_b_minus = get_WI(i+1,b_ij - 1);
					energy_t WI_bp_plus_j_minus = get_WI(bp_ij + 1,j-1);
					tmp = WI_i_plus_b_minus + WI_bp_plus_j_minus;
					if (tmp < min){
						min = tmp;
						best_row = 2;
					}
				}
				//case 3
				// Hosna April 9th, 2007
				// checking the borders as they may be negative
				if((tree.tree[i].parent->index) > 0 && (tree.tree[j].parent->index) > 0 && Bp_ij >= 0 && B_ij >= 0  && b_ij >= 0 && bp_ij>= 0){
					energy_t WI_i_plus_Bp_minus = get_WI(i+1,Bp_ij - 1);
					energy_t WI_B_plus_b_minus = get_WI(B_ij + 1,b_ij - 1);
					energy_t WI_bp_plus_j_minus = get_WI(bp_ij +1,j - 1);
					tmp = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
					if (tmp < min){
						min = tmp;
						best_row = 3;
					}
				}
				//case 4
				pair_type ptype_closingip1jm1 = pair[S_[i+1]][S_[j-1]];
				if(tree.tree[i+1].pair < 0 && tree.tree[j-1].pair < 0 && ptype_closingip1jm1 > 0){
					tmp = get_e_stP(i,j)+ get_VP(i+1,j-1);
					if (tmp < min){
						min = tmp;
						best_row = 4;
					}
				}
				
				cand_pos_t min_borders = std::min((cand_pos_tu) Bp_ij, (cand_pos_tu) b_ij);
				cand_pos_t edge_i = std::min(i+MAXLOOP+1,j-TURN-1);
				min_borders = std::min({min_borders,edge_i});
				for (cand_pos_t k = i+1; k < min_borders; ++k){
					// Hosna: April 20, 2007
					// i and ip and j and jp should be in the same arc
					// it should also be the case that [i+1,ip-1] && [jp+1,j-1] are empty regions
					if (tree.tree[k].pair < -1 && (tree.up[(k)-1] >= ((k)-(i)-1))){
						// Hosna, April 9th, 2007
						// whenever we use get_borders we have to check for the correct values
						cand_pos_t max_borders = std::max(bp_ij,B_ij)+1;
						cand_pos_t edge_j = k+j-i-MAXLOOP-2;
						max_borders = std::max({max_borders,edge_j});
						for (cand_pos_t l = j-1; l > max_borders ; --l){
							pair_type ptype_closingkj = pair[S_[k]][S_[l]];
							if (tree.tree[l].pair < -1 && ptype_closingkj>0 && (tree.up[(j)-1] >= ((j)-(l)-1))){
								// Hosna: April 20, 2007
								// i and ip and j and jp should be in the same arc
								tmp = get_e_intP(i,k,l,j) + get_VP(k,l);
								if (tmp < min){
									min = tmp;
									best_row = 5;
									best_ip = k;
									best_jp = l;
								}
							}
						}
					}
				}

				cand_pos_t min_Bp_j = std::min((cand_pos_tu) tree.b(i,j), (cand_pos_tu) tree.Bp(i,j));
				cand_pos_t max_i_bp = std::max(tree.B(i,j),tree.bp(i,j));
				
				for (cand_pos_t k = i+1; k < min_Bp_j; ++k){
					tmp = get_WIP(i+1,k-1) + get_VP(k,j-1) + ap_penalty + 2* bp_penalty;
					if (tmp < min){
						min = tmp;
						best_row = 6;
						best_r = k;
					}
					
				}

				for (cand_pos_t k = max_i_bp+1; k < j; ++k){
					tmp = get_VP(i+1,k) + get_WIP(k+1,j-1) + ap_penalty + 2* bp_penalty;
					if (tmp < min){
						min = tmp;
						best_row = 7;
						best_r = k;
					}
					
				}


				for (cand_pos_t k = i+1; k < min_Bp_j; ++k){
					tmp = get_WIP(i+1,k-1) + get_VPR(k,j-1) + ap_penalty + 2* bp_penalty;
					if (tmp < min){
						min = tmp;
						best_row = 8;
						best_r = k;
					}
					
				}

				
				for (cand_pos_t k = max_i_bp+1; k < j; ++k){
					tmp = get_VPL(i+1,k) + get_WIP(k+1,j-1) + ap_penalty + 2* bp_penalty;
					if (tmp < min){
						min = tmp;
						best_row = 9;
						best_r = k;
					}
					
				}

				switch (best_row)
				{
					case 1:
						if (i+1 <= Bp_ij-1){
							insert_node(i+1,Bp_ij-1,P_WI);
						}
						if ( B_ij+1 <= j-1){
							insert_node( B_ij+1,j-1,P_WI);
						}
						break;
					case 2:
						if (i+1 <= b_ij-1){
							insert_node(i+1,b_ij-1,P_WI);
						}
						if (bp_ij+1 <= j-1){
							insert_node(bp_ij+1,j-1,P_WI);
						}
						break;
					case 3:
						if (i+1 <= Bp_ij-1){
							insert_node(i+1,Bp_ij-1,P_WI);
						}
						if (B_ij+1 <= b_ij-1){
							insert_node(B_ij+1,b_ij-1,P_WI);
						}
						if (bp_ij+1 <= j-1){
							insert_node(bp_ij+1,j-1,P_WI);
						}
						break;
					case 4:
						if (i+1 <= j-1){
							insert_node(i+1,j-1,P_VP);
						}
						break;
					case 5:
						if (best_ip <= best_jp){
							insert_node(best_ip,best_jp,P_VP);
						}
						break;
					case 6:
						if (i+1 <= best_r-1){
							insert_node(i+1,best_r-1,P_WIP);
						}
						if (best_r <= j-1){
							insert_node(best_r,j-1,P_VP);
						}
						break;
					case 7:
						if (i+1 <= best_r){
							insert_node(i+1,best_r,P_VP);
						}
						if (best_r+1 <= j-1){
							insert_node(best_r+1,j-1,P_WIP);
						}
						break;
					case 8:
						if (i+1 <= best_r-1){
							insert_node(i+1,best_r-1,P_WIP);
						}
						if (best_r <= j-1){
							insert_node(best_r,j-1,P_VPR);
						}
						break;
					case 9:
						if (i+1 <= best_r){
							insert_node(i+1,best_r,P_VPL);
						}
						if (best_r+1 <= j-1){
							insert_node(best_r+1,j-1,P_WIP);
						}
						break;
				}
			}
				break;

			case P_VPL:
			{
				cand_pos_t i = cur_interval->i;
				cand_pos_t j = cur_interval->j;
				if (i >= j) return;
				energy_t min = INF, tmp = INF;
				cand_pos_t best_k = -1;

				cand_pos_t min_Bp_j = std::min((cand_pos_tu) tree.b(i,j), (cand_pos_tu) tree.Bp(i,j));
				for(cand_pos_t k = i+1; k<min_Bp_j; ++k){
					bool can_pair = tree.up[k-1] >= (k-i);
					if(can_pair) tmp = static_cast<energy_t>((k-i)*cp_penalty) + get_VP(k,j);
					if(tmp < min){
						best_k = k;
						min = tmp;
					}
				}

				if (best_k != -1){
					insert_node(best_k,j,P_VP);
				}


			}
				break;

			case P_VPR:
			{
				cand_pos_t i = cur_interval->i;
				cand_pos_t j = cur_interval->j;
				if (i >= j) return;
				energy_t min = INF, tmp = INF;
				cand_pos_t best_k = INF, best_row = -1;

				cand_pos_t max_i_bp = std::max(tree.B(i,j),tree.bp(i,j));

				for(cand_pos_t k = max_i_bp+1; k<j; ++k){
					tmp = get_VP(i,k) + get_WIP(k+1,j);
					if(tmp < min){
						best_k = k;
						best_row = 1;
						min = tmp;
					}
				}

				for(cand_pos_t k = max_i_bp+1; k<j; ++k){
					energy_t VP_energy = get_VP(i,k);
					bool can_pair = tree.up[j-1] >= (j-k);
					if(can_pair) tmp = VP_energy + static_cast<energy_t>((j-k)*cp_penalty);
					if(tmp < min){
						best_k = k;
						best_row = 2;
						min = tmp;
					}
				}

				switch (best_row)
				{
					case 1:
						if (best_k != -1){
							insert_node(i,best_k,P_VP);
						}
						if (best_k != -1){
							insert_node(best_k+1,j,P_WIP);
						}
						break;
					case 2:
						if (best_k != -1){
							insert_node(i,best_k,P_VP);
						}
						break;
				
				}

			}
				break;
		case P_WI:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = cur_interval->j;
			if (i >= j){
				return;
			}
			energy_t min = INF, tmp = INF;
			cand_pos_t best_row = -1, best_t= -1;

				tmp = V->get_energy(i,j-1) + PPS_penalty;
				if(tmp<min){
					min = tmp;
					best_row = 1;
				} 

				tmp = get_WMB(i,j) + PSP_penalty + PPS_penalty;
				if(tmp<min){
					min = tmp;
					best_row = 2;
				}

				for (cand_pos_t t = i+1; t< j; t++){
					energy_t wi_1 = get_WI(i,t-1);
					energy_t v_energy = wi_1 + V->get_energy(t,j) + PPS_penalty;
					if(v_energy < min){
						min = v_energy;
						best_row = 3;
						best_t = t;
					}
				}

				for (cand_pos_t t = i+1; t< j; t++){
					energy_t wi_1 = get_WI(i,t-1);
					energy_t wmb_energy = wi_1 + get_WMB(t,j) + PSP_penalty + PPS_penalty;
					if(wmb_energy < min){
						min = wmb_energy;
						best_row = 4;
						best_t = t;
					}
				}
				if (tree.tree[j].pair < 0){
					tmp = get_WI(i,j-1) + PUP_penalty; 
					if(tmp<min){
						min = tmp;
						best_row = 5;
					} 
				} 
			
			switch (best_row)
			{
				case 1:
					if (i < j){
						insert_node(i,j,LOOP);
					}
					break;

				case 2:
					if (i < j){
						insert_node(i,j-1,P_WMB);
					}
					break;
				case 3:
					if (best_t != -1){
						if (i <= best_t-1){
							insert_node(i,best_t-1,P_WI);
						}
						if (best_t < j){
							insert_node(best_t,j,LOOP);
						}
					}
					break;
				case 4:
					if (best_t != -1){
						if (i <= best_t-1){
							insert_node(i,best_t-1,P_WI);
						}
						if (best_t < j){
							insert_node(best_t,j,P_WMB);
						}
					}
					break;
				case 5:
					if (i < j){
						insert_node(i,j-1,P_WI);
					}
					break;
			}
		}
			break;
		case P_BE:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = tree.tree[i].pair;
			cand_pos_t ip = cur_interval->j;
			cand_pos_t jp = tree.tree[ip].pair;
			if (i > ip || i > j || ip > jp || jp > j){
				return;
			}

			f[i].pair = j;
			f[j].pair = i;
			this->structure[i] = '(';
			this->structure[j] = ')';
			f[i].type = P_BE;
			f[j].type = P_BE;
			f[ip].pair = jp;
			f[jp].pair = ip;
			this->structure[ip] = '(';
			this->structure[jp] = ')';
			f[ip].type = P_BE;
			f[jp].type = P_BE;

			energy_t min = INF, tmp = INF;
			cand_pos_t best_row = -1, best_l = INF;
			//case 1
			if (tree.tree[i+1].pair == j-1){
				tmp = get_e_stP(i,j) + get_BE(i+1,j-1,ip,jp,tree);
				if(tmp < min){
					min = tmp;
					best_row = 1;
				}
			}
			for (cand_pos_t l = i+1; l<= ip ; l++){
				if (tree.tree[l].pair >= 0 && jp <= tree.tree[l].pair && tree.tree[l].pair < j){
				cand_pos_t lp = tree.tree[l].pair;

				bool empty_region_il = (tree.up[(l)-1] >= l-i-1); //empty between i+1 and lp-1
				bool empty_region_lpj = (tree.up[(j)-1] >= j-lp-1); // empty between l+1 and ip-1
				bool weakly_closed_il = tree.weakly_closed(i+1,l-1); // weakly closed between i+1 and lp-1
				bool weakly_closed_lpj = tree.weakly_closed(lp+1,j-1); // weakly closed between l+1 and ip-1

				if (empty_region_il && empty_region_lpj){
					tmp = get_e_intP(i,l,lp,j)+ get_BE(l,lp,ip,jp,tree);
					if (min > tmp){
						min = tmp;
						best_row = 2;
						best_l = l;
					}
				}

				// case 3
				if (weakly_closed_il && weakly_closed_lpj){
					tmp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp,tree) + get_WIP(lp+1,j-1);
					if (min > tmp){
						min = tmp;
						best_row = 3;
						best_l = l;
					}
				}

				// case 4
				if (weakly_closed_il && empty_region_lpj){
					// Hosna: July 5th, 2007
					// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
					// in a multiloop that spans a band
					tmp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp,tree) + c_penalty * (j-lp+1) + ap_penalty + 2* bp_penalty;
					if (min > tmp){
						min = tmp;
						best_row = 4;
						best_l = l;
					}
				}

				// case 5
				if (empty_region_il && weakly_closed_lpj){
					// Hosna: July 5th, 2007
					// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
					// in a multiloop that spans a band
					tmp = ap_penalty + 2* bp_penalty+ c_penalty * (l-i+1) + get_BE(l,lp,ip,jp,tree) + get_WIP(lp+1,j-1);
					if (min > tmp){
						min = tmp;
						best_row = 5;
						best_l = l;
					}
				}
				}
			}
			switch(best_row)
			{
				case 1:
					if (i+1 <= ip){
						insert_node(i+1,ip,P_BE);
					}
					break;
				case 2:
					if (best_l != INF && best_l <= ip){
						insert_node(best_l,ip,P_BE);
					}
					break;
				case 3:
					if (best_l != INF){
						if (i+1 <= best_l-1){
							insert_node(i+1,best_l-1,P_WIP);
						}
						if (best_l <= ip){
							insert_node(best_l,ip,P_BE);
						}
						if (tree.tree[best_l].pair  <= j-1){
							insert_node(tree.tree[best_l].pair +1,j-1,P_WIP);
						}
					}
					break;
				case 4:
					if (best_l != INF){
						if (i+1 <= best_l-1){
							insert_node(i+1,best_l-1,P_WIP);
						}
						if (best_l <= ip){
							insert_node(best_l,ip,P_BE);
						}
					}
					break;
				case 5:
					if (best_l != INF){
						if (best_l <= ip){
							insert_node(best_l,ip,P_BE);
						}
						if (tree.tree[best_l].pair <= j-1){
							insert_node(tree.tree[best_l].pair +1,j-1,P_WIP);
						}
					}
					break;
			}

		}
			break;
		case P_WIP:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = cur_interval->j;
			if (i == j){
				return;
			}
			energy_t min = INF, tmp = INF;
			cand_pos_t best_row = -1, best_k = INF;

			tmp = V->get_energy(i,j) + bp_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 1;
			}

			tmp = get_WMB(i,j) + PSM_penalty + bp_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 2;
			}

			for (cand_pos_t k = i+1; k < j-TURN-1; ++k){
				bool can_pair = tree.up[k-1] >= (k-i);
				energy_t wi_1 = get_WIP(i,k-1);
				tmp = wi_1 + V->get_energy(k,j);
				if (tmp < min){
					min = tmp;
					best_row = 3;
					best_k = k;
				}
			}
			for (cand_pos_t k = i+1; k < j-TURN-1; ++k){
				energy_t wi_1 = get_WIP(i,k-1);
				tmp = wi_1 + get_WMB(k,j);
				if (tmp < min){
					min = tmp;
					best_row = 4;
					best_k = k;
				}
			}
			for (cand_pos_t k = i+1; k < j-TURN-1; ++k){
				bool can_pair = tree.up[k-1] >= (k-i);
				if(can_pair) tmp = static_cast<energy_t>((k-i)*cp_penalty) + V->get_energy(k,j);
				if (tmp < min){
					min = tmp;
					best_row = 5;
					best_k = k;
				}
			}
			for (cand_pos_t k = i+1; k < j-TURN-1; ++k){
				bool can_pair = tree.up[k-1] >= (k-i);
				if(can_pair) tmp = static_cast<energy_t>((k-i)*cp_penalty) + get_WMB(k,j);
				if (tmp < min){
					min = tmp;
					best_row = 6;
					best_k = k;
				}
			}
			//case 2
			if (tree.tree[j].pair < 0){
				tmp = get_WIP(i,j-1) + cp_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 7;
				}
			}

			switch(best_row)
			{
				case 1:
					if (i < j){
						insert_node(i,j,LOOP);
					}
					break;
				case 2:
					if (i < j){
						insert_node(i,j,P_WMB);
					}
					break;
				case 3:
					if (best_k != INF){
						if (i <= best_k-1){
							insert_node(i,best_k-1,P_WIP);
						}
						if (best_k < j){
							insert_node(best_k,j,LOOP);
						}
					}
					break;
				case 4:
					if (best_k != INF){
						if (i <= best_k-1){
							insert_node(i,best_k-1,P_WIP);
						}
						if (best_k <= j){
							insert_node(best_k,j,P_WMB);
						}
					}
					break;
				case 5:
					if (best_k != INF){
						if (best_k < j){
							insert_node(best_k,j,LOOP);
						}
					}
					break;
				case 6:
					if (best_k != INF){
						if (best_k <= j){
							insert_node(best_k,j,P_WMB);
						}
					}
					break;
				case 7:
					if (i <= j-1){
						insert_node(i,j-1,P_WIP);
					}
					break;
			}
		}
			break;
	}
}

void pseudo_loop::insert_node(int i, int j, char type)
{
	seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;

}

void pseudo_loop::set_stack_interval(seq_interval *stack_interval){
	this->stack_interval = stack_interval;
}
