#include "VM_final.h"
#include "externs.h"
#include "h_externs.h"

VM_final::VM_final(int *seq, int len)
{
	length = len;
	sequence = seq;
	this->v = NULL;
	this->wmb = NULL;

	index = new int[length];    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    int total_length = (length *(length+1))/2;
    index[0] = 0;
    int i;
    for (i=1; i < length; i++)
        index[i] = index[i-1]+length-i+1;

    WM = new int [total_length];
    if (WM == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (i=0; i < total_length; i++) WM[i] = INF;

    VM = new int [total_length];
    if (VM == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (i=0; i < total_length; i++) VM[i] = INF;

//    printf("an object of VM_final was successfully created! \n");
}

VM_final::~VM_final()
{
	delete [] index;
    delete [] WM;
    delete [] VM;
}

void VM_final::compute_energy(int i, int j, str_features *fres){
	// Hosna June 26, 2007
	// I have to figure out how to calculate the energy here

	// here comes the copied part from simfold with all dangling energies
	int min = INF, tmp, k;
    int iplus1k;
    int kplus1jminus1;
    int iplus2k;
    int kplus1jminus2;
    for (k = i+2; k <= j-3; k++)
    {
        iplus1k = index[i+1] + k -i-1;
        kplus1jminus1 = index[k+1] + j-1 -k-1;
        iplus2k = index[i+2] + k -i-2;
        kplus1jminus2 = index[k+1] + j-2 -k-1;


        tmp = WM[iplus1k] + WM[kplus1jminus1];
        if (tmp < min)
            min = tmp;

        if (fres[i+1].pair <= -1)
        {
            tmp = WM[iplus2k] + WM[kplus1jminus1] +
                dangle_top [sequence [i]]
                [sequence [j]]
                [sequence [i+1]] +
                misc.multi_free_base_penalty;
            if (tmp < min)
                min = tmp;
        }
        if (fres[j-1].pair <= -1)
        {
            tmp = WM[iplus1k] + WM[kplus1jminus2] +
                dangle_bot [sequence[i]]
                [sequence[j]]
                [sequence[j-1]] +
                misc.multi_free_base_penalty;
            if (tmp < min)
                min = tmp;
        }
        if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
        {
            tmp = WM[iplus2k] + WM[kplus1jminus2] +
                dangle_top [sequence [i]]
                [sequence [j]]
                [sequence [i+1]] +
                dangle_bot [sequence[i]]
                [sequence[j]]
                [sequence[j-1]] +
                2 * misc.multi_free_base_penalty;
            if (tmp < min)
            min = tmp;
        }
    }

    min += misc.multi_helix_penalty + misc.multi_offset +
           AU_penalty (sequence[i], sequence[j]);

	int wmb_energy = this->wmb->get_energy(i,j) + a_penalty + PSM_penalty;
	int ij = index[i]+j-i;
	VM[ij] = MIN(min,wmb_energy);
	//printf("VM[%d,%d] = %d \n",i,j,VM[ij]);

}

int VM_final::get_energy(int i, int j){
	int ij = index[i]+j-i;
	if (i >= j){
		return INF;
	}
	if (wmb->is_weakly_closed(i,j) == 1 || wmb->is_empty_region(i,j) == 1){
		return VM[ij];
	}
	return INF;
}

/*
int VM_final::get_energy_pk_only(int i, int j, str_features *fres){
	int ij = index[i]+j-i;
	if (i >= j || !(fres[i].pair== j && fres[j].pair==i)){
		return INF;
	}
	if (wmb->is_weakly_closed(i,j) == 1 || wmb->is_empty_region(i,j) == 1){
		return VM[ij];
	}
	return INF;
}
*/
 
/**
 *  PRE: simfold's WM matrix has been filled for i and j
 *  and now we need to fill in the WM matrix that hfold needs
 *
 */
void VM_final::WM_compute_energy(int i, int j){

	int s_wm = s_vm->get_energy_WM(i,j);

	// Hosna: July 5th, 2007
	// add a b_penalty to this case to match the previous cases
	int wmb_energy = wmb->get_energy(i,j)+PSM_penalty+b_penalty;
	int min = MIN(s_wm, wmb_energy);
	int ij = index[i]+j-i;
	this->WM[ij] = min;
//	printf("hfold's WM min = %d \n",min);
}



int VM_final::get_energy_WM(int i, int j){
	if (i >= j || wmb->is_weakly_closed(i,j) != 1 ){
		return INF;
	}
	int ij = index[i]+j-i;
//	printf("hfold's WM(%d,%d) = %d \n", i,j,WM[ij]);
	return this->WM[ij];

}
