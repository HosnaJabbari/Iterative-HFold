
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "V_final.h"
#include "h_struct.h"

V_final::V_final(){
	index = new int [MAXSLEN];
    int total_length = (MAXSLEN *(MAXSLEN+1))/2;
    index[0] = 0;
    int i;
    for (int i=1; i < MAXSLEN; i++)
        index[i] = index[i-1]+MAXSLEN-i+1;
    
    type = new int[total_length];
    if (type == NULL) giveup ("Cannot allocate memory", "V_final");
    for (i = 0; i < total_length; i++) type[i] = -1;
//	printf("an object of V_final was successfully created! \n");

	
}

V_final::~V_final(){
	delete [] index;
	delete [] type;	
}

void V_final::setloops(s_energy_matrix *v, VM_final *vm){
	this->v = v;
	this->vm = vm;

//	printf("V_final loops were successfully set! \n");
}

int V_final::get_energy(int i, int j){
	// Hosna: June 28th, 2007
	if (i >= j || (fres[i].pair > -1 && fres[i].pair != j) || (fres[j].pair > -1 && fres[j].pair != i)){
		return INF;
	}	
	
	int v_energy = v->get_energy(i,j);
	/*if (i>=13 && j<=39 && v_energy<INF){
	printf("V_final: v_energy(%d,%d) = %d and type = %c\n", i,j,v_energy, get_type(i,j));
	}*/
	int vm_energy = vm->get_energy(i,j);
//	printf("V_final: vm_energy(%d,%d) = %d \n", i,j,vm_energy);
	int ij = index[i]+j-i;
	if (v_energy < vm_energy){
		type[ij] = 0;
	}else{
		type[ij] = 1;
	}
	return MIN(v_energy,vm_energy);	
}

char V_final::get_type(int i, int j){
	int ij = index[i]+j-i;
	if (type[ij] == 0) // comes from v
	{
		return v->get_type(i,j);
	}
	return MULTI;
}

