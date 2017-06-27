#ifndef V_FINAL_H_
#define V_FINAL_H_

#include "h_struct.h"
#include "s_energy_matrix.h"
#include "VM_final.h"

class VM_final;
class V_final{
	public: 
	// constructor
	V_final();
	~V_final();
	void setloops(s_energy_matrix *v, VM_final *vm);
	int get_energy(int i, int j);
	
	char get_type (int i, int j);
    // return the type at v_final(i,j)
    
    void set_features(str_features *f){fres = f;}
		
	protected:
	s_energy_matrix *v;
	VM_final *vm;
	int *index;
	int *type;
	str_features *fres;
	
	
};
#endif /*V_FINAL_H_*/
