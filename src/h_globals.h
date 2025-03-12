#ifndef H_GLOBALS_H_
#define H_GLOBALS_H_

// Hosna, Nov. 1st, 2011 changed the parameter values based on HOtKNots v2
double PS_penalty = -138; //960; 		//exterior pseudoloop initiation penalty (9.6 Kcal/mol)
double PSM_penalty = 1007; //1500;		//penalty for introducing pseudoknot inside a multiloop (15 Kcal/mol)
double PSP_penalty = 1500;		//penalty for introducing pseudoknot inside a pseudoloop (15 Kcal/mol)
double PB_penalty = 246; //20;		//band penalty (0.2 Kcal/mol)
double PUP_penalty	= 6; //10;		//penalty for an un-paired base in a pseudoloop or a band (0.1 Kcal/mol)
double PPS_penalty = 96;//10;		//penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band(0.1 Kcal/mol)

double e_stP_penalty = 0.89; //0.83		// e_stP = 0.83 * e_s
double e_intP_penalty = 0.74; //0.83;		// e_intP = 0.83 * e_int

// Hosna, Nov. 2nd, 2011
// changed the multiloop penalties to be the same as the simfold's values
double a_penalty = 339; //The newest value is from DP09 //misc.multi_offset;//340;		//penalty for introducing a multiloop (3.4 Kcal/mol)
double b_penalty = 3;  //The newest value is from DP09 //misc.multi_helix_penalty; //40;			//penalty for base pair in a multiloop (0.4 Kcal/mol)
double c_penalty = 2;  //The newest value is from DP09 //misc.multi_free_base_penalty; //0;			//penalty for un-paired base in a multi-loop

double ap_penalty = 341; //340;			//penalty for introducing a multiloop that spans a band (3.4 Kcal/mol)
double bp_penalty = 56; //40;			//base pair penalty for a multiloop that spans a band (0.4 Kcal/mol)
double cp_penalty = 12; //0;			//penalty for unpaired base in a multiloop that spans a band

// Hosna November 16, 2015
double start_hybrid_penalty = 310.0;

#endif /*H_GLOBALS_H_*/
