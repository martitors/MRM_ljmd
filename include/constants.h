#ifndef CONSTANT_H
#ifdef __cpluscplus

extern "C" {
	#endif
	#define CONSTANT_H

    	/* a few physical constants */

    		#define BLEN 200
    		extern const double kboltz;  /* boltzman constant in kcal/mol/K */
    		extern const double mvsq2e; /* m*v^2 in kcal/mol */
	#ifdef __cpluscplus
	}
	#endif
#endif
