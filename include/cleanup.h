#ifndef CLEANUP_H
#ifdef __cplusplus
extern "C" {	
	#endif

	#define CLEANUP_H
		
		#include <stdlib.h>
		#include <stdio.h>
		#include "types.h"
		void cleanup(mdsys_t sys, FILE *erg, FILE *traj);
	
	#ifdef  __cplusplus
	}
	#endif
#endif
