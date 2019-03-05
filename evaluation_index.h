#ifndef LILY_EVALUATION_INDEX_H
#define LILY_EVALUATION_INDEX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "caea.h"
#include "read_network.h"

void uModularity(individual_t *, int,char * );
void wModularity(individual_t *, int,char * );
void uoModularity(individual_t *pop, int seq,char *mod_name);
void gNMI(individual_t *,char *,int,char *);
void label2comm(char *name,char *input_path);

#endif
