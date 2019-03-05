#ifndef JRIVER_MEASN_COMMUNITY_H
#define JRIVER_MEASN_COMMUNITY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "clique_network.h"

typedef struct community {
    links_t head,rear; //use rear is easy to connect
    double Pin, Pout;
    int size;
} *community_t;

community_t new_community(int);

//double between_p(community_t,community_t);
double between(community_t,community_t);

double merge(community_t,community_t,double);

double tightness(community_t);

void   free_community(community_t);

double tightness_inc(community_t,community_t,double);
int label_to_community(int*,community_t*);
community_t find_community(int,community_t*,int);
void   community_to_label(community_t*,int*,int);
int joint_size(community_t c1,community_t c2);
void joint(community_t c1,community_t c2);

#endif
