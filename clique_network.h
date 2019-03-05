#ifndef LILY_CLIQUE_NETWORK_INCLUDED_H
#define LILY_CLIQUE_NETWORK_INCLUDE_H

#include <string.h>
#include "maxcliques.h"
#include "read_network.h"


//the function for finding k-clique network
int construct_clique_network();
void output_clique_network(char *name);
int num_of_clique_node(int i);
int * node_of_clique(int i);

//the function for community detection
void parameter_of_clique_network();
int  network_size();
int is_network_ok();
double weight(int,int);
double link_sum(int);
links_t neighbor(int);
void  free_network();
double degree(int);
int is_network_ok();
void free_net_clique();

#endif