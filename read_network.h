#ifndef LILY_READ_NETWORK_INCLUDED_H
#define LILY_READ_NETWORK_INCLUDED_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef struct link{
    int to;
    struct link* next;
} *links_t;

char * strcatstr(char *,char *);
int orig_network_size();
int orig_network_nEdge();
int orig_network_degree(int);
double orig_network_AdjMatrix(int,int);
links_t orig_network_Neighbor(int);
void orig_neighbor(int,int);

void read_pajek_unweighted_social(const char* );
void read_pajek_unweighted(const char* );
void free_orig_network();

void free_links(links_t );

#endif