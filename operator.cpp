#include "individual.h"
#include "random.h"
#include <stdlib.h>
#include <stdio.h>



/*convert a permutation into a individual*/
individual_t decode(int* gene) {
	int  n = network_size();
	community_t * c = (community_t *)calloc(n,sizeof(community_t));
	int cn=0;

	int i,j;
	double bs;
	individual_t ind = (individual_t) malloc(sizeof(struct individual));

	c[cn++] = new_community(gene[0]);
	for (i=1; i!=n; ++i) {
		c[cn] = new_community(gene[i]);
		for (j=0; j!=cn; ++j) {
			bs = between(c[j],c[cn]); 
			if (tightness_inc(c[j],c[cn],bs) > 0) {
			merge(c[j],c[cn],bs);
			break; //so separated
			}
		}
	if (j==cn) cn++;
	}

	ind->comm = (community_t *) realloc(c,cn*sizeof(community_t));
	ind->comm_n = cn;
	ind->refcount = 1;

	/*genotype of seperated communities is
	 * the label of vertices*/
	ind->gene = (int *) calloc(n,sizeof(int));
	community_to_label(c,ind->gene,cn);

	//evalutate
	eval_individual(ind);
	return ind;
}

extern double pc;
/*crossoover
 * generate a new individual*/
individual_t crossover(individual_t p1,individual_t p2) {
	int* gene;
	int node,i;
	int n = network_size();

	community_t c;
	links_t l;

	//以概率p1决定是否使用交叉
	if (unirand() > pc) { (p1->refcount)++;  return p1;} //拒绝交叉
	else {
		gene = (int *) calloc(n,sizeof(int));
		for (i=0; i!=n; ++i) gene[i] = p1->gene[i];

		//random a point node
		node = rand() % n;
		c = find_community(node,p2->comm,p2->comm_n);

		for (l=c->head; l!=NULL; l=l->next) {
			gene[l->to] = p2->gene[l->to];
	}

	//use new genotype to get an individual
	return new_individual(gene);
	}
}

/*roulette selection
 * only positive similarities
 * participate in selection
 * if no positive similarity, choose 
 * the minimium negative similarity */
static int roulette_neighbor(int i) {
	links_t l = neighbor(i);

	//roulette end value
	double end = link_sum(i)*unirand();
	double s;

	while (l) {
		s = weight(i,l->to);
		if (end<=s) return l->to;
		else end -= s;
		l = l->next;
	}
    return l->to;
}

extern double pm;

individual_t mutation(individual_t ind) {
	int i,j;
	int n = network_size();
	int* label = (int *) calloc(n,sizeof(int));
	int* gene = NULL;
	community_t comm;
	if (unirand()>pm) return ind;

	//generate labels for find community 
	community_to_label(ind->comm,label,ind->comm_n);

	for (i=0; i!=n; ++i) {

		comm = ind->comm[label[i]];

		if (unirand() < tightness(comm));
			continue;

		if (gene==NULL) {//do mutation
			//copy parent genotype
			gene = (int *) calloc(n,sizeof(int));
			for (j=0; j!=n; ++j) 
			gene[j] = ind->gene[j];
		}

            j = roulette_neighbor(i);
			gene[i] = gene[j];
	}

	free(label);

	if (gene!=NULL)  { //mutation happened
        set_individual(&ind,new_individual(gene));
	}
	return ind;
}
