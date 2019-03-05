#include "community.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

extern double alpha;

/*make a community with a vertice i
 * community_t use a linklist to store its vertices
 * */
community_t new_community(int i) {

	community_t c = (community_t) malloc(sizeof(struct community));
	/*insert node to community*/
	links_t t = (links_t) malloc(sizeof(struct link));
	t->to = i;
	t->next = NULL;
	c->head = c->rear = t;

	/*init pin,pout,nin,nout as only 1 vertices in community,there is no edges (both + and -)in the community, 
	* and the edges intercommunity is all the edge of the vertice */
	c->Pin = 0.0;
	c->Pout = link_sum(i);
	c->size = 1;
	return c;
}

/*calculate the sum of similarity between
 * two communities c1 and c2
 * between_p returns the sum of all the
 * positive similarities between c1,c2
 * the value is useful in merge two communities*/

double between(community_t c1,community_t c2) {
	double b=0;
	int i,j;
	links_t p1,p2;
	for (p1=c1->head; p1!=NULL; p1=p1->next){
		for (p2=c2->head; p2!=NULL; p2=p2->next) {
			i = p1->to; j = p2->to;
			b += weight(i,j);
		}
	}
	return b;
}

/*merge two communities c1,c2
 * the sum similarity between c1,c2
 * is often calculated before, so
 * we use it as bp(+) bn(-) directly
 * */
double merge(community_t c1,community_t c2,double bs) {
	/*conect linklist*/
	c1->rear->next = c2->head;
	c1->rear = c2->rear;

	/*update pin,pout,nin,nout*/
	c1->Pin += c2->Pin + 2*bs;
	c1->Pout += c2->Pout - 2*bs;

	c1->size += c2->size;

	/*make c2 a empty community*/
	c2->head = c2->rear = NULL;
	free(c2);
	return 0.0;
}

/*calculate tightness of community*/
double tightness(community_t c) {
	return c->Pin / pow(c->Pin + c->Pout,alpha);
}

/*release a community*/
void free_community(community_t c) {
	links_t t;
	while (c->head) {
	t = c->head;
	c->head = t->next;
	free(t);
	}
	c->rear = NULL;
	free(c);
}

/*calculate Tightness(c1 | c2) - Tightess(c1)
 *it is used before merge
 *as bp,bn is calculated before
 *judging merge or not,
 * and merge operation is *fast.
 */
double tightness_inc(community_t c1,community_t c2, double bs) {

	double t1 = tightness(c1);

	/*pin,pout,nin,nout if merged*/
	double Sin = c1->Pin + c2->Pin + 2 * bs;
	double Sout = c1->Pout + c2->Pout - 2*bs;

	return Sin/pow(Sin+Sout,alpha) - t1;
}

int community_size(community_t c) {
	return c->size;
}

/*number of joint vertices in two communities
 * c1 and c2. It is used in overlapping
 * community detection*/
int joint_size(community_t c1,community_t c2) {
	char* label = (char *) calloc(network_size(),sizeof(char));
	links_t p;
	int n=0;
	//mark label if a vertice in c1
	for (p=c1->head; p!=NULL; p=p->next)
	label[p->to] = 1;

	//check label of vertices in c2
	//if labeled, it is a joint vertices
	for (p=c2->head; p!=NULL; p=p->next) {
	if (label[p->to]) n++;
	}
	free(label);
	return n;
}

void joint(community_t c1,community_t c2) {
	char* label = (char *) calloc(network_size(),sizeof(char));
	double bs;
	int i;
	community_t t;
	links_t p;

	//mark label
	for (p=c1->head; p!=NULL; p=p->next)
	label[p->to] = 1;

	while (c2->head) {
	p = c2->head;
	c2->head = p->next;
	i = p->to;
	if (label[i]==0) { //vertice i not in c1

		/*add vertice i to community c1.
		 *using new_community here makes
		 a vertice into a single node community
		 so we can use between and merge
		 to merge 2 communities with
		 pin,pout,nin,nout updating*/

		t = new_community(i);
		bs = between(c1,t);
		merge(c1,t,bs);
	}
	free(p);
	/*after walk through c2, c2 will be empty*/
	}
	free(c2);
	free(label);
}

/*convert a label array p into communities c
 *and return the number of communitis */
int label_to_community(int* p,community_t* c) {

	int i,label,cn = 0;//number of communities
	int n = network_size();
	double bs;
	community_t comm;
	/*int* l records the communities 
	 * the vertice with specific label used to belong*/
	int* l = (int *) calloc(n,sizeof(int));
	memset(l,0,n*sizeof(int));

	for (i=0; i!=n; ++i) {
		label = p[i];

		c[cn] = new_community(i);

		if (l[label]==0) l[label] = ++cn;
		else {
			comm = c[l[label]-1];
			bs = between(c[cn],comm);
			merge(comm,c[cn],bs);
		}
	}

	free(l);
	return cn;
}

/*find a community contains the vertice node
 * in cn communities c.
 * for overlapping communitis, the return community
 * is the first community contains the vertice*/
community_t find_community(int node,community_t* c,int cn) {
	int i;
	links_t p;
	for (i=0; i!=cn; ++i) {
		for (p=c[i]->head; p!=NULL; p=p->next)
			if (p->to == node) return c[i];
	}
	return NULL;
}

/*convert community structure into labels
 *it is only used in separated community detection*/
void community_to_label(community_t* c,int* l,int n) {
	links_t p;
	int i;
	for (i=0; i!=n; ++i) {
		for (p=c[i]->head; p!=NULL; p=p->next)
			l[p->to] = i;
	}
}
