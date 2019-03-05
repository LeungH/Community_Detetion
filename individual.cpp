#include "individual.h"
#include <stdlib.h>
#include <stdio.h>

extern int chidu;
extern int verbose;
extern bool obj_diff_flag;
extern double obj_diff;
extern double minx,maxx,miny,maxy;
using namespace std;

/*get f_pos_in and f_neg_out of individual*/
void eval_individual(individual_t ind) {

	community_t* comm = ind->comm;
	int n = ind->comm_n;

	int i,m;
	double pin,pout;
	int N = network_size();
	double temp = 0;
	ind->obj[0] = 0.0;
	ind->obj[1] = 0.0;

	for (i=0; i!=n; ++i) {
		pin = comm[i]->Pin;
		pout = comm[i]->Pout;
		m = comm[i]->size;
		temp += (pin/m);
		ind->obj[1] += (pout/m);
	}
    //ind->obj[0] = 2*(N-n) - temp;
    ind->obj[0] = (2*(N-n) - temp)/10;
    if (obj_diff_flag){
//        ind->obj[0]/=obj_diff;
        ind->obj[0]=(ind->obj[0]-minx)/(maxx-minx);
        ind->obj[1]=(ind->obj[1]-miny)/(maxy-miny);
        cout << "<eval" << ind->obj[0]<<","
            <<ind->obj[1] << ">" << endl;
    }
    
}

void min_max_normalization(individual_t ind,double minx,double maxx,double miny,double maxy){
    ind->obj[0]=(ind->obj[0] - minx)/(maxx-minx);
    ind->obj[1]=(ind->obj[1] - miny)/(maxy-miny);
}

/*new a individual with genotype geno
 * actually it is only used in separate
 * as the genotype here is the label
 * the genotype array will NOT be copied */
individual_t new_individual(int* gene) {
	int i;
	links_t j;
	individual_t ind = (individual_t) malloc(sizeof(struct individual));
	ind->gene = gene;
	ind->comm = (community_t *) calloc(network_size(),sizeof(community_t));
	//convert label to communities
	ind->comm_n = label_to_community(gene,ind->comm);
	ind->comm = (community_t *)realloc(ind->comm,ind->comm_n * sizeof(community_t));

	//evaluate individual
	eval_individual(ind);
	
	ind->refcount = 1;
	return ind;
}

/*calc_sectorial_index */
void calc_sectorial_index(individual_t ind,vector<double> observer_point,int objetive_num,int n,vector<vector<double>>anchor_point){
	double *f_obj = (double*)malloc(objetive_num*sizeof(double));
	double denominator = 0.0;

	for (int i = 0; i < objetive_num; i++){
		if (chidu == 0)
			f_obj[i] = ind->obj[i] - observer_point[i];
		else{
			if (anchor_point[abs(i - 1)][i] - anchor_point[i][i] == 0){
				f_obj[i] = ind->obj[i] - observer_point[i];
				//n0++;
			}
			else{
				//n1++;
				f_obj[i] = (ind->obj[i] - anchor_point[i][i]) / (anchor_point[abs(i - 1)][i] - anchor_point[i][i]);
			}
		}
		denominator += f_obj[i];
	}
	if (denominator == 0)
		ind->sectorial_angle = 0;
	else
		ind->sectorial_angle = f_obj[0] / denominator;

	free(f_obj);
	ind->sectorial_index = (int)floor((n - 1)*ind->sectorial_angle + 1 / 2);
}


/*release individual
 * first we decrease the refcount of an individual
 * if refcount==0 means no one use this individual
 * do real free*/
void free_individual(individual_t* iptr) {
	int i;
	individual_t ind = *iptr;

	if (iptr==NULL) return;

	if (ind==NULL) return;

	if (--(ind->refcount)) return;

	//real free
	free(ind->gene);
	for (i=0; i!=ind->comm_n; ++i) 
	free_community(ind->comm[i]);
	free(ind->comm);
	free(ind);
	*iptr = NULL;
}

/*set individual(d=s)
 * decrease refcount of d and increase refcount of s
 * if no one uses d do real free
 * so there is only one copy of the same individual
 * in memory, and copy an individual is fast*/
void set_individual(individual_t* d,individual_t s) {
	free_individual(d);
	(*d) = s;
	(s->refcount)++;
}



/*tchebycheff function
 * max (lambda_pos * |f_pos_in-f_pos_in_ref|,
 *      lambda_neg  * |f_neg_in-f_neg_out_ref|) */
double tchebycheff(individual_t ind, double* ref,double lambda) {

	double lp = lambda; //lambda_pos
	double ln = 1-lambda; //lambda_neg as lp+ln===1
	double fp = ind->obj[0]; //f_pos_in
	double fn = ind->obj[1]; //f_neg_out
	double rp = ref[0]; //f_pos_in_ref
	double rn = ref[1]; //f_neg_out_ref

	/*as minmium problem, f* is always less than f
	 * so we need not use abs */
	double t1 = lp * (fp-rp);
	double t2 = ln * (fn-rn);

	if (t1>t2) return t1;
	else return t2;
}
