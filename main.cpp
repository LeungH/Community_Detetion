#include <time.h>
#include <string.h>
#include <stdio.h>
#include "individual.h"
#include "read_network.h"
#include "maxcliques.h"
#include "clique_network.h"
#include "evaluation_index.h"

#include "caea.h"

void get_parameters(int argc,char* argv[]);


extern int generation;
extern int popsize;
extern int verbose;
extern char input_path[1000];
extern char output_path[1000];
extern char label_real_or_generated[50];


#define MAX_LENGTH_FOR_STR 1000


int main(int argc,char* argv[]){

	int i, j;
	double t;

	// get parameters from terminal
	get_parameters(argc,argv);  //set all parameters for EA containing generation,popsize...

	char input_network[MAX_LENGTH_FOR_STR],input_community[MAX_LENGTH_FOR_STR],output_community[MAX_LENGTH_FOR_STR];
	char output_clique[MAX_LENGTH_FOR_STR],output_index[MAX_LENGTH_FOR_STR];
	char output_MOD[MAX_LENGTH_FOR_STR],output_T[MAX_LENGTH_FOR_STR];

	// The ground truth of generated network is known, but the ground truth of real-world network is unknown
	// So there is some differences between these two networks
	// if label_real_or_generated is 'generated', indicate the input network is synthetic networks
	if(label_real_or_generated[0] == 'g'){
		sprintf(input_network,"%s.network",input_path);
		sprintf(input_community,"%s.comm",input_path);
		sprintf(output_community,"%s",output_path);
		sprintf(output_clique,"%s",output_community);
		sprintf(output_MOD,"%sModularity",output_path);
		sprintf(output_index,"%sgNMI",output_path);
		sprintf(output_T,"%sesclipe_t.dat",output_clique);
	}
	// if label_real_or_generated is 'real', indicate the input network is real-world networks
	else if(label_real_or_generated[0] == 'r'){
		//sprintf(input_community,"%s.comm",input_path);
		//sprintf(output_index,"%sgNMI",output_path);


		sprintf(input_network,"%s.network",input_path);
		sprintf(output_community,"%s",output_path);
		sprintf(output_clique,"%s",output_community);
		sprintf(output_MOD,"%sModularity",output_path);
		sprintf(output_T,"%sesclipe_t.dat",output_clique);
	}
	else
		printf("please check the parameters\n");


	//test
	printf("input_network: %s\n",input_network);
	printf("input_community: %s\n",input_community);
	printf("output_community: %s\n",output_community);
	printf("output_clique: %s\n",output_clique);
	printf("output_index: %s\n",output_index);
	printf("output_T: %s\n",output_T);
	
	// //--------------write community--------------------
	// label2comm(input_community,input_path);
	// return 0;
	
	FILE *fp;
	time_t l_time,m_time,n_time,clique_time;
//    double l_time,m_time,n_time,clique_time;
	fp = fopen(output_T,"w");
	l_time = time(NULL);

	//read data from network file
	if(label_real_or_generated[0] == 'r'){
		read_pajek_unweighted_social(input_network);
		printf("read data from social network unweighted\n");
	}

	else if(label_real_or_generated[0] == 'g') {
		read_pajek_unweighted(input_network);
		printf("read data from LFM unweighted\n");
	}
	else
		printf("please check parameters\n");

	//finding all complete subgraph, and attributing node to the biggest clique which belong to
	printf("find maximum cliques for every node\n");
	maxclique(output_clique);

	//convert node network to clique network and output clique network
	printf("------------constructing clique network----------------");
	construct_clique_network();
	output_clique_network(output_clique);

	//get parameters of maximal-clique graph
	printf("----------calculating parameters of clique network--------------\n");
	parameter_of_clique_network();

	// struct population *pop=NULL;
	m_time=time(NULL);
	clique_time = m_time-l_time;
	fprintf(fp,"%lf\n",(double)clique_time);
	//printf ("%lf\t%lf\t%lf\n",m_time , l_time ,clique_time);
	

	//community detection by using CAEA
	CAEA caea;
	for(j=0;j<10;j++){
		m_time=time(NULL);

		if (verbose>=1) printf("initalizing population\n");

		// pop = init_population();

		individual_t *pop = caea.init_population();
		caea.evolution(pop);

//        for (i=1; i<=generation; ++i) {
//            if (verbose>=1) printf("---------------" "%3d/%-3d""---------------\n",i,generation);
//            // evolve_population(pop);
//        }

		//recoding the esclipe time
		n_time = time(NULL);
		t=n_time - m_time;
		fprintf(fp,"%lf\n",(double)t);

		caea.dump_population(output_community,pop,j);

		//if(u==0 && on==0 && om==0)  uoModularity( pop,j,output_MOD);   // for social network which without known ground truth
		uoModularity( pop,j,output_MOD);
		if(label_real_or_generated[0] == 'g') gNMI(pop,input_community,j,output_index);


		caea.free_population(pop);
	}

	free_network();
	free_cliques();
	free_net_clique();
	free_orig_network();

	fclose(fp);

	return 0;
}
