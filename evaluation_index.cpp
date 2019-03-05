#include "evaluation_index.h"
#include "caea.h"


#define NUM_OF_COMMUNITY 1000  //the maximal number of community
#define MAX_NUM_COMMUNITY_REAL 1000 //the maximal number of nodes in every communtiy
#define MAX_NUM_COMMUNITY_EVAL 1000
//#define RAND_MAX 500000
int * mem_num;
extern int popsize;
CAEA temp;

struct partition{
	int size;
	int * comm;
};

struct partition *realcomm;
struct partition *evalcomm;
int lengthofrealcomm; //the number of "real" communities


static void calc_mem_num(int n,int * comm){
	int i;
	for(i=0;i<n;i++)
		mem_num[comm[i]]++;
}

int they_are_mate(int q,int r,int sizeofnode,int *comm){
	int i,flag1=0,flag2=0,flag=0;

	for(i=0;i<sizeofnode;i++){
		if(q == comm[i]) flag1=1;
		if(r == comm[i]) flag2=1;
	}
	if(flag1==1 && flag2==1) flag=1;
	return flag;
}


// the evaluation index Modularity, which is used for unweighted and undirected network
// pop: the results of community detection;
// seq: to get a stable result, all networks will be run multiple times, this parameter means the 'seq'th time
// mod_name: used to store calculted value of Modularity
// evaluation for overlapping structures  EQ=1/2m(sum(sum(1/OvOw[Avw-kvkw/2m])))
void uoModularity(individual_t *pop, int seq,char *mod_name){
	int i,j,k,m,n;
	int maxlengthofcomm = 0, numofnode = 0;
	int *node_overlapping_membership;
	double *Mod;
	FILE *fp;
	char *output_for_mod = (char *) calloc(strlen(mod_name)+10,sizeof(char));
	int N = orig_network_size();
	double A;
	sprintf(output_for_mod,"%s_%d.dat",mod_name,seq);
	fp = fopen(output_for_mod,"w");

	//printf("#######################%d######\n",popsize);
	//system("pause");

	Mod = (double *)malloc(popsize * sizeof(double));
	if(Mod == NULL) printf("Mod_for_unweighted_network - alloc memory failed");
	memset(Mod,0,popsize * sizeof(double));

	node_overlapping_membership = (int *)malloc(N * sizeof(double));
	if(node_overlapping_membership  == NULL) printf("Mod_for_unweighted_network - alloc memory failed");
	memset(node_overlapping_membership ,0,N * sizeof(int));

	for(i=0;i<popsize;i++){
		Mod[i] = 0;
		if(maxlengthofcomm<pop[i]->comm_n)
			maxlengthofcomm = pop[i]->comm_n;
	}

	evalcomm = (struct partition *) calloc(maxlengthofcomm,sizeof(struct partition));
	for(i=0;i<maxlengthofcomm;i++){
		evalcomm[i].comm = (int *)calloc(MAX_NUM_COMMUNITY_EVAL,sizeof(int));
		if(evalcomm[i].comm == NULL) printf("evalcomm - alloc memory failed");
		memset(evalcomm[i].comm,0,MAX_NUM_COMMUNITY_EVAL * sizeof(int));
		evalcomm[i].size = 0;
	}

	for(i=0;i<popsize;i++){
		printf("------------calculating Modularity %d / %d--------------\n",i,popsize);
		for(j=0;j<N;j++) node_overlapping_membership[j] = 0;

		// convert clique network to original network
		for(j=0;j<pop[i]->comm_n;j++){
			numofnode = temp.convert_clique_to_node(pop[i]->comm[j],evalcomm[j].comm);
			evalcomm[j].size = numofnode;
			for(k=0;k<evalcomm[j].size;k++){
			//	printf("%d  ",evalcomm[j].comm[k]);
			}
			//printf("\n");
		}
		// calculate Ov and Ow
		for(k=0;k<pop[i]->comm_n;k++){
			for(j=0;j<evalcomm[k].size;j++){
				node_overlapping_membership[evalcomm[k].comm[j]]++;
			}
		}
		//for(k=0;k<N;k++) printf("%d\t%d\n",k,node_overlapping_membership[k]);
		//printf("\n");

	//	for(k=0;k<N;k++) printf("%d\t%d\n",k,network_degree(k));
	//	printf("\n");
		// calculate Avw
		A = 0;
		for(k=0;k<pop[i]->comm_n;k++){
			for(m=0;m<evalcomm[k].size;m++){
			   // printf("--------%d--------\n",evalcomm[k].comm[m]);
				for(n=0;n<evalcomm[k].size;n++){
		//		    printf("Adj between %d and %d is : %lf\n", evalcomm[k].comm[m],evalcomm[k].comm[n],
		  //              network_AdjMatrix(evalcomm[k].comm[m],evalcomm[k].comm[n]));
					A = 0;
					if(m == n) continue;
					if(orig_network_AdjMatrix(evalcomm[k].comm[m],evalcomm[k].comm[n]) > 0) {
						double temp1 = orig_network_degree(evalcomm[k].comm[m]) * orig_network_degree(evalcomm[k].comm[n]);
						double temp2 = node_overlapping_membership[evalcomm[k].comm[m]] * node_overlapping_membership[evalcomm[k].comm[n]];
						A = (orig_network_AdjMatrix(evalcomm[k].comm[m],evalcomm[k].comm[n]) - temp1/(2*orig_network_nEdge()))/temp2;
					//	printf("%lf  ",A);
					}
					Mod[i] += A;
					//printf("Modularity %d is : %lf\n",i,Mod[i]);
				}
			}
		}
		fprintf(fp,"%lf\n",Mod[i]/(2 * orig_network_nEdge()));
	}

	if(Mod != NULL){
		free(Mod); Mod = NULL;
	}

	if(node_overlapping_membership != NULL){
		free(node_overlapping_membership); node_overlapping_membership = NULL;
	}

	for(i=0;i<maxlengthofcomm;i++){
		if(evalcomm[i].comm != NULL) {
			free(evalcomm[i].comm);
			evalcomm[i].comm = NULL;
		}
	}
	free(evalcomm);
	fclose(fp);
	return;
}


double HX_f(int k){

	int N = orig_network_size();
	double int_to_double = (double) k;
	double p = int_to_double/N;
	if(p == 0 || p == 1) return 0;
	else return (-p * log(p) - (1-p) *log(1-p));
}


int intersection_size(struct partition real, struct partition eval){
	int i,j,temp;
	links_t a;
	int num = 0;
	int m = real.size;
	int n = eval.size;
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			//if(real.comm[i]<eval.comm[j]) break;
			if(real.comm[i] == eval.comm[j]) {num++;break;}
		}
	}

	return num;
}

//the index for overlapping community detection
void gNMI(individual_t *pop,char *name,int seq,char *gNMI_name){ //name:the name for community; seq:times of repetition for EA
    int i,j,c1,c2,max_lengthofevalcomm = 0;
	int *set_node,numofnode = 0;
	double P00 = 0, P10 = 0, P01 = 0, P11 = 0;
	double H00 = 0 ,H10 = 0, H01 = 0 ,H11 = 0;
	double HX1Y2 = 0, HY2X1 = 0;
	double HXY = 0, HYX = 0;
	double *HX, *HY, *minHXY, *minHYX, *gNMI,ave_gNMI = 0;
	int *nbestc1, *nbestc2;
	int num = 0,num1 = 0,num2 = 0;
	int N = orig_network_size();
	int k=0; //k:the number of communtity of original partition C1;
	FILE *f,*fp;
	char *output_for_gnmi = (char *) calloc(strlen(gNMI_name)+10,sizeof(char));
	sprintf(output_for_gnmi,"%s_%d.dat",gNMI_name,seq);
	fp = fopen(output_for_gnmi,"w");

	//--------------read community--------------------
	if((f= fopen(name,"r"))==NULL){
			printf("open real community failed\n");
			return;
	}
	realcomm = (struct partition *) calloc(NUM_OF_COMMUNITY,sizeof(struct partition));
	for(i=0;i<NUM_OF_COMMUNITY;i++){
		realcomm[i].comm = (int *)calloc(MAX_NUM_COMMUNITY_REAL,sizeof(int));
		if(realcomm[i].comm == NULL) printf("realcomm - alloc memory failed");
		memset(realcomm[i].comm,0,MAX_NUM_COMMUNITY_REAL * sizeof(int));
		realcomm[i].size = 0;
	}

	lengthofrealcomm = 0;
	while (fscanf(f,"%5d%5d",&i,&j)==2){
		int num = realcomm[j-1].size;
		if(j>lengthofrealcomm) lengthofrealcomm=j;
		realcomm[j-1].comm[num] = i -1;
		realcomm[j-1].size++;
	}
	fclose(f);
	//------------------------------------------
	// calculate the maximum of community in all population
	for(i=0;i<popsize;i++){
		if(max_lengthofevalcomm<pop[i]->comm_n)
			max_lengthofevalcomm = pop[i]->comm_n;
	}

	//initial our evaluations and store in evalcomm
	evalcomm = (struct partition *) calloc(max_lengthofevalcomm,sizeof(struct partition));
	for(i=0;i<max_lengthofevalcomm;i++){
		evalcomm[i].comm = (int *)calloc(MAX_NUM_COMMUNITY_EVAL,sizeof(int));
		if(evalcomm[i].comm == NULL) printf("evalcomm - alloc memory failed");
		memset(evalcomm[i].comm,0,MAX_NUM_COMMUNITY_EVAL * sizeof(int));
		evalcomm[i].size = 0;
	}
	//initial the parameters for calculating NMI
	HX = (double *) calloc(max_lengthofevalcomm,sizeof(double));
	if(HX == NULL) printf("HX - alloc memory failed");
	memset(HX,0,max_lengthofevalcomm * sizeof(double));

	HY = (double *) calloc(lengthofrealcomm,sizeof(double));
	if(HY == NULL) printf("HY - alloc memory failed");
	memset(HY,0,lengthofrealcomm * sizeof(double));

	nbestc1 = (int *) calloc(max_lengthofevalcomm,sizeof(int));
	if(nbestc1 == NULL) printf("nbestc1 - alloc memory failed");
	memset(nbestc1,0,max_lengthofevalcomm * sizeof(int));

	nbestc2 = (int *) calloc(lengthofrealcomm,sizeof(double));
	if(nbestc2 == NULL) printf("nbestc2 - alloc memory failed");
	memset(nbestc2,0,lengthofrealcomm * sizeof(int));

	minHXY = (double *)calloc(max_lengthofevalcomm,sizeof(double));
	if(minHXY == NULL) printf("minHXY - alloc memory failed");
	memset(minHXY,0,max_lengthofevalcomm * sizeof(double));

	minHYX = (double *)calloc(lengthofrealcomm,sizeof(double));
	if(minHYX == NULL) printf("minHYX - alloc memory failed");
	memset(minHYX,0,lengthofrealcomm * sizeof(double));

	gNMI = (double *)malloc(popsize * sizeof(double));
	if(gNMI == NULL) printf("gNMI - alloc memory failed");
	memset(gNMI,0,popsize * sizeof(double));

	//calculate HY
	//printf("\ndisplay HY\n");
	for(j=0;j<lengthofrealcomm;j++){
		HY[j] = HX_f(realcomm[j].size);
	//	printf("%lf  ",HY[j]);
	}

	for(j=0;j<max_lengthofevalcomm;j++)  minHXY[j] = RAND_MAX;
	for(j=0;j<lengthofrealcomm;j++)  minHYX[j] = RAND_MAX;

	//calculate HX and HY

	for(i=0;i<popsize;i++){
		//double NXY = 0;
		printf("------------calculating gNMI %d / %d--------------\n",i,popsize);
		//initilization
		if(i>1){
			for(c1=0;c1<pop[i-1]->comm_n;c1++){
				for(j=0;j<evalcomm[c1].size;j++) {
					evalcomm[c1].comm[j] = 0;
				}
				evalcomm[c1].size = 0;
			}
		}
		if(i>1){
			for(c1=0;c1<pop[i-1]->comm_n;c1++){
				HX[c1] = 0;
				minHXY[c1] = RAND_MAX;
				nbestc1[c1] = 0;
			}
			for(c2=0;c2<lengthofrealcomm;c2++){
				minHYX[c2] = RAND_MAX;
				nbestc2[c2] = 0;
			}
		}
		//convert clique network to node network
        
		for(c1=0;c1<pop[i]->comm_n;c1++){
            numofnode = temp.convert_clique_to_node(pop[i]->comm[c1],evalcomm[c1].comm);
			evalcomm[c1].size = numofnode;
		}
		//calculate the HX
		//printf("display HX\n");
		for(c1=0;c1<pop[i]->comm_n;c1++){
			HX[c1] = HX_f(evalcomm[c1].size);
		//	printf("%lf  ",HX[c1]);
		}

		//calculate H(Xk,Yl),HXkY,HYkX
		for(c1=0;c1<pop[i]->comm_n;c1++){
			//numofnode = convert_clique_to_node(pop[i].ind->comm[c2],set_node);
			//num1 = numofnode;
			num1 = evalcomm[c1].size;
			for(c2=0;c2<lengthofrealcomm;c2++){
				P00 = 0; P01 = 0; P10 = 0; P11 = 0;
				H00 = 0; H01 = 0; H10 = 0; H11 = 0;
				HXY = 0; HY2X1 = 0;
				num2 = realcomm[c2].size;
				num = intersection_size(evalcomm[c1],realcomm[c2]);
				P11 = (num)/(double)N;
				P10 = (num1 - num)/(double)N;
				P01 = (num2 - num)/(double)N;
				P00 = (N - (num1 + num2 -num))/(double)N;

				if(P11 > 0) H11 = -P11 * log(P11);	else H11 = 0;
				if(P10 > 0) H10 = -P10 * log(P10);	else H10 = 0;
				if(P01 > 0) H01 = -P01 * log(P01);	else H01 = 0;
				if(P00 > 0) H00 = -P00 * log(P00);	else H00 = 0;

				HX1Y2 = (H11+H10+H01+H00) - HY[c2];
				HY2X1 = (H11+H10+H01+H00) - HX[c1];

				if( (H11+H00) > (H01+H10)){
					nbestc1[c1]++;
					if (HX1Y2 < minHXY[c1])
						minHXY[c1] = HX1Y2;
					nbestc2[c2]++;
					if (HY2X1 < minHYX[c2])
						minHYX[c2] = HY2X1;
				}
			}
		}

	//calculate H(X|Y) and H(Y|X)
		for(c1=0;c1<k;c1++){
			if(nbestc1[c1]>0 && HX[c1]>0)
				HXY = HXY + minHXY[c1]/HX[c1];
			else HXY = HXY +1;
		}
		HXY = HXY / (pop[i]->comm_n);

		for(c2=0;c2<lengthofrealcomm;c2++){
			if(nbestc2[c2]>0 && HY[c2]>0)
				HYX = HYX + minHYX[c2]/HY[c2];
			else HYX = HYX +1;
		}
		HYX = HYX /lengthofrealcomm;
		gNMI[i] = 1 - (HXY + HYX)/2;
		fprintf(fp,"%lf\n",gNMI[i]);
	}
	//calculate the average gNMI of all population
	for(i=0;i<popsize;i++){
		ave_gNMI += gNMI[i];
	}
	ave_gNMI = ave_gNMI / popsize;
	fprintf(fp,"%lf\n",ave_gNMI);
	fclose(fp);
	//free variables
	if(HX != NULL) {
		free(HX); HX = NULL;
	}
	if(HY != NULL) {
		free(HY); HY = NULL;
	}
	if(nbestc1 != NULL){
		free(nbestc1); nbestc1 = NULL;
	}
	if(nbestc2 != NULL){
		free(nbestc2); nbestc2 = NULL;
	}
	if(minHXY != NULL){
		free(minHXY); minHXY = NULL;
	}
	if(minHYX != NULL){
		free(minHYX); minHYX = NULL;
	}
	if(gNMI != NULL){
		free(gNMI); gNMI = NULL;
	}
		//free realcomm and evalcomm
	for(i=0;i<NUM_OF_COMMUNITY;i++){
		if(realcomm[i].comm != NULL) {
			free(realcomm[i].comm);
			realcomm[i].comm = NULL;
		}
	}
	free(realcomm);
	for(i=0;i<max_lengthofevalcomm;i++){
		if(evalcomm[i].comm != NULL) {
			free(evalcomm[i].comm);
			evalcomm[i].comm = NULL;
		}
	}
	free(evalcomm);

}

void label2comm(char *name,char *input_path){
    
    FILE *f,*fp;
    int i,j;
    char output_for_comm[500];
    //--------------read community--------------------
    if((f= fopen(name,"r"))==NULL){
        printf("open real community failed\n");
        return;
    }
    realcomm = (struct partition *) calloc(NUM_OF_COMMUNITY,sizeof(struct partition));
    for(i=0;i<NUM_OF_COMMUNITY;i++){
        realcomm[i].comm = (int *)calloc(MAX_NUM_COMMUNITY_REAL,sizeof(int));
        if(realcomm[i].comm == NULL) printf("realcomm - alloc memory failed");
        memset(realcomm[i].comm,0,MAX_NUM_COMMUNITY_REAL * sizeof(int));
        realcomm[i].size = 0;
    }
    
    lengthofrealcomm = 0;
    while (fscanf(f,"%5d%5d",&i,&j)==2){
        int num = realcomm[j-1].size;
        if(j>lengthofrealcomm) lengthofrealcomm=j;
        realcomm[j-1].comm[num] = i -1;
        realcomm[j-1].size++;
    }
    fclose(f);
    
    //--------------write community--------------------
    sprintf(output_for_comm,"%s.dat",input_path);
    fp = fopen(output_for_comm,"w");
    
    for(i=0;i<NUM_OF_COMMUNITY;i++){
        if(realcomm[i].size!=0)
            fprintf(fp,"#community %d size %d\n",i,realcomm[i].size);
    }
    fclose(fp);
}
