#include "read_network.h"

struct NetworkOrig{
	int n;  
	int nEdge;
	double *AdjMatrix;
	int *degree;
	links_t *Neighbor;
};

static struct NetworkOrig *net_orig;
extern int verbose;
static double weight_sum = 0;

char * strcatstr(char *str1,char *str2){

	char *temp = (char *) calloc(strlen(str1)+strlen(str2)+10,sizeof(char));
	if(temp == NULL) printf("temp - alloc memory failed");
	memset(temp,0,(strlen(str1)+strlen(str2)+2)*sizeof(char));

	strcpy(temp,str1);
	strcat(temp,str2);
	return temp;

	if(temp != NULL){
		free(temp);
		temp = NULL;
	}

}

int orig_network_size(){
	return net_orig->n;
}

int orig_network_nEdge(){
	return net_orig->nEdge;
}

int orig_network_degree(int i){
	return net_orig->degree[i];
}

double orig_network_AdjMatrix(int i,int j){
	return net_orig->AdjMatrix[i*net_orig->n+j];
}

links_t orig_network_Neighbor(int i){
	return net_orig->Neighbor[i];
}

void orig_neighbor(int i,int j) {
	links_t t = (links_t) malloc(sizeof(struct link));
	t->to = j;
	t->next = net_orig->Neighbor[i];
	net_orig->Neighbor[i] = t;
}


void read_pajek_unweighted(const char* name) {
	int i,j,n;
	int en=0; //number of edges
	FILE *f;
	net_orig = (struct NetworkOrig *) malloc(sizeof(struct NetworkOrig));

	if((f= fopen(name,"r"))==NULL){
		printf("open network file failed");
		return;
	}
	if(verbose>=1)
		printf("reading network : %s\n",name);

	fscanf(f,"%5d",&n);//ignore '*Vertices'

	if(verbose>=2) printf("vertices : %d\n",n);
	net_orig->n = n;

	net_orig->AdjMatrix = (double *) calloc(n*n,sizeof(double));
	net_orig->degree = (int *) calloc(n,sizeof(int));
	net_orig->Neighbor =(links_t *)calloc(n,sizeof(struct link));

	while (fscanf(f,"%d%d",&i,&j)==2){
		en ++;
		orig_neighbor(i-1,j-1);
		net_orig->AdjMatrix[(i-1)*n+(j-1)] = 1;
		net_orig->degree[i-1]++;
	}
	net_orig->nEdge = en;
	if(verbose>=2) printf("edges : %d\n",en);
	printf("\n");
}

// read social network
void read_pajek_unweighted_social(const char* name) {
	int i,j,n;
	int en=0; //number of edges
	FILE *f,*fp;
	net_orig = (struct NetworkOrig *) malloc(sizeof(struct NetworkOrig));
	
	// add for degree info of the network
	int d_max=0,d_min=1000;
	double d_sum=0;
	char *output_degree_info = (char *) calloc(strlen(name)+10,sizeof(char));
	sprintf(output_degree_info,"%s.info",name);
	fp = fopen(output_degree_info,"w");
	
	if((f= fopen(name,"r"))==NULL){
			printf("open network file failed\n");
			return;
	}

	if(verbose>=1)
		printf("reading network : %s\n",name);

	fscanf(f,"%d",&n);//ignore '*Vertices'

	if(verbose>=2) printf("vertices : %d\n",n);
	net_orig->n = n;

	net_orig->AdjMatrix = (double *) calloc(n*n,sizeof(double));
	net_orig->degree = (int *) calloc(n,sizeof(int));
	net_orig->Neighbor =(links_t *)calloc(n,sizeof(struct link));

	while (fscanf(f,"%d%d",&i,&j)==2) {
		en ++;
		orig_neighbor(i-1,j-1);
		orig_neighbor(j-1,i-1);
		net_orig->AdjMatrix[(j-1)*n+(i-1)] = 1;
		net_orig->AdjMatrix[(i-1)*n+(j-1)] = 1;
		net_orig->degree[i-1]++;
		net_orig->degree[j-1]++;
	}
	net_orig->nEdge = en;
	if(verbose>=2) printf("edges : %d\n",en);
	if (verbose>=1) printf("\n");
	for(int i=0;i<n;++i){
		if (net_orig->degree[i]>d_max) d_max=net_orig->degree[i];
		if (net_orig->degree[i]<d_min) d_min=net_orig->degree[i];
		d_sum+=net_orig->degree[i];
	}
	fprintf(fp,"%s & %f & %d & %d &\n",name,d_sum/n,d_min,d_max);
	fclose(fp);
}

void free_links(links_t l) {
	links_t t;
	while (l) {
	t = l;
	l = l->next;
	free(t);
	}
}
// free original network
void free_orig_network(){
	int i;
	free(net_orig->AdjMatrix);
	free(net_orig->degree);
	for (i=0; i!=net_orig->n; ++i)
		free_links(net_orig->Neighbor[i]);
	free(net_orig->Neighbor);
	free(net_orig);
	net_orig = NULL;
}
