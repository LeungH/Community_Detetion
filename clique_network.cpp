#include "clique_network.h"


struct clique_network{
	int size_n; //the number of clique nodes
	struct clique_node_info *clique;
	double *Adj;  //adjacent matrix of maximal-clique graph
};

struct clique_node_info{
	int k;  //the number of nodes of maximal clique
	int *node; // the nodes of maximal clique
};

struct network{
	int n; // the number of nodes of the maximal-clique graph
	double *w,*ls; // w: the link strength between different clique nodes; 
				  //ls: the sum of link strength for each node
	double *degree; // the degree for each node
	links_t *e;  // the neighbors for each node
};

struct clique_network *net_clique;
/*global network data*/
static struct network* net;

extern int CliqueCount;//the number of max clique
double para_node = 0.33333; // the parameter for overlapping nodes between two communities
double para_edge_overlapping = 0.33333;// the parameter for overlapping edges between two communities
double para_edge_joint = 0.33333; // the parameter for joint edges between two communities
static double weight_sum_orig;  // the sum of edges in orginal network
int *label;
double clique_ave_weight;  // the average weight of edges for clique network

extern int verbose;

int num_of_clique_node(int i){
	return net_clique->clique[i].k;
}

int * node_of_clique(int i){
	return net_clique->clique[i].node;
}

// calculate the sum of link strength of original network
double clique_weight(){
	int i,j;
	double temp = 0;
	for(i=0;i<orig_network_size();i++)
		for(j=i+1;j<orig_network_size();j++)
			temp += orig_network_AdjMatrix(i,j);
	return temp;
}
// calculate the overlapping nodes between different maximal cliques
int clique_Node(int * a,int size_a,int * b,int size_b,int *node){
	int num=0;
	int i,j;
	for(i=0;i<size_a;i++)
		for(j=0;j<size_b;j++){
			if(a[i] == b[j]){
				node[num] = a[i];
				label[a[i]] = 1;
				num++; break;
			}
		}
	return num;
}

// calculate the overlapping edges and joint edges between different maximal cliques
double clique_Edge(int *a,int size_a,int *b,int size_b, int inter_num,int *node){
	int i,j;
	int flag = 0;
	double edge1 = 0,edge2 =0;
	/*printf("\n");
	for(i=0;i<inter_num;i++) printf("%d  ",node[i]);*/

	// divide into two conditions: 
	// 1. overlapping edges 
	// 2. joint edges
	for(i=0;i<size_a;i++){
		if(label[a[i]] == 1) continue;
		for(j=0;j<size_b;j++){
			if(label[b[j]] == 1) continue;
			flag = 0;
			edge1 += orig_network_AdjMatrix(a[i],b[j]);
		}
	}
	//
	for(i=0;i<inter_num;i++)
		for(j=i+1;j<inter_num;j++)
			edge2 += orig_network_AdjMatrix(node[i],node[j]);

	return para_edge_joint * edge1 + para_edge_overlapping * edge2;

}

void free_net_clique(){

	free(net_clique->Adj);
	free(net_clique->clique);
	free(net_clique);
	net_clique = NULL;
}

// construct the maximal-clique graph from the original graph
int construct_clique_network(){
	int i,j,n=0,num,p,en_original = 0;
	double node_rate,edge_rate;

	int *intersection = (int *)malloc(orig_network_size()*sizeof(int));
	if(intersection == NULL) printf("intersection - alloc memory failed");
	memset(intersection,0,orig_network_size()*sizeof(int));

	label = (int *)malloc(orig_network_size()*sizeof(int));
	if(label == NULL) printf("label - alloc memory failed");
	memset(label,0,orig_network_size()*sizeof(int));

	net_clique = (struct clique_network *) malloc (sizeof(struct clique_network));

	net_clique->size_n = 0;
	for(i=0;i<CliqueCount;i++) {
		net_clique->size_n += MaxClique_size(i);
	}

	num = net_clique->size_n;  // the sum number of maximum cliques

	weight_sum_orig = clique_weight();

	net_clique->Adj = (double *)malloc(num * num *sizeof(double));
	if(net_clique->Adj == NULL) printf("net_clique->Adj alloc memory failed\n");
	memset(net_clique->Adj,0,num * num *sizeof(double));

	net_clique->clique = (struct clique_node_info *) malloc(num * sizeof(struct clique_node_info));

	// determine the clique nodes for maximal-clique graph
	for(i=0;i<CliqueCount;i++){
		if(MaxClique_size(i) == 0) continue;
		for(j=0;j<MaxClique_size(i);j++){
			net_clique->clique[n].k = i+1;
			net_clique->clique[n].node = MaxClique_Clique(i,j);
			n++;
		}
	}

	// evaluate the link strength between different maximal clique graph
	for(i = 0;i<num;i++){
		int intersection_node = 0;
		double intersection_edge = 0;
		if(net_clique->clique[i].k == 1) continue;
		for(j = i+1;j<num;j++){
			for(p=0;p<orig_network_size();p++) { intersection[p] = 0; label[p] = 0;}
			// the ratio of overlapping nodes
			intersection_node = clique_Node(net_clique->clique[i].node, net_clique->clique[i].k,
											net_clique->clique[j].node, net_clique->clique[j].k,intersection);
			node_rate = para_node * intersection_node /(net_clique->clique[i].k + net_clique->clique[j].k-intersection_node);
			// the ratios of overlapping edges and joint edges
			intersection_edge = clique_Edge(net_clique->clique[i].node, net_clique->clique[i].k,
											net_clique->clique[j].node, net_clique->clique[j].k,
											intersection_node,intersection)/weight_sum_orig;
			edge_rate = intersection_edge;
			net_clique->Adj[i*num + j] = node_rate + edge_rate;
		}
	}

	clique_ave_weight = 0;
	//calculate the average of network
	for(i=0;i<net_clique->size_n;i++)
		for(j=i+1;j<net_clique->size_n;j++){
			double w = net_clique->Adj[i * net_clique->size_n + j];
			if(w>0){
				en_original++;
				clique_ave_weight += w;
			}
	}
	clique_ave_weight = clique_ave_weight/en_original;

	if(intersection != NULL){
		free(intersection);
		intersection = NULL;
	}

	if(label != NULL){
		free(label);
		label = NULL;
	}

	return 0;
}


/*nubmer of vertices*/
int network_size() {
	return net->n;
}

/*neighbors, return a link list*/
links_t neighbor(int i) {
	return net->e[i];
}

/*get W(i,j) weight of edge from i to j*/
double weight(int i,int j) {
	return net->w[i*(net->n) + j];
}

// the sum of links for each clique nodes
double link_sum(int i) {
	return net->ls[i];
}


double degree(int i) {
	return net->degree[i];
}

/*free a linklist*/

/*free all the network*/
void free_network() {
	int i;
	free(net->ls);
	free(net->w);
	free(net->degree);
	for (i=0; i!=net->n; ++i)
		free_links(net->e[i]);
		//free(net->e[i]);
	free(net->e);
	free(net);
	net = NULL;
}


int is_network_ok() {
	return net!=NULL;
}

//add a edge from i to j with weight w to network, this function only add a directed link from i to j, but not j to i.
// the network is undirected, add_link(i,j,w), add_link(j,i,w);
static void add_link(int i,int j,double w) {
	links_t t = (links_t) malloc(sizeof(struct link));
	t->to = j;
	t->next = net->e[i];
	net->e[i] = t;
	net->w[i*(net->n) + j] = w;
}

// calculate the sum of links for each clique node
static double* sum_link() {
	int n = net_clique->size_n;
	int i,j;
	double x;
	links_t l;
	double *s = (double *) calloc(n,sizeof(double));
	if(s == NULL) printf("s - alloc memory failed\n");
	memset(s,0,n*sizeof(double));
	for (i=0; i!=n; ++i)
		for (l=net->e[i]; l!=NULL; l=l->next) {
			j = l->to;
			x = weight(i,j);
			s[i] += x;
		}
	return s;
}


/* read a pajek file */
void parameter_of_clique_network() {
	int n;
	int i,j;
	double w;
	int en_threshold = 0; //number of edges
	net = (struct network *) malloc(sizeof(struct network));
	n = net_clique->size_n;
	if(verbose>=2) printf("vertices : %d\n",n);

	net->n = n;
	net->w = (double *) calloc(n*n,sizeof(double));
	net->degree = (double *) calloc(n,sizeof(double));
	net->e = (links_t *) calloc(n,sizeof(links_t));

	for(i=0;i<net_clique->size_n;i++)
		for(j=i+1;j<net_clique->size_n;j++){
			w = net_clique->Adj[i * net_clique->size_n + j];
			if(w < clique_ave_weight) net_clique->Adj[i * net_clique->size_n + j] = 0;

			else{
				add_link(i,j,w);
				add_link(j,i,w);
				net->degree[i]++;
				net->degree[j]++;
				en_threshold++;
			}
	}

	if(verbose>=2) printf("edges : %d\n",en_threshold);
	if(verbose>=2) printf("ave_weight : %lf\n",clique_ave_weight);
	net->ls = sum_link();

}

// output the information for maximal-clique graph 
void output_clique_network(char *name){
	int i,j;
	FILE *fp1,*fp2;
//	FILE *fp3;
	int n = net_clique->size_n;

	char of_1[1000],of_2[1000];

	sprintf(of_1,"%sclique_network_ori.dat",name);
	sprintf(of_2,"%sclique_network_cut.dat",name);

	fp1 = fopen(of_1,"w");
	fp2 = fopen(of_2,"w");

	printf("\n--------------saving network---------------\n");
	for(i=0;i<net_clique->size_n;i++){

/*
		fprintf(fp3,"%d: ",i+1);
		for(q=net_clique->clique[i].node;q!=NULL;q=q->next){
			fprintf(fp3,"%d  ",(q->to)+1);
		}
		fprintf(fp3,"\n");
*/
		for(j=0;j<net_clique->size_n;j++){
			if(net_clique->Adj[i * n + j] > 0)
				fprintf(fp1,"%d  %d  %lf\n",i+1,j+1,net_clique->Adj[i * n + j]);
			if(net_clique->Adj[i * n + j] > clique_ave_weight)
				fprintf(fp2,"%d  %d  %lf\n",i+1,j+1,net_clique->Adj[i * n + j]);
		}

	}
	fclose(fp1);fclose(fp2);
	//fclose(fp3);
}
