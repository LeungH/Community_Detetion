#include <string.h>
#include "maxcliques.h"

#define MAX_COUNT_CLIQUE_EVERY_SIZE 50000
#define MAX_CLIQUE_NUM 2000

typedef int * pint;

typedef struct clique{
	int k;    //the size of clique is k
	int num_of_cliques; //the number of k-cliques
	int **MaxClique; //
}*cliques;

static struct clique *CLP;
static int *flag_node; //if 1，the node has been assigned to maximal clique;
					  // if 0, the node has not been assigned to maximal clique;
static int *flag_for_whole_clique;   // if 1, delete this node
int CliqueCount;
int *setA;
int *setB;
int **label_avoid_same_clique;


int MaxClique_size(int k){
	return CLP[k].num_of_cliques;
}

int * MaxClique_Clique(int k,int m){
	return CLP[k].MaxClique[m];
}

int * BubbleSort(int * node_degree_sequence){
	int i,j,temp;
	for(i=0;i<orig_network_size();i++)
		for(j=0;j<orig_network_size()-i-1;j++){
			if(orig_network_degree(node_degree_sequence[j]) > orig_network_degree(node_degree_sequence[j+1])){
				temp = node_degree_sequence[j];
				node_degree_sequence[j]=node_degree_sequence[j+1];
				node_degree_sequence[j+1]=temp;
			}
		}
		return node_degree_sequence;
}

void oneclique(int x){
	//calculate CLP
	int m = CLP[0].num_of_cliques;
	CLP[0].MaxClique[m][0] = x;
	CLP[0].num_of_cliques++;
	//label the node to 1, which indicate the node has been assigned to maximal cliques
	flag_for_whole_clique[x] = 1;
}

int intersection(int *setA,int n,int num,int *setB,int * num_set,int num_setB){
	int i,flag;
	links_t j;
	if(n == -1) return 0;
	setA[num_set[0]] = n;
	num_set[0]++;
	setB[num] = -1;
	num_set[1]--;
	if(num_set[1]==0) return 0;
	for(i = 0;i<num_setB;i++){
		int x = setB[i];
		flag=0;
		if(setB[i] == -1) continue;
		for(j = orig_network_Neighbor(n);j!=NULL;j=j->next){
			if(j->to == x) {
				flag = 1;
				break;
			}
			if(j->to < x) break;
		}
		if(flag==0) {
			setB[i] = -1;
			num_set[1]--;
		}
	}

	return 1;
}

// this function is used to find whether the node could consitute the k-clique 
int Clique(int node,int sizeofclique,int *flag_node_lower_std){ 
	int i,num1 = 0,num2 = 0; 
	links_t m,n,q;
	// A basic condition for one node to consitute k-clique 
	// is that the number of its neighbors whose degree more than 
	// k must more than k-1
	int k = orig_network_degree(node) - sizeofclique + 1; 

	// cut 1: less than (k-1) adjacent nodes have a degree no smaller than (k−1)
	for(m = orig_network_Neighbor(node);m != NULL;m = m->next){
		if(num1>k)  {return 0;}
		if(flag_node[m->to]==1) {
			num1++; continue;
		}
		if(orig_network_degree(m->to) < sizeofclique-1) {
			flag_node_lower_std[m->to] = 1;  // if 1, this node is not suited this clique
			num1++;
			if(num1>k) { return 0; }
		}
	}

	for(m= orig_network_Neighbor(node);m!=NULL;m=m->next){
		int x = m->to;
		int num_intersection_node = 0;
		if(flag_node_lower_std[m->to]==1 || flag_node[m->to]==1)  {num2++; continue;}
		for(n = orig_network_Neighbor(x);n!=NULL;n=n->next){
			if(flag_node_lower_std[n->to]==1 || flag_node[n->to]==1 || n->to == node)  continue;
			for(q = orig_network_Neighbor(node);q!=NULL;q=q->next){
				if(flag_node_lower_std[q->to]==1 || flag_node[q->to]==1)  continue;
				if(n->to==q->to) {
					num_intersection_node++;break;
				}
				if(n->to>q->to) break;
			}
		}
		if(num_intersection_node<sizeofclique-2) {num2++;flag_node_lower_std[x]=1;}
		if(num2>k) return 0;
	}
	/*
	printf("\n");
	for(i=0;i<network_size();i++){
		printf("the label of node %d is: %d  \n",i,flag_node_lower_std[i]);
	}*/

	//if the above consition is satified, find the k-clique for node
	// divide into two conditions:
	// 1. the node could only consitute one k-clique
	// 2. the node could consitute multiple k-cliques

	// the first condition
	if(k==num2){
		int i,flag=0,num=0,temp;
		links_t m,n,a;
		int Count1, Count2;  
		for(m = orig_network_Neighbor(node);m != NULL;m = m->next){
			if(flag_node_lower_std[m->to] == 1 || flag_node[m->to]==1) continue;
			for(n = orig_network_Neighbor(node);n!=NULL;n=n->next){
				if(flag_node_lower_std[n->to]==1 || flag_node[n->to]==1 ||m->to==n->to) continue;
				//printf("\nthe adjacient of node %d and node %d is %lf",m->to,n->to,network_AdjMatrix(m->to,n->to));
				if(orig_network_AdjMatrix(m->to,n->to)==0) return 0;
			}
		}

		//printf("the number of %d-clique is: %d",sizeofclique,CLP[sizeofclique-1].num_of_cliques);
		Count1 = 1;  
		Count2 = CLP[sizeofclique-1].num_of_cliques;
		CLP[sizeofclique-1].MaxClique[Count2][0] = node;
	//	printf("\n--------------insert test-----------\n");
	//	printf("%d  ",CLP[sizeofclique-1].MaxClique[Count2][0]+1);
		for(m = orig_network_Neighbor(node);m!=NULL;m=m->next){
			if(flag_node[m->to]==1 || flag_node_lower_std[m->to]==1) continue;
			CLP[sizeofclique-1].MaxClique[Count2][Count1] = m->to;
		//	printf("%d  ",CLP[sizeofclique-1].MaxClique[Count2][Count1]+1);
			Count1++;
		}
		CLP[sizeofclique-1].num_of_cliques++;
		flag_node[node] = 1;
		return 1;
	}

	// the second condition
	else{
			int i,j,num_exist_clique = 0,num_node_for_set[2];
			links_t m,n,q;

			//intilize setA and setB
			for(m = orig_network_Neighbor(node);m != NULL;m=m->next){
				int num = 0,p=0,flag = 0,flag_avoid_same_clique=0,total_num_setB=0;
				if(flag_node_lower_std[m->to] == 1 || flag_node[m->to] == 1) continue;
				
				num_node_for_set[0] = 2;
				setA[0] = node;
				setA[1] = m->to;
				
				for(n = orig_network_Neighbor(node);n!=NULL;n=n->next){
					if( flag_node_lower_std[n->to]==1 || n->to==m->to ||flag_node[n->to] == 1) continue;
					for(q=orig_network_Neighbor(m->to);q!=NULL;q=q->next){
						if(n->to>q->to) break;
						if(n->to==q->to){ setB[p] = n->to;p++;break;}
					}
				}
				num_node_for_set[1] = p; //the number of nodes in setA
				total_num_setB = p;

				while(num_node_for_set[0] <sizeofclique && num_node_for_set[1] >0){
					intersection (setA,setB[num],num,setB,num_node_for_set,total_num_setB);
					num++;
				}

				if(num_node_for_set[0] == sizeofclique){
					// detect whether there is overlaps between different cliques
					for(i = 0;i<num_exist_clique;i++){
						int flag_clique = 0;
						for(j = 0;j<sizeofclique;j++)
							flag_clique += label_avoid_same_clique[i][setA[j]];
						if(flag_clique == sizeofclique) {flag_avoid_same_clique = 1; break;}
					}

					//construct new clique
					if(flag_avoid_same_clique == 0){
						//************ test ************************
						int pp,qq;
						int Count1, Count2;
						/*
						for(qq = 0;qq<MAX_COUNT_CLIQUE_EVERY_SIZE;qq++){
							printf("CLP %d: ",qq);
							for(pp = 0;pp<sizeofclique; pp++){
								printf("%d ",CLP[sizeofclique-1].MaxClique[qq][pp]);
							}
							printf("\n");
						}*/
						//*******************************************
						//将新的派系加入到CLP和label_avoid_same_clique中
						Count1 = 1;
						Count2 = CLP[sizeofclique-1].num_of_cliques;
						CLP[sizeofclique-1].MaxClique[Count2][0] = node;
					   // printf("\n--------------insert test:  %d----------------\n",Count2);
					   // printf("%d  ",CLP[sizeofclique-1].MaxClique[Count2][0]+1);
						label_avoid_same_clique[num_exist_clique][node] =1;
						for(i = 1;i<sizeofclique;i++){
							CLP[sizeofclique-1].MaxClique[Count2][Count1] = setA[i];
						//    printf("%d  ",CLP[sizeofclique-1].MaxClique[Count2][Count1]+1);
							Count1++;
							label_avoid_same_clique[num_exist_clique][setA[i]] = 1;
						}
						CLP[sizeofclique-1].num_of_cliques++;
						num_exist_clique++;
						label_avoid_same_clique[num_exist_clique] = (int *)calloc(orig_network_size(),sizeof(int));
					}
				}
			}
			flag_node[node] = 1;
			return 1;
	}
}


int maxclique(char *name){
	int i,j,k,m,n,p,q;
	links_t a;
	int num_a = 0;
		char output_clique[1000];
	FILE *fp;
	int *node_degree_sequence, *flag_node_lower_std;

	node_degree_sequence = (int *)calloc(orig_network_size(),sizeof(int));
	if(node_degree_sequence == NULL) printf("node_degree sequence - alloc memory failed\n");
	memset(node_degree_sequence,0,orig_network_size()*sizeof(int));

	flag_node_lower_std = (int *)calloc(orig_network_size(),sizeof(int));
	if(flag_node_lower_std == NULL) printf("flag_node_lower_std - alloc memory failed\n");
	memset(flag_node_lower_std,0,orig_network_size()*sizeof(int));

	flag_for_whole_clique = (int *)calloc(orig_network_size(),sizeof(int));
	if(flag_for_whole_clique == NULL) printf("flag_for_whole_clique - alloc memory failed\n");
	memset(flag_for_whole_clique,0,orig_network_size()*sizeof(int));

	flag_node = (int *)calloc(orig_network_size(),sizeof(int));
	if(flag_node == NULL) printf("flag_node - alloc memory failed\n");
	memset(flag_node,0,orig_network_size()*sizeof(int));

	setA = (int *)malloc(orig_network_size() * sizeof(int));  //存储目标派系
	if(setA == NULL) printf("setA - alloc memory failed\n");
	memset(setA,0,orig_network_size() * sizeof(int));

	setB = (int *)malloc(orig_network_size() * sizeof(int));  //存储该节点的邻居节点
	if(setB == NULL) printf("setB - alloc memory failed\n");
	memset(setB,0,orig_network_size() * sizeof(int));

	label_avoid_same_clique = (int **)malloc(MAX_CLIQUE_NUM * sizeof(int *));
	for(i=0;i<MAX_CLIQUE_NUM;i++){
		label_avoid_same_clique[i] = (int *)calloc(orig_network_size(),sizeof(int));
		if(label_avoid_same_clique[i]==NULL) printf("label_avoid_same_clique - alloc memory failed\n");
		memset(label_avoid_same_clique[i],0,orig_network_size()*sizeof(int));
	}

	//sort node ascending by degree
	for(i=0;i<orig_network_size();i++){
		node_degree_sequence[i] =i;
	}
	node_degree_sequence = BubbleSort(node_degree_sequence);

	//initial CLP
	CliqueCount = orig_network_degree(node_degree_sequence[orig_network_size()-1])+1;

	CLP = (struct clique *)calloc(CliqueCount,sizeof(struct clique));
	if(CLP == NULL) printf("CLP - alloc memory failed\n");
	memset(CLP,0,sizeof(struct clique)*CliqueCount);

	for(i=0;i<CliqueCount;i++){
		CLP[i].k = i + 1;
		CLP[i].num_of_cliques = 0;
		CLP[i].MaxClique = (int **)malloc(MAX_COUNT_CLIQUE_EVERY_SIZE * sizeof(int *));
		for(j=0;j<MAX_COUNT_CLIQUE_EVERY_SIZE;j++){
			CLP[i].MaxClique[j]= (int *)calloc(i+1,sizeof(int));
			if(CLP[i].MaxClique[j] == NULL) printf("CLP[%d].MaxClique[%d] - alloc memory failed\n",i,j);
			memset(CLP[i].MaxClique[j], 0, (i+1) * sizeof(int));
		}
	}
/*
	for(i=0;i<CliqueCount;i=i+10){
		printf("clique %d : %d\n",CLP[i].k,CLP[i].num_of_cliques);
		for(j=0;j<MAX_COUNT_CLIQUE_EVERY_SIZE;j++){
			printf("the %d clique: ",j);
			for(k=0;k<i+1;k++){
				printf("%d ",CLP[i].MaxClique[j][k]);
			}
			printf("\n");
		}
	}
*/
    
//    printf("orig_network_degree... \n");
//    for(i=0;i<orig_network_size();i++){
//        printf("the degree of %d is %d\n",i,orig_network_degree(i));
//    }

	printf("after sorting ... \n");
	for(i=0;i<orig_network_size();i++){
		printf("the degree of %d is %d\n",node_degree_sequence[i],orig_network_degree(node_degree_sequence[i]));
	}

	//find 1-clique; if the degree of the node less than 1, then this node belong to
	m = 1;
	while(m>0){
		int x = orig_network_degree(node_degree_sequence[m-1]);
		printf("--------%d    %d\n",x,node_degree_sequence[m-1]);
		if(x > 0)  break;  //if the degree of nodes surpass 0, no 1-clique exist
		else if(x==0)  oneclique(node_degree_sequence[m-1]);  //construct 1-clique
		m++;
	}
	printf("\n");
	//寻找派系

	for(k=orig_network_degree(node_degree_sequence[orig_network_size()-1])+1;k>1;k--){
		printf("-------------%d-clique---------------\n",k);
		//set <flag_node> to zero
		for(p=0;p<orig_network_size();p++)
			flag_node[p]=0;
		//Detect cliques
		for(i = orig_network_size()-1; i >= 0; i--){
			int x = node_degree_sequence[i];
			if(orig_network_degree(x) < k-1) break;   // first condition: the degree of the node must more than k-1
			if(flag_for_whole_clique[x] == 1) continue; // second condition: the node does not found maximum cliques
			else{
				// label_avoid_same_clique is used to avoid two cliques contrain same nodes, but different sequence
				// flag_node_lower_std: used to label the node which is not suited for clique
				for(p=0;p<MAX_CLIQUE_NUM;p++)
					for(q=0;q<orig_network_size();q++)
						label_avoid_same_clique[p][q] = 0;
				for(p=0;p<orig_network_size();p++) {
					flag_node_lower_std[p] = 0;
					setA[p] = 0; setB[0] = 0;
					}
				//find the largest maximal clique of each node belongs to
				Clique(node_degree_sequence[i],k,flag_node_lower_std); 
			}
		}
		// if the node has been assigned to one maximal clique, it is labeled 
		// since only the largest maximal clique is considered here
		if(CLP[k-1].num_of_cliques>0){
			for(m=0;m<CLP[k-1].num_of_cliques;m++){
				printf("the %d clique: ",m+1);
				for(n=0;n<k;n++){
					flag_for_whole_clique[CLP[k-1].MaxClique[m][n]] = 1;
					printf("%d   ",CLP[k-1].MaxClique[m][n]);
				}
				printf("\n");
			}
		}
	}
	//ouput all cliques
	num_a = 1;
	sprintf(output_clique,"%sclique_node.dat",name);
	printf("output_clique = %s\n",output_clique);
	fp = fopen(output_clique,"w");
	for(k=0;k<orig_network_degree(node_degree_sequence[orig_network_size()-1]);k++){
		printf("the %d-clique-------number of clique: %d\n",k+1,CLP[k].num_of_cliques);
		for(m=0;m<CLP[k].num_of_cliques;m++){
			fprintf(fp,"%d: ",num_a);
			num_a ++ ;
			for(n=0;n<k+1;n++){
				fprintf(fp,"%d  ",CLP[k].MaxClique[m][n]+1);
				printf("%d  ",CLP[k].MaxClique[m][n]+1);
			}
			fprintf(fp,"\n");
			printf("\n");
		}
	}
	fclose(fp);
	if(node_degree_sequence!=NULL){
		free(node_degree_sequence);
		node_degree_sequence = NULL;
	}

	if(flag_for_whole_clique!=NULL){
		free(flag_for_whole_clique);
		flag_for_whole_clique = NULL;
	}

	if(flag_node!=NULL){
		free(flag_node);
		flag_node = NULL;
	}

	if(flag_node_lower_std!=NULL){
		free(flag_node_lower_std);
		flag_node_lower_std = NULL;
	}

	if(setA!=NULL){
		free(setA); setA=NULL;
	}

	if(setB!=NULL){
		free(setB); setB=NULL;
	}

	for(i=0;i<MAX_CLIQUE_NUM;i++){
		if(label_avoid_same_clique[i] != NULL){
			free(label_avoid_same_clique[i]);
			label_avoid_same_clique[i] = NULL;
		}
	}

	return 1;
}

int free_cliques(){
	int i,j;
	if(CLP == NULL) return 0;

	for(i=0;i<CliqueCount;i++){
		for(j=0;j<MAX_COUNT_CLIQUE_EVERY_SIZE;j++){
			if(CLP[i].MaxClique[j] != NULL) {
				free(CLP[i].MaxClique[j]);
				CLP[i].MaxClique[j] = NULL;
			}
		}
	}
	free(CLP);
	CLP = NULL;
	return 1;
}
