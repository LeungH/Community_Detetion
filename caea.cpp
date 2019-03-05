#include<iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "caea.h"
#include "random.h"
#include "individual.h"


using namespace std;

int chidu = 1;
int objetive_num=2;
long rnd_uni_init = -(long)(time(NULL) + 23) % 1377;
double minx=100000.0,maxx=0.0,miny=100000.0,maxy=0.0;
bool obj_diff_flag=false;
double obj_diff = 9999999.0;

extern int generation;
extern int verbose;
extern int popsize;

individual_t* CAEA::init_population(){
	vector<bool> sector_pop_flag;//set initial population as false(null)
	vector<individual_t> initial_pop;	//initial population
    obj_diff_flag=false;

	int n = network_size();
	int i,j;
    double diff_temp;
    
//    struct obj_diff {
//        int diff; //genotype
//        int count;
//    };
//    obj_diff *obj_diff_point = (obj_diff *) calloc(5, sizeof(obj_diff));
//    vector<obj_diff> obj_diff_vec;
//    obj_diff obj_diff_value;
//    for (i=0;i!=5;++i){
//        obj_diff_value->diff=10^i;
//        obj_diff_value->cout
//    }

	individual_t *pop = (individual_t *)calloc(popsize, sizeof(individual_t));

	int *perm = (int *)calloc(n, sizeof(int));
	for (i = 0; i!= n; ++i) perm[i] = i;

	//initial population
	for (int i = 0; i < popsize; i++){
		shuffle(perm, n);
		sector_pop_flag.push_back(false);
		initial_pop.push_back(decode(perm));
        
        if(minx>initial_pop[i]->obj[0]) minx=initial_pop[i]->obj[0];
        if(maxx<initial_pop[i]->obj[0]) maxx=initial_pop[i]->obj[0];
        if(miny>initial_pop[i]->obj[1]) miny=initial_pop[i]->obj[1];
        if(maxy<initial_pop[i]->obj[1]) maxy=initial_pop[i]->obj[1];
        
//        diff_temp=initial_pop[i]->obj[0]/initial_pop[i]->obj[1];
//        if(diff_temp<obj_diff) obj_diff=diff_temp;
	}
	free(perm);
    
    for (int i = 0; i < popsize; i++){
        initial_pop[i]->obj[0]=(initial_pop[i]->obj[0]-minx)/(maxx-minx);
        initial_pop[i]->obj[1]=(initial_pop[i]->obj[1]-miny)/(maxy-miny);
        
//        cout << "<first" << initial_pop[i]->obj[0]<<","
//                    << initial_pop[i]->obj[1] << ">" << endl;
    }
    obj_diff_flag=true;
    
//    //Min-Max Normalization
//    for (int i = 0; i < popsize; i++){
//        min_max_normalization(initial_pop[i],minx,maxx,miny,maxy);
//    }
    
    

	//ck1104
//    for (int i = 0; i != 20; i++) {
//        cout << "<" << initial_pop[i].obj[0]<<","
//            << initial_pop[i].obj[1] << ">" << endl;
//    }
//    for (int i = 0; i != 20; i++) {
//        cout << "< " << initial_pop[0].gene[i]<<" >"<<endl;
//    }

	//≥ı ºªØπ€≤Ïµ„°¢√™µ„°¢◊Ó¥Û÷µµ„°¢≤Œøºµ„
    
    vector<double> obj_value;
    for (int i = 0; i < objetive_num; i++) obj_value.push_back(initial_pop[0]->obj[i]);
	for (int i = 0; i < objetive_num; i++){
		anchor_point.push_back(obj_value);

        true_nadir_point.push_back(initial_pop[0]->obj[i]);//’‚¿Ô–≈œ¢÷ÿ∏¥¥Ê¥¢
		
		observer_point.push_back(anchor_point[i][i]);
	}

	//original
	if (chidu == 0){
		for (int i = 0; i < objetive_num; i++){
			reference_point.push_back(true_nadir_point[i] + 1e3 * (true_nadir_point[i] - observer_point[i]));//≤Œøºµ„£®Ω”Ω¸Œﬁ«Ó¥Ûµƒµ„£¨º∆À„±ﬂ‘µ◊”«¯”Ú“™”√µΩ£©
		}
	}
	//normalizing
	else{
		for (int i = 0; i < objetive_num; i++)	{
			if (anchor_point[abs(i - 1)][i] - anchor_point[i][i]==0){
				reference_point.push_back(true_nadir_point[i] + 1e3 * (true_nadir_point[i] - observer_point[i]));
			}
			else
				reference_point.push_back(true_nadir_point[i] + 1e3 * 
					(true_nadir_point[i] - anchor_point[i][i]) / (anchor_point[abs(i - 1)][i] - anchor_point[i][i]));//≤Œøºµ„£®Ω”Ω¸Œﬁ«Ó¥Ûµƒµ„£¨º∆À„±ﬂ‘µ◊”«¯”Ú“™”√µΩ£©
		}
	}
	for (int i = 0; i < popsize; i++){
		update_extreme_point(initial_pop[i]);
	}

	//first round: select the solutions
	//select the solutions which have the small conical area as the initialize the population
	//some solution maybe will not be selected because of the not small enough conical area, and 
	//they will be save in vector initial_remain_pop
	vector<individual_t> initial_remain_pop;
	for (int n = 0; n < popsize; n++){
		calc_sectorial_index(initial_pop[n],observer_point,objetive_num,popsize,anchor_point);
		if (sector_pop_flag[initial_pop[n]->sectorial_index] == false){
			sector_pop_flag[initial_pop[n]->sectorial_index] = true;
			pop[initial_pop[n]->sectorial_index] = initial_pop[n];
		}
		else{
			bool bReplaced;
			individual_t remain_ind = initial_pop[n];
			individual_t rest_ind = pop[initial_pop[n]->sectorial_index];
			bReplaced = pareto_hypervolume_compare_sectorial_grid(remain_ind, rest_ind);
			if (bReplaced){
				initial_remain_pop.push_back(rest_ind);
				pop[initial_pop[n]->sectorial_index] = initial_pop[n];
			}
			else{
				initial_remain_pop.push_back(initial_pop[n]);
			}
		}	
	}
	
	//second round: select the solutions
	//select the remain solutions
	for (int i = 0; i < popsize; i++){
		if (sector_pop_flag[i] == false){
			double b_min_angle, b_angle;
			int flag_index = 0;
			individual_t temp;
			int size = (int)initial_remain_pop.size();
			
			b_min_angle = initial_remain_pop[0]->sectorial_angle - double(i) / (popsize - 1);
			for (int j = 1; j < size; j++){
				b_angle = initial_remain_pop[j]->sectorial_angle - double(i) / (popsize - 1);
				if (b_angle < b_min_angle){
					b_min_angle = b_angle;
					flag_index = j;
				}
			}
			//select
			pop[i] = initial_remain_pop[flag_index];
			sector_pop_flag[i] = true;
			
			//pop back the selected solution
			temp = initial_remain_pop[flag_index];
			initial_remain_pop[flag_index] = initial_remain_pop[size - 1];
			initial_remain_pop[size - 1] = temp;
			initial_remain_pop.pop_back();
		}

		//ck1104
		/*cout << "ind " << i << "\t";		
		cout << "<" << pop[i].obj[0] << ","
			<< pop[i].obj[1] << ">" << endl;

		for (int ii = 0; ii < pop[i].comm_n; ii++)
			cout << pop[i].comm[ii]->size << "\t";
		cout << endl << "comm num " << pop[i].comm_n << endl;
		cout << "id£∫" << pop[i].id << endl;
		system("pause");*/
		/*
		cout << "sectorial_index£∫" << pop[i].sectorial_index << endl;
		cout << "sectorial_angle£∫" << pop[i].sectorial_angle << endl;
		cout << "gene:" << endl;
		for (int ii = 0; ii < net.network_size; ii++)
		{
			cout << pop[i].gene[ii] << "\t";
		}
		cout << endl;
		cout << "comm:" << endl;*/
		
	}

	//clear
	initial_remain_pop.clear();
	vector<individual_t>(initial_remain_pop).swap(initial_remain_pop);
	sector_pop_flag.clear();
	vector<bool>(sector_pop_flag).swap(sector_pop_flag);
//    for (int i = 0; i != 20; i++) {
//        cout << "< 78 " << pop[78].gene[i]<<" >"<<endl;
//    }
	return pop;
}

void CAEA::evolution(individual_t *pop){

	int parent_index1, parent_index2;

//    individual_t *child  =(individual_t*)calloc(1,sizeof(individual_t));
//    individual_t child2 =(individual_t)calloc(1,sizeof(individual_t));
//    child.gene = (int* )calloc(network_size(), sizeof(int));
//    child2->gene = (int* )calloc(network_size(), sizeof(int));
    int i,j;
    individual_t child;
    
	int iL = 0,len = ceil(popsize / 10.0);
	//population evolution
	for (i = 0; i < generation; ++i){
        cout<<"---------"<<i<<setw(3)<<"/"<<generation<<"--------\n";

		//choose two parents
        iL = iL % len;
        if (iL == 0)
            parent_index1 = 0;
        else if (iL == 1)
            parent_index1 = popsize - 1;
        else
            parent_index1 = tour_selection_hv(pop);
        do parent_index2 = rand() % popsize; while (parent_index2 == parent_index1);
        iL++;
//        parent_index2 = tour_selection_hv(pop);
        cout<<"parent_index1,parent_index2 "<<parent_index1<<","<<parent_index2<<endl;
//        for (int i=0; i<20; i++)
//        {
//            cout<<(*pop)[parent_index1].gene[i]<<"   "<<(*pop)[parent_index2].gene[i]<<endl;
//        }

        child = crossover(pop[parent_index1],pop[parent_index2]);
//        min_max_normalization(child,minx,maxx,miny,maxy);
		child = mutation(child);
//        min_max_normalization(child,minx,maxx,miny,maxy);

		//recalculate the index of
		//individual in population
		if (update_extreme_point(child)){
			for (int i = 0; i < popsize; i++)
				calc_sectorial_index(pop[i],observer_point,
					objetive_num, popsize, anchor_point);
		}

		calc_sectorial_index(child,observer_point,
			objetive_num, popsize,anchor_point);

		update_pop(pop, child);

		free_individual(&child);
	}
}

/*dump all pops in population*/
void CAEA::dump_population(char *network_file,individual_t * pop,int m) {
	int i,*set_for_output_community;
	char of[1000];

	/*comm_of = (char *)calloc(strlen(community_file)+6*sizeof(int)+10,sizeof(char));
	if(comm_of == NULL) printf("comm_of - alloc memory failed");
	memset(comm_of,0,strlen(community_file)+6*sizeof(int)+10);*/

	set_for_output_community = (int *)calloc(orig_network_size(),sizeof(int));
	if(set_for_output_community == NULL) printf("set_for_output_community alloc memory failed");
	memset(set_for_output_community,0,orig_network_size()*sizeof(int));

	for (i=0; i!=popsize; ++i) {
		if (verbose>=1) population_info(pop[i]);
		printf("-------------saving  %d  population-----------\n",i+1);
		//generate output filename
		sprintf(of,"%s%03d_%02d.dat",network_file,i,m);
		//sprintf(comm_size,"%s_comm_size_%04d_%02d.dat",network_file,i,m);
		dump_individual(of,pop[i],set_for_output_community);
	}

	if(set_for_output_community != NULL){
		free(set_for_output_community); set_for_output_community = NULL;
	}

}

/*show a infomation of a pop
* format:
* [id](lambda) f=(f1,f2) #(number of communities)
* */
void CAEA::population_info(individual_t p) {
    printf("[%04d]",p->id);
//    printf("(%.5lf) ",p->lambda);
    printf("f=(%.2lf,%.2lf) ",p->obj[0],p->obj[1]);
    printf("#%-5d\n",p->comm_n);
}

/*dump a ind into a file*/
void CAEA::dump_individual(const char* name, individual_t ind,int *set) {
    int i,j;
    FILE* f = fopen(name,"w");
    //FILE *f_comm = fopen(comm_size,"w");
    fprintf(f,"#id      %d\n",ind->id);
    //    fprintf(f,"#lambda  %lf\n",ind->lambda);
    fprintf(f,"#comm_in  %lf\n",ind->obj[0]);
    fprintf(f,"#comm_out %lf\n",ind->obj[1]);
    fprintf(f,"#number  %d\n",ind->comm_n);
    //int i;
    //links_t p;
    for (i=0; i!=ind->comm_n; ++i) {
        int numofnode = 0;
        for(j=0;j<orig_network_size();j++) set[j] = 0;
        numofnode = convert_clique_to_node(ind->comm[i],set);
        //fprintf(f_comm,"%d\n",numofnode);
        fprintf(f,"#community %d size %d\n",i,numofnode);
        for(j = 0;j<numofnode;j++)    fprintf(f,"%d\n",set[j]+1);
    }
    fclose(f);
    //fclose(f_comm);
}

/*release population*/
void CAEA::free_population(individual_t * pop) {
	int i;
	for (i=0; i!=popsize; ++i) {
	free_individual(&(pop[i]));
	}
	free(pop);
}

//convert clique network to node network
int CAEA::convert_clique_to_node(community_t comm,int *set){
    int i,j,num=0;
    links_t n;
    
    for(i=0;i<num_of_clique_node(comm->head->to);i++){
        set[num] = node_of_clique(comm->head->to)[i];
        num++;
    }
    
    for(n = comm->head->next;n!=NULL;n=n->next){
        int temp = 0;
        for(i=0;i<num_of_clique_node(n->to);i++){
            int flag = 0;
            for(j = 0;j < num; j++){
                if(node_of_clique(n->to)[i] == set[j]) {flag = 1; break;}
            }
            if(flag == 1) continue;
            else {
                set[num+temp] = node_of_clique(n->to)[i]; temp++;
            }
        }
        num = num + temp;
    }
    //printf("\n转换之后的值： ");
    //for(i=0;i<num;i++) printf("%d  ",set[i]);
    return num;
}


//void CAEA::save_hv(char fileName[1024],vector<double> hv){
//	std::fstream fout;
//	fout.open(fileName, std::ios::out);
//	for (int i = 0; i < hv.size(); i++){
//		fout << i << "   " << hv[i] << "\n";
//	}
//	fout.close();
//}


bool CAEA::pareto_hypervolume_compare_sectorial_grid(individual_t ind1,individual_t ind2){
	double contribute1,contribute2;
	//individual_t ind2 = pop[ind1.sectorial_index];//the individual ind2 lies in this sectorial before individual ind1
    if (ind1->sectorial_index == ind2->sectorial_index){
        contribute1 = fastigiate_hyper_volume(ind1, ind1->sectorial_index, reference_point);
		contribute2 = fastigiate_hyper_volume(ind2, ind2->sectorial_index, reference_point);
		if (contribute1 < contribute2)	return true;
	}
	else {
		return true;
	}
	return false;
}

//use
bool CAEA::update_extreme_point(individual_t& ind){
	bool anchor_update = false;
	bool true_nadir_update = false;
	for (int i = 0; i < objetive_num; i++){
		if (ind->obj[i] < anchor_point[i][i] || ((ind->obj[i] == anchor_point[i][i])
			&& ind->obj[abs(i - 1)] < anchor_point[i][abs(i - 1)])){
			anchor_update = true;
			//anchor_point[i] = ind->obj;
            for (int j=0;j < objetive_num; j++) anchor_point[i][j] = ind->obj[j];
		}
		if (ind->obj[i] > true_nadir_point[i]){
			true_nadir_update = true;
			true_nadir_point[i] = ind->obj[i];
		}
	}
	for (int j = 0; j < objetive_num; j++){
		if (anchor_update){
			observer_point[j] = anchor_point[j][j];
		}
		if (anchor_update || true_nadir_update){
			if (chidu == 0){
				reference_point[j] = (true_nadir_point[j] + 1e3 * (true_nadir_point[j] - observer_point[j]));
			}
			else{
				if (anchor_point[abs(j - 1)][j] - anchor_point[j][j] == 0)
					reference_point[j] = (true_nadir_point[j] + 1e3 * (true_nadir_point[j] - observer_point[j]));

				else
					reference_point[j] = true_nadir_point[j] + 1e3 * (true_nadir_point[j] - 
						anchor_point[j][j]) / (anchor_point[abs(j - 1)][j] - anchor_point[j][j]);
			}		
		}
	}	
	return anchor_update;
}

void CAEA::reset_angle(individual_t* pop){

	for (int i = 0; i < popsize; i++)
		calc_sectorial_index(pop[i],observer_point,objetive_num,popsize,anchor_point);
}

int CAEA::tour_selection_hv(individual_t* population)
{
//    int p1 = int(rnd_uni(&rnd_uni_init)*popsize);
//    int p2 = int(rnd_uni(&rnd_uni_init)*(popsize - 1));
    
    int p1 = int(unirand()*popsize);
    int p2 = int(unirand()*(popsize - 1));
	if (p2 >= p1) p2++;

	double hv1 = tour_selection_difference(p1, population);
	double hv3 = tour_selection_difference(p2, population);
	if (hv1 >= hv3)
		return p1;
	else
		return p2;
}

double CAEA::tour_selection_difference(int p, individual_t* population)
{
	int num = 0;
	double hv1 = fastigiate_hyper_volume(population[p], population[p]->sectorial_index, reference_point);
	double hv3, hv4, hv_difference = 0.0;
	if (p - 1 >= 0)
	{
		hv3 = fastigiate_hyper_volume(population[p - 1], population[p - 1]->sectorial_index, reference_point);
		hv_difference += hv3 - hv1;
		num++;
	}
	if (p + 1 <= popsize - 1){
		hv4 = fastigiate_hyper_volume(population[p + 1], population[p + 1]->sectorial_index, reference_point);
		hv_difference += hv4 - hv1;
		num++;
	}
	hv_difference = hv_difference / num;

	return hv_difference;
}

double CAEA::fastigiate_hyper_volume(individual_t ind, int index, vector<double> &reference_point){
    double fastigVolume;
    double *normalizedf = (double* )calloc(objetive_num,sizeof(double));
    if (chidu == 0){
        for (int i = 0; i < objetive_num; i++)
            normalizedf[i] = ind->obj[i] - observer_point[i];
    }
    else{
        for(int i = 0; i < objetive_num; i++){
            double temp = anchor_point[abs(i - 1)][i] - anchor_point[i][i];
            if ( temp == 0){
                normalizedf[i] = ind->obj[i] - observer_point[i];
            }
            else{
                normalizedf[i] = (ind->obj[i] - anchor_point[i][i]) / temp;
            }
        }
    }
    
    if (index == 0)
        fastigVolume = (0.5*(index + 0.5) / (popsize - index - 1.5))*normalizedf[1] * normalizedf[1]
        + (reference_point[1] - normalizedf[1])*normalizedf[0];
    else if (index == popsize - 1)
        fastigVolume = (0.5*(popsize - index - 0.5) / (index - 0.5))*normalizedf[0] * normalizedf[0]
        + (reference_point[0] - normalizedf[0] * normalizedf[1]);
    else
        fastigVolume = (0.5*(popsize - index - 0.5) / (index - 0.5))*normalizedf[0] * normalizedf[0]
        + (0.5*(index + 0.5) / (popsize - index - 1.5))*normalizedf[1] * normalizedf[1]
        - normalizedf[0] * normalizedf[1];
    
    free(normalizedf);
    return fastigVolume;
}


/*convert a permutation into a individual*/
individual_t CAEA::decode(int* gene) {
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


//void CAEA::save_pop_front(vector<vector<double>> pop_front, char *saveFilename)
//{
//    std::fstream fout;
//    fout.open(saveFilename, std::ios::out);
//    for (int n = 0; n<pop_front.size(); n++)
//    {
//        if (_3f0 == 1)
//        {
//            fout << (1.0 - pop_front[n][0]) / 3.0 << "  ";
//        }
//        else
//        {
//            fout << (1.0 - pop_front[n][0])<< "  ";
//        }
//        fout <<  pop_front[n][1] << "  ";
//        fout << "\n";
//    }
//    fout.close();
//}


void CAEA::update_pop(individual_t *pop, individual_t child){

	if (pareto_hypervolume_compare_sectorial_grid(child, pop[child->sectorial_index])) {
		child->id = pop[child->sectorial_index]->id;
		child->refcount++;
		pop[child->sectorial_index]->refcount = 0;
		pop[child->sectorial_index] = child;
		//cout << "pop[child->sectorial_index] ref  " << pop[child->sectorial_index].ref << endl;
		population_info(pop[child->sectorial_index]);
	}
}

////what
//int CAEA::comparechild(individual_t child1, individual_t child2)
//{
//    int num1 = 0, num2 = 0, num3 = 0;
//    for (int i = 0; i < objetive_num; i++){
//        if (child1.obj[i]<child2.obj[i]){
//            num1++;
//        }
//        if (child1.obj[i] > child2.obj[i]){
//            num2++;
//        }
//        if (child1.obj[i] == child2.obj[i]){
//            num3++;
//        }
//    }
//    if (num3 == objetive_num){
//        return 2;
//    }
//    if (num1 + num3 == objetive_num){
//        return 1;
//    }
//    if (num2 + num3 == objetive_num){
//        return -1;
//    }
//    else{
//        return 0;
//    }
//}

//what
void CAEA::set_value(int &a, int &a1, int &a2)
{
	if (a == 0)
	{
		a1 = 2;
	}
	else if (a == popsize - 1)
	{
		a2 = popsize - 3;
	}
	else
	{
		a1 = a - 1;
		a2 = a + 1;
	}
}

