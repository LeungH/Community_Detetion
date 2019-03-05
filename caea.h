#pragma once
#include "random.h"
#include "individual.h"
#include <vector>
#include <ctime>
#include <iomanip>
#include <stdlib.h>

using std::vector;
using namespace std;

class CAEA{
public:
    
    individual_t decode(int* gene);

	//initial population
	individual_t* init_population();
    
    individual_t* new_individual(int*);
    

	void evolution(individual_t* );

    void free_population(individual_t ** pop);
    
    
    void update_pop(individual_t *pop, individual_t child);

    void dump_population(char *network_file,individual_t * pop,int m);
    void dump_individual(const char* name, individual_t ind,int *set);
    int  convert_clique_to_node(community_t comm,int *set);
    void population_info(individual_t p);
    
    int    tour_selection_hv(individual_t* population);
    double tour_selection_difference(int p, individual_t* population);
    double fastigiate_hyper_volume(individual_t  ind, int ind_index, vector <double> &reference_point);
    
    bool pareto_hypervolume_compare_sectorial_grid(individual_t ,individual_t );
    
	void save_population(vector<individual_t> &mypopulation, char saveFilename[1024],vector<individual_t>&);  // save the pareto front into files
	
	vector <individual_t>  ps;

	int best_ind_index;

   /* anchor_point   //
	* true_nadir_point	  max of objective function
	* observer_point
	* reference_point     an approximate infinity point dominated by all feasible solutions
	*/
	vector<vector <double>> anchor_point;    
	vector <double> true_nadir_point;     //useful
	vector <double> observer_point;
	vector <double> reference_point;	  //useful	

	
	bool update_extreme_point(individual_t& ind);
	void reset_angle(individual_t* );
	
	//return modularity_best_id_value
	void save_best_ind(individual_t** , int);

	//save complete information of population 
	void save_pop_all(vector<individual_t> mypopulation);

	//modularity-based tournament selection
	int  tour_selection_modularity(individual_t** );

	

	void population2front(vector <individual_t>  &mypopulation, vector <vector<double>> &population_front, int);//vector <CSNGAInd>

	int roulette_neighbor(int i);

	void save_pop_front(vector< vector<double> >,char *);

	void set_individual(individual_t *d, individual_t *s);

	


	void free_population(individual_t* pop);

	

	int comparechild(individual_t,individual_t);
	void set_value(int &,int &,int &);
	
};
