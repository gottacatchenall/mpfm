
#ifndef POPULATION_H
#define POPULATION_H

#include "include.h"

class Population{
    private:
        static int id_counter;
        int id;
        double x;
        double y;
        double K;
        MPFM* mpfm;
        std::vector<Individual*> migrant_queue;
        std::vector<Individual*>* individuals;
        std::vector<Individual*>* next_gen;
        std::vector<double> env_factors;
        double exp_num_off_this_population;
    public:
        Population(MPFM* run, double x, double y, double K);
        double get_x();
        double get_y();
        double get_K();
        int get_id();
        int get_size();
        void add_individual(Individual* indiv);
        void remove_individual(Individual* indiv);
        std::vector<Individual*> get_all_individuals();
        std::vector<double> get_env_factors();
        double get_exp_num_off_this_pop();


        void selection();
        double beverton_holt_prob(int n, double k_prime);
        std::vector<std::vector<Individual*>> split_by_sex();

        void add_to_migrant_queue(Individual* indiv);
        void add_migrants_to_population();
        std::vector<Individual*> pick_n_random_indivs(int n);
        Individual* get_random_individual();

        void add_to_next_gen(Individual* indiv);
        void replace_current_gen();

};

#endif
