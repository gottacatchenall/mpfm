
#ifndef MPFM_H
#define MPFM_H

#include "include.h"

class MPFM{
    private:
        // Initialization
        void init_random_generators();
        void init_env_factors();
        void init_populations();
        std::vector<std::vector<double>> init_kernel(double decay);
        void init_dispersal_kernel();
        void init_mating_kernel();
        void init_genome_dict();
        void init_individuals();
        void init_genomes();

        std::vector<int> gen_random_perm_of_unique_ints(int max_val, int length);
        void update_progress_bar(int gen);

        // Core
        void burn_in();
        void fragmentation();

        void run_generation(int gen);

        void dispersal(int gen);
        void selection();
        void reproduction();


        // Genome Dict
        std::vector<std::vector<int>> loci_per_ef;
        std::vector<double> selection_strengths;
        std::vector<int> chromo_map;

        // Kernels
        std::vector<std::vector<double>> dispersal_kernel;
        std::vector<std::vector<double>> mating_kernel;

        DataWrangler* data_wrangler;
    public:
        MPFM(std::unordered_map<std::string, float> params);

        // External
        void init();
        void start();

        // Parameters
        int N_POPULATIONS;
        int N_INDIVIDUALS;
        double BASE_MIGRATION_RATE;
        double DISPERSAL_DECAY_BURN_IN;
        double DISPERSAL_DECAY_FRAGMENTATION;
        double MATING_DECAY;
        int BURN_IN_LENGTH;
        int FRAGMENTATION_LENGTH;
        double INIT_NUM_ALLELES_MEAN;
        double PROP_FRAGMENTED;
        double MIN_RESISTANCE_AMOUNT;
        double MAX_RESISTANCE_AMOUNT;
        double EF_H_VALUE;
        std::string EF_STOCHASTIC_TYPE;
        int EF_NUMBER;
        int EF_GRANULARITY;
        double MUTATION_RATE;
        double GENOME_LENGTH;
        int N_NEUTRAL_LOCI;
        int N_FITNESS_LOCI_PER_EF;
        int N_LOCI;
        int N_CHROMOSOMES;
        double LOCUS_WEIGHT_MIN;
        double LOCUS_WEIGHT_MAX;
        double AVG_NUM_OFFSPRING_PER_INDIV;
        double SELECTION_STRENGTH;
        double SIGMA_S;
        int LD_LOG_FREQ;
        int CENSUS_FREQ;
        int FRAG_LOG_FREQ;
        int RS_EF;
        int RS_PATCHES;
        int RS_GENOME;
        int RS_MAIN;

        // Other (not canonical) params
        int CUSTOM_POPULATION_LOCATIONS; 
        int HARD_SELECTION;
        bool EF_STOCHASTIC;
        int REPRODUCTION_MODE;

        int FILL_UP_TO_K;

        // Random generators
        std::mt19937* main_gen;
        std::mt19937* ef_gen;
        std::mt19937* genome_gen;
        std::mt19937* patch_gen;

        // Getters
        int get_fitness_loci(int ef, int num);
        double get_selection_strength(int locus);
        int get_first_locus_of_chromo(int chromo_num);
        int get_total_pop_size();

        // Distributions
        int uniform_int(std::mt19937* gen, int a, int b);
        double uniform_real(std::mt19937* gen, double a, double b);
        int poisson(std::mt19937* gen, double rate);
        double normal(std::mt19937* gen, double mu, double sigma);
        int binomial(std::mt19937* gen, double p, double n);

        // EFs
        std::vector<EnvFactor*> env_factors;
        std::vector<Population*> populations;

};

#endif
