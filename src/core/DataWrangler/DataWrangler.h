
#ifndef DATA_WRANGLER_H
#define DATA_WRANGLER_H

#include "include.h"

typedef struct dependent_allele{
    dependent_allele(int l, double val){
        locus = l;
        allele_val = val;
        n_total = 0;
    };
    int locus;
    double allele_val;
    int n_total;
    std::unordered_map<int,int> freq_map;
} dependent_allele;

typedef struct allele{
    allele(int l, double val, int n_loci){
        locus = l;
        allele_val = val;
        n_total = 0;

        for (int i = 0; i < n_loci; i++){
            loci.push_back(std::vector<dependent_allele*>());
        }
    };


    int locus;
    double allele_val;
    int n_total;
    std::unordered_map<int,int> freq_map;
    std::vector<std::vector<dependent_allele*>> loci;
} allele;

class DataWrangler{
    private:
        std::ofstream* demography_file;
        std::ofstream* dispersal_file;
        std::ofstream* pairwise_ld_file;
        std::ofstream* genome_file;
        std::ofstream* populations_file;
        std::ofstream* between_loci_ld_file;

        std::ofstream* pop_dist_file;

        MPFM* mpfm;
        std::vector<int> pop_ct;
        std::vector<double> fitness;

        std::vector<std::vector<int>> attempted_migration_matrix;
        std::vector<std::vector<allele*>> allele_map;
    public:
        DataWrangler(MPFM* run);
        std::ofstream* init_file(std::string file_path, std::string header);

        void census(int gen);
        void get_ld(int gen);


        void clear_dispersal_data();
        void note_attempted_migration(int from, int to, int count);

        void log_populations(int gen, int pop, double x, double y, double w_mean, double prop_of_k, double eff_mig, std::vector<double> efs);
        void log_demography(int gen, int pop, double mean_num_loci_fixed, double mean_alleles_per_locus);
        void log_dispersal(int gen, int pop_from, int pop_to, int ct, double att_mig);
        void log_pairwise_population_data(int gen, int pop1, int pop2, double d, double f_st);
        void log_ld_between_loci(int gen, int l1, int l2, double d);
        void log_genome(int locus, double locus_weight, int chromo);
        void log_pop_dist(int pop1, int pop2, double dist);

        // LD nonsense s
        double get_population_pairwise_ld(int gen, int pop1, int pop2);
        double get_population_pairwise_fst(int gen, int pop1, int pop2);
        void get_global_ld();
        double calc_ld(double p_ab, double p_a, double p_b);
        void get_ld_between_loci(int gen, int l1, int l2);

        // table
        void construct_allele_table();
        void update_tracker(int locus, double allele_val, int patch_id);
        void add_dependent_allele(allele* primary_allele, int dependent_locus, double dependent_allele_val, int patch_id);
        void note_allele_seen(dependent_allele* al, int patch_id);
        void note_allele_seen(allele* al, int patch_id);
        allele* find_allele(int locus, double allele_val);

};

#endif
