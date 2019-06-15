
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "include.h"

class Individual{
    private:
        static int id_counter;
        int id;
        Population* patch_born_in;
        Population* patch;
        double* haplotype0;
        double* haplotype1;
        double w;
        MPFM* mpfm;

        int sex;
    public:
        bool has_migrated;
        double parent_was_migrant;
        double exp_num_off;

        Individual(MPFM* run, Population* patch, bool parent_was_migrant);
        ~Individual();
        Population* get_patch();
        Population* get_patch_born_in();
        int get_id();
        double get_fitness();
        void set_locus(int locus, int haplotype, double val);
        double get_locus(int locus, int haplotype);

        void migrate(Population* patch);

        void set_exp_num_off(double val);
        double get_exp_num_off();

        void calc_fitness();
        double calc_pref(Population* patch_i);
        void gen_haplotype(Individual* parent, int offspring_haplo);
        std::vector<double> get_crossing_over_points();

        int get_sex();

};

#endif
