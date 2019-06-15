#include "Individual.h"
#include "Population.h"
//#include "GenomeDict.h"
//#include "MigrationTracker.h"

int Individual::id_counter = 0;

Individual::Individual(MPFM* run, Population* patch, bool parent_was_migrant){
    this->mpfm = run;
    int genome_size = this->mpfm->N_LOCI;

    this->id = this->id_counter++;
    this->patch = patch;
    this->patch_born_in = patch;
    this->has_migrated = false;
    this->haplotype0 = new double[genome_size];
    this->haplotype1 = new double[genome_size];
    this->parent_was_migrant = parent_was_migrant;

}

int Individual::get_sex(){
    return this->sex;
}

Individual::~Individual(){
    delete this->haplotype0;
    delete this->haplotype1;
}

int Individual::get_id(){
    return this->id;
}

Population* Individual::get_patch(){
    return this->patch;
}
Population* Individual::get_patch_born_in(){
    return this->patch_born_in;
}

void Individual::set_locus(int locus, int haplotype, double val){
    if (haplotype == 0){
        this->haplotype0[locus] = val;
        return;
    }
    else if (haplotype == 1){
        this->haplotype1[locus] = val;
        return;
    }
    assert(0 && "invalid haplo!");
}

double Individual::get_locus(int locus, int haplotype){
    double al = -1;
    if (haplotype == 0){
        al =  this->haplotype0[locus];
    }
    else if (haplotype == 1){
        al = this->haplotype1[locus];
    }
    /*
    if (al > 1 || al < 0){
        assert(0 && "alleles still not working...");
    }*/
    return al;
    assert(0 && "invalid haplotype");
}


// =========================================================================
// Migration
// =========================================================================

void Individual::migrate(Population* patch){
    this->patch = patch;
    //migration_tracker->note_attempted_migration(this->patch_born_in, patch);
    if (patch != this->patch_born_in){
        this->has_migrated = true;
    }
}


// =========================================================================
// Selection
// =========================================================================
void Individual::calc_fitness(){
    int n_ef = this->mpfm->EF_NUMBER;
    int n_fit_loci_per_ef = this->mpfm->N_FITNESS_LOCI_PER_EF;
    double sigma_s = this->mpfm->SIGMA_S;
    double sel_strength = this->mpfm->SELECTION_STRENGTH;

    int locus;
    double theta_i, x_i, s_i, w_i, diff, gaussian;

    double w = 1.0;
    double s_max_i = -1*sel_strength;

    double locus_weight;

    std::vector<double> theta = this->patch->get_env_factors();
    for (int i = 0; i < n_ef; i++){
        theta_i = theta[i];
        for (int j = 0; j < n_fit_loci_per_ef; j++ ){
            locus = this->mpfm->get_fitness_loci(i, j);
            locus_weight = this->mpfm->get_selection_strength(locus);
            x_i = this->get_locus(locus, 0);
            diff = x_i - theta_i;

            gaussian = exp( (-1*(pow(diff,2))/sigma_s ) );
            s_i = s_max_i * locus_weight * gaussian;
            w_i = 1.0 + s_i;
            w = w * w_i;

            x_i = this->get_locus(locus, 1);
            diff = x_i - theta_i;

            gaussian = exp( (-1*(pow(diff,2))/sigma_s ) );

            s_i = s_max_i * locus_weight * gaussian;
            w_i = 1.0 + s_i;
            w = w * w_i;
        }
    }

    this->w = w;
}

double Individual::get_fitness(){
    return this->w;
}

void Individual::set_exp_num_off(double val){
    this->exp_num_off = val;
}

double Individual::get_exp_num_off(){
    return this->exp_num_off;
}

void Individual::gen_haplotype(Individual* parent, int offspring_haplo){

    std::vector<double> crossing_over_points = get_crossing_over_points();

    int curr_crossover = 0;
    int curr_chromo = 0;
    int current_haplo = this->mpfm->uniform_int(this->mpfm->main_gen, 0,1);
    int n_loci = this->mpfm->N_LOCI;

    int n_chromo = this->mpfm->N_CHROMOSOMES;

   // pull n mutations from a binomial
    int n_mutations = this->mpfm->binomial(this->mpfm->main_gen, this->mpfm->MUTATION_RATE, n_loci);

    std::vector<int> mutation_sites;

    for (int i = 0; i < n_mutations; i++){
        int random_site = this->mpfm->uniform_int(this->mpfm->main_gen, 0, n_loci-1);
        mutation_sites.push_back(random_site);
    }

    int mutation_site_cter = 0;


    for (int locus = 0; locus < n_loci; locus++){
        // Check if at start of a new chromosome
        if (curr_chromo < n_chromo){
            if (locus == this->mpfm->get_first_locus_of_chromo(curr_chromo+1)){
                curr_chromo++;
                current_haplo = this->mpfm->uniform_int(this->mpfm->main_gen, 0,1);
            }
        }

        if (crossing_over_points.size() > 0){
            if (locus >= crossing_over_points[curr_crossover]){
                curr_crossover++;
                current_haplo = !(current_haplo);
            }
        }

        bool sets_something = false;
        bool mutation_this_site = false;

        if (n_mutations > 0){
            if (mutation_sites[mutation_site_cter] == locus){
                mutation_this_site = true;
            }
        }

        if (mutation_this_site){
            double new_al = this->mpfm->uniform_real(this->mpfm->main_gen, 0.0, 1.0);
            assert(new_al > 0 && new_al < 1);
            sets_something = true;
            this->set_locus(locus, offspring_haplo, new_al);
        }
        else{
            double new_al =  parent->get_locus(locus, current_haplo);
            assert(new_al >= 0 && new_al <= 1);
            sets_something = true;
            this->set_locus(locus, offspring_haplo, new_al);
        }

        assert(sets_something);
    }
}

std::vector<double> Individual::get_crossing_over_points(){
    double cM = this->mpfm->GENOME_LENGTH;

    double exp_num_co = cM / 100.0;
    int n_crossover_events = this->mpfm->poisson(this->mpfm->main_gen, exp_num_co);

    std::vector<double> crossing_over_points;

    int n_loci = this->mpfm->N_LOCI;

    // Generate a random permutation of unique loci to be the crossover points
    int r;
    bool exists;
    for(int i = 0; i < n_crossover_events; i++){
        r = this->mpfm->uniform_real(this->mpfm->main_gen, 0,n_loci-1);
        exists = (std::find(crossing_over_points.begin(), crossing_over_points.end(), r) != crossing_over_points.end());

        if (!exists){
            crossing_over_points.push_back(r);
        }
    }

    // Sort crossing_over_points
/*
    printf("crossovers: ");

    for (int co : crossing_over_points){
        printf(" %d", co);
    }
    printf("\n");*/




    std::sort(crossing_over_points.begin(), crossing_over_points.end());
    return crossing_over_points;
}
