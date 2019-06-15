#include "MPFM.h"
#include "EnvFactor.h"
#include "Population.h"
#include "Individual.h"
#include "DataWrangler.h"

MPFM::MPFM(std::unordered_map<std::string, float> params){
    this->N_POPULATIONS = params["N_POPULATIONS"];
    this->N_INDIVIDUALS = params["N_INDIVIDUALS"];
    this->BASE_MIGRATION_RATE = params["BASE_MIGRATION_RATE"];
    this->DISPERSAL_DECAY_BURN_IN = params["DISPERSAL_DECAY_BURN_IN"];
    this->DISPERSAL_DECAY_FRAGMENTATION = params["DISPERSAL_DECAY_FRAGMENTATION"];
    this->MATING_DECAY = params["MATING_DECAY"];
    this->BURN_IN_LENGTH = params["BURN_IN_LENGTH"];
    this->FRAGMENTATION_LENGTH = params["FRAGMENTATION_LENGTH"];
    this->INIT_NUM_ALLELES_MEAN = params["INIT_NUM_ALLELES_MEAN"];
    this->PROP_FRAGMENTED = params["PROP_FRAGMENTED"];
    this->MIN_RESISTANCE_AMOUNT = params["MIN_RESISTANCE_AMOUNT"];
    this->MAX_RESISTANCE_AMOUNT = params["MAX_RESISTANCE_AMOUNT"];
    this->EF_H_VALUE = params["EF_H_VALUE"];
    this->EF_NUMBER = params["EF_NUMBER"];
    this->EF_GRANULARITY = params["EF_GRANULARITY"];
    this->MUTATION_RATE = params["MUTATION_RATE"];
    this->GENOME_LENGTH = params["GENOME_LENGTH"];
    this->N_NEUTRAL_LOCI = params["N_NEUTRAL_LOCI"];
    this->N_FITNESS_LOCI_PER_EF = params["N_FITNESS_LOCI_PER_EF"];
    this->N_CHROMOSOMES = params["N_CHROMOSOMES"];
    this->AVG_NUM_OFFSPRING_PER_INDIV = params["AVG_NUM_OFFSPRING_PER_INDIV"];
    this->SELECTION_STRENGTH = params["SELECTION_STRENGTH"];
    this->SIGMA_S = params["GAUSSIAN_SELECTION_STRENGTH"];
    this->LD_LOG_FREQ = params["LD_LOG_FREQ"];
    this->CENSUS_FREQ = params["CENSUS_FREQ"];
    this->FRAG_LOG_FREQ = params["FRAG_LOG_FREQ"];
    this->RS_EF = params["RS_EF"];
    this->RS_PATCHES = params["RS_PATCHES"];
    this->RS_GENOME = params["RS_GENOME"];
    this->RS_MAIN = params["RS_MAIN"];
    this->LOCUS_WEIGHT_MIN = params["LOCUS_WEIGHT_MIN"];
    this->LOCUS_WEIGHT_MAX = params["LOCUS_WEIGHT_MAX"];

    this->N_LOCI = this->N_NEUTRAL_LOCI + this->EF_NUMBER* this->N_FITNESS_LOCI_PER_EF;


    this->EF_STOCHASTIC = params["EF_STOCHASTIC"];
    this->HARD_SELECTION = params["HARD_SELECTION"];
    this->REPRODUCTION_MODE = params["REPRODUCTION_MODE"];
    this->FILL_UP_TO_K = params["FILL_UP_TO_K"];
    this->data_wrangler = new DataWrangler(this);

}

// ======================================================
// Core Functions
// ======================================================

void MPFM::init(){
    this->init_random_generators();
    this->init_env_factors();
    this->init_populations();
    this->init_genome_dict();
    this->init_mating_kernel();
    this->init_dispersal_kernel();
    this->init_individuals();
    this->init_genomes();
}

void MPFM::start(){
    this->burn_in();
    this->fragmentation();
}


void MPFM::update_progress_bar(int gen){
    int n_gen = this->BURN_IN_LENGTH;

    double progress = double(gen)/double(n_gen);
    int barWidth = 50;
       std::cout << "\t[";
       int pos = barWidth * progress;
       for (int i = 0; i < barWidth; ++i) {
           if (i <= pos) std::cout << "=";
           else std::cout << " ";
       }
       std::cout << "] [" << gen << " / " << n_gen << "] \r";
       std::cout.flush();
}

void MPFM::burn_in(){
    int burn_in_length = this->BURN_IN_LENGTH;

    for (int gen = 0; gen <= burn_in_length; gen++){
        this->run_generation(gen);
    }
}

void MPFM::fragmentation(){

    int frag_length = this->FRAGMENTATION_LENGTH;

    // Change the dispersal kernel to the "fragmented" state using a new dispersal decay value
    this->dispersal_kernel = this->init_kernel(this->DISPERSAL_DECAY_FRAGMENTATION);

    // Change Log Params
    this->LD_LOG_FREQ = this->FRAG_LOG_FREQ;
    this->CENSUS_FREQ = this->FRAG_LOG_FREQ;

    // Run Generations
    int burn_in = this->BURN_IN_LENGTH;

    for (int gen = burn_in; gen <= burn_in + frag_length; gen++){
        this->run_generation(gen);
    }
}

void MPFM::run_generation(int gen){

    this->dispersal(gen);
    this->selection();
    this->reproduction();

    // Logging
    if (gen % this->CENSUS_FREQ == 0){
        this->data_wrangler->census(gen);
        this->data_wrangler->clear_dispersal_data();
    }
    if (gen % this->LD_LOG_FREQ == 0){
        this->data_wrangler->get_ld(gen);
    }

    if (gen % 10 == 0){
        update_progress_bar(gen);
    }
}

void MPFM::dispersal(int gen){
    double base_disp = this->BASE_MIGRATION_RATE;
    int n_pops = this->N_POPULATIONS;
    std::vector<Individual*> indivs;
    Population* pop1;
    Population* pop2;
    double p1_to_p2_before_rounding, p1_to_p2_kernel, remainder;
    int num_from_p1_to_p2, p1_to_p2_after_rounding, additional_migrant, pop1_size;

    for (int p1 = 0; p1 < n_pops; p1++){
        pop1 = this->populations[p1];
        pop1_size = pop1->get_size();
        for (int p2 = 0; p2 < n_pops; p2++){
            // Calc dist
            pop1 = this->populations[p1];
            if (p1 != p2){
                pop2 = this->populations[p2];
                p1_to_p2_kernel = this->dispersal_kernel[p1][p2];
                p1_to_p2_before_rounding = (base_disp * p1_to_p2_kernel * double(pop1_size));


                p1_to_p2_after_rounding = int(floor(p1_to_p2_before_rounding));
                remainder = p1_to_p2_before_rounding - p1_to_p2_after_rounding;

                if (this->uniform_real(this->main_gen, 0.0, 1.0) < remainder){
                    additional_migrant = 1;
                }
                else{
                    additional_migrant = 0;
                }

                num_from_p1_to_p2 = p1_to_p2_after_rounding + additional_migrant;
                num_from_p1_to_p2 = p1_to_p2_after_rounding;
                //printf("dist: %f, p1_to_p2_kernel: %f, base_disp: %f, pop_size: %d, before_rounding: %f, num: %d\n", dist, p1_to_p2_kernel, base_disp, pop1_size, p1_to_p2_before_rounding, num_from_p1_to_p2);

                if (gen % this->CENSUS_FREQ == 0){
                    this->data_wrangler->note_attempted_migration(p1, p2, num_from_p1_to_p2);
                }

                // REMEMBER pick_n_random_indivs removes them from the patch
                indivs = pop1->pick_n_random_indivs(num_from_p1_to_p2);
                for (Individual* indiv_i: indivs){
                    pop2->add_to_migrant_queue(indiv_i);
                    indiv_i->migrate(pop2);
                    //pop1->remove_individual(indiv_i);
                }
            }
        }
    }

    for (Population* pop: populations){
        pop->add_migrants_to_population();
    }

}

void MPFM::selection(){
    for (Population* pop: populations){
        std::vector<Individual*> indivs = pop->get_all_individuals();
        for (Individual* indiv: indivs){
            indiv->calc_fitness();
        }
    }

    for (Population* pop: populations){
        pop->selection();
    }
}

void MPFM::reproduction(){
    bool parent_migrated;

    Population* pop1;

    Individual* random1;
    Individual* random2;
    Individual* offspring;
    std::vector<Individual*> indivs;

    int n_pops = this->N_POPULATIONS;

    for (int p1 = 0; p1 < n_pops; p1++){
        pop1 = this->populations[p1];
        int n = pop1->get_size();

        if (n > 1){
            indivs = pop1->get_all_individuals();

            double exp_num_off = pop1->get_exp_num_off_this_pop();
            int n_off = int(exp_num_off);

            for (int i = 0; i < n_off; i++ ){
                random1 = pop1->get_random_individual();
                random2 = random1;

                while (random2 == random1){
                    random2 = pop1->get_random_individual();
                }

                parent_migrated = (random1->has_migrated || random2->has_migrated);
                offspring = new Individual(this, pop1, parent_migrated);
                offspring->gen_haplotype(random1, 0);
                offspring->gen_haplotype(random2, 1);
                pop1->add_to_next_gen(offspring);
            }
        }

        pop1->replace_current_gen();
    }
}


// ======================================================
// Initialization
// ======================================================

void MPFM::init_random_generators(){
    this->main_gen = new std::mt19937(this->RS_MAIN);
    this->ef_gen =  new std::mt19937(this->RS_EF);
    this->genome_gen = new std::mt19937(this->RS_GENOME);
    this->patch_gen = new std::mt19937(this->RS_PATCHES);
}

void MPFM::init_env_factors(){
    EnvFactor* tmp;
    for (int ef = 0; ef < this->EF_NUMBER; ef++){
        tmp = new EnvFactor(this, this->EF_GRANULARITY, this->EF_H_VALUE);
        this->env_factors.push_back(tmp);
    }
}

void MPFM::init_populations(){
    int n_pops = this->N_POPULATIONS;
    int n_indiv = this->N_INDIVIDUALS;

    double k_mean = double(n_indiv)/double(n_pops);

    for (int i = 0; i < n_pops; i++){
        double k = k_mean;
        double x = this->uniform_real(this->patch_gen, 0.0, 1.0);
        double y = this->uniform_real(this->patch_gen, 0.0, 1.0);
        //printf("pop: %d, x: %f, y: %f\n", i, x, y);
        Population *tmp = new Population(this, x, y, k);
        this->populations.push_back(tmp);
    }

    Population* pop1;
    Population* pop2;
    double delta_x, delta_y, dist;

    for (int p1 = 0; p1 < n_pops; p1++){
        pop1 = this->populations[p1];
        for (int p2 = 0; p2 < n_pops; p2++){
            pop2 = this->populations[p2];
            delta_x = pop1->get_x() - pop2->get_x();
            delta_y = pop1->get_y() - pop2->get_y();

            dist = sqrt((delta_x*delta_x) + (delta_y*delta_y));
            this->data_wrangler->log_pop_dist(p1, p2, dist);
        }
    }
}

std::vector<std::vector<double>> MPFM::init_kernel(double decay){
    int n_pops = this->N_POPULATIONS;

    // special case where n pops is 2. migration between both demes is equal to the base migration rate
    if (n_pops == 2){
        std::vector<std::vector<double>> kern = {{1,1}, {1,1}};
        return kern;
    }

    Population* pop1;
    Population* pop2;

    double deltaX, deltaY, dist, kern_val;
    //double s = 0.0;
    std::vector<std::vector<double>> kernel;

    for (int p1 = 0; p1 < n_pops; p1++){
        kernel.push_back(std::vector<double>());
    //    s = 0.0;
        pop1 = this->populations[p1];
        for (int p2 = 0; p2 < n_pops; p2++){
                pop2 = this->populations[p2];

                deltaX = pop1->get_x() - pop2->get_x();
                deltaY = pop1->get_y() - pop2->get_y();

                dist = sqrt((deltaX*deltaX) + (deltaY*deltaY));

                if (p1 != p2){
                    kern_val = exp(-1*decay*dist);
                }
                else {
                    kern_val = 0;
                }

                kernel[p1].push_back(kern_val);
        }
    }

    return kernel;
}

void MPFM::init_dispersal_kernel(){
    this->dispersal_kernel = this->init_kernel(this->DISPERSAL_DECAY_BURN_IN);
}

void MPFM::init_mating_kernel(){
    this->mating_kernel = this->init_kernel(this->MATING_DECAY);
}

void MPFM::init_genome_dict(){
    int nlef = this->N_FITNESS_LOCI_PER_EF;
    int n_fitness_loci = nlef*this->EF_NUMBER;
    int n_loci = this->N_LOCI;
    int max_locus_index = n_loci-1;

    for (int l = 0; l < n_loci; l++){
        this->selection_strengths.push_back(0.0);
    }


    std::vector<int> random_perm = this->gen_random_perm_of_unique_ints(max_locus_index, n_fitness_loci);
    int ct = 0;
    for (int ef = 0; ef < this->EF_NUMBER; ef++){
        std::vector<int> this_ef;
        for (int l = 0; l < nlef; l++){
            int locus = random_perm[ct];
            this_ef.push_back(random_perm[ct]);
            double lw = this->uniform_real(this->genome_gen, this->LOCUS_WEIGHT_MIN, this->LOCUS_WEIGHT_MAX);
            this->selection_strengths[locus] = lw;
            ct++;
        }
        this->loci_per_ef.push_back(this_ef);
    }


    // selection_strengths
    // Chromosome Map
    int n_chromo = this->N_CHROMOSOMES;
    this->chromo_map = this->gen_random_perm_of_unique_ints(max_locus_index, n_chromo);

    std::sort(this->chromo_map.begin(), this->chromo_map.end());

    int chromo_ct = 0;

    for (int l = 0; l < n_loci; l++){
        if (l >= this->chromo_map[chromo_ct]){
            chromo_ct++;
        }
        this->data_wrangler->log_genome(l, this->selection_strengths[l], chromo_ct);
    }
}

void MPFM::init_individuals(){
    int n_pops = this->N_POPULATIONS;
    int k;
    Population* pop;
    Individual* tmp;

    for (int p = 0; p < n_pops; p++){
        pop = populations[p];
        k = int(pop->get_K());
        for (int indiv = 0; indiv < k; indiv++){
            tmp = new Individual(this, pop, false);
            pop->add_individual(tmp);
        }
    }
}

void MPFM::init_genomes(){
    int n_loci = this->N_LOCI;
    int n_alleles;
    double al_val;

    std::vector<std::vector<double>> alleles;
    for (int l = 0; l < n_loci; l++){
        n_alleles = this->poisson(this->genome_gen, this->INIT_NUM_ALLELES_MEAN);
        if (n_alleles == 0){
            n_alleles = 1;
        }
        alleles.push_back(std::vector<double>());
        for (int al = 0;  al < n_alleles; al++){
            al_val = this->uniform_real(this->genome_gen, 0.0, 1.0);
            alleles[l].push_back(al_val);
        }
    }

    std::vector<Individual*> indivs;
    int allele_ct, rand_allele_index;
    double val;

    for (Population* p : populations){
        indivs = p->get_all_individuals();
        for (Individual* indiv : indivs){
            for (int l = 0; l < n_loci; l++){
                allele_ct = alleles[l].size();
                rand_allele_index = this->uniform_int(this->genome_gen, 0, allele_ct-1);
                val = alleles[l][rand_allele_index];
                indiv->set_locus(l, 0, val);

                rand_allele_index = this->uniform_int(this->genome_gen, 0, allele_ct-1);
                val = alleles[l][rand_allele_index];
                indiv->set_locus(l, 1, val);
            }
        }
    }
}

std::vector<int> MPFM::gen_random_perm_of_unique_ints(int max_val, int length){
    std::vector<int> result;
    int r;
    int exists;
    for(int i = 0; i < length; i++){
        do{
            exists = 0;
            r = this->uniform_int(this->genome_gen, 0, max_val);
            for (int j = 0; i < j; j++){
                if (result[j] == r){
                    exists = 1;
                    break;
                }
            }
        } while(exists);
        result.push_back(r);
    }
    return result;
}

// ======================================================
// Getters
// ======================================================

int MPFM::get_fitness_loci(int ef, int num){
    return this->loci_per_ef[ef][num];
}

double MPFM::get_selection_strength(int locus){
    return this->selection_strengths[locus];
}

int MPFM::get_first_locus_of_chromo(int chromo_num){
    return this->chromo_map[chromo_num];
}

int MPFM::get_total_pop_size(){
    int sum = 0;
    for (Population* pop : this->populations){
        sum += pop->get_size();
    }
    return sum;
}

// ======================================================
// Random Generators
// ======================================================

int MPFM::uniform_int(std::mt19937* gen, int a, int b){
    std::uniform_int_distribution<int> dis(a,b);
    return dis(*gen);
}

double MPFM::uniform_real(std::mt19937* gen, double a, double b){
    std::uniform_real_distribution<double> dis(a,b);
    return dis(*gen);
}

double MPFM::normal(std::mt19937* gen, double mu, double sigma){
    std::normal_distribution<double> norm_dis(mu,sigma);
    double val = norm_dis(*gen);
    return val;
}

int MPFM::poisson(std::mt19937* gen, double rate){
    std::poisson_distribution<int> p_dis(rate);
    return p_dis(*gen);
}

int MPFM::binomial(std::mt19937* gen, double p, double n){
    std::binomial_distribution<int> b_dis(n, p);
    return b_dis(*gen);
}
