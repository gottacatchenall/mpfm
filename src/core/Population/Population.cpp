#include "include.h"
#include "Population.h"
#include "Individual.h"
#include "EnvFactor.h"

int Population::id_counter = 0;

Population::Population(MPFM* run, double x, double y, double K){
    this->x = x;
    this->y = y;
    this->K = K;
    this->id = id_counter++;
    this->mpfm = run;
    this->individuals = new std::vector<Individual*>;
    this->next_gen = new std::vector<Individual*>;
}

double Population::get_x(){
    return this->x;
}

double Population::get_y(){
    return this->y;
}

double Population::get_K(){
    return this->K;
}

double Population::get_exp_num_off_this_pop(){
    return this->exp_num_off_this_population;
}

void Population::add_individual(Individual* indiv){
    this->individuals->push_back(indiv);
}

void Population::remove_individual(Individual* indiv){
    int index = std::distance(this->individuals->begin(), std::find(this->individuals->begin(), this->individuals->end(), indiv));
    this->individuals->erase(this->individuals->begin()+index);
}

int Population::get_size(){
    return this->individuals->size();
}

int Population::get_id(){
    return this->id;
}

std::vector<Individual*> Population::get_all_individuals(){
   return (*this->individuals);
}

std::vector<double> Population::get_env_factors(){
    double x = this->x;
    double y = this->y;

    std::vector<double> ef_vector;

    for (EnvFactor* ef_i : (mpfm->env_factors)){
        ef_vector.push_back(ef_i->get_env_factor_value(x,y));
    }

    return ef_vector;
}

void Population::selection(){
    std::vector<Individual*> indivs = this->get_all_individuals();
    bool surv;

    int K = this->K;
    int n = indivs.size();

    int b = mpfm->AVG_NUM_OFFSPRING_PER_INDIV;

    // abs fitness could be very important

    // just do beverton holt w/ abs fitness
    // calc realized K based on absolute fitness of pop
    // all parents equally likely after BH

    // try just BH at first, then see what happens when you weight K by mean abs fitness across pop

    double this_pop_abs_fitness_sum = 0;
    int this_pop_ct = 0;

    for (Individual* indiv : indivs){
        double abs_fitness = indiv->get_fitness();
        double k_prime = double(K * abs_fitness);
        double prob = this->beverton_holt_prob(n, k_prime);

        surv = false;
        double num = mpfm->uniform_real(mpfm->main_gen, 0.0, 1.0);
        if (num < prob){
            surv =  true;
        }

        if (!surv){
            remove_individual(indiv);
            delete indiv;
        }
        else {
            this_pop_abs_fitness_sum += abs_fitness;
            this_pop_ct++;
            indiv->set_exp_num_off(b);
        }
    }

    double this_pop_mean_abs_fitness = double(this_pop_abs_fitness_sum)/double(this_pop_ct);

    if (this->mpfm->FILL_UP_TO_K){
        this->exp_num_off_this_population = K;
    }
    else{
        this->exp_num_off_this_population = this_pop_mean_abs_fitness*K;
    }
}

double Population::beverton_holt_prob(int n, double k_prime){
    int b = mpfm->AVG_NUM_OFFSPRING_PER_INDIV;

    double prop_full = double(n)/double(k_prime);
    double prob = double(1.0) / double(1 + (double(b)/double(2) - 1)*(prop_full));
    return prob;
}

void Population::add_to_migrant_queue(Individual* indiv){
    this->migrant_queue.push_back(indiv);
}

void Population::add_migrants_to_population(){
    for (Individual* indiv: this->migrant_queue){
        this->add_individual(indiv);
    }
    this->migrant_queue.clear();
}

Individual* Population::get_random_individual(){
    int size = this->individuals->size();
    int index = this->mpfm->uniform_int(mpfm->main_gen, 0,size-1);

    Individual* tmp = (*this->individuals)[index];
    return tmp;
}

std::vector<Individual*> Population::pick_n_random_indivs(int n){
    int index;
    int size;
    Individual* indiv;
    std::vector<Individual*> indivs;

    int ct = 0;
    while (ct < n){
        //printf("ct: [%d, %d]\n", ct, n);
        size = this->individuals->size();
        /*printf("size: %d, n= :%d\n", size, n);
        for (Individual* i : (*this->individuals)){
            if(!i->has_migrated){
                mig_ct++;
            }
        }

        if ((size-mig_ct) >= ct){
            printf("All have migr!\n");
        }*/
        index = mpfm->uniform_int(mpfm->main_gen, 0,size-1);
        indiv = (*this->individuals)[index];
        if (!indiv->has_migrated){
            indiv->has_migrated = true;
            indivs.push_back(indiv);
            ct++;
            this->remove_individual(indiv);
        }
    }

    return indivs;
}


void Population::add_to_next_gen(Individual* indiv){
    this->next_gen->push_back(indiv);
}

void Population::replace_current_gen(){
    for (Individual* indiv: *(this->individuals)){
        delete indiv;
    }
    delete this->individuals;
    this->individuals = this->next_gen;
    this->next_gen = new std::vector<Individual*>;
}


std::vector<std::vector<Individual*>> Population::split_by_sex(){
    std::vector<Individual*> males;
    std::vector<Individual*> females;

    for (Individual* indiv: this->get_all_individuals()){
        if (indiv->get_sex()== 1){
           males.push_back(indiv);
        }
        else if (indiv->get_sex() == 0){
           females.push_back(indiv);
        }
    }
    std::vector<std::vector<Individual*>> res;
    res.push_back(females);
    res.push_back(males);
    return res;
}
