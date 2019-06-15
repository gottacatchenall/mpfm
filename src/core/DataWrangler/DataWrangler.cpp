
#include "DataWrangler.h"
#include "Population.h"
#include "Individual.h"

DataWrangler::DataWrangler(MPFM* run){

    this->mpfm = run;

    for (int l = 0; l < this->mpfm->N_LOCI; l++){
        this->allele_map.push_back(std::vector<allele*>());
    }

    int n_pops = this->mpfm->N_POPULATIONS;
    for (int p1 = 0; p1 < n_pops; p1++){
        this->attempted_migration_matrix.push_back(std::vector<int>());
        for (int p2 = 0; p2 < n_pops; p2++){
            this->attempted_migration_matrix[p1].push_back(0);
        }
    }

    std::string populations_header = "generation,population,x,y,w_mean,prop_of_k,effective_migration,";
    std::string demography_header = "generation,population,prop_loci_fixed,mean_alleles_per_locus\n";
    std::string dispersal_header = "generation,populationA,populationB,ct,prop_of_popB_from_popA\n";
    std::string pairwise_ld_header = "generation,population1,population2,ld,f_st\n";
    std::string genome_header = "locus,weight,chromosome\n";

    std::string pop_dist_header = "pop1,pop2,dist\n";
    std::string between_loci_ld_header = "generation,loci1,loci2,ld\n";

    int n_ef = this->mpfm->EF_NUMBER;

    for (int i = 0; i < n_ef; i++){
        populations_header +=  "ef" + std::to_string(i) + ",";
    }
    populations_header += "\n";

    this->pop_dist_file = this->init_file("pop_distance.csv", pop_dist_header);
    this->between_loci_ld_file = this->init_file("between_loci_ld.csv", between_loci_ld_header);
    this->demography_file = this->init_file("demography.csv", demography_header);
    this->dispersal_file = this->init_file("dispersal.csv", dispersal_header);
    this->pairwise_ld_file = this->init_file("pairwise_ld.csv", pairwise_ld_header);
    this->genome_file = this->init_file("genomes.csv", genome_header);
    this->populations_file = this->init_file("populations.csv", populations_header);
}

void DataWrangler::clear_dispersal_data(){
    int n_pops = this->mpfm->N_POPULATIONS;

    for (int p1 = 0; p1 < n_pops; p1++){
        for (int p2 = 0; p2 < n_pops; p2++){
            this->attempted_migration_matrix[p1][p2] = 0;
        }
    }
}

void DataWrangler::note_attempted_migration(int from, int to, int count){
    this->attempted_migration_matrix[from][to] += count;
}

void DataWrangler::census(int gen){
    std::vector<Population*> pops = this->mpfm->populations;
    int n_pops = this->mpfm->N_POPULATIONS;

    std::vector<Individual*> indivs;

    // Populations
    for (int p = 0; p < n_pops; p++){
        indivs = (this->mpfm->populations)[p]->get_all_individuals();
        double K = (this->mpfm->populations)[p]->get_K();
        int eff_migrants = 0;
        int ct = indivs.size();
        double fitness_sum = 0;

        for (Individual* indiv: indivs){
            double w =  indiv->get_fitness();
            //printf("w: %f\n", w);
            fitness_sum += w;

            if (indiv->parent_was_migrant){
                eff_migrants++;
            }
        }

        double eff_mig = double(eff_migrants)/double(ct);
        double mean_w = double(fitness_sum)/double(ct);
        double prop_of_k = double(ct)/double(K);
        double x = (this->mpfm->populations)[p]->get_x();
        double y = (this->mpfm->populations)[p]->get_y();
        std::vector<double> efs = (this->mpfm->populations)[p]->get_env_factors();

        this->log_populations(gen, p, x, y, mean_w, prop_of_k, eff_mig, efs);

    }

    // Dispersal
    for (int p1 = 0; p1 < n_pops; p1++){
        for (int p2 = 0; p2 < n_pops; p2++){
            if (p1 != p2){
                int p2_size = (this->mpfm->populations)[p2]->get_size();
                double att_mig = double(this->attempted_migration_matrix[p1][p2])/double(p2_size);

                if (att_mig > 0){
                    log_dispersal(gen, p1, p2, this->attempted_migration_matrix[p1][p2], att_mig);
                    // dispersal logging!
                    //printf("gen: %d, p1: %d, p2: %d, ct: %d, att_mig: %f\n", gen, p1, p2, this->attempted_migration_matrix[p1][p2], att_mig);
                }
            }
        }
    }

    this->populations_file->flush();
    this->dispersal_file->flush();
}

std::ofstream* DataWrangler::init_file(std::string file_path, std::string header){
    std::ofstream* file = new std::ofstream;
    file->open(file_path.c_str(), std::ios::app);
    (*file) << header;
    file->flush();
    return file;
}

void DataWrangler::log_dispersal(int gen, int pop_from, int pop_to, int ct, double att_mig){
    (*this->dispersal_file) << gen << "," << pop_from << "," << pop_to << "," << ct << "," << att_mig << "\n";
}

void DataWrangler::log_populations(int gen, int pop, double x, double y, double w_mean, double prop_of_k, double eff_mig, std::vector<double> efs){
    (*this->populations_file) << gen << "," << pop << "," << x << "," << y << "," << w_mean << "," << prop_of_k << "," << eff_mig << ",";

    for (double val : efs){
        (*this->populations_file) << val << ",";
    }
    (*this->populations_file) << "\n";
}

void DataWrangler::log_demography(int gen, int pop, double mean_num_loci_fixed, double mean_alleles_per_locus){
    (*this->demography_file) << gen << "," << pop << "," << mean_num_loci_fixed << "," << mean_alleles_per_locus << "\n";
}

void DataWrangler::log_pairwise_population_data(int gen, int pop1, int pop2, double d, double f_st){
    (*this->pairwise_ld_file) <<  gen << "," << pop1 << "," << pop2 << "," << d << "," << f_st << "\n";
}

void DataWrangler::log_ld_between_loci(int gen, int l1, int l2, double d){
    (*this->between_loci_ld_file) <<  gen << "," << l1 << "," << l2 << "," << d << "\n";
}

void DataWrangler::log_genome(int locus, double locus_weight, int chromo){
    (*this->genome_file) << locus << "," << locus_weight << "," << chromo << "\n";
    this->genome_file->flush();
}

void DataWrangler::log_pop_dist(int pop1, int pop2, double dist){
    (*this->pop_dist_file) << pop1 << "," << pop2 << "," << dist << "\n";
    this->pop_dist_file->flush();
}



void DataWrangler::get_ld(int gen){
    int n_loci = this->mpfm->N_LOCI;
    this->construct_allele_table();



    // Logging Pairwise Population Things
    int n_pops = this->mpfm->N_POPULATIONS;
    for (int p1 = 0; p1 < n_pops; p1++){
        for (int p2 = p1+1; p2 < n_pops; p2++){
            double ld = this->get_population_pairwise_ld(gen, p1, p2);
            double f_st = this->get_population_pairwise_fst(gen, p1, p2);
            this->log_pairwise_population_data(gen, p1, p2, ld, f_st);
        }
    }


    std::vector<allele*> alleles_l1;

    // Get freqencies and polymorphism ct
    for (int p = 0; p < n_pops; p++){
        int num_loci_fixed = 0;
        int polymorphism_ct = 0;
        for (int l1 = 0; l1 < n_loci; l1++){
            int polymorphism_ct_this_locus = 0;
            alleles_l1 = this->allele_map[l1];
            for (allele* al : alleles_l1){
                if (al->freq_map[p] > 0){
                    polymorphism_ct_this_locus++;
                }
            }

            polymorphism_ct += polymorphism_ct_this_locus;

            if (polymorphism_ct_this_locus == 1){
                num_loci_fixed++;
            }
        }

        double avg_num_polymorphisms = double(polymorphism_ct)/double(n_loci);

        this->log_demography(gen, p, num_loci_fixed, avg_num_polymorphisms);
    }

    // ==============================================
    // Get linkage network
    // ==============================================
    for (int l1 = 0; l1 < n_loci; l1++){
        for (int l2 = l1+1; l2 < n_loci; l2++){
            this->get_ld_between_loci(gen, l1, l2);
        }
    }

    // ==============================================
    // Free up memory
    // ==============================================
    for (int l = 0; l < n_loci; l++){
        for (allele* al : (this->allele_map[l])){
            std::vector<std::vector<dependent_allele*>> dep_als = al->loci;
            for (int l2 = 0; l2 < n_loci; l2++){
                for (dependent_allele* d_al : dep_als[l2]){
                    delete d_al;
                }
            }
            delete al;
        }
        this->allele_map[l].clear();
    }

    // Flush buffers
    this->demography_file->flush();
    this->pairwise_ld_file->flush();

}

void DataWrangler::get_ld_between_loci(int gen, int l1, int l2){
    double sum = 0;
    int ct = 0;

    std::vector<allele*> alleles_l1;
    std::vector<dependent_allele*> alleles_l2;

    alleles_l1 = this->allele_map[l1];

    // get total pop size * 2
    double eff_pop_size = 2*this->mpfm->get_total_pop_size();

    for (allele* al1: alleles_l1){
        double f_al1 = double(al1->n_total)/double(eff_pop_size);
        alleles_l2 = al1->loci[l2];
        for (dependent_allele* al2 : alleles_l2){
            double f_both = double(al2->n_total)/double(eff_pop_size);

            allele* al2_s = this->find_allele(al2->locus, al2->allele_val);
            double f_al2 = double(al2_s->n_total)/double(eff_pop_size);

            double ld = this->calc_ld(f_both, f_al1, f_al2);
            //printf("ld: %f\n", ld);
            ct++;
            sum += ld;
        }
    }


    double this_pair_avg_ld = double(sum)/double(ct);
//    printf("l1,l2: %d %d, sum: %f, ct: %d, this_pair_avg_ld: %f\n", l1, l2, sum, ct, this_pair_avg_ld);

    this->log_ld_between_loci(gen, l1, l2, this_pair_avg_ld);
}

double DataWrangler::get_population_pairwise_ld(int gen, int pop1, int pop2){
    int n_loci = this->mpfm->N_LOCI;
    Population* p1 = this->mpfm->populations[pop1];
    Population* p2 = this->mpfm->populations[pop2];

    int n1 = p1->get_size();
    int n2 = p2->get_size();

    int eff_pop_size = 2*(n1+n2);

    std::vector<allele*> alleles_l1;
    std::vector<dependent_allele*> alleles_l2;


    double total_sum = 0;
    int total_ct = 0;

    for (int l1 = 0; l1 < n_loci; l1++){
        alleles_l1 = this->allele_map[l1];
        for (int l2 = l1+1; l2 < n_loci; l2++){
            double sum = 0.0;
            int ct = 0;
            for (allele* al1: alleles_l1){
                double f_al1 = double(al1->freq_map[pop1] + al1->freq_map[pop2])/double(eff_pop_size);

                // this allele exists in either pop1 or pop2
                // if al1->freq_map is 0, the length of al2 at in that patch should also be 0;
                //if (al1->freq_map[pop1] > 0 || al1->freq_map[pop2] > 0){
                alleles_l2 = al1->loci[l2];
                // THIS ONLY CONSIDERS ALLELES AT AL1 that have been seen with AL2 somewhere globally!
                // if there exists an allele at l2 that is present in the pop, but never seen with l1, what happens?
                    // if allele exists at l2, it must be seen with something at l1, so it will be covered
                for (dependent_allele* al2 : alleles_l2){
                    double f_both = double(al2->freq_map[pop1] + al2->freq_map[pop2])/double(eff_pop_size);

                    allele* al2_s = this->find_allele(al2->locus, al2->allele_val);
                    double f_al2 = double(al2_s->freq_map[pop1] + al2_s->freq_map[pop2])/double(eff_pop_size);

                    double ld = this->calc_ld(f_both, f_al1, f_al2);

                    ct++;
                    sum += ld;

                    //double this_pair_avg_ld = double(sum)/double(ct);
                    total_ct++;
                    total_sum += ld;
                }

                //}
            }

            /*double this_pair_avg_ld = double(sum)/double(ct);

            total_ct++;
            total_sum += this_pair_avg_ld;*/
        }
    }

    // difference between heterozygosity and purely frequency based fst

    double avg_ld = double(total_sum)/double(total_ct);
    // print the dist
    // double dist = sqrt(pow(p1->get_x() - p2->get_x(), 2) + pow(p1->get_y() - p2->get_y(), 2));
    // printf("dist: %f\tld: %f \n", dist, avg_ld );
    return avg_ld;
}

double DataWrangler::get_population_pairwise_fst(int gen, int pop1, int pop2){
    int n_loci = this->mpfm->N_LOCI;
    Population* p1 = this->mpfm->populations[pop1];
    Population* p2 = this->mpfm->populations[pop2];

    int n1 = p1->get_size();
    int n2 = p2->get_size();
    
    std::vector<allele*> alleles_this_locus;
    double f_st_sum = 0;


    for (int l = 0; l < n_loci; l++){
        alleles_this_locus = this->allele_map[l];

        double J_1 = 0.0;
        double J_2 = 0.0;
        double J_12 = 0.0;
        double J_T = 0.0;

        // Let f_xy be the frequency of the yth allele in xth pop
        for (allele* al_i : alleles_this_locus){

            double f_1i = al_i->freq_map[pop1] / double(n1);
            double f_2i = al_i->freq_map[pop2] / double(n2);

            double w_1 = double(n1)/double(n1+n2);
            double w_2 = 1.0 - w_1;

            J_1 += (f_1i * f_1i);
            J_2 += (f_2i * f_2i);

            J_12 += (f_1i * f_2i);

            double sum_internal =  w_1*f_1i + w_2*f_2i;
            J_T += sum_internal * sum_internal;
        }

        double D_ST = 0.25 * (0.5*(J_1 + J_2) - J_12);
        double F_ST = D_ST / J_T;
        f_st_sum += F_ST;
    }

    double mean_fst = f_st_sum / double(n_loci);
    return mean_fst;
}


double DataWrangler::calc_ld(double p_ab, double p_a, double p_b){
    double D;

    D = p_ab - p_a*p_b;

    if (D == 0){
        return 0;
    }

    double denom = p_a*(1.0-p_a)*p_b*(1.0-p_b);
    double d2 = pow(D, 2);

    return d2/denom;
}

void DataWrangler::get_global_ld(){

}


void DataWrangler::construct_allele_table(){
    // construct allele table
    int n_loci = this->mpfm->N_LOCI;
    allele* al1_1;
    allele* al1_2;

    for (Population* pop: this->mpfm->populations){
        int patch_id = pop->get_id();
        std::vector<Individual*> indivs = pop->get_all_individuals();
        for (Individual* indiv: indivs){
            for (int l1 = 0; l1 < n_loci; l1++){
                this->update_tracker(l1, indiv->get_locus(l1, 0), patch_id);
                this->update_tracker(l1, indiv->get_locus(l1, 1), patch_id);
            }
            for (int l1 = 0; l1 < n_loci; l1++){
                al1_1 = find_allele(l1, indiv->get_locus(l1, 0));
                al1_2 = find_allele(l1, indiv->get_locus(l1, 1));
                for (int l2 = l1+1; l2 < n_loci; l2++){
                    this->add_dependent_allele(al1_1, l2, indiv->get_locus(l2, 0), patch_id);
                    this->add_dependent_allele(al1_2, l2, indiv->get_locus(l2, 1), patch_id);
                }
            }
        }
    }


    // do some checks to make sure this isn't stupid
    // (its still going to be stupid)

    int n_pop = this->mpfm->N_POPULATIONS;
    std::vector<allele*> alleles_l1;
    std::vector<dependent_allele*> alleles_l2;


    for (int p = 0; p < n_pop; p++){
        for (int l1 = 0; l1 < n_loci; l1++){
            alleles_l1 = this->allele_map[l1];
            for (int l2 = l1+1; l2 < n_loci; l2++){
                int al1_ct = 0;
                for (allele* al1 : alleles_l1){
                    al1_ct++;
                    alleles_l2 = al1->loci[l2];
                }
            }
        }
    }
}

void DataWrangler::update_tracker(int locus, double allele_val, int patch_id){
    int n_loci = this->mpfm->N_LOCI;
    allele* allele_struct = find_allele(locus, allele_val);
    if (allele_struct != NULL){
        note_allele_seen(allele_struct, patch_id);
    }
    else{
        allele* new_allele = new allele(locus, allele_val, n_loci);
        note_allele_seen(new_allele, patch_id);
        this->allele_map[locus].push_back(new_allele);
    }
}

void DataWrangler::add_dependent_allele(allele* primary_allele, int dependent_locus, double dependent_allele_val, int patch_id){
    assert(primary_allele != NULL);

    for (dependent_allele* al2 : primary_allele->loci[dependent_locus]){
        if (al2->allele_val == dependent_allele_val){
           note_allele_seen(al2, patch_id);
           return;
        }
    }


    dependent_allele* dep_al = new dependent_allele(dependent_locus, dependent_allele_val);
    primary_allele->loci[dependent_locus].push_back(dep_al);

    note_allele_seen(dep_al, patch_id);
}

void DataWrangler::note_allele_seen(dependent_allele* al, int patch_id){
    al->n_total++;
    al->freq_map[patch_id]++;
}

void DataWrangler::note_allele_seen(allele* al, int patch_id){
    al->n_total++;
    al->freq_map[patch_id]++;
}

allele* DataWrangler::find_allele(int locus, double allele_val){
    std::vector<allele*> allele_vec = this->allele_map[locus];
    for(allele* al: allele_vec) {
        if (al->allele_val == allele_val){
            return al;
        }
    }
    return NULL;
}
