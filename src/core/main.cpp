#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <assert.h>
#include <cmath>
#include <map>
#include <set>
#include <unordered_map>
#include <unistd.h>

#include "MPFM.h"

std::unordered_map<std::string, float> read_params_file(){
    std::ifstream infile("params.ini");
    std::string line;

    std::unordered_map<std::string, float> params;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::vector<std::string> record;
        while (iss){
          std::string s;
          if (!getline(iss, s, ',' )) break;
          record.push_back( s );
        }

        std::string name = record[0];
        double val = atof(record[1].c_str());

        params.insert(std::pair<std::string, double>(name, val));
    }

    return params;
}

int main(int argc, char* argv[]){
    if (argv[1]){
        chdir(argv[1]);
    }

    std::unordered_map<std::string, float> params = read_params_file();

    MPFM* mpfm_run = new MPFM(params);
    mpfm_run->init();
    mpfm_run->start();


    return 0;
}
