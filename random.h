//
// Created by eurek on 6/5/2019.
//

#ifndef COMPASS_V2_RANDOM_H
#define COMPASS_V2_RANDOM_H


#include <random>

class Random {
public:
    Random();
    static Random *get_instance();
    double random_number();
private:
    static Random* instance;
    std::mt19937_64 mt;
    std::uniform_real_distribution<> dist;
};


#endif //COMPASS_V2_RANDOM_H
