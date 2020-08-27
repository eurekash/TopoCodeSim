#include "random.h"

Random* Random::instance = nullptr;

Random* Random::get_instance() {
    if (instance == nullptr) {
        instance = new Random();
    }

    return instance;
}

Random::Random(): mt(std::random_device()()), dist(0,1) { }

double Random::random_number() {
    return dist(mt);
}

