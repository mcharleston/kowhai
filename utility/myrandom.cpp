/**
 * Bunch of random number generation things
 * MAC, Dec 21, 2022
 */

#include <chrono>
#include <cmath>
#include <random>
#include <string>

#include "myrandom.h"

using namespace std;
//using namespace kowhai;

bool _silent(false);

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

std::random_device rd;
std::mt19937 generator{rd()}; // or std::default_random_engine e{rd()};
//std::default_random_engine generator(seed);
std::uniform_real_distribution<float> runif(0.0, 1.0);
std::uniform_real_distribution<double> dunif(0.0, 1.0);

float fran() {
	return dunif(generator);
}

double dran(double mult) {
	return mult * dunif(generator);
}

int iran(int max) {
	return floor(max * 1.0 * fran());
}
uint plran(float l, float u, float r) {
	float y(fran());
	float ex(r+1.0);
	return static_cast<unsigned int>( std::pow( std::pow(u,ex) - std::pow(l, ex)*y + std::pow(l, ex), 1.0/ex ) );
}

