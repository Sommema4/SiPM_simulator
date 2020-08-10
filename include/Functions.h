#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <iostream>

struct light
{
    int index;
    double time;
    double wavelenght; // saves wavelenght of light or index of microcell (look at sipm.simulate()
    std::string origin;
};

struct coordinates
{
    std::vector<double> x; // wavelenght
    std::vector<double> y; // QE value || emission probability || absorption probability
};

struct sipm_par
{
    double C_eq;
    double a1;
    double a2;
    double T_i;
    double T_d;
    double a_m_1;
    double a_m_2;
    double T_i2;
    double T_d2;
};

// prototypes
bool compare_light_time(const light &, const light &);
void add_pulse(int, double, std::string, std::vector<light>&);
double linear_interpolation(double, double, double, double, double);
double getClosest(double, double, double);
int get_nearest_index(std::vector<double>, double);
double get_rand_uni_dist(void);
int get_rand_poiss_dist(double);
void print_seed(void);
void sort_vector(std::vector<light>&);
void print_primaries(std::vector<light>, int);
double integrate(std::vector<double>);
double bisection_method(double, double, double, double, std::vector<double>, std::vector<double>, double (*f)(double, double, std::vector<double>, std::vector<double>));
double decay_function(double, double, std::vector<double>, std::vector<double>);
bool is_file_exist(const std::string&);
std::vector<double> linspace(double, double, int);

#endif // FUNCTIONS_H
