#include "Functions.h"

#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <ctime>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include "json.hpp"

time_t current_time = time(NULL);
std::mt19937 generator(current_time);
std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);

void parse_settings_json(settings *set, std::string filename)
{
    std::ifstream i(filename);
    nlohmann::json init;
    i >> init;

    /* SIMULATION VARIABLES */
    set->sim_name = init["simulation"]["name"]; // the name of your simulation
    set->out_dir = init["simulation"]["parent directory"]; // the name of the directory where results will be saved - relative to the binary file
    set->sim_type = init["simulation"]["calculation type"]; // the type of the simulation
    set->out_aux = init["simulation"]["auxiliary output"]; // decides if auxiliary output should be written to file [boolean]
    set->overvoltage = init["simulation"]["overvoltage"];
    set->geometry_factor = init["simulation"]["geometry factor"];
    set->R_l = init["simulation"]["load resistor"]; // load resistor
    set->I_th = init["simulation"]["threshold current"]; // threshold level of current to sustain the avalanche = 0.0001 A

    /* GLOBAL VARIABLES */
    set->repetition = init["global variables"]["repetition"]; // the number of simuation repetitons for each deposited energy
    set->timestep = init["global variables"]["timestep"]; // the lenght of timestep in the simulation [seconds]
    set->pulse_lenght = init["global variables"]["pulse lenght"]; // the maximum lenght of the calculation [seconds]
    set->deposited_energy = init["global variables"]["deposited energy"].get<std::vector<double>>();; // the list of energies of the simulation [MeV]

    /* SCINTILLATOR VARIABLES */
    set->scint_name = init["scintillator"]["name"]; // the full name of scintillator
    set->scint_lib = init["scintillator"]["library"]; // the scintillator library name
    set->particle = init["scintillator"]["particle type"]; // the abbreviation of particle and its full name
    set->photon_list = init["scintillator"]["photon list"]; // decides if list of all generated photons should be written to auxiliary output [boolean]

    /* SiPM VARIABLES */
    set->sipm_name = init["sipm"]["name"]; // the full name of sipm
    set->sipm_lib = init["sipm"]["library"]; // the sipm library name
    set->crosstalk = init["sipm"]["crosstalk"]; // decides if crosstalk should be taken into account during calculation [boolean]
    set->afterpulse = init["sipm"]["afterpulse"]; // decides if afterpulse should be taken into account during calculation [boolean]
    set->dark_current = init["sipm"]["dark current"]; // decides if dark current should be taken into account during calculation [boolean]

    /* ALGORITHM PARAMETERS */
    set->interpolation = init["algorithm parameters"]["interpolation"]; // chooses the type of interpolation technique
    set->seed = init["algorithm parameters"]["seed"]; // the seed number, if seed=0 the number is chosen based on the time
}

void print_seed(void)
{
    std::cout << "seed: " << current_time << std::endl;
}

double get_rand_uni_dist(void)
{
    return uniform_distribution(generator);
}

int get_rand_poiss_dist(double mean)
{
    std::poisson_distribution<int> poisson_dist(mean);
    return poisson_dist(generator);
}

bool compare_light_time(const light &a, const light &b)
{
    return a.time < b.time;
}

double linear_interpolation(double x0, double x1, double y0, double y1, double target)
{
    double a = (y1 - y0) / (x1 - x0);
    double y = a * (target - x0) + y0;
    return y;
}

double getClosest(double val1, double val2, double target)
{
    if (target - val1 >= val2 - target)
        return val2;
    else
        return val1;
}

int get_nearest_index(std::vector<double> arr, double target)
{
    double difference = 1e+30;
    double temp;
    int output = 0;
    for (int i=0;i<arr.size();i++){
        temp = std::abs(arr[i] - target);
        //std::cout << temp << " " << arr[i] << " " << target << std::endl;
        if (temp <= difference){
            difference = temp;
            output = i;
        }
    }
    //std::cout << output << std::endl;
    return output;
}

void add_pulse(int index, double time, std::string origin, std::vector<light>& buffer)
{
    light temp;
    temp.index = index;
    temp.time = time;
    temp.origin = origin;
    buffer.push_back(temp);
    std::sort(buffer.begin(), buffer.end(), compare_light_time);
}

void sort_vector(std::vector<light>& buffer)
{
    std::sort(buffer.begin(), buffer.end(), compare_light_time);
}

void print_primaries(std::vector<light> buffer, int n=5)
{
    /* SHOWS EXAMPLE OF GENERATED LIGHT PRIMARIES */
    if (n == 0) // if n ==0 show nothing
        return;

    int size_vec = buffer.size();
    int temp = std::min(n, size_vec);
    std::cout << "----PRIMARY LIGHT EXAMPLE----" << std::endl;
    for (int i=0;i<temp;i++)
        std::cout << buffer[i].index << " " << buffer[i].time << " " << buffer[i].wavelenght << " " << buffer[i].origin << std::endl;
    std::cout << "-----------------------------" << std::endl;
}

double integrate(std::vector<double> a)
{
    double sum_of_elems = std::accumulate(a.begin(), a.end(), 0.0);
    return sum_of_elems;
}

double bisection_method(double a, double b, double y, double threshold, std::vector<double> decay_weight, std::vector<double> decay_component, double (*f)(double, double, std::vector<double>, std::vector<double>))
{
    double c;
    int kk = 0;
    //std::cout << y << std::endl;
    if(f(a, y, decay_weight, decay_component) * f(b, y, decay_weight, decay_component) >= 0.0)
        std::cout << "Incorrect a and b" << std::endl;

    while ((b - a) >= threshold){
        kk++;
        c = (a + b)/ 2.0;
        //std::cout << a << " " << b << " " << y << std::endl;
        //std::cout << f(a, y, decay_weight, decay_component) << " " << f(b, y, decay_weight, decay_component) << std::endl;

        if (f(c, y, decay_weight, decay_component) == 0.0)
            break;
        else if (f(c, y, decay_weight, decay_component) * f(a, y, decay_weight, decay_component) < 0.0)
            b = c;
        else
            a = c;
    }
    //std::cout << c << " " << kk << std::endl;
    return c;
}

double decay_function(double time, double y, std::vector<double> decay_weight, std::vector<double> decay_component)
{

    double integration_constant = 0;
    for (int i=0;i<decay_weight.size();i++)
        integration_constant += decay_weight[i] * decay_component[i];

    //std::cout<< "integration constant: " << integration_constant << std::endl;

    double output = 0;
    for (int i=0;i<decay_component.size();i++){
        //std::cout << time << " " << decay_weight[i] << " " << decay_component[i] << std::endl;
        output += decay_weight[i] * decay_component[i] * (1 - exp(-time / decay_component[i]));
    }

    //return output - y;
    return (output / integration_constant) - y;
}

bool is_file_exist(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}

std::vector<double> linspace(double a, double b, int N)
{
    double h = (b - a) / (N-1);
    std::vector<double> xs(N);
    typename std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}
