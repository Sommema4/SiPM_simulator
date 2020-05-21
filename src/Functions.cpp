#include "Functions.h"

#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <ctime>
#include <iostream>
#include <cmath>
#include <sys/stat.h>

time_t current_time = time(NULL);
std::mt19937 generator(current_time);
std::uniform_real_distribution<double> uniform_distribution(0.0,1.0);

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
