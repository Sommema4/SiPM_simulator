#include "Scintillator.h"
#include "Functions.h"
#include "spline.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <ctime>


Scintillator::Scintillator(std::string nm, double ly, std::vector<double> h_c, std::vector<double> h_w)
{
    std::ifstream i(lib);
    json init;
    i >> init;

    std::string scint_name[2] = init[name]["name"]; // the abbreviation of scintillator and its full name
    std::string scint_lib = init["scintillator"]["library"]; // the scintillator library name
    std::string particle[2] = init["scintillator"]["particle_type"]; // the abbreviation of particle and its full name
    std::bool photon_list = init["scintillator"]["photon list"]; // decides if list of all generated photons should be written to auxiliary output [boolean]

    assert(h_c.size() == h_w.size()); // test the size of component vectors
    name = nm;
    light_yield = ly;
    halftime_component = h_c;
    halftime_weight = h_w;
    for (int i=0;i<halftime_component.size();i++){
        halftime_amplitude.push_back(halftime_weight[i] / halftime_component[i]);
        //std::cout << halftime_amplitude[i] << std::endl;
    }
    //ctor
}

Scintillator::~Scintillator()
{
    //dtor
}

void Scintillator::read_emission_spectrum(std::string name)
{
    // spectrum must be in cumulative form!!!!!!!!!!
    // .txt file must be tab delimited

    std::ifstream file;
    file.open(name);
    if (file.is_open()){
        //std::cout << name << std::endl;
        double a, b;
        while(file >> a >> b)
        {
            emission_spectrum.x.push_back(a);
            emission_spectrum.y.push_back(b);
        }
        file.close();
    }
    else
        std::cout << "The emission spectrum file (" << name << ") could not be open." << std::endl;
}

std::vector<light> Scintillator::generate_light(double deposited_energy)
{
    // spectrum must be in cumulative form!!!!!!!!!!


    photons = int((deposited_energy * light_yield) + 0.5); // round to integer
    std::vector<light> output;
    spline spectrum;
    spectrum.set_points(emission_spectrum.y, emission_spectrum.x, true);

    for (int i=0;i<photons;i++){
        double dice = get_rand_uni_dist();
        double time = bisection_method(0.0, 0.001, dice, 1e-15, halftime_amplitude, halftime_component, decay_function);
        dice = get_rand_uni_dist();
        double wavelenght = spectrum(dice);
        light temp;
        temp.index = 0;
        temp.time = time;
        temp.wavelenght = wavelenght;
        temp.origin = "scintillator_light";
        //std::cout << temp.time << " " << temp.wavelenght << " " << temp.origin << std::endl;
        output.push_back(temp);
    }
    sort_vector(output);

    return output;
}

void Scintillator::print_scintillator(void)
{
    std::cout << std::endl;
    std::cout << "----------SCINTILLATOR PARAMETERS---------- " << std::endl;
    std::cout << "name: " << name << std::endl;
    std::cout << "light_yield: " << light_yield << " photons/MeV" << std::endl;
    std::cout << "halftime components: ";
    for (int i=0; i<halftime_component.size();i++)
        std::cout << halftime_component[i] << " s, ";
    std::cout << std::endl;
    std::cout << "halftime component weights: ";
    for (int i=0; i<halftime_weight.size();i++)
        std::cout << halftime_weight[i] << " ";
    std::cout << std::endl;
    std::cout << "------------------------------------------- " << std::endl;
    std::cout << std::endl;
}

int Scintillator::get_photons(void)
{
    return photons;
}

std::string Scintillator::get_name(void)
{
    return name;
}
