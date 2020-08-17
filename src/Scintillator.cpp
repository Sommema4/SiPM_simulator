#include "Scintillator.h"
#include "Functions.h"
#include "spline.h"
#include "json.hpp"

#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <ctime>


Scintillator::Scintillator(std::string scint_name, std::string scint_lib, std::string par_type)
{
    std::ifstream i(scint_lib);
    nlohmann::json scin;
    i >> scin;
    assert(scin.find(scint_name) !=  scin.end());
    name = scint_name;
    light_yield = scin[scint_name]["light_yield"]["value"];
    particle_type = par_type;
    particle_abbrev = scin[scint_name]["particle"][par_type]["abbreviation"].get<std::string>();
    halftime_component = scin[scint_name]["particle"][particle_type]["decay component"].get<std::vector<double>>();
    halftime_weight = scin[scint_name]["particle"][particle_type]["decay intensity"].get<std::vector<double>>();
    assert(halftime_component.size() == halftime_weight.size()); // test the size of component vectors
    for(int i=0;i<halftime_component.size();i++)
        halftime_amplitude.push_back(halftime_weight[i] / halftime_component[i]);
    emission_spectrum.x = scin[scint_name]["emission spectrum"]["wavelenght"].get<std::vector<double>>();
    emission_spectrum.y = scin[scint_name]["emission spectrum"]["cumulative probability"].get<std::vector<double>>();
    //ctor
}

Scintillator::~Scintillator()
{
    //dtor
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

std::string Scintillator::get_extension(void)
{
    return particle_abbrev;
}
