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

void read_scintillator_lib(std::string str)
{
    // spectrum must be in cumulative form!!!!!!!!!!
    /*
    std::string name = "./scintillator_lib/" + str;
    std::ifstream file;
    file.open(name);
    if (file.is_open()){
        std::vector<std::string> result;
        while (std::getline(file, line)){
            boost::split(result, line, boost::is_any_of(": "));
            if (line.find("light yield [#/MeV]") != string::npos)
                light_yield = std::stod(result[1]);
            else if (line.find("decay component [ns, %, type]:") != string::npos)

            else if (line.find("") != string::npos)

            else if (line.find("") != string::npos)

            else if (line.find("") != string::npos)

            else if (line.find("") != string::npos)
        }
        file.close();
    }
    else
        std::cout << "The emission spectrum file (" << name << ") could not be open." << std::endl;
    */
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
