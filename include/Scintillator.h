#ifndef SCINTILLATOR_H
#define SCINTILLATOR_H

#include "Functions.h"

class Scintillator
{
    public:
        Scintillator(std::string, std::string, std::string);
        virtual ~Scintillator();
        std::vector<light> generate_light(double);
        void print_scintillator(void);
        int get_photons(void);
        std::string get_name(void);
        std::string get_extension(void);

    protected:

    private:
        std::string name;
        double light_yield;
        int photons;
        std::string particle_type;
        std::string particle_abbrev;
        std::vector<double> halftime_component;
        std::vector<double> halftime_weight;
        std::vector<double> halftime_amplitude;
        coordinates emission_spectrum;
};

#endif // SCINTILLATOR_H
