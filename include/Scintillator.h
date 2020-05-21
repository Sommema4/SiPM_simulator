#ifndef SCINTILLATOR_H
#define SCINTILLATOR_H

#include "Functions.h"

class Scintillator
{
    public:
        Scintillator(std::string, double, std::vector<double>, std::vector<double>);
        virtual ~Scintillator();
        void read_scintillator_lib(std::string);
        void read_emission_spectrum(std::string);
        std::vector<light> generate_light(double);
        void print_scintillator(void);
        int get_photons(void);
        std::string get_name(void);

    protected:

    private:
        std::string name;
        double light_yield;
        int photons;
        std::vector<double> halftime_component;
        std::vector<double> halftime_weight;
        std::vector<double> halftime_amplitude;
        coordinates emission_spectrum;
};

#endif // SCINTILLATOR_H
