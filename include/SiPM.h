#ifndef SIPM_H
#define SIPM_H

#include "Microcell.h"
#include "Scintillator.h"
#include "spline.h"
#include "Functions.h"

#include <vector>
#include <map>
#include <random>
#include <ctime>

class SiPM
{
    public:
        friend class Microcell;
        SiPM(std::string, std::string, double, double, double, double, double, double);
        virtual ~SiPM();
        void simulate(std::vector<light>);
        void map_light(std::vector<light>&);
        void dark_current(std::vector<light>&);
        void reset();
        void read_absorption_spectrum(std::string);
        void read_crosstalk_probability(std::string);
        void read_afterpulse_probability(std::string);
        void print_sipm(void);
        void tally_pulse_shape(std::string);
        friend std::ostream& operator<<(std::ostream&, const SiPM *);
        void tally_integration(std::string, double, int);
        void generate_pwl_file(std::string);

        void linearity_simulation(std::string, Scintillator); // simulates the linerarity of SiPM (.csv file)
        void statistical_simulation(std::string, Scintillator, double, int); // simulates lots of pulses and saves their waveforms (.csv file)
        void wide_simulation(std::string, Scintillator, std::vector<double>, int); // simulates lots of pulses with different energies and saves their waveforms in various files (.csv)
        void spice_simulation(std::string, Scintillator, double); // simulate one pulse and saves output so it can be read by LTSpice (.pwl file)

    protected:

    private:
        std::string name;
        double geometry_factor;
        double overvoltage; // overvoltage
        double R_q; // quenching resistor
        double C_q; // capacitance of quenching resistor
        double C_d; //
        double R_d; //
        double R_l; // load resistor
        double C_m; // parasitic capacitance of the grid
        double I_th; // threshold level of current to sustain the avalanche
        int N_c; // number of microcells
        double C_eq; //
        double a1;
        double a2;
        double T_q;
        double T_i;
        double T_d;
        double a_m_1;
        double a_m_2;
        double T_i2;
        double T_d2;
        std::vector<double> afterpulse_weights;
        std::vector<double> afterpulse_components;
        std::vector<double> afterpulse_amplitude;
        double dark_count_rate;
        std::vector<class Microcell> SiPM_microcells;
        std::vector<class Microcell> Microcell_discharge;
        std::vector<class Microcell> Microcell_charge;
        std::map<double, coordinates> absorption_spectrum;
        spline crosstalk_function;
        spline afterpulse_function;
        double timestep;
        double pulse_lenght;
        int output_size;
        std::vector<double> arr_bins;
        std::vector<double> arr_light;
        std::vector<double> arr_afterpulse;
        std::vector<double> arr_crosstalk;
        std::vector<double> arr_dark;
        sipm_par* LUT;
};

#endif // SIPM_H
