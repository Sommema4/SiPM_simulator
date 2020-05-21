#include "SiPM.h"
#include "Functions.h"
#include "Scintillator.h"

#include <iostream>

int main()
{
    /* DEFINE SEMI-GLOBAL VARIABLES */
    int repetition = 100;
    double timestep = 10*1e-12; // 0.01 ns timestep
    double pulse_lenght = 2*1e-06;

    std::vector<double> deposited_energy = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0}; // [MeV]
    //std::vector<double> deposited_energy = {0.1, 0.2}; // [MeV]
    double light_yield = 8600; // number of photons per 1 MeV

    /* STILBENE */
    /*
    std::vector<double> decay_component {5.21 * 1e-09, 21.33 * 1e-09, 134.77 * 1e-09}; // GAMMA
    std::vector<double> decay_weight {0.95, 0.03, 0.02}; // GAMMA
    //std::vector<double> decay_component {5.01 * 1e-09, 27.7 * 1e-09, 253.19 * 1e-09}; // NEUTRON
    //std::vector<double> decay_weight {0.95, 0.04, 0.01}; // NEUTRON
    */

    /* NE-213 */
    /*
    //std::vector<double> decay_component {4.1 * 1e-09, 32 * 1e-09, 160 * 1e-09, 1870 * 1e-09}; // GAMMA
    //std::vector<double> decay_weight {0.89, 0.07, 0.03, 0.01}; // GAMMA
    //std::vector<double> decay_component {5.2 * 1e-09, 32 * 1e-09, 130 * 1e-09, 510 * 1e-09}; // NEUTRON
    //std::vector<double> decay_weight {0.51, 0.23, 0.17, 0.09}; // NEUTRON
    */

    /* PLASTIC SCINTILLATOR EJ-276 */

    std::vector<double> decay_component {4.0 * 1e-09, 16.0 * 1e-09, 98.0 * 1e-09, 690.0 * 1e-09}; // GAMMA
    std::vector<double> decay_weight {0.71, 0.12, 0.08, 0.09}; // GAMMA
    //std::vector<double> decay_component {3.9 * 1e-09, 18.0 * 1e-09, 106.0 * 1e-09, 800.0 * 1e-09}; // NEUTRON
    //std::vector<double> decay_weight {0.47, 0.13, 0.13, 0.27}; // NEUTRON


    print_seed(); // prints seed number which is based on current time

    /* DEFINE SCINTILLATOR AND READ ITS EMISSION SPECTRUM */
    Scintillator crystal("EJ276", light_yield, decay_component, decay_weight);
    // ARGUMENTS: deposited energy, light yield, halflife components, halflife component weights
    crystal.read_emission_spectrum("EJ-276.txt"); // read emission spectrum - spectrum must be in cumulative form
    crystal.print_scintillator();

    /* DEFINE SiPM AND SAVE ITS ABSORPTION SPECTRUM */
    std::vector<double> afterpulse_weights = {0.5, 0.5};
    std::vector<double> afterpulse_components = {10*1e-09, 100*1e-09};
    SiPM sipm("MicroFC60035", 1.0, 2.5, 633000, 20 * 1e-15, 170.1 * 1e-15, 300, 10.0, 3 * 1e-15, 18980, afterpulse_weights, afterpulse_components, 10 * 1e+06, timestep, pulse_lenght);
    // ARGUMENTS: geometry factor, overvoltage, R_q, C_q, C_d, R_d, R_load, C_m, number of microcells, af_ts, af_tl, dark count rate, timestep, pulse lenght
    sipm.read_absorption_spectrum("MicroFC_60035.txt"); // read absorption spectrum - spectrum must be in its natural form
    sipm.read_afterpulse_probability("Afterpulse.txt"); // read absorption spectrum - spectrum must be in its natural form
    sipm.read_crosstalk_probability("Crosstalk.txt"); // read absorption spectrum - spectrum must be in its natural form
    sipm.print_sipm();

    /* RUN SIMULATION ROUTINE IN FOR CYCLE */
    //sipm.linearity_simulation("Integration.csv", crystal);
    //sipm.statistical_simulation("Photon_pulse.csv", crystal, deposited_energy[0], repetition);
    sipm.wide_simulation("_p", crystal, deposited_energy, repetition);
    //sipm.spice_simulation("Pulse.pwl", crystal, deposited_energy[0]);

    return 0;
}
