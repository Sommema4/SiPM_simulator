#include "SiPM.h"
#include "Functions.h"
#include "Scintillator.h"
#include "json.hpp"

#include <iostream>
#include <fstream>

int main()
{
    using json = nlohmann::json;
    std::ifstream i("setting.json");
    json init;
    i >> init;

    /* SIMULATION VARIABLES */
    std::string sim_name = init["simulation"]["name"]; // the name of your simulation
    std::string out_dir = init["simulation"]["directory"]; // the name of the directory where results will be saved - relative to the binary file
    std::string sim_type = init["simulation"]["calculation_type"]; // the type of the simulation
    bool out_aux = init["simulation"]["aux_output"]; // decides if auxiliary output should be written to file [boolean]

    /* GLOBAL VARIABLES */
    int repetition = init["global variables"]["repetition"]; // the number of simuation repetitons for each deposited energy
    double timestep = init["global variables"]["timestep"]; // the lenght of timestep in the simulation [seconds]
    double pulse_lenght = init["global variables"]["pulse_lenght"]; // the maximum lenght of the calculation [seconds]
    std::vector<double> deposited_energy = init["global variables"]["deposited_energy"]; // the list of energies of the simulation [MeV]

    /* SCINTILLATOR VARIABLES */
    std::string scint_name = init["scintillator"]["name"]; // the full name of scintillator
    std::string scint_lib = init["scintillator"]["library"]; // the scintillator library name
    std::string particle = init["scintillator"]["particle_type"]; // the abbreviation of particle and its full name
    bool photon_list = init["scintillator"]["photon list"]; // decides if list of all generated photons should be written to auxiliary output [boolean]

    /* SiPM VARIABLES */
    std::string sipm_name = init["sipm"]["name"]; // the full name of sipm
    std::string sipm_lib = init["sipm"]["library"]; // the sipm library name
    bool crosstalk = init["sipm"]["crosstalk"]; // decides if crosstalk should be taken into account during calculation [boolean]
    bool afterpulse = init["sipm"]["afterpulse"]; // decides if afterpulse should be taken into account during calculation [boolean]
    bool dark_current = init["sipm"]["dark_current"]; // decides if dark current should be taken into account during calculation [boolean]

    R_l = sipm[sipm_name]["name"]; // load resistor
    I_th = sipm[sipm_name]["name"]; // threshold level of current to sustain the avalanche = 0.0001 A

    /* SIMULATION PARAMETERS */
    std::string interpolation = init["simulation parameters"]["interpolation"]; // chooses the type of interpolation technique
    long double seed = init["simulation parameters"]["seed"]; // the seed number, if seed=0 the number is chosen based on the time

    print_seed(); // prints seed number which is based on current time

    /* CONSTRUCTS SCINTILLATOR AND READS ATA FROM SCINTILLATOR LIBRARY */
    Scintillator crystal(scint_name, scint_lib, particle);
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
    //sipm.wide_simulation("_p", crystal, deposited_energy, repetition);
    //sipm.spice_simulation("Pulse.pwl", crystal, deposited_energy[0]);

    return 0;
}
