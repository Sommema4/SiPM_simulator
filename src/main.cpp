#include "SiPM.h"
#include "Functions.h"
#include "Scintillator.h"
#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

#include <iostream>

int main()
{
    settings set; // struct defined in Functions.h - serve only for storing a setting data
    parse_settings_json(&set, "settings.json"); // reads settings file and saves data to settings struct

    // Create a directory
    std::string path = set.out_dir + set.sim_name;
    if (mkdir(path.c_str(), 0777) == -1)
        std::cerr << strerror(errno) << std::endl;
    else
        std::cout << "Directory created: " << path << std::endl;

    print_seed(); // prints seed number which is based on current time

    /* CONSTRUCTS SCINTILLATOR AND READS DATA FROM SCINTILLATOR LIBRARY */
    Scintillator crystal(set.scint_name, set.scint_lib, set.particle);
    crystal.print_scintillator();

    /* CONSTRUCTS SiPM AND READS DATA FROM SiPM LIBRARY */
    SiPM sipm(set.sipm_name, set.sipm_lib, set.overvoltage, set.geometry_factor, set.R_l, set.I_th, set.timestep, set.pulse_lenght, set.afterpulse, set.crosstalk, set.dark_current);
    sipm.print_sipm();

    /* RUN SIMULATION ROUTINE */
    if (set.sim_type == "waveform")
        sipm.waveform_simulation(path, crystal, set.deposited_energy, set.repetition);
    else if(set.sim_type == "spice")
        sipm.spice_simulation("Pulse.pwl", crystal, set.deposited_energy[0]);

    return 0;
}
