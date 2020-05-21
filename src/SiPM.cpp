#include "SiPM.h"
#include "Microcell.h"
#include "Functions.h"
#include "spline.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <iomanip>

SiPM::SiPM(std::string nm, double gf, double ov, double rq, double cq, double cd, double rd, double rl, double cm, int nc, std::vector<double> af_w, std::vector<double> af_c, double dcr, double ts, double pl)
{
    name = nm;
    geometry_factor = gf; // approximation of geometry factor
    overvoltage = ov; // maximum overvoltage
    R_q = rq; // quenching resistor
    C_q = cq; // capacitance of quenching resistor
    C_d = cd; //
    R_d = rd; //
    R_l = rl; // load resistor
    C_m = cm; // parasitic capacitance of the grid
    I_th = 0.0001; // threshold level of current to sustain the avalanche
    N_c = nc; // number of microcells
    C_eq = (N_c - 1) * ((C_d * C_q) / (C_d + C_q)) + N_c * C_m;
    a1 = (R_d * C_d * (R_q + R_l) + R_q * C_q * (R_d + R_l) + C_eq * R_l * (R_q + R_d)) / (R_q + R_d + R_l);
    a2 = ((R_d * R_q * R_l) / (R_q + R_d + R_l)) * (C_eq * (C_d + C_q) + C_d * C_q);
    T_q = R_q * C_q;
    T_i = a2 / a1;
    T_d = a1;
    a_m_1 = R_q * (C_d + C_q) + R_l * (C_eq + C_d);
    a_m_2 = R_q * R_l * (C_eq * (C_d + C_q) + C_d * C_q);
    T_i2 = (2 * a_m_2) / (a_m_1 + sqrt(pow(a_m_1, 2) - 4 * a_m_2));
    T_d2 = (2 * a_m_2) / (a_m_1 - sqrt(pow(a_m_1, 2) - 4 * a_m_2));
    afterpulse_weights = af_w;
    afterpulse_components = af_c;
    dark_count_rate = dcr;
    for (int i=0;i<nc;i++){
        Microcell temp(i, overvoltage);
        SiPM_microcells.push_back(temp);
    };
    timestep = ts;
    pulse_lenght = pl;
    output_size = int(pulse_lenght / timestep);
    spline afterpulse_function;
    spline crosstalk_function;
    for(int k=0;k<output_size;k++){
        arr_bins.push_back(k * timestep);
        arr_light.push_back(0);
        arr_afterpulse.push_back(0);
        arr_crosstalk.push_back(0);
        arr_dark.push_back(0);
    }
    for (int i=0;i<afterpulse_weights.size();i++)
        afterpulse_amplitude.push_back(afterpulse_weights[i] / afterpulse_components[i]);
    //ctor
}

SiPM::~SiPM()
{
    //dtor
}

void SiPM::dark_current(std::vector<light>& buffer)
{
    double mean_dark_count = dark_count_rate * arr_bins.back();
    int number = get_rand_poiss_dist(mean_dark_count);

    for (int i=0;i<number;i++)
    {
        double time = get_rand_uni_dist() * arr_bins.back();
        double wavelenght = 0;
        light temp;
        temp.index = 0;
        temp.time = time;
        temp.wavelenght = 0;
        temp.origin = "dark_current";
        //std::cout << temp.time << " " << temp.origin << std::endl;
        buffer.push_back(temp);
    }
    //std::cout << buffer.size() << std::endl;
    std::sort(buffer.begin(), buffer.end(), compare_light_time);
}

void SiPM::map_light(std::vector<light>& buffer)
{
    int ind;
    for (int i=0;i<buffer.size();i++){
        ind = rand() % N_c;
        //std::cout << ind << std::endl;
        buffer[i].index = ind;
    }
}

void SiPM::simulate(std::vector<light> buffer)
{
    while (buffer.size() != 0)
    {
        SiPM_microcells[buffer[0].index].discharge(buffer[0], buffer, this);
        buffer.erase(buffer.begin());
    }

    for (int i=0;i<SiPM_microcells.size();i++){
        //std::cout << arr_bins[-1] << std::endl;
        SiPM_microcells[i].load_discharge(arr_bins[output_size-1], this);
    }
}

void SiPM::reset()
{
    for (int i=0;i<SiPM_microcells.size();i++)
        SiPM_microcells[i].reset(this);

    for(int i=0;i<output_size;i++){
        arr_light[i] = 0;
        arr_afterpulse[i] = 0;
        arr_crosstalk[i] = 0;
        arr_dark[i] = 0;
    }
}

void SiPM::read_absorption_spectrum(std::string name)
{
    std::vector<double> x, y;
    double last_c = -10;

    std::ifstream file;
    file.open(name);
    if (file.is_open()){
        while(!file.eof()){
            double a, b, c;
            file >> a >> b >> c; // extracts 2 floating point values separated by whitespace
            //std::cout << a << " " << b / 100.0 << " " << c << std::endl;
            x.push_back(a);
            y.push_back(b / 100.0);
            if ((c != last_c && last_c != -10) || file.eof()){
                //std::cout << "Here is Johnny! " << c << std::endl;
                x.pop_back();
                y.pop_back();
                std::vector<double> temp_x, temp_y;
                coordinates temp;
                temp.x = x;
                temp.y = y;
                absorption_spectrum.insert({last_c, temp});
                x.erase(x.begin(),x.end());
                y.erase(y.begin(),y.end());
                x.push_back(a);
                y.push_back(b);
            }
            last_c = c;
        }
    }
    else
        std::cout << "The emission spectrum file (" << name << ") could not be open." << std::endl;

    /*
    for (auto itr = absorption_spectrum.begin(); itr != absorption_spectrum.end(); ++itr) {
        std::cout << itr->first << std::endl;
        std::vector<double> temp1 = (itr->second).x;
        std::vector<double> temp2 = (itr->second).y;
        for (int i=0;i<temp1.size();i++)
            std::cout << temp1[i] << " " << temp2[i] << std::endl;
    }
    */
}

void SiPM::read_crosstalk_probability(std::string name)
{
    std::vector<double> x, y;
    std::ifstream file;
    file.open(name);
    if (file.is_open()){
        while(!file.eof()){
            double a, b;
            file >> a >> b; // extracts 2 floating point values separated by whitespace
            //std::cout << a << " " << b << std::endl;
            x.push_back(a);
            y.push_back(b);
        }
    }
    else
        std::cout << "The emission spectrum file (" << name << ") could not be open." << std::endl;

    //spline crosstalk_function;
    crosstalk_function.set_points(y, x);
}

void SiPM::read_afterpulse_probability(std::string name)
{
    std::vector<double> x, y;
    std::ifstream file;
    file.open(name);

    if (file.is_open()){
        while(!file.eof()){
            double a, b;
            file >> a >> b; // extracts 2 floating point values separated by whitespace
            //std::cout << a << " " << b << std::endl;
            x.push_back(a);
            y.push_back(b);
        }
    }
    else
        std::cout << "The emission spectrum file (" << name << ") could not be open." << std::endl;

    //spline afterpulse_function;
    afterpulse_function.set_points(x, y);
}

void SiPM::print_sipm(void)
{
    std::cout << std::endl;
    std::cout << "----------SiPM PARAMETERS----------" << std::endl;
    std::cout << "name: " << name << std::endl;
    std::cout << "geometry_factor: " << geometry_factor << std::endl;
    std::cout << "overvoltage: " << overvoltage << " V" << std::endl;
    std::cout << "R_q: " << R_q << " ohm" << std::endl;
    std::cout << "C_q: " << C_q << " F" << std::endl;
    std::cout << "C_d: " << C_d << " F" << std::endl;
    std::cout << "R_d: " << R_d << " ohm" << std::endl;
    std::cout << "R_l: " << R_l << " ohm" << std::endl;
    std::cout << "C_m: " << C_m << " F" << std::endl;
    std::cout << "I_th: " << I_th << " A" << std::endl;
    std::cout << "N_c: " << N_c << " #" << std::endl;
    std::cout << "C_eq: " << C_eq << " F" << std::endl;
    std::cout << "a1: " << a1 << " ?" << std::endl;
    std::cout << "a2: " << a2 << " ?" << std::endl;
    std::cout << "T_q: " << T_q << " s" << std::endl;
    std::cout << "T_i: " << T_i << " s" << std::endl;
    std::cout << "T_d: " << T_d << " s" << std::endl;
    std::cout << "a_m_1: " << a_m_1 << " ?" << std::endl;
    std::cout << "a_m_2: " << a_m_2 << " ?" << std::endl;
    std::cout << "T_i2: " << T_i2 << " s" << std::endl;
    std::cout << "T_d2: " << T_d2 << " s" << std::endl;
    std::cout << "Microcell vector size: " << SiPM_microcells.size() << " #" << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << std::endl;
}

std::ostream& operator<<(std::ostream& os, const SiPM* en)
{
    for (int i=0;i<en->arr_bins.size();i++)
        os << en->arr_bins[i] << "," << en->arr_light[i] << "," << en->arr_afterpulse[i] << "," << en->arr_crosstalk[i] << "," << en->arr_dark[i] << std::endl;
    return os;
}

void SiPM::tally_pulse_shape(std::string name)
{
    std::ofstream myfile;
    if (is_file_exist(name)){
        myfile.open(name, std::ios_base::app);
        for (int i=0;i<arr_light.size();i++)
            myfile << arr_light[i] << " ";
        myfile << std::endl;
        for (int i=0;i<arr_afterpulse.size();i++)
            myfile << arr_afterpulse[i] << " ";
        myfile << std::endl;
        for (int i=0;i<arr_crosstalk.size();i++)
            myfile << arr_crosstalk[i] << " ";
        myfile << std::endl;
        for (int i=0;i<arr_dark.size();i++)
            myfile << arr_dark[i] << " ";
        myfile << std::endl;
    }
    else{
        myfile.open(name);
        for (int i=0;i<arr_bins.size();i++)
            myfile << arr_bins[i] << " ";
        myfile << std::endl;
        for (int i=0;i<arr_light.size();i++){
            myfile << arr_light[i] << " ";
            //std::cout << arr_light[i] << std::endl;
        }

        myfile << std::endl;
        for (int i=0;i<arr_afterpulse.size();i++)
            myfile << arr_afterpulse[i] << " ";
        myfile << std::endl;
        for (int i=0;i<arr_crosstalk.size();i++)
            myfile << arr_crosstalk[i] << " ";
        myfile << std::endl;
        for (int i=0;i<arr_dark.size();i++)
            myfile << arr_dark[i] << " ";
        myfile << std::endl;
    }

    //myfile << this;
    myfile.close();
}

void SiPM::tally_integration(std::string filename, double energy, int photons)
{
    double integ_light = integrate(arr_light);
    double integ_crosstalk = integrate(arr_crosstalk);
    double integ_afterpulse = integrate(arr_afterpulse);
    double integ_dark_current = integrate(arr_dark);

    std::ofstream myfile;
    myfile.open(filename, std::ios_base::app);
    myfile << energy << "," << photons << "," << integ_light << "," << integ_crosstalk << "," << integ_afterpulse << "," << integ_dark_current << std::endl;
    std::cout << energy << "," << photons << "," << integ_light << "," << integ_crosstalk << "," << integ_afterpulse << "," << integ_dark_current << std::endl;
    myfile.close();
}

void SiPM::generate_pwl_file(std::string name)
{
    std::ofstream myfile;
    myfile.open(name);
    for (int i=0;i<arr_bins.size();i++)
        myfile << std::fixed << std::setprecision(12) << arr_bins[i] << " " << arr_light[i] << std::endl;
    myfile.close();
}

void SiPM::linearity_simulation(std::string filename, Scintillator crystal)
{
    std::vector<double> energy_range = linspace(0.1, 50, 40);
    for (int i=0;i<energy_range.size();i++)
    {
        std::cout << "Energy: " << energy_range[i] << std::endl;
        reset(); // reset each microcell of sipm before simulation
        std::vector<light> primaries = crystal.generate_light(energy_range[i]); // generate primaries, primaries is array of light
        int photons = crystal.get_photons();
        dark_current(primaries); // add dark current to primaries
        map_light(primaries); // maps primary photons to individual microcells
        simulate(primaries); // output is array of tally
        tally_integration(filename, energy_range[i], photons); // prints integral of the pulse
    }
}

void SiPM::statistical_simulation(std::string filename, Scintillator crystal, double deposited_energy, int n)
{
    for (int i=0;i<n;i++){
        std::cout << i << "/" << n << std::endl;
        reset(); // reset each microcell of sipm before simulation
        std::vector<light> primaries = crystal.generate_light(deposited_energy); // generate primaries, primaries is array of light
        dark_current(primaries); // add dark current to primaries
        map_light(primaries); // maps primary photons to individual microcells
        print_primaries(primaries, 0); // show first 5 photons in container
        simulate(primaries); // output is array of tally
        tally_pulse_shape(filename); // tally the matrix of pulse shape including total, afterpulse, crosstalk and dark current
    }
}

void SiPM::wide_simulation(std::string extension, Scintillator crystal, std::vector<double> deposited_energy, int rep)
{
    for (int i=0;i<deposited_energy.size();i++)
    {
        std::cout << deposited_energy[i] << " MeV" << std::endl;
        std::string sc = crystal.get_name();
        std::string filename = name + "_" + sc + "_" + std::to_string(deposited_energy[i]) + extension + ".csv";
        for (int j=0;j<rep;j++){
            std::cout << j << "/" << rep << std::endl;
            reset(); // reset each microcell of sipm before simulation
            std::vector<light> primaries = crystal.generate_light(deposited_energy[i]); // generate primaries, primaries is array of light
            dark_current(primaries); // add dark current to primaries
            map_light(primaries); // maps primary photons to individual microcells
            //print_primaries(primaries, 0); // show first 5 photons in container
            simulate(primaries); // output is array of tally
            tally_pulse_shape(filename); // tally the matrix of pulse shape including total, afterpulse, crosstalk and dark current
        }
    }
}

void SiPM::spice_simulation(std::string filename, Scintillator crystal, double deposited_energy)
{
    reset(); // reset each microcell of sipm before simulation
    std::vector<light> primaries = crystal.generate_light(deposited_energy); // generate primaries, primaries is array of light
    dark_current(primaries); // add dark current to primaries
    map_light(primaries); // maps primary photons to individual microcells
    simulate(primaries); // output is array of tally
    generate_pwl_file(filename);
}
