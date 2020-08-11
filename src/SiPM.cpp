#include "SiPM.h"
#include "Microcell.h"
#include "Functions.h"
#include "spline.h"
#include "json.hpp"

#include <assert.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <iomanip>

#ifndef NDEBUG
#   define M_Assert(Expr, Msg) \
    __M_Assert(#Expr, Expr, __FILE__, __LINE__, Msg)
#else
#   define M_Assert(Expr, Msg) ;
#endif

void __M_Assert(const char* expr_str, bool expr, const char* file, int line, const char* msg)
{
    if (!expr)
    {
        std::cerr << "Assert failed:\t" << msg << "\n"
            << "Expected:\t" << expr_str << "\n"
            << "Source:\t\t" << file << ", line " << line << "\n";
        abort();
    }
}

SiPM::SiPM(std::string sipm_name, std::string sipm_lib, double ov, double geo_factor, double load_resistor, double current_threshold, double time_step, double pls_lenght)
{
    /* Open up a .json file */
    std::ifstream i(sipm_lib);
    nlohmann::json sipm;
    i >> sipm;

    /* Save inner parameters of SiPM to SiPM object */
    assert(sipm.find(sipm_name) !=  sipm.end());
    name = sipm_name;
    R_q = sipm[sipm_name]["inner parameters"]["R_q"]; // quenching resistor
    C_q = sipm[sipm_name]["inner parameters"]["C_q"]; // capacitance of quenching resistor
    C_d = sipm[sipm_name]["inner parameters"]["C_d"];
    R_d = sipm[sipm_name]["inner parameters"]["R_d"];
    C_m = sipm[sipm_name]["inner parameters"]["C_m"]; // parasitic capacitance of the grid
    N_c = sipm[sipm_name]["inner parameters"]["N_c"]; // number of microcells
    T_q = R_q * C_q;

    /* Parse global simulation parameters to SiPM object */
    geometry_factor = geo_factor;
    overvoltage = ov;
    I_th = current_threshold; // threshold current to sustain an avalanche
    R_l = load_resistor;

    /* Save parameters of noise effects to SiPM object and use spline interpolation */
    double a;
    M_Assert(a = sipm[sipm_name]["noise effects"]["dark current"]["dark acount rate"], "asasd");
    //M_Assert(("A must be equal to B", a == sipm[sipm_name]["noise effects"]["dark current"]["dark acount rate"]));
    //assert(a = sipm[sipm_name]["noise effects"]["dark current"]["dark count rate"]);
    std::cout << "spatny a: " << a << std::endl;
    afterpulse_weights = sipm[sipm_name]["noise effects"]["afterpulse"]["afterpulse intensity"].get<std::vector<double>>();
    afterpulse_components = sipm[sipm_name]["noise effects"]["afterpulse"]["afterpulse component"].get<std::vector<double>>();
    for (int i=0;i<afterpulse_weights.size();i++)
        afterpulse_amplitude.push_back(afterpulse_weights[i] / afterpulse_components[i]);
    std::vector<double> x_afterpulse = sipm[sipm_name]["noise effects"]["afterpulse"]["function"]["overvoltage"].get<std::vector<double>>();
    std::vector<double> y_afterpulse = sipm[sipm_name]["noise effects"]["afterpulse"]["function"]["probability"].get<std::vector<double>>();
    std::vector<double> x_crosstalk = sipm[sipm_name]["noise effects"]["crosstalk"]["function"]["kacer"].get<std::vector<double>>();
    std::vector<double> y_crosstalk = sipm[sipm_name]["noise effects"]["crosstalk"]["function"]["probability"].get<std::vector<double>>();
    spline afterpulse_function;
    afterpulse_function.set_points(x_afterpulse, y_afterpulse);
    spline crosstalk_function;
    crosstalk_function.set_points(y_crosstalk, x_crosstalk);
    dark_count_rate = sipm[sipm_name]["noise effects"]["dark current"]["dark count rate"].get<double>();

    /* Save absorption spectrum to SiPM object */
    for (auto& el : sipm[sipm_name]["absorption spectrum"].items()){
        if(el.key() != "source"){
            coordinates temp;
            double voltage = sipm[sipm_name]["absorption spectrum"][el.key()]["value"].get<double>();
            //std::cout << voltage << std::endl;
            temp.x = sipm[sipm_name]["absorption spectrum"][el.key()]["wavelenght"].get<std::vector<double>>();
            temp.y = sipm[sipm_name]["absorption spectrum"][el.key()]["probability"].get<std::vector<double>>();
            absorption_spectrum.insert({voltage, temp});
        }
    }

    /* Create and save a look-up table to SiPM object */
    LUT = new sipm_par[N_c+1];
    for (int i=0;i<N_c+1;i++){
        LUT[i].C_eq = (N_c - i) * ((C_d * C_q) / (C_d + C_q)) + N_c * C_m;
        LUT[i].a1 = (R_d * C_d * (R_q + i * R_l) + R_q * C_q * (R_d + i * R_l) + LUT[i].C_eq * R_l * (R_q + R_d)) / (R_q + R_d + i * R_l);
        LUT[i].a2 = ((R_d * R_q * R_l) / (R_q + R_d + i * R_l)) * (LUT[i].C_eq * (C_d + C_q) + i * C_d * C_q);
        LUT[i].T_i = (2 * LUT[i].a2) / LUT[i].a1 + sqrt(pow(LUT[i].a1, 2) - 4 * LUT[i].a2);
        LUT[i].T_d = (2 * LUT[i].a2) / LUT[i].a1 - sqrt(pow(LUT[i].a1, 2) - 4 * LUT[i].a2);
        LUT[i].a_m_1 = R_q * (C_d + C_q) + R_l * (LUT[i].C_eq + i * C_d);
        LUT[i].a_m_2 = R_q * R_l * (LUT[i].C_eq * (C_d + C_q) + i * C_d * C_q);
        LUT[i].T_i2 = (2 * LUT[i].a_m_2) / (LUT[i].a_m_1 + sqrt(pow(LUT[i].a_m_1, 2) - 4 * LUT[i].a_m_2));
        LUT[i].T_d2 = (2 * LUT[i].a_m_2) / (LUT[i].a_m_1 - sqrt(pow(LUT[i].a_m_1, 2) - 4 * LUT[i].a_m_2));
    }

    for (int i=0;i<N_c;i++){
        Microcell temp(i, overvoltage);
        SiPM_microcells.push_back(temp);
    };

    timestep = time_step;
    pulse_lenght = pls_lenght;
    output_size = int(pulse_lenght / timestep);

    for(int k=0;k<output_size;k++){
        arr_bins.push_back(k * timestep);
        arr_light.push_back(0);
        arr_afterpulse.push_back(0);
        arr_crosstalk.push_back(0);
        arr_dark.push_back(0);
    }
    //ctor
}

SiPM::~SiPM()
{
    if(LUT) // True if LUT is not a null pointer
        delete[] LUT;
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
        buffer.push_back(temp);
    }
    std::sort(buffer.begin(), buffer.end(), compare_light_time);
}

void SiPM::map_light(std::vector<light>& buffer)
{
    int ind;
    for (int i=0;i<buffer.size();i++){
        ind = rand() % N_c;
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
