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

SiPM::SiPM(std::string sipm_name, std::string sipm_lib, double ov, double geo_factor, double load_resistor, double current_threshold, double time_step, double pls_lenght, bool ap, bool ct, bool dc)
{
    /* Open up a .json file */
    std::ifstream i(sipm_lib);
    nlohmann::json sipm;
    i >> sipm;

    /* Check the existence of parameters in .json file */
    std::cout << sipm_name << std::endl;
    M_Assert(sipm.find(sipm_name) != sipm.end(), "SiPM was not found in the .json file");
    M_Assert(sipm[sipm_name].find("inner parameters") != sipm[sipm_name].end(), "'inner parameters' was not found in the .json file");
    M_Assert(sipm[sipm_name]["inner parameters"].find("R_q") != sipm[sipm_name]["inner parameters"].end(), "'R_q' was not found in the .json file");
    M_Assert(sipm[sipm_name]["inner parameters"].find("C_q") != sipm[sipm_name]["inner parameters"].end(), "'C_q' was not found in the .json file");
    M_Assert(sipm[sipm_name]["inner parameters"].find("C_d") != sipm[sipm_name]["inner parameters"].end(), "'C_d' was not found in the .json file");
    M_Assert(sipm[sipm_name]["inner parameters"].find("R_d") != sipm[sipm_name]["inner parameters"].end(), "'R_d' was not found in the .json file");
    M_Assert(sipm[sipm_name]["inner parameters"].find("C_m") != sipm[sipm_name]["inner parameters"].end(), "'C_m' was not found in the .json file");
    M_Assert(sipm[sipm_name]["inner parameters"].find("N_c") != sipm[sipm_name]["inner parameters"].end(), "'N_c' was not found in the .json file");
    M_Assert(sipm[sipm_name].find("noise effects") != sipm[sipm_name].end(), "'noise effects' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"].find("afterpulse") != sipm[sipm_name]["noise effects"].end(), "'afterpulse' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"].find("crosstalk") != sipm[sipm_name]["noise effects"].end(), "'crosstalk' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"].find("dark current") != sipm[sipm_name]["noise effects"].end(), "'dark current' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"]["afterpulse"].find("function") != sipm[sipm_name]["noise effects"]["afterpulse"].end(), "'function' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"]["crosstalk"].find("function") != sipm[sipm_name]["noise effects"]["crosstalk"].end(), "'function' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"]["afterpulse"].find("afterpulse intensity") != sipm[sipm_name]["noise effects"]["afterpulse"].end(), "'afterpulse intensity' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"]["afterpulse"].find("afterpulse component") != sipm[sipm_name]["noise effects"]["afterpulse"].end(), "'afterpulse component' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"]["afterpulse"]["function"].find("overvoltage") != sipm[sipm_name]["noise effects"]["afterpulse"]["function"].end(), "'overvoltage' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"]["afterpulse"]["function"].find("probability") != sipm[sipm_name]["noise effects"]["afterpulse"]["function"].end(), "'probability' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"]["crosstalk"]["function"].find("qe") != sipm[sipm_name]["noise effects"]["crosstalk"]["function"].end(), "'overvoltage' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"]["crosstalk"]["function"].find("probability") != sipm[sipm_name]["noise effects"]["crosstalk"]["function"].end(), "'probability' was not found in the .json file");
    M_Assert(sipm[sipm_name]["noise effects"]["dark current"].find("dark count rate") != sipm[sipm_name]["noise effects"]["dark current"].end(), "'dark count rate' was not found in the .json file");
    M_Assert(sipm[sipm_name].find("absorption spectrum") != sipm[sipm_name].end(), "'absorption spectrum' was not found in the .json file");

    /* Save inner parameters of SiPM to SiPM object */
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
    afterpulse_swich = ap;
    crosstalk_swich = ct;
    dark_current_swich = dc;

    /* Save parameters of noise effects to SiPM object and use spline interpolation */
    afterpulse_weights = sipm[sipm_name]["noise effects"]["afterpulse"]["afterpulse intensity"].get<std::vector<double>>();
    afterpulse_components = sipm[sipm_name]["noise effects"]["afterpulse"]["afterpulse component"].get<std::vector<double>>();
    for (int i=0;i<afterpulse_weights.size();i++)
        afterpulse_amplitude.push_back(afterpulse_weights[i] / afterpulse_components[i]);
    std::vector<double> x_afterpulse = sipm[sipm_name]["noise effects"]["afterpulse"]["function"]["overvoltage"].get<std::vector<double>>();
    std::vector<double> y_afterpulse = sipm[sipm_name]["noise effects"]["afterpulse"]["function"]["probability"].get<std::vector<double>>();
    std::vector<double> x_crosstalk = sipm[sipm_name]["noise effects"]["crosstalk"]["function"]["qe"].get<std::vector<double>>();
    std::vector<double> y_crosstalk = sipm[sipm_name]["noise effects"]["crosstalk"]["function"]["probability"].get<std::vector<double>>();
    afterpulse_function.set_points(x_crosstalk, y_crosstalk);
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

    /* Create and save a look-up table into SiPM object */
    LUT = new sipm_par[N_c+1];
    for (int i=0;i<N_c+1;i++){
        LUT[i].C_eq = (N_c - i) * ((C_d * C_q) / (C_d + C_q)) + N_c * C_m;
        LUT[i].a1 = (R_d * C_d * (R_q + i * R_l) + R_q * C_q * (R_d + i * R_l) + LUT[i].C_eq * R_l * (R_q + R_d)) / (R_q + R_d + i * R_l);
        LUT[i].a2 = ((R_d * R_q * R_l) / (R_q + R_d + i * R_l)) * (LUT[i].C_eq * (C_d + C_q) + i * C_d * C_q);
        LUT[i].T_i = (2 * LUT[i].a2) / (LUT[i].a1 + sqrt(pow(LUT[i].a1, 2) - 4 * LUT[i].a2));
        LUT[i].T_d = (2 * LUT[i].a2) / (LUT[i].a1 - sqrt(pow(LUT[i].a1, 2) - 4 * LUT[i].a2));
        LUT[i].a_m_1 = R_q * (C_d + C_q) + R_l * (LUT[i].C_eq + i * C_d);
        LUT[i].a_m_2 = R_q * R_l * (LUT[i].C_eq * (C_d + C_q) + i * C_d * C_q);
        LUT[i].T_i2 = (2 * LUT[i].a_m_2) / (LUT[i].a_m_1 + sqrt(pow(LUT[i].a_m_1, 2) - 4 * LUT[i].a_m_2));
        LUT[i].T_d2 = (2 * LUT[i].a_m_2) / (LUT[i].a_m_1 - sqrt(pow(LUT[i].a_m_1, 2) - 4 * LUT[i].a_m_2));
    }

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
    /* RANDOMLY GENERATES PHOTONS CREATED DUE TO DARK CURRENT */
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
    /* RANDOMLY DECIDES WHAT MICROCELL WAS HIT BY PHOTON */
    int ind;
    for (int i=0;i<buffer.size();i++){
        ind = rand() % N_c;
        buffer[i].index = ind;
    }
}

void SiPM::simulate(std::vector<light> buffer)
{
    for (int i=0;i<arr_bins.size();i++){ // loop over all discrete times of the simulation
        if (buffer.size() == 0 && Microcell_discharge.size() == 0 && Microcell_charge.size() == 0)
            break;

        /* While statement fills discharge and charge container with fired microcells */
        while (std::abs(buffer[0].time - arr_bins[i]) < (timestep / 2.0)){ // if true the photon hits a microcell
            auto it = find_if(Microcell_discharge.begin(), Microcell_discharge.end(), [&buffer] (Microcell& obj) {return obj.get_ID() == buffer[0].index;}); // find iterator to object with the same index in discharge container
            if (it != Microcell_discharge.end()){ // the microcell is in the discharge container already
                buffer.erase(buffer.begin()); // photon is deleted - the microcell cannot be discharged two times at the same time
            }
            it = find_if(Microcell_charge.begin(), Microcell_charge.end(), [&buffer] (Microcell& obj) {return obj.get_ID() == buffer[0].index;}); // find iterator to object with the same index in charge container
            if (it != Microcell_charge.end()){ // the microcell is in the charge container already
                if (it->discharge_init(buffer[0], arr_bins[i], this, buffer)){ // check whether photon fires the microcell
                    //std::cout << "asdasd" << std::endl;
                    Microcell_discharge.push_back(*it);
                    Microcell_charge.erase(it); // delete the microcell from charge container
                }
                else
                    buffer.erase(buffer.begin()); // photon did not fire the microcell therefore the photon is deleted
            }
            else { // the microcell is not present therefore it must be created
                Microcell temp(buffer[0].index, overvoltage); // call a constructor of Microcell
                if (temp.discharge_init(buffer[0], arr_bins[i], this, buffer)){ // check whether photon fires the microcell
                    Microcell_discharge.push_back(temp); // and push the object to container
                }
                else
                    buffer.erase(buffer.begin()); // photon did not fire the microcell therefore the photon is deleted
            }
        }

        /* This block of code controls discharging of microcells */
        if (Microcell_discharge.size() != 0 || Microcell_charge.size() != 0){
            int idx = Microcell_discharge.size();
            sipm_par par = LUT[idx]; // find correct parameters of SiPM based on the number of the active microcells

            /* DISCHARGE */
            for (int j=0;j<Microcell_discharge.size();j++)
                Microcell_discharge[j].discharge(arr_bins[i], i, this, par); // discharge Microcell

            /* This block of code manages the correct operation of discharge containers - sorts container and pops Microcells which are discharged and moves them to charge container */
            sort(Microcell_discharge.begin(), Microcell_discharge.end(), [](const Microcell& obj1, const Microcell& obj2) {return obj1.discharged < obj2.discharged;}); // sort Microcell_discharge
            for (auto it =  Microcell_discharge.rbegin(); it != Microcell_discharge.rend(); ++it){ // backward iteration through discharge container
                if (it->discharged == 1){ // if Microcell is discharged
                    Microcell_charge.push_back(*it); // move it to charge container
                    Microcell_discharge.pop_back(); // pop it from discharge container
                }
            }

            /* CHARGE */
            for (int k=0;k<Microcell_charge.size();k++){
                Microcell_charge[k].charge(arr_bins[i], i, this, par); // charge Microcell
            }

            /* This block of code manages the correct operation of charge containers - sorts container and pops Microcells which are charged */
            sort(Microcell_charge.begin(), Microcell_charge.end(), [](const Microcell& obj1, const Microcell& obj2) {return obj1.charged < obj2.charged;}); // sort Microcell_charge
            for (auto it =  Microcell_charge.rbegin(); it != Microcell_charge.rend(); ++it){
                if (it->charged == 1){ // if Microcell is charged
                    Microcell_charge.pop_back(); // pop it from charge container
                }
            }
        }
        else
            continue;
    }
}

void SiPM::reset()
{
    Microcell_discharge = {};
    Microcell_charge = {};
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
    std::cout << "------------------------------------" << std::endl;
    std::cout << std::endl;
}

std::ostream& operator<<(std::ostream& os, const SiPM* en)
{
    for (int i=0;i<en->arr_bins.size();i++)
        os << en->arr_bins[i] << "," << en->arr_light[i] << "," << en->arr_afterpulse[i] << "," << en->arr_crosstalk[i] << "," << en->arr_dark[i] << std::endl;
    return os;
}

void SiPM::tally_waveform(std::string name)
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

void SiPM::waveform_simulation(std::string path, Scintillator crystal, std::vector<double> deposited_energy, int rep)
{
    for (int i=0;i<deposited_energy.size();i++)
    {
        std::cout << deposited_energy[i] << " MeV" << std::endl;
        std::string sc = crystal.get_name();
        std::string extension = crystal.get_extension();
        std::string filename = name + "_" + sc + "_" + std::to_string(deposited_energy[i]) + extension + ".csv";
        filename = path + "/" + filename;
        for (int j=0;j<rep;j++){
            std::cout << j << "/" << rep << std::endl;
            reset(); // reset each microcell of sipm before simulation
            std::vector<light> primaries = crystal.generate_light(deposited_energy[i]); // generate primaries, primaries is array of light
            if (dark_current_swich == true)
                dark_current(primaries); // add dark current to primaries
            map_light(primaries); // maps primary photons to individual microcells
            simulate(primaries); // output is array of tally
            tally_waveform(filename); // tally the matrix of pulse shape including total, afterpulse, crosstalk and dark current
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
