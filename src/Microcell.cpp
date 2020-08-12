#include "Microcell.h"
#include "SiPM.h"
#include "spline.h"
#include "Functions.h"

#include <cmath>
#include <string>
#include <iostream>
#include <random>
#include <ctime>
#include <cmath>

Microcell::Microcell(int i, double ov)
{
    id = i;
    V_ov = ov;
    last_trigger = 0;
    last_quench = 0;
    //ctor
}

Microcell::~Microcell()
{
    //dtor
}

bool Microcell::discharge_init(light in, double time)
{
    if (in.origin == "scintillator_light"){
        condition = sipm->geometry_factor * qe; // calculates PDE
        qe = calculate_QE(in.wavelenght, sipm);
    }
    else{
        condition = 1; // in case of dark current, afterpulses or optical crosstalk
    }

    double dice = get_rand_uni_dist();
    if (dice <= condition){
        last_trigger = time;
        return true;
    }
    else
        return false;
}

void Microcell::discharge_stop(double time)
{
    last_quench = time;
    discharge = false;
    charge = true;
}

void Microcell::load_discharge(double time, SiPM *sipm)
{
    if (last_quench == 0.0 && last_trigger == 0.0){ // in case there has been no hit of the microcell, do nothing
        //std::cout << "no hit" << std::endl;
        return;
    }
    //std::cout << "again" << std::endl;

    //std::cout << "load_discharge: " << this->id << " " << last_quench << " " << last_trigger << " " << time << std::endl;
    int index_trigger = get_nearest_index(sipm->arr_bins, last_trigger); // find the closest time of last trigger in arr_bins
    int index_quench = get_nearest_index(sipm->arr_bins, last_quench); // find the closest time of last quench in arr_bins
    int index_avalanche = get_nearest_index(sipm->arr_bins, time); // find the closest time of current trigger in arr_bins

    //std::cout << last_trigger << " " << last_quench << " " << time << std::endl;
    //std::cout << index_trigger << " " << index_quench << " " << index_avalanche << std::endl;

    double out = 0;

    for (int i=index_trigger; i<=index_quench; i++){ // discharge phase
        //std::cout << "i " << i << std::endl;
        //std::cout << "exp " << (1 - ((sipm->T_q - sipm->T_i) / (sipm->T_d - sipm->T_i)) * exp(-(sipm->arr_bins[i] - sipm->arr_bins[index_trigger]) / sipm->T_i) + ((sipm->T_q - sipm->T_d) / (sipm->T_d - sipm->T_i)) * exp(-(sipm->arr_bins[i] - sipm->arr_bins[index_trigger]) / sipm->T_d)) << std::endl;
        out = I_f * (1 - ((sipm->T_q - sipm->T_i) / (sipm->T_d - sipm->T_i)) * exp(-(sipm->arr_bins[i] - sipm->arr_bins[index_trigger]) / sipm->T_i) + ((sipm->T_q - sipm->T_d) / (sipm->T_d - sipm->T_i)) * exp(-(sipm->arr_bins[i] - sipm->arr_bins[index_trigger]) / sipm->T_d));

        if (last_origin == "scintillator_light")
            sipm->arr_light[i] += out;
        //std::cout << "origin " << last_origin << std::endl;
        //std::cout << "arr_light[" << i << "] " << sipm->arr_light[i] << std::endl;
        if (last_origin == "afterpulse"){
            //std::cout << "afterpulse " << out << std::endl;
            sipm->arr_afterpulse[i] += out;
        }
        if (last_origin == "dark_current"){
            //std::cout << "dark_current " << out << std::endl;
            sipm->arr_dark[i] += out;
        }
        if (last_origin == "crosstalk"){
            //std::cout << "crosstalk " << out << std::endl;
            sipm->arr_crosstalk[i] += out;
        }

    }

    //std::cout << "kacer1 " << std::endl;

    double y_T = out;
    double y_M = calculate_microcell_current(sipm->arr_bins[index_quench], sipm);
    double Q_cd = (V_ov - sipm->R_d * y_M) * sipm->C_d;
    double Q_cq = (sipm->R_l * y_T + sipm->R_d * y_M - V_ov) * sipm->C_q;
    double Q_ceq = (sipm->R_l * y_T) * sipm->C_eq;
    double T_q2 = (sipm->T_q * Q_cd + sipm->R_q * sipm->C_d * Q_cq + sipm->R_q * (sipm->C_d + sipm->C_q) * Q_ceq) / (Q_ceq + Q_cd);

    for (int i=index_quench+1; i<=index_avalanche; i++){ // recharge phase
        out = ((Q_ceq + Q_cd) / (sipm->T_d2 - sipm->T_i2)) * (((T_q2 - sipm->T_i2) / sipm->T_i2) * exp(-(sipm->arr_bins[i] - sipm->arr_bins[index_quench]) / sipm->T_i2) - ((T_q2 - sipm->T_d2) / sipm->T_d2) * exp(-(sipm->arr_bins[i] - sipm->arr_bins[index_quench]) / sipm->T_d2));
        if (std::isnan(out)){
            std::cout << "asdasdasd" << std::endl;
            abort();
        }

        if (last_origin == "scintillator_light")
            sipm->arr_light[i] += out;
        if (last_origin == "afterpulse"){
            //std::cout << out << std::endl;
            sipm->arr_afterpulse[i] += out;
        }

        if (last_origin == "dark_current")
            sipm->arr_dark[i] += out;
        if (last_origin == "crosstalk")
            sipm->arr_crosstalk[i] += out;
    }
    //std::cout << "kacer3 " << std::endl;
}

double Microcell::calculate_overvoltage(double time, SiPM *sipm)
{
    if (last_trigger == 0.0){
        return sipm->overvoltage;
    }
    else{
        return sipm->overvoltage * (1 - exp(-time / sipm->T_d2));
    }
}

void Microcell::calculate_I_f(SiPM *sipm)
{
    //std::cout << "V_ov " << V_ov << std::endl;
    I_f = V_ov / (sipm->R_q + sipm->R_d + sipm->R_l);
}

double Microcell::calculate_microcell_current(double time, SiPM *sipm)
{
    //std::cout << sipm->T_d2 << " " << sipm->T_i << " " << sipm->T_i2 << " " << sipm->T_d << " " << time << " " << I_f << std::endl;
    return I_f * (1 + ((((sipm->T_d2 - sipm->T_i) * (sipm->T_i2 - sipm->T_i)) / (sipm->T_i * (sipm->T_d - sipm->T_i))) * exp(-time / sipm->T_i) - (((sipm->T_d2 - sipm->T_d) * (sipm->T_i2 - sipm->T_d)) / (sipm->T_d * (sipm->T_d - sipm->T_i))) * exp(-time / sipm->T_d)));
}

void Microcell::discharge(light in, std::vector<light>& buffer, SiPM *sipm)
{
    double condition, ov, qe;
    //std::cout << "ahoj " << last_trigger << " " << last_quench << std::endl;

    ov = calculate_overvoltage(in.time, sipm);
    if (in.origin == "scintillator_light"){
        condition = sipm->geometry_factor * qe; // calculates PDE
        qe = calculate_QE(in.wavelenght, sipm);
    }
    else{
        //std::cout << "in.origin " << in.origin << std::endl;
        condition = 1; // in case of dark current, afterpulses or optical crosstalk
    }


    //std::cout << in.time << " " << last_trigger << " + " << last_quench << std::endl;
    if (in.time <= last_quench){ // case in which the light hits microcell during quenching the previous avalanche (no effect)
        //std::cout << "stop " << std::endl;
        condition = 0;
    }

    double T = 0;
    bool crosstalk = false;
    double dice = get_rand_uni_dist();
    //std::cout << "condition " << condition << std::endl;
    if (dice <= condition)
    {
        load_discharge(in.time, sipm);
        V_ov = ov;
        QE = qe;
        calculate_I_f(sipm); // calculate I_f
        int i = 0;
        while (1){ // numerically finds the quenching time T
            double t = i * sipm->timestep;
            double I = calculate_microcell_current(t, sipm);
            if (I < sipm->I_th){
                if (i == 0)
                    T = sipm->timestep;
                else{
                    double t_1 = (i-1) * sipm->timestep;
                    double I_1 = calculate_microcell_current(t_1, sipm);
                    T = linear_interpolation(I_1, I, t_1, t, sipm->I_th);
                }
                break;
            }
            i += 1;
        }
        //std::cout << "T " << T << std::endl;
        if (in.origin == "scintillator_light"){
            afterpulse(in.index, in.time, sipm,  buffer);
            optical_crosstalk(in.time, sipm, buffer);
        }

        last_trigger = in.time;
        last_quench = in.time + T;
        last_origin = in.origin;
    }
}

double Microcell::calculate_QE(double wavelenght, SiPM *sipm)
{
    std::vector<double> x, y;

    for (auto itr = sipm->absorption_spectrum.begin(); itr != sipm->absorption_spectrum.end(); ++itr) {
        //std::cout << itr->first << " " << itr->second.x.size() << std::endl;
        x.push_back(itr->first);

        spline QE_function;
        QE_function.set_points((itr->second).x, (itr->second).y);
        double prob = QE_function(wavelenght);
        //std::cout << "prob= " << prob << " " << wavelenght << std::endl;
        y.push_back(prob);
    }

    //std::cout << "V_ov= " << V_ov << std::endl;
    //std::cout << x[0] << " " << x[1] << " " << y[0] << " " << y[1] << " " << V_ov << std::endl;
    if (x.size() == 2)
        QE = linear_interpolation(x[0], x[1], y[0], y[1], V_ov);
    else if (x.size() > 2){
        spline QE_overvoltage;
        QE_overvoltage.set_points(x, y);
        QE = QE_overvoltage(V_ov);
    }
    //std::cout << "QE " << QE << std::endl;

    if (QE > 1.0 || QE < 0.0)
        std::cout << "Error: QE = " << QE << std::endl;
    return QE;
}

void Microcell::calculate_crosstalk_probability(SiPM *sipm)
{
    //std::cout << "QE a crosstalk probability1 " << QE << std::endl;
    crosstalk_probability = sipm->crosstalk_function(QE);
    //std::cout << "QE a crosstalk probability2 " << QE << " " << crosstalk_probability << std::endl;
}

void Microcell::calculate_afterpulse_probability(SiPM *sipm)
{
    afterpulse_probability = sipm->afterpulse_function(V_ov);
}

void Microcell::afterpulse(int id, double time, SiPM *sipm, std::vector<light>& buffer)
{
    calculate_afterpulse_probability(sipm);
    int lenght = get_rand_poiss_dist(afterpulse_probability);

    if (lenght > 0){
        for (int i=0;i<lenght;i++){
            double dice = get_rand_uni_dist();
            double T = bisection_method(0.0, 0.001, dice, 1e-15, sipm->afterpulse_amplitude, sipm->afterpulse_components, decay_function);
            //std::cout << lenght << " " << id << " " << T << " " << "afterpulse" << " " << std::endl;
            add_pulse(id, time + T, "afterpulse", buffer);
        }
    }
}

void Microcell::optical_crosstalk(double time, SiPM *sipm, std::vector<light>& buffer)
{
    double dice = get_rand_uni_dist();
    calculate_crosstalk_probability(sipm);
    int lenght = get_rand_poiss_dist(crosstalk_probability);
    for (int i=0;i<lenght;i++){
        int id = rand() % sipm->N_c;
        //std::cout << lenght << " " << id << " " << time + 1e-09 << " " << "optical_crosstalk" << " " << std::endl;
        add_pulse(id, time + 1e-09, "crosstalk", buffer);
    }
}

void Microcell::reset(SiPM *sipm)
{
    last_trigger = 0;
    last_quench = 0;
    V_ov = sipm->overvoltage;
}
int Microcell::get_ID()
{
    return id;
}

bool Microcell::operator<(const Microcell& right)const
{
      return (this->id < right.id);
}

bool Microcell::operator==(const Microcell& right)const
{
      return (this->id == right.id);
}
