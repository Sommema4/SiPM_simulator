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
    //ctor
}

Microcell::~Microcell()
{
    //dtor
}

bool Microcell::discharge_init(light in, double time, SiPM *sipm, std::vector<light>& buffer)
{
    if (in.origin == "scintillator_light"){
        PDE = sipm->geometry_factor * calculate_QE(in.wavelenght, sipm); // calculates PDE
        double dice = get_rand_uni_dist(); // Alea iacta est
        if (dice <= PDE){
            calculate_I_f(sipm);
            last_trigger = time;
            last_origin = in.origin;
            if (sipm->afterpulse_swich == true)
                afterpulse(id, time, sipm, buffer);
            if (sipm->crosstalk_swich == true)
                optical_crosstalk(time, sipm, buffer);
            return true;
        }
    }
    else{
        calculate_I_f(sipm);
        last_trigger = time;
        last_origin = in.origin;
        return true;
    }
    return false;
}

void Microcell::calculate_I_f(SiPM *sipm)
{
    //std::cout << "V_ov " << V_ov << std::endl;
    I_f = V_ov / (sipm->R_q + sipm->R_d + sipm->R_l);
}

double Microcell::calculate_overvoltage(double time, SiPM *sipm, sipm_par par)
{
    //std::cout << "par.T_d2: " << par.T_d2 << std::endl;
    return sipm->overvoltage * (1 - exp(-(time) / par.T_d2));
}

double Microcell::calculate_QE(double wavelenght, SiPM *sipm)
{
    std::vector<double> x, y;
    for (auto itr = sipm->absorption_spectrum.begin(); itr != sipm->absorption_spectrum.end(); ++itr) {
        x.push_back(itr->first);

        spline QE_function;
        QE_function.set_points((itr->second).x, (itr->second).y);
        double prob = QE_function(wavelenght);
        y.push_back(prob);
    }

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

void Microcell::afterpulse(int id, double time, SiPM *sipm, std::vector<light>& buffer)
{
    calculate_afterpulse_probability(sipm);
    int lenght = get_rand_poiss_dist(afterpulse_probability);

    if (lenght > 0){
        for (int i=0;i<lenght;i++){
            double dice = get_rand_uni_dist();
            double T = bisection_method(0.0, 0.001, dice, 1e-15, sipm->afterpulse_amplitude, sipm->afterpulse_components, decay_function);
            //std::cout << "afterpulse added" << std::endl;
            add_pulse(id, time + T, "afterpulse", buffer);
        }
    }
}

void Microcell::calculate_afterpulse_probability(SiPM *sipm)
{
    afterpulse_probability = sipm->afterpulse_function(V_ov);
}

void Microcell::optical_crosstalk(double time, SiPM *sipm, std::vector<light>& buffer)
{
    double dice = get_rand_uni_dist();
    calculate_crosstalk_probability(sipm);
    int lenght = get_rand_poiss_dist(crosstalk_probability);
    for (int i=0;i<lenght;i++){
        int id = rand() % sipm->N_c;
        //std::cout << "crosstalk added" << std::endl;
        add_pulse(id, time + 1e-09, "crosstalk", buffer);
    }
}

void Microcell::calculate_crosstalk_probability(SiPM *sipm)
{
    crosstalk_probability = sipm->crosstalk_function(QE);
}

double Microcell::calculate_microcell_current(double time, sipm_par par)
{
    return I_f * (1 + ((((par.T_d2 - par.T_i) * (par.T_i2 - par.T_i)) / (par.T_i * (par.T_d - par.T_i))) * exp(-time / par.T_i) - (((par.T_d2 - par.T_d) * (par.T_i2 - par.T_d)) / (par.T_d * (par.T_d - par.T_i))) * exp(-time / par.T_d)));
}

double Microcell::discharge(double time, int index, SiPM *sipm, sipm_par par)
{
    charged = 0;
    double rel_time = time - last_trigger;
    double out = I_f * (1 - ((sipm->T_q - par.T_i) / (par.T_d - par.T_i)) * exp(-(rel_time) / par.T_i) + ((sipm->T_q - par.T_d) / (par.T_d - par.T_i)) * exp(-(rel_time) / par.T_d));

    if (last_origin == "scintillator_light")
        sipm->arr_light[index] += out;
    if (last_origin == "afterpulse")
        sipm->arr_afterpulse[index] += out;
    if (last_origin == "dark_current")
        sipm->arr_dark[index] += out;
    if (last_origin == "crosstalk")
          sipm->arr_crosstalk[index] += out;

    double I = calculate_microcell_current(rel_time, par);
    if (I < sipm->I_th){
        last_current = out;
        last_quench = time;
        y_M = calculate_microcell_current(last_quench, par);
        Q_cd = (V_ov - sipm->R_d * y_M) * sipm->C_d;
        Q_cq = (sipm->R_l * last_current + sipm->R_d * y_M - V_ov) * sipm->C_q;
        discharged = 1;
    }
    return out;
}

double Microcell::charge(double time, int index, SiPM *sipm, sipm_par par)
{
    discharged = 0;
    double rel_time = time - last_quench;
    Q_ceq = (sipm->R_l * last_current) * par.C_eq;
    T_q2 = (sipm->T_q * Q_cd + sipm->R_q * sipm->C_d * Q_cq + sipm->R_q * (sipm->C_d + sipm->C_q) * Q_ceq) / (Q_ceq + Q_cd);
    double out = ((Q_ceq + Q_cd) / (par.T_d2 - par.T_i2)) * (((T_q2 - par.T_i2) / par.T_i2) * exp(-(rel_time) / par.T_i2) - ((T_q2 - par.T_d2) / par.T_d2) * exp(-(rel_time) / par.T_d2));

    if (last_origin == "scintillator_light")
        sipm->arr_light[index] += out;
    if (last_origin == "afterpulse")
        sipm->arr_afterpulse[index] += out;
    if (last_origin == "dark_current")
        sipm->arr_dark[index] += out;
    if (last_origin == "crosstalk")
        sipm->arr_crosstalk[index] += out;


    V_ov = calculate_overvoltage(rel_time, sipm, par);
    if (V_ov >= 0.97 * sipm->overvoltage){
        V_ov = sipm->overvoltage;
        charged = 1;
    }
    return out;
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
