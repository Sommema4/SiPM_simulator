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
    //std::cout << "It is gonna crash 3" << std::endl;
    if (in.origin == "scintillator_light"){
        //std::cout << "It is gonna crash 4" << std::endl;
        PDE = sipm->geometry_factor * calculate_QE(in.wavelenght, sipm); // calculates PDE
        //std::cout << "It is gonna crash 5" << std::endl;
        double dice = get_rand_uni_dist(); // Alea iacta est
        //std::cout << "dice: " << dice << ", PDE: " << PDE << std::endl;
        if (dice <= PDE){
            //std::cout << "It is gonna crash 6" << std::endl;
            calculate_I_f(sipm);
            last_trigger = time;
            last_origin = in.origin;
            //std::cout << "It is gonna crash 7" << std::endl;
            //afterpulse(id, time, sipm, buffer);
            //std::cout << "It is gonna crash 8" << std::endl;
            //optical_crosstalk(time, sipm, buffer);
            //std::cout << "It is gonna crash 9" << std::endl;
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
    //std::cout << "out: " << I_f << ", " << sipm->T_q << ", " << par.T_i << ", " << par.T_d << std::endl;
    double out = I_f * (1 - ((sipm->T_q - par.T_i) / (par.T_d - par.T_i)) * exp(-(rel_time) / par.T_i) + ((sipm->T_q - par.T_d) / (par.T_d - par.T_i)) * exp(-(rel_time) / par.T_d));
    //std::cout << "time= " << time << ", discharge out = " << out << std::endl;
    if (last_origin == "scintillator_light")
      sipm->arr_light[index] += out;
    if (last_origin == "afterpulse")
        sipm->arr_afterpulse[index] += out;
    if (last_origin == "dark_current")
        sipm->arr_dark[index] += out;
    if (last_origin == "crosstalk")
          sipm->arr_crosstalk[index] += out;

    std::cout << sipm->arr_light[index] << std::endl;

    double I = calculate_microcell_current(rel_time, par);
    //std::cout << "I = " << I << std::endl;
    if (I < sipm->I_th){
        //std::cout << "I_th reached" << std::endl;
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
    std::cout << "time= " << time << ", charge out = " << out << std::endl;

    if (last_origin == "scintillator_light"){
        sipm->arr_light[index] += out;
    }
    if (last_origin == "afterpulse")
        sipm->arr_afterpulse[index] += out;
    if (last_origin == "dark_current")
        sipm->arr_dark[index] += out;
    if (last_origin == "crosstalk")
        sipm->arr_crosstalk[index] += out;

    std::cout << sipm->arr_light[index] << std::endl;

    V_ov = calculate_overvoltage(rel_time, sipm, par);
    //std::cout << "V_ov: " << V_ov << std::endl;
    if (V_ov >= 0.97 * sipm->overvoltage){
        std::cout << "que?" << std::endl;
        V_ov = sipm->overvoltage;
        charged = 1;
        //return time - last_quench;
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
