#ifndef MICROCELL_H
#define MICROCELL_H

#include "Functions.h"
#include "SiPM.h"

class Microcell
{
    public:
        friend class SiPM;
        Microcell(int, double);
        virtual ~Microcell();

        bool discharge_init(light, double, class SiPM *, std::vector<light>&); // flags microcell for discharge, decides whether it undergoes afterpulse or optical crosstalk

        double discharge(double, int, SiPM *, sipm_par); // controls discharge of microcell
        double charge(double, int, SiPM *, sipm_par); // manages the saving data to tallies

        double calculate_overvoltage(double, SiPM *, sipm_par); // calculates overvoltage
        double calculate_QE(double, SiPM *); // calculates quantum efficiency
        double calculate_microcell_current(double time, sipm_par); // calculates current going through microcell
        void calculate_I_f(SiPM *);

        void afterpulse(int, double, SiPM *, std::vector<light>&); // manages the afterpulse including adding new pulses to the buffer
        void calculate_afterpulse_probability(SiPM *); // calculates the overall afterpulse probability

        void optical_crosstalk(double, SiPM *, std::vector<light>&); // manages the optical crosstalk including adding new pulses to the buffer
        void calculate_crosstalk_probability(SiPM *); // calculates optical crosstalk prabability

        int get_ID(); // returns id of the object
        bool operator<(const Microcell&)const;
        bool operator==(const Microcell&)const;
    protected:

    private:
        int id;
        double V_ov;
        double last_trigger;
        double last_quench;
        double last_current;
        std::string last_origin;
        int discharged = 0;
        int charged = 0;
        double QE;
        double PDE;
        double I_f;
        double afterpulse_probability;
        double crosstalk_probability;
        double y_T;
        double y_M;
        double Q_cd;
        double Q_cq;
        double Q_ceq;
        double T_q2;
};

#endif // MICROCELL_H
