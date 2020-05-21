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

        void discharge(light, std::vector<light>&, class SiPM *); // controls discharge of microcell
        void load_discharge(double, SiPM *); // manages the saving data to tallies
        double calculate_overvoltage(double, SiPM *); // calculates overvoltage
        double calculate_QE(double, SiPM *); // calculates quantum efficiency
        double calculate_microcell_current(double time, SiPM *); // calculates current going through microcell
        void calculate_I_f(SiPM *);

        void afterpulse(int, double, SiPM *, std::vector<light>&); // manages the afterpulse including adding new pulses to the buffer
        void calculate_afterpulse_probability(SiPM *); // calculates the overall afterpulse probability
        double afterpulse_function(double, double); // calculates the time decay function

        void optical_crosstalk(double, SiPM *, std::vector<light>&); // manages the optical crosstalk including adding new pulses to the buffer
        void calculate_crosstalk_probability(SiPM *); // calculates optical crosstalk prabability

        void reset(SiPM *);
    protected:

    private:
        int id;
        double V_ov;
        double last_trigger;
        double last_quench;
        std::string last_origin;
        double QE;
        double I_f;
        double afterpulse_probability;
        double crosstalk_probability;
};

#endif // MICROCELL_H
