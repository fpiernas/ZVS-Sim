/**---------------------------------------------------------------------------------------
ZVS Mazzilli Driver Simulator Program V1.1
@author Piernas Diaz, Fran
University of Granada, Spain

MIT License

Copyright (c) 2019 Francisco Piernas Diaz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Changelog 1.1:
The program now simulates the circuit as in V1.0, but if a convergence error happens, it will
retry the simulation with a bigger L1 inductor or smaller time step.

***The program:
This program simulates a ZVS Mazzilli Driver in time domain, it will help you know the currents
and voltages involved when the circuit reaches a steady state, if you know the parameters of your circuit.

The program solves the mesh currents I1 (mesh current across the voltage source and half of the primary,
I2 (same as I1 but across the other half of the primary), I3 (mesh current across the capacitor) and I4
(current in secondary) using a simple second order Euler algorithm.

Important: to avoid current convergence errors, the inductance of L1 should be large, like 0.1, it does not
matter if your circuit uses a different value for this inductor.
The value of this inductor has no effect on the currents or voltages when in steady state, it will
only make the simulation slower on reaching the steady state, you can easily compensate this by simulating
more time.

***Parameters:
L1 is the inductance of the inductor placed after the voltage source (if convergence error happens, the program will
   retry the simulation with a bigger L1 until the simulation gives no errors)
L2 is the inductance of half of the primary. Usually Mazzilli drivers have a 5+5 turns primary,
   here it is asked the inductance of those 5 turns. If you had a 3+3 primary, then L2 should be the
   inductance of 3 turns.
L3 is not asked, since it has the same value as L2
L4 is the inductance of the secondary winding
V is the voltage of the voltage source
C is the value of the capacitor or capacitor bank
delta_t is the time step of the simulation (usually 1e-9 is enough)
t_total is the total simulation time (usually 0.1 sec is enough to reach a steady state)
R is the maximum value of the resistor that pretends to be a Mosfet
slope_R is the slope of the turning on and off process of the resistors pretending to be Mosfets
        (recommended value = 0.0001)
R_Sec is the load, attached to the secondary
last_perc is the last percentage of data you want to save (usually just 1% is enough)
last_points is the number of points to save, equal to t_total*(1.0-last_perc/100.0)/delta_t

-------------------------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

class classR //This is a variable resistor class with a value periodic in time, it simulates to be a Mosfet
{
    public:
        long double p1, p2, p3, p4, T, max_value; //p1,2,3 and 4 are: Mosfet turning off, Mosfet off, Mosfet turning on and Mosfet on zones, max_value is the maximum value of the resistor that pretends to be a Mosfet

        void configure(long double T_,
                       long double slope_,
                       long double value_) //slope is the slope of the Mosfet or how fast it turns on and off, T_ is the resonant period and value_ is tha maximum vlue of the resistor that pretends to be the Mosfet
        {
            slope_/=100.0;
            p1=T_*slope_;
            p2=T_/2.0-T_*slope_/2.0;
            p3=T_/2.0+T_*slope_/2.0;
            p4=T_;
            T=T_;
            max_value=value_;
            return;
        }

        long double value(long double t) //t is the value of time that is being simulated
        {
            long double int_part;
            t=1.0*T*modf(1.0*t/T,&int_part);
            if(t<=p1) return max_value*t/p1; //Mosfet turning off zone
            if(t>p1&&t<=p2) return max_value; //Mosfet turned off zone
            if(t>p2&&t<=p3) return max_value*(1.0-((t-p2)/(p3-p2))); //Mosfet turning on
            if(t>p3) return 0; //Mosfet turned on
            else return 0;
        }
};

class zvs
{
    public:
        long double L1, L2, V, C, L4, delta_t, t_total, slope_R, T, f, R, R_Sec, last_perc;
        classR resist;
        void configure(long double L1_,
                       long double L2_,
                       long double L4_,
                       long double C_,
                       long double V_,
                       long double R_,
                       long double delta_t_,
                       long double t_total_,
                       long double slope_,
                       long double R_Sec_,
                       long double last_perc_)
        {
            L1=L1_;
            L2=L2_;
            L4=L4_;
            C=C_;
            V=V_;
            R=R_;
            delta_t=delta_t_;
            t_total=t_total_;
            slope_R=slope_;
            T=2.0*M_PI*sqrt(4.0*L2*C);
            f=1.0/T;
            R_Sec=R_Sec_;
            resist.configure(T,slope_R,R);
            last_perc=last_perc_;
            return;
        }

        void save_parameters(void)
        {
            ofstream parameters;
            parameters.open("parameters.dat");
            parameters<<"L1:                    "<<L1<<endl;
            parameters<<"L2:                    "<<L2<<endl;
            parameters<<"L4:                    "<<L4<<endl;
            parameters<<"V:                     "<<V<<endl;
            parameters<<"C:                     "<<C<<endl;
            parameters<<"Time step:             "<<delta_t<<endl;
            parameters<<"Total simulation time: "<<t_total<<endl;
            parameters<<"Resistor slope:        "<<slope_R<<endl;
            parameters<<"Max resistance:        "<<R<<endl;
            parameters<<"Secondary load:        "<<R_Sec<<endl;
            return;
        }

        bool simulate(void)
        {
            ofstream Vsec_dat, VC_dat, IL2_dat, ISource_dat, IC_dat;
            Vsec_dat.open("Vsec.dat");
            VC_dat.open("VC.dat");
            IL2_dat.open("IL2.dat");
            ISource_dat.open("ISource.dat");
            IC_dat.open("IC.dat");

            Vsec_dat.precision(15);
            VC_dat.precision(15);
            IL2_dat.precision(15);
            ISource_dat.precision(15);
            IC_dat.precision(15);

            long double I1=0, I2=0, I3=0, I4=0; //mesh currents
            long double I1d=0, I2d=0, I3d=0, I4d=0; //the d means time derivative value, at time t
            long double I1da=0, I2da=0, I3da=0, I4da=0; //the da means time derivative value at time t - delta_t
            long double int_I3=0; //integral of I3
            long double M24, R1, R2, t=0;
            long double progress=0;
            long double progress_int=0, progress_int_a=0;
            bool recording=false, convergence_error=false;
            M24=sqrt(L2*L4);


            for(t=0;t<t_total;t+=delta_t)
            {
                //Beginning of algorithm
                R1=resist.value(t);
                R2=R-R1;

                I1d=(V+L1*I2d-L2*I2d+2.0*L2*I3d+M24*I4d-R2*I1)/(L1+L2);
                I1+=I1d*delta_t+0.5*delta_t*(I1d-I1da);
                I1da=I1d;

                I2d=-1.0*(V+R1*I2+L2*(I1d-2.0*I3d)-L1*I1d-M24*I4d)/(L2+L1);
                I2+=I2d*delta_t+0.5*delta_t*(I2d-I2da);
                I2da=I2d;

                I3d=-1.0*(int_I3/C-2.0*L2*(I1d+I2d)+2.0*M24*I4d)/(4.0*L2);
                I3+=I3d*delta_t+0.5*delta_t*(I3d-I3da);
                I3da=I3d;

                I4d=-1.0*(R_Sec*I4+M24*(2.0*I3d-I1d-I2d))/(L4);
                I4+=I4d*delta_t+0.5*delta_t*(I4d-I4da);
                I4da=I4d;

                int_I3+=I3*delta_t;
                //End of algorithm

                //Check if it is time to start recording data
                if(t>t_total*(1.0-last_perc/100.0))
                {
                    if(!recording) {cout<<"Start recording data."<<endl; recording=true;}
                    Vsec_dat<<t<<" "<<I4*R_Sec<<endl;
                    VC_dat<<t<<" "<<int_I3/C<<endl;
                    IL2_dat<<t<<" "<<I3-I2<<endl;
                    ISource_dat<<t<<" "<<I1-I2<<endl;
                    IC_dat<<t<<" "<<I3<<endl;
                }

                //Check the progress and check if current is infinite->error and finish simulation
                progress=100.0*t/t_total;
                modf(progress,&progress_int);
                if(progress_int>progress_int_a) cout<<progress_int<<"%"<<endl;
                if(!isfinite(I1)||abs(I1)>1e10)
                {
                    convergence_error=true;
                    break;
                }
                progress_int_a=progress_int;
            }
            return convergence_error;
        }

        void read_parameters(void)
        {
            long double L1_input, L2_input, V_input, C_input, L4_input, delta_t_input, t_total_input, slopeR_input, R_input, R_Sec_input, last_points;
            cout<<"Set L1 value (recommended value = 0.1): ";
            cin>>L1_input;
            cout<<"Set L2 value: ";
            cin>>L2_input;
            cout<<"Set L4 value: ";
            cin>>L4_input;
            cout<<"Set V value: ";
            cin>>V_input;
            cout<<"Set C value: ";
            cin>>C_input;
            cout<<"Set time step value (recommended value = 1e-9): ";
            cin>>delta_t_input;
            cout<<"Set total simulation time value (usually 0.1 seconds is enough): ";
            cin>>t_total_input;
            cout<<"Set Mosfet slope value (recommended value = 0.0001): ";
            cin>>slopeR_input;
            cout<<"Set Mosfet max resistance value (recommended value = 100e6): ";
            cin>>R_input;
            cout<<"Set Secondary resistance value: ";
            cin>>R_Sec_input;
            cout<<"Set last number of points of data saved (100e3 points to plot is good): ";
            cin>>last_points;
            configure(L1_input,
                      L2_input,
                      L4_input,
                      C_input,
                      V_input,
                      R_input,
                      delta_t_input,
                      t_total_input,
                      slopeR_input,
                      R_Sec_input,
                      100.0*last_points/(t_total_input/delta_t_input));
            return;
        }
};


int main(void)
{
    zvs circuit1;
    cout<<"ZVS Simulator Program."<<endl;
    circuit1.read_parameters();
    /*
    read_parameters will ask the user to write the parameters
    You can alternatively call circuit1.configure(L1,L2,L4,C,V,R,delta_t,t_total,slopeR,R_Sec,last_perc)
    instead, if you already have the values of the parameters
    */
    while(circuit1.simulate()) //while there is an error
    {
        if(circuit1.L1<1.0)
        {
            circuit1.L1*=2.0;
            cout<<"Convergence error, readjusting L1 to "<<circuit1.L1<<endl;
        }
        if(circuit1.L1>=1.0&&circuit1.delta_t<=100e-9)
        {
            circuit1.L1+=2.0;
            cout<<"Convergence error, readjusting L1 to "<<circuit1.L1<<endl;
        }
        if(circuit1.L1>=1.0&&circuit1.delta_t>100e-9)
        {
            circuit1.delta_t/=2.0;
            cout<<"Convergence error, readjusting L1 to "<<circuit1.L1<<" and time step to "<<circuit1.delta_t<<endl;
        }
        if(circuit1.L1>20.0&&circuit1.delta_t<0.01e-9)
        {
            cout<<"Convergence error could not be solved."<<endl;
            break;
        }
    }
    circuit1.save_parameters();
    return 0;
}
