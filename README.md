# ZVS-Sim

**About the program**

This is a simple C++ implementation of the second order Taylor series for solving the circuit of a Mazzilli ZVS driver.
The simulated circuit is not the original ZVS driver but a model that uses two resistors instead of Mosfets. The value of these resistors vary in time like a pulse train with values from 0 to 100M (recommended value).

**Known bugs**

Sometimes the simulation will fail to converge due to the inductor that every Mazzilli drivers have next to the power source. If the value of this inductor is too low, the simulation will probably fail. If this happens, the program will retry the simulation with a bigger value for this inductor. 

It doesn't matter that you are using a different inductance, this inductance has very little effect on the voltages and currents in the resonant mesh.

**How to compile**

Compile using any compiler, there are no external libraries used.

**How to use**

The program will ask the user to enter these parameters:

-L1 inductance (this is the inductor next to the power source)

-L2 inductance: the inductance of half the primary (Mazzilli drivers usually have a 5+5 primary, here it's asked the inductance of those 5 turns)

-L4 inductance: the inductance of the secondary

-V: voltage of the power source

-C: capacitor used

-Time step: recommended 1e-9

-Total time to simulate: 0.2s is usually enough

-Mosfet slope value (how fast the Mosfets switch): I use 0.0001

-Mosfet max resistance: use 100e6

-Secondary resistance: the load that is connected to the secondary

-Points saved: use 100e3 and the program will save the last 100k points of data.

**What the program saves**

The program will save the data of the capacitor current (IC.dat), the primary current (IL2.dat), the current from the power source (ISource.dat), the voltage on the capacitor (VC.dat) and the voltage on secondary (Vsec.dat)

**Example of parameters**

L1=0.1

L2=3.45e-6

L4=20507e-6

V=48

C=136e-9

Time step=1e-9

Time to simulate=0.1

Mosfet slope=0.0001

Mosfet max resistance=100e6

Secondary resistance=200e3

And this will give you the simulation of a ZVS that has 3+3 turns around a 1.5mm gapped ETD59 ferrite core with
232 turns of secondary and a resonant working frequency of around 110kHz. The simulation gives similar values for primary and secondary voltage & current
simulatd with OrCAD.
