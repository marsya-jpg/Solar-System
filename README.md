# Solar-System
Simulation of a solar system with user selectable simulation parameters, incorporating the velocity Verlet algorithm to write an XYZ compatible trajectory file for the simulation.

Input textfiles:
1. p_setup.txt - Initial conditions of solar bodies
                (mass (kg), position (AU), velocity(AU/day))
2. sim_parameter.txt - Parameters of solar system 
                (G (AU3/kg Msun2), dt (day), numstep)
                
Output file:
1. track.dat - Positions of solar bodies in XYZ file format

*Only one satellite/moon (of planet 3)is allowed and written at the end of p_setup.txt
