"""
VELOCITY VERLET: velocity Verlet time integration of
a body moving in orbit

Produces plots of the position of the particle
and its energy, both as function of time. Also
saves both to file.

The potential is V(x) = a*x^4 - b*x^2, where
a and b are hard-coded in the main() method
and passed to the functions that
calculate force and potential energy.
"""

import sys
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3Dsos import Particle 
from scipy.signal import find_peaks, peak_prominences

def get_separation(p_list, N):
    """
    Method that calculates all necessary separations to return as an array
    
    :param p_list : list
    :return: N x N array
    """
    mod_array = np.zeros([N,N])
    separations = np.zeros([N,N,3])
    
    for i in range(N):
        for j in range(N):
            separations[i,j] = p_list[i].pos - p_list[j].pos
            mod_array[i,j] = np.linalg.norm(separations[i,j])
            
    return mod_array, separations    
    

def get_force(p_list, G, N):
    """
    Returns array of 3D vector of force from all pair interactions
    """

    sys_force = np.zeros([N,3])
    #Addresses tuple from get_separation
    sep, sep_vector = get_separation(p_list, N)
    force = np.zeros([N,N,3])
    
    for i in range(N):
        for j in range(N):
            #Gets the force numerator
            if i!=j:      
                force[i,j] = -G*p_list[i].mass*p_list[j].mass/(sep[i,j])**3 * sep_vector[i,j]
                
            else:
                continue
            
    sys_force = force.sum(axis=1)
    
    return sys_force

def update_vel(p_list, force, dt, N):
    """
    Velocity updater for a list of P3Ds
    """
    
    sys_velocity = np.zeros([N,3])
    for i in range(N):
        f = force[i]
        sys_velocity[i] = p_list[i].leap_vel(dt, f)
        
    return sys_velocity    
        
def update_pos(p_list, force, dt, N):
    """
    2nd order position updater for a list of P3Ds
    Uses dynamic methods for velocity and position as the velocity and ]
    position changes
    """
    sys_position = np.zeros([N,3])
    #Opens [N,N] array to hold pair separation values in scalar
    for i in range(N):
        f = force[i]
        sys_position[i] = p_list[i].leap_pos2nd(dt, f)
        
    
    return sys_position
    
def get_pos(p_list, N):
    """
    Get p_list[i].pos in an array
    """
    sys_position = np.zeros([N,3])
    sys_mod = np.zeros([N])
    for i in range(N):
        sys_position[i] = p_list[i].pos
        sys_mod[i] = np.linalg.norm(p_list[i].pos)
        
    return sys_position, sys_mod
    
def sys_pot_energy(p_list, G, dt, N):
    """
    Returns potential energy of system as a single value {scalar}
    Sum of potential from every interaction between bodies
    """
    potential = np.zeros([N,N])
    sep = get_separation(p_list, N)[0]
    for i in range(N):
        for j in range(N):
            if i!=j:
                potential[i,j] = -G*p_list[i].mass*p_list[j].mass/sep[i,j]
            else:
                continue
            
    sys_potential= np.triu(potential).sum()
    return sys_potential

     

# Begin main code
def main():
    # Read name of output file from command line
    if len(sys.argv)!=2:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file>")
        sys.exit()
    else:
        outfile_name = sys.argv[1]

    # Open output file
    outfile = open(outfile_name, "w")
    #Opens o2_vibration in symplectic and verlet
    
    #Reads in simulation parameters 
    with open("sim_parameter.txt","r") as my_input:
        #Reads first line only
        solar_data = my_input.readline().split(" ")
        G = float(solar_data[0])
        dt = float(solar_data[1])
        numstep = int(solar_data[2])        
        
    #Set up celestial bodies from input file
    f = open("new_setup.txt")
    p_list = []
    
    
   #Assigns new particle from input particle file
    count = sum(1 for line in open('new_setup.txt'))
    print(count)
    for i in range(count):
        val = f.readline()
        part = Particle.new_particle(val)
        p_list.append(part)      
   
    f.close()
    
    N = len(p_list)

    #COM  correction
    total_mass, com_vel = Particle.com_velocity(p_list)
    for i in range(N):
        p_list[i].vel = p_list[i].vel-com_vel
    
    #Setup initial condition
    time = 0

    # Write out initial conditions
    energy = sys_pot_energy(p_list, G, dt, N) + Particle.sys_kinetic(p_list)
    force = get_force(p_list, G, N)
    KE = Particle.sys_kinetic(p_list)
    time_list = [time]
    sys_pos, sys_mod = get_pos(p_list, N)
        
    # Initialise observables initial conditions and data list
    mod_array, seps = get_separation(p_list, N)
    planet_sep = np.zeros([numstep, N-2])
    #remove sun and moon
    aphelion = []
    perihelion = []
    
    moon_sep = np.zeros([numstep])
    apogee = []
    perigee = []
    
    energy_list = [energy]
    ke_list = [KE]
    force_list = [np.linalg.norm(force[3])]
    
    #Start write XYZ file 
    #Makes a list of particle names
    particle_label = [p_list[i].label for i in range(len(p_list))]
        
    def xyz(p_list, num):
        outfile.write("%s\n"%len(p_list))
        outfile.write("%s\n"%num)
        for x, particle in zip(sys_pos, particle_label):
            outfile.write("%s %.18g %.18g %.18g\n"%(particle, x[0], x[1], x[2]))      
    
    xyz(p_list, 0)
    
    
    # Start the time integration loop

    for num in range(numstep):
        # Update particle position
        update_pos(p_list, force, dt, N)
        
        # Update force
        force_new = get_force(p_list, G, N)
        
        # Update particle velocity by averaging
        # current and new forces
        update_vel(p_list, 0.5*(force+force_new), dt, N)
        
        sys_pos, sys_mod = get_pos(p_list, N)
        
        #Assign separations for apsis and periapsis
        mod_array, seps = get_separation(p_list, N)
        planet_sep[num] = mod_array[1:-1,0]
        moon_sep[num] = mod_array[3,2]
        #index of moon: earth
        
        # Re-define force value
        force = force_new

        KE = Particle.sys_kinetic(p_list)

        # Output particle information
        energy = sys_pot_energy(p_list, G, dt, N) + Particle.sys_kinetic(p_list)
        
        # Increase time
        time += dt

                                                                     
        #Writing XYZ file (track.dat)
        xyz(p_list, num+1)
        
        # Append information to data lists
        time_list.append(time)
        #pos_list[num+1] = np.linalg.norm(get_pos(p_list, N)[i])
                            
        
        energy_list.append(energy)
        ke_list.append(KE)
        
        force_list.append(np.linalg.norm(force[3]))
    #print("ps", planet_sep.shape)   (366,10)

    # Post-simulation:
    # Close output file
    outfile.close()
    
    aphelion.append(np.amax(planet_sep,axis=0))
    perihelion.append(np.amin(planet_sep, axis=0))
    print("aphelion of planets",aphelion)
    print("perihelion of planets", perihelion)
    
    apogee.append(np.amax(moon_sep))
    perigee.append(np.amin(moon_sep))
    print("apogee of moon", apogee)
    print("perigee of moon", perigee)  
    
    
    #Calculation of orbital period
    #Finding peaks
    planet_period = []
    print(planet_sep.shape)
    #Range of 0 to len(p_list)-1
    for i in range(0,2):
    #index 0: index last planet+1    
        peaks = find_peaks(planet_sep[:,i], height=0)  #Peaks of position
        peak_time = [time_list[i] for i in peaks[0]] 
        #print(peak_time)
        val = peak_time[1]-peak_time[0]
        
        planet_period.append(val)
    
    print("period",planet_period)
    
    #peaks = find_peaks(moon_sep, height=0)  #Peaks of position
    #peak_time = [time_list[i] for i in peaks[0]] 
    #moon_period = peak_time[1]-peak_time[0]
    
    #print("moon period",moon_period)
    
    # Plot particle energy to screen
    #pyplot.title('Velocity Verlet: total energy vs time')
    pyplot.title('Moon: Force vs time')
    #pyplot.xlim([0,40])
    pyplot.xlabel('Time')
    pyplot.ylabel('Force')
    pyplot.plot(time_list, force_list)
    
    pyplot.show()


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()

