                    """
PARTICLE3D: Class to describe N bodies in solar system 

 An instance describes a particle in Euclidean 3D space: 
 velocity and position are [3] arrays

 Includes time integrator methods +...

author: Marsya Irdina Mohammad Karim s1938504

"""
import sys
import math
import numpy as np
    

class Particle(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label: name of the body
    mass: mass of the body
    pos: position of the body
    vel: velocity of the body

        Methods:
    __init__
    __str__
    kinetic_e  - computes the kinetic energy
    momentum - computes the linear momentum
    leap_pos2nd - updates the position to 2nd order
    leap_vel - updates the velocity

        Static Methods:
    new_particle - initializes a P3D instance from a file handle
    sys_kinetic - computes total K.E. of a p3d list
    com_velocity - computes total mass and CoM velocity of a p3d list
    """
    def __init__(self, label, mass, pos, vel):
        """
        Initialises a particle in 3D space

        :param label: String w/ the name of the particle
        :param mass: float, mass of the particle
        :param position: [3] float array w/ position
        :param velocity: [3] float array w/ velocity
        """
        
        self.label = label
        self.mass = mass
        self.vel = vel
        self.pos = pos
        

    @staticmethod
    def new_particle(line):
        #Splits every element in the line into list called data
        data = line.split(" ")
        label = str(data[0])
        mass = float(data[1])
        pos = np.array(data[2:5],float)
        vel = np.array(data[5:],float)
                  
        return Particle(label, mass, pos, vel)
       
    def __str__(self):
        """
        Define output format.
        Prints out mass, D, re, alpha
        For particle p=(2.0, 0.5, 1.0) this will print as
        "x = 2.0, v = 0.5, m = 1.0"
        """
        
        return "x = " + str(self.pos) + ", v = " + str(self.vel) + ", m = " + str(self.mass)
    
    #INSTANCE
    def KE(self):
        """
        Returns the kinetic energy of a body

        :return ke: float, 1/2 m v**2
        """
        
        return 0.5*self.mass*np.linalg.norm(self.vel)**2
    
    def momentum(self):
        """
        Returns the linear momentum of a Particle3D instance
        :return p: (3) float np.array, m*v
        """
        return self.mass*self.velocity

    def leap_pos2nd(self, dt, force):
        """
        2nd order position update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        """  
        self.pos += dt*self.vel + (dt**2)*force/(2*self.mass)
        return self.pos 
        
    def leap_vel(self, dt, force):
        """
        Velocity update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        """    #Adds a factor to the current value of self.velocity
        self.vel += dt*force/self.mass
        
        return self.vel
    
    @staticmethod
    def sys_kinetic(p_list):
        """
        Returns the kinetic energy of the whole system in a list

        :param p3d_list: list in which each item is a P3D instance

        :return sys_ke: \sum 1/2 m_i v_i^2 
        """
        sys_ke = 0
        for particle in p_list:
            sys_ke += particle.KE()
            
        return sys_ke


    @staticmethod
    def com_velocity(p_list):
        """
        Computes the total mass and CoM velocity of a list of P3D's

        :param p3d_list: list in which each item is a P3D instance
        :return total_mass: The total mass of the system in
                            array of size N x 1
        :return com_vel: Centre-of-mass velocity in
        """
        N = len(p_list)
        sum_mom = np.zeros([3])
        total_mass = 0
        
        
        for particle in p_list:
            total_mass += particle.mass
            sum_mom += particle.mass*particle.vel
            
        com_vel = sum_mom/total_mass    
        return total_mass, com_vel
    
    
  






