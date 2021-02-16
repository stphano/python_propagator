# General imports
import numpy as np

# Tudatpy imports
import tudatpy
from tudatpy.kernel import constants
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.astro import frames
from tudatpy.kernel.math import interpolators

class ThrustGuidance:
    """
    Class that defines and updates the thrust guidance of the deep space propagator problem at each time step.
    Attributes
    ----------
    vehicle_body
    initial_time
    parameter_vector
    Methods
    -------
    get_current_thrust_magnitude(time)
    get_current_thrust_direction(time)
    """

    def __init__(self,
                 vehicle_body: str,
                 initial_time: float,
                 parameter_vector: list):
        """
        Constructor of the ThrustGuidance class.
        Parameters
        ----------
        vehicle_body : str
            Name of the vehicle to apply the thrust guidance to.
        initial_time: float
            Initial time of the simulation [s].
        parameter_vector : list
            List of thrust parameters to retrieve the thrust guidance from.
        Returns
        -------
        none
        """
        print('Initializing guidance...')
        # Set arguments as attributes
        self.vehicle_body = vehicle_body
        self.initial_time = initial_time
        self.parameter_vector = parameter_vector
        self.thrust_magnitude = parameter_vector[0]
        self.time_interval = parameter_vector[1]
        # Prepare dictionary for thrust angles
        self.thrust_angle_dict = {} # this defines the in-plane angle
        self.thrust_angle_2_dict = {} # this defines the out-of-plane angle
        # Initialize time
        current_time = initial_time
        # Loop over nodes
        for i in range(4): #((len(parameter_vector) - 2)/2): # qui mi da non risolto perche e' ancora basato sulla lunghezza precedente
            # Store time as key, thrust angle as value
            self.thrust_angle_dict[current_time] = parameter_vector[i + 2]
            self.thrust_angle_2_dict[current_time] = parameter_vector[i + 2 + 5]#(len(parameter_vector) - 2)/2]
            # Increase time
            current_time += self.time_interval
        # Create interpolator settings
        interpolator_settings = interpolators.linear_interpolation(
            boundary_interpolation=interpolators.use_boundary_value)
        # Create the interpolator between nodes and set it as attribute
        self.thrust_angle_interpolator = interpolators.create_one_dimensional_interpolator(self.thrust_angle_dict,
                                                                                           interpolator_settings)
        self.thrust_angle_2_interpolator = interpolators.create_one_dimensional_interpolator(self.thrust_angle_2_dict,
                                                                                           interpolator_settings)

    def get_current_thrust_direction(self,
                                     time: float) -> np.array:
        """
        Retrieves the direction of the thrust in the inertial frame.
        This function is needed to get the thrust direction at each time step of the propagation, based on the thrust
        parameters provided.
        Parameters
        ----------
        time : float
            Current time.
        Returns
        -------
        thrust_inertial_frame : np.array
            Thrust direction expressed in the inertial frame.
        """
        # Interpolate with time
        angle_1 = self.thrust_angle_interpolator.interpolate(time) #in-plane
        angle_2 = self.thrust_angle_2_interpolator.interpolate(time) #out-of-plane
        # Set thrust in vertical frame and transpose it
        # thrust_direction_vertical_frame = np.array([[0, np.sin(angle), - np.cos(angle)]]).T
        #thrust_direction_body_frame = np.array([[np.cos(angle_1)*np.cos(angle_2), np.sin(angle_1)*np.cos(angle_2), np.sin(angle_2)]]).T
        thrust_inertial_frame = np.array(
            [[np.cos(angle_1) * np.cos(angle_2), np.sin(angle_1) * np.cos(angle_2), np.sin(angle_2)]]).T
        # Retrieve spacecraft state to compute the rotation matrix
        #spacecraft_state = self.vehicle_body.get_state(time)
        #radial = np.array(spacecraft_state[0], spacecraft_state[1], spacecraft_state[2])
        #velocity = np.array(spacecraft_state[3], spacecraft_state[4], spacecraft_state[5])
        #x_axis = -radial / np.linalg.norm(radial)
        #z_axis = np.cross(velocity/np.linalg.norm(radial), x_axis)
        #y_axis = np.cross(z_axis, x_axis)
        #body_to_inertial_frame = np.matrix([[x_axis[0], y_axis[0], z_axis[0]], [x_axis[1], y_axis[1], z_axis[1]], [x_axis[2], y_axis[2], z_axis[2]]])
        # Compute the thrust in the inertial frame
        #thrust_inertial_frame = np.dot(body_to_inertial_frame,
         #                              thrust_direction_body_frame)
        return thrust_inertial_frame

    def get_current_thrust_magnitude(self,
                                     time: float) -> float:
        """
        Retrieves the magnitude of the thrust.
        This function is needed to get the thrust magnitude at each time step of the propagation, based on the thrust
        parameters provided. Currently, the thrust magnitude is constant.
        Parameters
        ----------
        time : float
            Current time.
        Returns
        -------
        float
            Thrust magnitude.
        """
        return self.thrust_magnitude



#########
def get_thrust_acceleration_model_from_parameters(thrust_parameters: list,
                                                  bodies: tudatpy.kernel.simulation.environment_setup.SystemOfBodies,
                                                  initial_time: float,
                                                  specific_impulse: float) -> \
        tudatpy.kernel.simulation.propagation_setup.acceleration.ThrustAccelerationSettings:
    """
    Creates the thrust acceleration models from the LunarAscentThrustGuidance class and sets it in the propagator.
    Parameters
    ----------
    thrust_parameters : list of floats
        List of thrust parameters.
    bodies : tudatpy.kernel.simulation.environment_setup.SystemOfBodies
        System of bodies present in the simulation.
    initial_time : float
        The start time of the simulation in seconds.
    specific_impulse : float
        Specific impulse of the vehicle.
    Returns
    -------
    tudatpy.kernel.simulation.propagation_setup.acceleration.ThrustAccelerationSettings
        Thrust acceleration settings object.
    """

######
    # Create Thrust Guidance object
    thrust_guidance = ThrustGuidance(bodies.get_body('CubeSat'),
                                                initial_time,
                                                thrust_parameters)
    # Retrieves thrust functions
    thrust_direction_function = thrust_guidance.get_current_thrust_direction
    thrust_magnitude_function = thrust_guidance.get_current_thrust_magnitude
    # Set thrust functions in the acceleration model
    thrust_direction_settings = propagation_setup.acceleration.custom_thrust_direction(thrust_direction_function)
    thrust_magnitude_settings = propagation_setup.acceleration.custom_thrust_magnitude(thrust_magnitude_function,
                                                                                       lambda time: specific_impulse)
    # Create and return thrust acceleration settings
    return propagation_setup.acceleration.ThrustAccelerationSettings(thrust_direction_settings,
                                                                     thrust_magnitude_settings)
