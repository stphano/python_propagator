import numpy as np
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.astro import conversion
from tudatpy.io import save2txt
from matplotlib import pyplot as plt

# Problem-specific imports
from defining_thrust import defining_thrust
from defining_thrust import get_thrust_acceleration_model_from_parameters
import defining_thrust as Util

# # define initial and end epoch
simulation_start_epoch = 10000.0
simulation_end_epoch = constants.JULIAN_DAY

# define thrust parameters
thrust_parameters = [0.5,
                     24 / 5,
                     -0.03344538412056863,
                     -0.06456210720352829,
                     0.3943447499535977,
                     0.5358478897251189,
                     -0.8607350478880107,
                     -0.03344538412056863,
                     -0.06456210720352829,
                     0.3943447499535977,
                     0.5358478897251189,
                     -0.8607350478880107,
                     ]
# # load spice ephemerides
spice_interface.load_standard_kernels()
# # define settings for the celestial bodies and propagation
bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus", "Mercury", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create,
    simulation_start_epoch,
    simulation_end_epoch,
    "SSB", "ECLIPJ2000")
bodies = environment_setup.create_system_of_bodies(body_settings)

# #  Create veichle
bodies.create_empty_body("3UCubeSat")
# set mass
bodies.get_body("3UCubeSat").set_constant_mass(4.0)
# set radiation coefficients
reference_area_radiation = 0.4
radiation_pressure_coefficient = 1.2
constant_specific_impulse = 311.0  # s
occulting_bodies = []
radiation_pressure_settings = environment_setup.radiation_pressure.cannonball(
    "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies
)
environment_setup.add_radiation_pressure_interface(
    bodies, "3UCubeSat", radiation_pressure_settings)

# # create accelerations
bodies_to_propagate = ["3UCubeSat"]
central_bodies = ["Sun"]
accelerations_settings_3u_cubesat = dict(
    bodies_to_propagate = [get_thrust_acceleration_model_from_parameters(thrust_parameters,
                                                              bodies,
                                                              simulation_start_epoch,
                                                              constant_specific_impulse)],
   # "3U-CubeSat" :
   # [
   #     propagation_setup.acceleration.thrust_acceleration()
   # ],
    Sun=
    [
        propagation_setup.acceleration.cannonball_radiation_pressure(),
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Earth=
    [
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Moon=
    [
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Mars=
    [
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Venus=
    [
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Mercury =
    [
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Jupiter =
    [
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Saturn =
    [
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Uranus=
    [
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Neptune=
    [
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Pluto=
    [
        propagation_setup.acceleration.point_mass_gravity()
    ]
)

acceleration_settings = {"3UCubeSat": accelerations_settings_3u_cubesat}
acceleration_models = propagation_setup.create_acceleration_models(
    bodies,
    acceleration_settings,
    bodies_to_propagate,
    central_bodies)

# # Create propagation settings
# define initial conditions
sun_gravitational_parameter = bodies.get_body("Sun").gravitational_parameter
initial_state = conversion.keplerian_to_cartesian(
    gravitational_parameter=sun_gravitational_parameter,
    semi_major_axis=1.3*149597870700.0,
    eccentricity=0.1,
    inclination=np.deg2rad(1.0),
    argument_of_periapsis=np.deg2rad(235.7),
    longitude_of_ascending_node=np.deg2rad(23.4),
    true_anomaly=np.deg2rad(139.87)
)
# define the variables we need
dependent_variables_to_save = [
    propagation_setup.dependent_variable.total_acceleration("3UCubeSat"),
    propagation_setup.dependent_variable.keplerian_state("3UCubeSat", "Sun"),
]
# create propagation settings
propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models,
    bodies_to_propagate,
    initial_state,
    simulation_end_epoch,
    output_variables=dependent_variables_to_save
)
# numerical integrator settings
fixed_step_size = 10.0

integrator_settings = propagation_setup.integrator.runge_kutta_4(
    simulation_start_epoch,
    fixed_step_size
)
# # propagate orbit
dynamics_simulator = propagation_setup.SingleArcDynamicsSimulator(
    bodies, integrator_settings, propagator_settings)

states = dynamics_simulator.state_history
dependent_variables = dynamics_simulator.dependent_variable_history

save2txt(solution=states,
         filename="PropagationHistory_3U_CubeSat_thrust.dat",
         directory="./",  # default = "./"  # default = None
         )
