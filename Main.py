# Calculator used to aid design choices for Aere 261 Semester Project

import math
import time
from ambiance import Atmosphere
from matplotlib import pyplot as plt
import numpy as np


def get_float_input(prompt):
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Invalid input. Please enter a number.")


def get_menu_input(prompt, num_selections):
    while True:
        try:
            user_selection = int(input(prompt))
            if user_selection in range(1, num_selections + 1):
                return user_selection
            else:
                print(f"Invalid input. Please enter a number, 1 to {num_selections}.")
        except ValueError:
            print(f"Invalid input. Please enter a number, 1 to {num_selections}.")


def get_yes_no_selection(prompt):  # returns boolean
    while True:
        user_response = input(prompt).lower().strip()
        if user_response == 'y':
            return True
        elif user_response == 'n':
            return False
        else:
            print("Invalid response. Please enter 'y' or 'n'.")


def get_value(name, value):
    # no value
    if value is None:
        return None
    # single value
    if isinstance(value, (int, float)):
        return value

    # multiple values
    if isinstance(value, list):
        print(f"\nMultiple values found for {name}:")
        for i, v in enumerate(value, start=1):
            print(f"{i}: {v}")
        choice = get_menu_input(f'\nChoose a value for {name}: ', len(value))
        return value[choice - 1]


def reynolds(vel, char_len, rho, mu):
    return (rho * vel * char_len) / mu


def mach(vel, temp):
    specific_gas = 287  # metric: J/kgK
    return vel / math.sqrt(1.4 * specific_gas * temp)


def skin_friction(re, mach_num):
    if re > 100000:
        print("\nReynolds is greater than 100,000. Turbulent Flow has been selected")
        return 0.455 / ((math.log10(re) ** 2.58) * ((1 + 0.144 * mach_num ** 2) ** 0.65))
    elif re <= 100000:
        print("\nReynolds is less than 100,000. Laminar flow has been selected\n")
        return 1.328 / math.sqrt(re)


def ff_wing(xc_m, tc, sweep_angle, mach_num):
    return (((1 + ((0.6 / xc_m) * tc) + (100 * (tc ** 4))) *
             (1.34 * mach_num ** 0.18) * ((math.cos(math.radians(sweep_angle))) ** 0.28)))


def ff_fuse(l, A_max):
    f = l / (math.sqrt((4 / math.pi) * A_max))
    return 0.9 + (5 / (f ** 1.5)) + (f / 400)


def ff_nacelle(l, A_max):
    f = l / (math.sqrt((4 / math.pi) * A_max))
    return 1 + (0.35 / f)


def s_wet_wing(s_expo, t_c):
    if t_c < 0.05:
        return 2.003 * s_expo
    else:
        return (1.977 + (0.52 * t_c)) * s_expo


def s_wet_fuse(a_top, a_side):
    return 3.4 * ((a_top + a_side) / 2)


# functions for basic equations of standard atmosphere
def alt_vari(alt):
    # define region variables. could add the other regions and metric vs english units, but not necessary for project
    regions = [
        (0, 11000, {'lapse_rate': -6.5e-3, 't_base': 288.16, 'h_base': 0,
                    'p_base': 101325, 'density_base': 1.225, 'g': 9.81, 'R': 287}
         )
    ]
    for alt_lower_limit, alt_upper_limit, region_variables in regions:
        if alt_lower_limit <= alt <= alt_upper_limit:
            return region_variables


def temp(alt):
    region_variables = alt_vari(alt)  # gather region variables

    return region_variables['t_base'] + (region_variables['lapse_rate'] * (alt - region_variables['h_base']))


def pres(alt):
    region_variables = alt_vari(alt)  # gather region variables

    return (region_variables['p_base'] * ((temp(alt) / region_variables['t_base'])
                                          ** ((-1 * region_variables['g']) / (
                    region_variables['lapse_rate'] * region_variables['R']))))


def density(alt):
    region_variables = alt_vari(alt)  # gather region variables

    return pres(alt) / (region_variables['R'] * temp(alt))


def v_rc(alt, K, Cd_0, W, S):
    return ((2 / density(alt)) * ((K / (3 * Cd_0)) ** (1 / 2)) * (W / S)) ** (1 / 2)


def t_r(alt, K, Cd_0, W, S):
    return (((1 / 2) * density(alt) * ((v_rc(alt, K, Cd_0, W, S)) ** 2) * S * Cd_0)
            + ((2 * K * (W ** 2)) / (density(alt) * (v_rc(alt, K, Cd_0, W, S) ** 2) * S)))


def t_r_const(alt, K, Cd_0, W, S, v):
    return (((1 / 2) * density(alt) * (v ** 2) * S * Cd_0)
            + ((2 * K * (W ** 2)) / (density(alt) * v ** 2 * S)))


# CALCULATOR MADE DURING PROJECT PART B
def cd_0_calculator():
    # Altitude Selection
    user_altitude = get_float_input("Enter your altitude: ")
    alt = Atmosphere(user_altitude)

    # C_f calculation
    print("\nNow gathering information for Reynolds Number: \n")
    velocity = get_float_input("Enter the free stream velocity: ")
    char_length = get_float_input("Enter the characteristic length, ex. chord length or fuselage length: ")

    R = reynolds(velocity, char_length, alt.density, alt.dynamic_viscosity)
    M = mach(velocity, alt.temperature)

    print(f"\nThe Reynold's Number, R, is {R}")
    print(f"The Mach Number, M, is {M}")

    # Calculate Skin Friction
    C_f = skin_friction(R, M)

    print(f'\nThe value for Skin Friction, C_f, is {C_f}')

    # FF_c calculation
    print(f'\nForm Factors:\n1: Wing, Tail, Strut, and Pylon\n2: Fuselage and Smooth Canopy\n3: '
          f'Nacelle and Smooth External Store')

    ff_selection = get_menu_input("\nSelect your Form Factor: ", 3)

    # Wing, Tail, Strut, and Pylon
    if ff_selection == 1:
        print("\nYou have selected Wing, Tail, Strut, and Pylon")
        max_thick = get_float_input("\nEnter Pt of Max Thickness, (x/c)_m: ")
        thick_ratio = get_float_input("\nEnter the Thickness Ratio, (t/c): ")
        sweepback_angle = get_float_input("\nEnter the Sweepback Angle (degrees), Lambda_m: ")

        FF_c = ff_wing(max_thick, thick_ratio, sweepback_angle, mach(velocity, alt.temperature))
        print(f'\nFF_c is: {FF_c}')

    # Fuselage and Smooth Canopy
    elif ff_selection == 2:
        print("\nYou have selected Fuselage and Smooth Canopy")
        length = get_float_input("\nEnter the length of the component, l: ")
        cross_sectional_area_max = get_float_input("\nEnter the Maximum Cross-Sectional Area, A_max: ")

        FF_c = ff_fuse(length, cross_sectional_area_max)
        print(f'\nThe value for FF_c is: {FF_c}')

    # Nacelle and Smooth External Store
    elif ff_selection == 3:
        print("\nYou have selected Nacelle and Smooth External Store")
        length = get_float_input("\nEnter the length of the component, l: ")
        cross_sectional_area_max = get_float_input("\nEnter the Maximum Cross-Sectional Area, A_max: ")

        FF_c = ff_nacelle(length, cross_sectional_area_max)
        print(f'\nFF_c is: {FF_c}')

    # Interference Factor
    Q_c = get_float_input('\nEnter your Interference Factor value, Q: ')

    print(f'\nThe value for the Component Interference Factor, Q, is {Q_c}')

    # Wetted Area Estimation
    # Allow the user to enter a custom wetted area gathered from somewhere like XFLR5
    print('\n1: Calculate Estimated Wetted Area\n2: Custom Wetted Area, ex. from XFLR5 ')
    wetted_area_selection = get_menu_input('\nSelect 1 or 2: ', 2)

    # Calculate an Estimated Wetted Area
    if wetted_area_selection == 1:
        print(f'\n1: Airfoil Shapes\n2: Fuselage')
        wetted_area_shape_selection = get_menu_input('\nSelect your shape: ', 2)

        if wetted_area_shape_selection == 1:
            S_exposed = get_float_input('\nEnter your EXPOSED planform area, S_exposed: ')
            S_wet = s_wet_wing(S_exposed, thick_ratio)
        elif wetted_area_shape_selection == 2:
            Area_top = get_float_input('\nEnter the Projected Area of the TOP of the fuselage, A_top: ')
            Area_side = get_float_input('\nEnter the Projected Area of the SIDE of the fuselage, A_side: ')
            S_wet = s_wet_fuse(Area_top, Area_side)

    elif wetted_area_selection == 2:
        S_wet = get_float_input('\nEnter your value for the Wetted Area, S_wet: ')

    # Calculate S_wet / S_wing
    S_wing = get_float_input('\nEnter a value for the wing planform area: ')
    wet_area_ratio = S_wet / S_wing

    print(f'\nThe value of S_wet / S_wing is: {wet_area_ratio}')

    # Calculate Cd_0
    return C_f * FF_c * Q_c * wet_area_ratio


# CALCULATOR MADE DURING PROJECT PART C
def cl_calculator():
    # Altitude Selection
    user_altitude = 2331.72#get_float_input("\nEnter your altitude: ")
    alt = Atmosphere(user_altitude)

    # Gather Aircraft Info
    gross_weight = 45 * 9.81#get_float_input('\nEnter the Gross Weight of the aircraft: ')
    velocity = 26.8224#get_float_input('\nEnter the Velocity of the aircraft: ')
    planform = 1.00649158#get_float_input('\nEnter the Planform Area of the aircraft: ')

    # Calculate Cl
    cl = gross_weight / (.5 * alt.density * (velocity ** 2) * planform)
    return cl


# CALCULATOR MADE DURING PROJECT PART C
def cd_calculator():
    print(f'\n1: Calculate a new Cd_0\n2: Use a stored value for Cd_0\n3: Enter a custom value for Cd_0')
    cd_0_choice = get_menu_input('\nSelect one: ', 3)
    if cd_0_choice == 1:
        cd_0 = cd_0_calculator()
    elif cd_0_choice == 2:
        cd_0 = get_value('Cd_0', stored_calculations.get('Cd_0'))
        if cd_0 is None:
            print("No stored value found, running Cd_0 calculator first: ")
            cd_0 = cd_0_calculator()
    elif cd_0_choice == 3:
        cd_0 = get_float_input('\nEnter a value for Cd_0: ')
    else:
        cd_0 = get_float_input('\nAn error occurred, please enter a value for Cd_0: ')

    print(f'\nThe value for Cd_0 is: {cd_0}')

    print(f'\n1: Calculate a new Cl\n2: Use a stored value for Cl\n3: Enter a custom value for Cl')
    cl_choice = get_menu_input('\nSelect one: ', 3)
    if cl_choice == 1:
        cl = cl_calculator()
    if cl_choice == 2:
        cl = get_value('C_l', stored_calculations.get('C_l'))
        if cl is None:
            print("No stored value found, running Cl calculator first: ")
            cl = cl_calculator()
    elif cl_choice == 3:
        cl = get_float_input('\nEnter a value for Cl: ')
    else:
        cl = get_float_input('\nAn error occurred, please enter a value for Cl: ')
    print(f'\nThe value for Cl is: {cl}')

    print('\nNow calculating Cd:')
    b = get_float_input('\nEnter your wingspan: ')
    s = get_float_input('\nEnter your planform area: ')

    span_eff_factor = get_float_input('\nEnter your span efficiency factor, e: ')

    aspect_ratio = b ** 2 / s

    k = 1 / (math.pi * span_eff_factor * aspect_ratio)
    cd = cd_0 + (k * (cl ** 2))
    return cd


# Calculator for Project Part E
def rate_of_climb_calculator():
    # aircraft specifications
    W = 294.3  # get_float_input("Enter the Weight, in Newtons: ")
    S = 1.0065  # get_float_input("Enter the Planform Area: ")
    Cl_max = 1.41  # get_float_input("Enter CL_max: ")
    K = 0.038976721  # get_float_input("Enter a value for K: ")
    Cd_0 = 0.0178  # get_float_input("Enter a value for Cd_0: ")
    n = 0.9  # get_float_input("Enter a value for Prop Efficiency, n: ")
    P_A = 2240  # get_float_input("Enter a value for Power Available in Watts: ")

    v_cruise = 26.8224  # get_float_input("Enter a value for Cruise Velocity: ")  # used for constant Vel calculations

    h = np.arange(11001)  # 0 to 11 km w/ 1 meter increments

    # Plotting temp, pressure, and density from 0 to 11 km w/ d_m = 1 meter

    # # TEMPERATURE PLOT
    # plt.figure()
    # plt.plot(np.array([temp(h_i) for h_i in h]), h)
    # plt.xlabel('Temperature (K)')
    # plt.ylabel('Altitude (m)')
    # plt.title('Altitude vs Temperature')
    # plt.grid()
    #
    # # PRESSURE PLOT
    # plt.figure()
    # plt.plot(np.array([pres(h_i) for h_i in h]), h)
    # plt.xlabel('Pressure (Pa)')
    # plt.ylabel('Altitude (m)')
    # plt.title('Altitude vs Pressure')
    # plt.grid()
    #
    # # DENSITY PLOT
    # plt.figure()
    # plt.plot(np.array([density(h_i) for h_i in h]), h)
    # plt.xlabel('Density (Kg/m^3)')
    # plt.ylabel('Altitude (m)')
    # plt.title('Altitude vs Density')
    # plt.grid()
    #
    # plt.show()

    # Plot Power Required, Power Available, Excess Power for a constant velocity
    plt.figure()
    plt.plot(h, np.array([(((t_r_const(h_i, K, Cd_0, W, S, v_cruise)) * v_cruise) / 1000) for h_i in h]),
             label='Power Required', color='orange')
    plt.plot(h, np.array([((P_A * n * (density(h_i) / density(0))) / 1000) for h_i in h]),
             label='Power Available', color='blue')
    plt.plot(h,
             np.array(
                 [(((P_A * n * (density(h_i) / density(0))) / 1000) -
                   ((t_r_const(h_i, K, Cd_0, W, S, v_cruise) * v_cruise) / 1000)) for h_i in h]),
             color='gray', label='Excess Power')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Power (kW)')
    plt.title('Constant Cruise Velocity Power Curves vs Altitude')
    plt.legend()
    plt.ylim(0)
    plt.grid()

    plt.show()

    # Plot Rate of Climb for a constant velocity. (R/C)Cruise
    plt.figure()
    plt.plot(h,
             np.array(
                 [(((P_A * n * (density(h_i) / density(0))) -
                    (t_r_const(h_i, K, Cd_0, W, S, v_cruise) * v_cruise)) / W) for h_i in h]),
             color='gray', label='(R/C)Cruise')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Velocity (m/s)')
    plt.title('Rate of Climb at Constant Cruise Velocity')
    plt.legend()
    plt.ylim(0)
    plt.grid()

    plt.show()

    # Plot V_stall and V_(r/c)max then overlay the plot

    # Standalone plot of V_stall vs altitude
    plt.figure()
    plt.plot(h, np.array([(((2 / density(h_i)) * (W / S) * (1 / Cl_max)) ** (1 / 2)) for h_i in h]))
    plt.xlabel('Altitude (m)')
    plt.ylabel('Velocity (m/s)')
    plt.title('V_stall vs Altitude')
    plt.grid()

    # Standalone plot of V_(r/c)max
    plt.figure()
    plt.plot(h, np.array([v_rc(h_i, K, Cd_0, W, S) for h_i in h]), color='orange', label='V_(r/c)max')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Velocity (m/s)')
    plt.title('V_(r/c)max vs Altitude')
    plt.legend()
    plt.grid()

    # Plot of V_stall and V_(r/c)max
    plt.figure()
    plt.plot(h, np.array([(((2 / density(h_i)) * (W / S) * (1 / Cl_max)) ** (1 / 2)) for h_i in h]), label='V_stall')
    plt.plot(h, np.array([v_rc(h_i, K, Cd_0, W, S) for h_i in h]), color='orange', label='V_(r/c)max')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Velocity (m/s)')
    plt.title('V_stall and V_(r/c)max vs Altitude')
    plt.legend()
    plt.grid()

    plt.show()

    # Plot P_r, then P_r and P_A on one plot, then plot excess power, P_r, and P_A on one plot

    # Plot of Power Required Only
    plt.figure()
    plt.plot(h, np.array([((t_r(h_i, K, Cd_0, W, S) * v_rc(h_i, K, Cd_0, W, S)) / 1000) for h_i in h]))
    plt.xlabel('Altitude (m)')
    plt.ylabel('Power Required (kW)')
    plt.title('Power Required vs Altitude')
    plt.legend()
    plt.grid()

    # Plot of Power Required and Power Available
    plt.figure()
    plt.plot(h, np.array([((t_r(h_i, K, Cd_0, W, S) * v_rc(h_i, K, Cd_0, W, S)) / 1000) for h_i in h]),
             color='orange', label='Power Required')
    plt.plot(h, np.array([((P_A * n * (density(h_i) / density(0))) / 1000) for h_i in h]),
             label='Power Available')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Power (kW)')
    plt.title('Power Curves vs Altitude')
    plt.legend()
    plt.ylim(0)
    plt.grid()

    # Plot of Power Required, Power Available, and Excess Power
    plt.figure()
    plt.plot(h, np.array([((t_r(h_i, K, Cd_0, W, S) * v_rc(h_i, K, Cd_0, W, S)) / 1000) for h_i in h]),
             color='orange', label='Power Required')
    plt.plot(h, np.array([((P_A * n * (density(h_i) / density(0))) / 1000) for h_i in h]),
             label='Power Available')
    plt.plot(h,
             np.array(
                 [(((P_A * n * (density(h_i) / density(0))) / 1000) -
                   ((t_r(h_i, K, Cd_0, W, S) * v_rc(h_i, K, Cd_0, W, S)) / 1000)) for h_i in h]),
             color='gray', label='Excess Power')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Power (kW)')
    plt.title('Variable Velocity Power Curves vs Altitude')
    plt.legend()
    plt.ylim(0)
    plt.grid()

    plt.show()

    # Plot (R/C)max from sea level to (R/C)max = 0 then identify the service and absolute ceilings
    RC_vals = np.array(([(((P_A * n * (density(h_i) / density(0))) -
                           (t_r(h_i, K, Cd_0, W, S) * v_rc(h_i, K, Cd_0, W, S))) / W) for h_i in h]))

    # determine the ceilings
    abs_ceiling = np.where(RC_vals <= 0)[0][0]
    serv_ceiling = np.where(RC_vals <= 0.508)[0][0]
    print(f'The service ceiling is {serv_ceiling}')
    print(f'The absolute ceiling is {abs_ceiling}')

    # Plot of (R/C)max

    plt.figure()
    plt.plot(h, RC_vals,
             color='blue', label='Rate of Climb Max')
    plt.axvline(x=serv_ceiling, color='gray', label=f'Service Ceiling')
    plt.axvline(x=abs_ceiling, color='orange', label='Absolute Ceiling')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Velocity (m/s)')
    plt.title('Rate of Climb Max vs Altitude')
    plt.legend()
    plt.ylim(0)
    plt.grid()
    plt.show()

    # Plot of (R/C)max and (R/C)Cruise
    plt.figure()

    plt.plot(h, RC_vals,
             color='blue', label='(R/C)Max')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Velocity (m/s)')

    plt.plot(h,
             np.array(
                 [(((P_A * n * (density(h_i) / density(0))) -
                    (t_r_const(h_i, K, Cd_0, W, S, v_cruise) * v_cruise)) / W) for h_i in h]),
             color='gray', label='(R/C)Cruise')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Velocity (m/s)')

    plt.title('Rate of Climb Max vs Altitude')
    plt.legend()
    plt.ylim(0)
    plt.grid()
    plt.show()

    # Plot of Pr_(R/C)max and Pr_Cruise
    plt.figure()
    plt.plot(h, np.array([((t_r(h_i, K, Cd_0, W, S) * v_rc(h_i, K, Cd_0, W, S)) / 1000) for h_i in h]),
             color='Green', label='Pr_(R/C)max')
    plt.plot(h, np.array([(((t_r_const(h_i, K, Cd_0, W, S, v_cruise)) * v_cruise) / 1000) for h_i in h]),
             label='Pr_Cruise', color='Red')
    plt.xlabel('Altitude (m)')
    plt.ylabel('Power (kW)')

    plt.title('Power Required Curves vs Altitude')
    plt.legend()
    plt.ylim(0)
    plt.grid()
    plt.show()

    return abs_ceiling, serv_ceiling


# Calculator for Project Part F
def maneuvers():
    # Altitude Selection
    user_altitude = 2331.72  # get_float_input("\nEnter your altitude: ")
    alt = Atmosphere(user_altitude)

    # Gather Aircraft Info
    gross_weight = 294.3  # get_float_input('\nEnter the Gross Weight of the aircraft: ')
    velocity = 26.8224  # get_float_input('\nEnter the Cruise Velocity of the aircraft: ')
    planform = 1.00649158  # get_float_input('\nEnter the Planform Area of the aircraft: ')
    n_lower = -0.872  # get_float_input('\nEnter the Lower Structural Load Factor: ')
    n_upper = 1.51  # get_float_input('\nEnter the Upper Structural Load Factor: ')

    cl_min = -0.698  # get_float_input('\nEnter the Min Coefficient of Lift: ')
    cl_max = 1.42  # get_float_input('\nEnter the Max Coefficient of Lift: ')

    lift_min = cl_min * .5 * alt.density * velocity ** 2 * planform
    lift_max = cl_max * .5 * alt.density * velocity ** 2 * planform

    n_aero_min = lift_min / gross_weight
    n_aero_max = lift_max / gross_weight

    print(f'\nThe Aerodynamic Load Factor, n, is: {n_aero_max}')

    # Calculate Pull-Up Radii
    # pullup_rad_aero_min = ((gross_weight/9.81) * velocity**2) / (lift_min - gross_weight)
    pullup_rad_aero_max = (velocity ** 2) / (9.81 * (n_aero_max - 1))
    # pullup_rad_struct_lower = velocity**2 / (9.81 * (n_lower - 1))
    pullup_rad_struct_upper = (velocity ** 2) / (9.81 * (n_upper - 1))

    # Calculate Level-Turn Radii
    # levelturn_rad_aero_min = velocity ** 2 / (9.81 * (((n_aero_min ** 2) - 1) ** (1 / 2)))
    levelturn_rad_aero_max = velocity ** 2 / (9.81 * (((n_aero_max ** 2) - 1) ** (1 / 2)))
    # levelturn_rad_struct_lower = velocity ** 2 / (9.81 * (((n_lower**2) - 1) ** (1/2)))
    levelturn_rad_struct_upper = velocity ** 2 / (9.81 * (((n_upper ** 2) - 1) ** (1 / 2)))

    pullup_limit = "Structural" if pullup_rad_struct_upper > pullup_rad_aero_max else "Aerodynamic"
    levelturn_limit = "Structural" if levelturn_rad_struct_upper > levelturn_rad_aero_max else "Aerodynamic"

    v_star = ((2 * gross_weight * n_upper) / (alt.density * planform * cl_max)) ** 0.5
    above_below = "above" if velocity > v_star else "below"

    print("\nPull-Up Radii:")
    print(f"  Structural Limit (Upper): {pullup_rad_struct_upper:} m")
    print(f"  Aerodynamic Limit (Upper): {pullup_rad_aero_max:} m")
    print(f"  Limiting Case: {pullup_limit}")

    print("\nLevel Turn Radii:")
    print(f"  Structural Limit (Upper): {levelturn_rad_struct_upper:} m")
    print(f"  Aerodynamic Limit (Upper): {levelturn_rad_aero_max:} m")
    print(f"  Limiting Case: {levelturn_limit}")

    print(f"\nManeuvering Speed (V*): {v_star:} m/s")
    print(f"Cruise Velocity ({velocity:} m/s) is {above_below} V*")


# PROJECT PART G
def takeoff_and_landing_calculator():
    from ambiance import Atmosphere
    import numpy as np
    import matplotlib.pyplot as plt

    user_altitude = 583.4  # get_float_input("\nEnter the Airport's altitude: ")
    alt = Atmosphere(user_altitude)

    weight = 294.3  # get_float_input("Enter the Aircraft's Weight: ")
    mass = weight / 9.81

    planform = 1.00649158  # get_float_input("Enter the Planform Area: ")
    CL_max = 1.42  # get_float_input("Enter CL_max: ")
    CL_rolling = .83853  # get_float_input("Enter CL_rolling: ")
    CD_0 = 0.017800  # get_float_input("Enter CD_0: ")
    k = 0.039  # get_float_input("Enter value of k: ")

    power_available = 3 * .9 * ((alt.density / 1.225))  # get_float_input("Enter the Power Available in kW: ")

    mu_rolling = .04  # get_float_input("Enter the Rolling Friction (typ 0.03 - 0.05): ")
    mu_braking = .4  # get_float_input("Enter the Braking Friction (typ 0.3 - 0.5): ")
    n_flare = 1.1  # get_float_input("Enter the Flare Load Factor: ")
    gamma = 4  # get_float_input("Enter the Climb Angle (deg): ")
    obstacle_height = 15.24  # get_float_input("Enter the Obstacle Height (15.24 m typ): ")
    approach_angle = 5.71  # get_float_input("Enter the Approach Angle (deg): ")

    rho = alt.density
    v_stall = np.sqrt((2 * weight) / (rho * planform * CL_max))
    v_liftoff = 1.2 * v_stall
    v_inf_takeoff = 0.7 * v_liftoff

    thrust = (power_available * 1000) / v_inf_takeoff
    print(f"\nThrust available: {thrust} N")

    drag_liftoff = 0.5 * rho * v_inf_takeoff ** 2 * planform * (CD_0 + k * CL_rolling ** 2)
    lift_liftoff = 0.5 * rho * v_inf_takeoff ** 2 * planform * CL_rolling

    dis_ground_takeoff = (1.44 * (weight / planform)) / (9.81 * rho * CL_max *
                                                         ((thrust / weight) - (drag_liftoff / weight) -
                                                          mu_rolling * (1 - lift_liftoff / weight)))

    R_transition = (v_liftoff ** 2) / ((n_flare - 1) * 9.81)
    h_transition = R_transition * (1 - np.cos(np.radians(gamma)))
    dis_transition_takeoff = R_transition * np.sin(np.radians(gamma))

    dis_climb_takeoff = max(0, (obstacle_height - h_transition) / np.tan(np.radians(gamma)))

    dis_takeoff = dis_ground_takeoff + dis_transition_takeoff + dis_climb_takeoff
    print(f"Minimum TAKEOFF distance: {dis_takeoff} m")

    v_touchdown = 1.3 * v_stall
    v_inf_landing = 0.7 * v_touchdown

    drag_landing = 0.5 * rho * v_inf_landing ** 2 * planform * (CD_0 + k * CL_rolling ** 2)
    lift_landing = 0.5 * rho * v_inf_landing ** 2 * planform * CL_rolling

    dis_ground_landing = (1.69 * (weight / planform)) / (9.81 * rho * CL_max *
                                                         ((drag_landing / weight) + mu_braking * (
                                                                     1 - lift_landing / weight)))

    R_flare = (v_touchdown ** 2) / ((n_flare - 1) * 9.81)
    h_flare = R_flare * (1 - np.cos(np.radians(approach_angle)))
    dis_flare_landing = R_flare * np.sin(np.radians(approach_angle))

    dis_approach_landing = max(0, (obstacle_height - h_flare) / np.tan(np.radians(approach_angle)))

    dis_landing = dis_ground_landing + dis_flare_landing + dis_approach_landing
    print(f"Minimum LANDING distance: {dis_landing} m")

    # grafin
    x_ground = np.linspace(0, dis_ground_takeoff, 200)
    y_ground = np.zeros_like(x_ground)

    theta = np.linspace(0, np.radians(gamma), 100)
    x_transition = dis_ground_takeoff + R_transition * np.sin(theta)
    y_transition = R_transition * (1 - np.cos(theta))

    x_climb = np.linspace(x_transition[-1],
                          x_transition[-1] + dis_climb_takeoff,
                          200)
    y_climb = y_transition[-1] + np.tan(np.radians(gamma)) * (x_climb - x_transition[-1])

    plt.figure(figsize=(10, 6))
    plt.plot(x_ground, y_ground, label="Ground Roll")
    plt.plot(x_transition, y_transition, label="Transition Arc")
    plt.plot(x_climb, y_climb, label="Climb")
    plt.axhline(obstacle_height, color='r', linestyle='--', label="Obstacle Height")
    plt.title("Takeoff Trajectory")
    plt.xlabel("Distance (m)")
    plt.ylabel("Height (m)")
    plt.grid()
    plt.legend()

    x_approach = np.linspace(0, dis_approach_landing, 200)
    y_approach = obstacle_height - np.tan(np.radians(approach_angle)) * x_approach

    theta2 = np.linspace(0, np.radians(approach_angle), 100)
    x_flare = dis_approach_landing + ((-R_flare) * np.sin(theta2)) + dis_flare_landing
    y_flare = R_flare * (1 - np.cos(theta2))

    x_landing_ground = np.linspace(dis_approach_landing + dis_flare_landing,
                                   dis_approach_landing + dis_flare_landing + dis_ground_landing, 200)
    y_landing_ground = np.zeros_like(x_landing_ground)

    plt.figure(figsize=(10, 6))
    plt.plot(x_approach, y_approach, label="Approach")
    plt.plot(x_flare, y_flare, label="Flare Arc")
    plt.plot(x_landing_ground, y_landing_ground, label="Ground Roll")
    plt.axhline(obstacle_height, linestyle="--", color="r", label="Obstacle Height")
    plt.title("Landing Trajectory")
    plt.xlabel("Distance (m)")
    plt.ylabel("Height (m)")
    plt.grid()
    plt.legend()
    plt.show()


calc_dict = {1: "Cd,0", 2: "Cl", 3: "Cd", 4: "Rate Of Climb Plots and Operating Ceilings", 5: "Maneuvers",
             6: "Takeoff and Landing", 7: "Range and Endurance"}

stored_calculations = {}

print("\n\nWARNING: ALL UNITS MUST BE METRIC\n\n")  # std atm package only gives metric outputs and i aint doing allat
time.sleep(1)


def range_and_endurance_calculator():
    user_altitude = 2331.72
    alt = Atmosphere(user_altitude)

    CD_0 = 0.017888
    CL_max = 1.41
    k = 0.038976721
    bsfc = 9.7612e-7  # Wf dot / P ((N/s)/Watt); for Lycoming: ((0.48[kg/hr] * 9.81[m^2/sec]) / 3600[sec]) / 1340[Watts]
    fuel_mass = 1  # [kg] or 3.56 US gallons
    fuel_weight = fuel_mass * 9.81
    gtow = 30 * 9.81
    eta = 0.9
    planform = 1.00649158
    vel_cruise = 26.8224

    CL_three_halves_to_CD = 0.25 * ((3 / (k * (CD_0 ** (1 / 3)))) ** (3 / 4))
    CL_to_CD = (1 / (4 * k * CD_0)) ** (1 / 2)

    v_stall = ((2 * gtow) / (alt.density * planform * CL_max)) ** (1 / 2)
    print(f'The Stall Velocity is {v_stall}\n')



    max_endurance = ((eta / bsfc) * ((2 * alt.density * planform) ** (1 / 2)) * CL_three_halves_to_CD
                     * (((gtow - fuel_weight) ** (-1 / 2)) - (gtow ** (-1 / 2))))
    vel_max_endurance = ((2 / alt.density) * (gtow / planform) * ((k / (3 * CD_0)) ** (1 / 2))) ** (1 / 2)

    print(f'The Max Endurance is {max_endurance / 3600} hours')
    print(f'The Velocity for Max Endurance is {vel_max_endurance}\n')

    thrust_req = 0.5 * alt.density * (vel_cruise ** 2) * planform * CD_0 + ((2 * k * (gtow ** 2)) / (alt.density * (vel_cruise ** 2) * planform))
    p_req = 0.5 * alt.density * (vel_cruise ** 3) * planform * CD_0 + ((2 * k * (gtow ** 2)) / (alt.density * vel_cruise * planform))

    print(f'The Thrust Required at Cruise Velocity is {thrust_req}\n')
    print(f'The Power Required at Cruise Velocity is {p_req}\n')

    max_range = (eta / bsfc) * CL_to_CD * math.log(gtow / (gtow - fuel_weight))
    vel_max_range = ((2 / alt.density) * (gtow / planform) * ((k / CD_0) ** (1 / 2))) ** (1 / 2)
    print(f'The Max Range is {max_range / 1000} km')
    print(f'The Velocity for Max Range is {vel_max_range}')


# main loop
while True:

    # Calculator Selection
    print('\nCurrent Calculators:')

    for calc in calc_dict:
        print(f'{calc}: {calc_dict[calc]}')

    calc_selection = get_menu_input("\nSelect your calculator: ", len(calc_dict))
    print(f'\nYou have selected the {calc_dict[calc_selection]} calculator')

    # Calculators:

    # Cd,0
    if calc_selection == 1:

        Swing = 1.0065
        Sbackwheel = 0.00098564
        Sfrontwheel = 0.00062650
        Sfrontstrut = 0.00111361
        Srearstrut = 0.00110969

        Cd_0backwheel = 0.25 * (Sbackwheel / Swing) * 1  # more than a diameter away => q = 1
        Cd_0frontwheel = 0.25 * (Sfrontwheel / Swing) * 1  # more than a diameter away => q = 1

        Cd_0frontstrut = 1 * (Sfrontstrut / Swing) * 1.05  # not streamlined, fork
        Cd_0rearstrut = 0.05 * (Srearstrut / Swing) * 1.05  # streamlined

        Cd_0Gear = Cd_0frontstrut + Cd_0rearstrut + Cd_0frontwheel + 2 * Cd_0backwheel

        print(f'The value for CD_0_Landing Gear is {Cd_0Gear}')

        cd_0_aircraft_dict = {}

        while True:
            component_name = input('\nEnter the name of the component: ')
            Cd_0_c = cd_0_calculator()

            print(f'\nThe value for Cd_0_{component_name} is {Cd_0_c}')

            cd_0_aircraft_dict[component_name] = float(Cd_0_c)

            if not get_yes_no_selection('\nWould you like to calculate Cd_0 for another component? (y/n): '):
                break
        while True:
            if get_yes_no_selection('\nDo you want to enter additional, known component values for Cd_0? (y/n): '):
                component_name = input('\nEnter the name of the component: ')
                Cd_0_c = get_float_input(f'\nEnter the value for Cd_0_{component_name}: ')

                cd_0_aircraft_dict[component_name] = float(Cd_0_c)
            else:
                break

        for key in cd_0_aircraft_dict.keys():
            print(f'\nCd_0_{key} = {cd_0_aircraft_dict[key]}')

        cd_0_aircraft = sum((cd_0_aircraft_dict.values()))
        print(f'\n\nThe value of Cd_0 for the entered components is: {cd_0_aircraft}')

        if "Cd_0" not in stored_calculations:  # no stored values yet
            stored_calculations["Cd_0"] = [float(cd_0_aircraft)]
        else:
            stored_calculations["Cd_0"].append(float(cd_0_aircraft))

    # Cl
    if calc_selection == 2:
        C_l = cl_calculator()
        print(f'\n\nThe value of C_l is: {C_l}')

        if "C_l" not in stored_calculations:  # no stored values yet
            stored_calculations["C_l"] = [float(C_l)]
        else:
            stored_calculations["C_l"].append(float(C_l))

    # Cd
    if calc_selection == 3:
        C_d = cd_calculator()
        print(f'\n\nThe value of C_d is: {C_d}')

        if "C_d" not in stored_calculations:  # no stored values yet
            stored_calculations["C_d"] = [float(C_d)]
        else:
            stored_calculations["C_d"].append(float(C_d))

    if calc_selection == 4:
        (abs_ceiling, serv_ceiling) = rate_of_climb_calculator()

    if calc_selection == 5:
        maneuvers()

    if calc_selection == 6:
        takeoff_and_landing_calculator()

    if calc_selection == 7:
        range_and_endurance_calculator()

    # ENDING PROGRAM PROCEDURE
    if get_yes_no_selection("\nWould you like to use another calculator? y/n: "):
        continue
    else:
        print('\nThe final values for all used calculators were as follows: \n')
        for val in stored_calculations:
            print(f'The value(s) of {val} were: {stored_calculations[val]}\n')
        break
