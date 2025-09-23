# Calculator used to aid design choices for Aere 261 Semester Project

import math
import time

from ambiance import Atmosphere


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


def reynolds(vel, char_len, rho, mu):
    return (rho * vel * char_len) / mu


def mach(vel, temp):
    specific_gas = 287  # metric: J/kgK
    return vel / math.sqrt(1.4 * specific_gas * temp)


def skin_friction(re, mach_num):
    if R > 100000:
        print("\nReynolds is greater than 100,000. Turbulent Flow has been selected")
        return 0.455 / ((math.log10(15645290) ** 2.58) * ((1 + 0.144 * 0.282 ** 2) ** 0.65))
    elif R <= 100000:
        print("\nReynolds is less than 100,000. Laminar flow has been selected\n")
        return 1.328 / math.sqrt(R)


def ff_wing(xc_m, tc, sweep_angle):
    return ((1 + ((0.6 / xc_m) * tc) + (100 * (tc ** 4))) *
            ((1.34 * ((mach(velocity, alt.temperature)) ** 0.18)) * ((math.cos(math.radians(sweep_angle))) ** 0.28)))


def ff_fuse(l, A_max):
    f = l / (math.sqrt((4 / math.pi) * A_max))
    return 0.9 + (5 / (f ** 1.5)) + (f / 400)


def ff_nacelle(l, A_max):
    f = l / (math.sqrt((4 / math.pi) * A_max))
    return 1 + (0.35 / f)


def s_wet_wing(s_expo):
    if thick_ratio < 0.05:
        return 2.003 * S_exposed
    else:
        return (1.977 + (0.52 * thick_ratio)) * S_exposed


def s_wet_fuse(a_top, a_side):
    return 3.4 * ((a_top + a_side) / 2)


calc_dict = {1: "Cd,0", 2: "Cl"}

print("\n\nWARNING: ALL UNITS MUST BE METRIC\n\n")  # std atm package only gives metric outputs and i aint doing allat
time.sleep(1)

# User Selection loops

# while True:
#     try:
#         raw_unit_selection = input("Enter 1 for Metric Units or 2 for English Units: ")
#         int(raw_unit_selection)
#     except ValueError:
#         print('Invalid Selection')
#     else:
#         unit_selection = int(raw_unit_selection)
#         if unit_selection == 1:
#             print("Metric Units are Selected\n")
#         else:
#             print("English Units are Selected\n")
#         break

# Calculator Selection

print('\nCurrent Calculators:')

for calc in calc_dict:
    print(f'{calc}: {calc_dict[calc]}')

calc_selection = get_menu_input("Select your calculator: ", len(calc_dict))
print(f'You have selected the {calc_dict[calc_selection]} calculator\n')

# Altitude Selection

user_altitude = get_float_input("Enter you altitude: ")
alt = Atmosphere(user_altitude)

# Calculators:

# Cd,0
if calc_selection == 1:

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

        FF_c = ff_wing(max_thick, thick_ratio, sweepback_angle)
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
            S_wet = s_wet_wing(S_exposed)
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
    Cd_0 = C_f * FF_c * Q_c * wet_area_ratio

    print(f'\n\nThe value for Cd_0 is: {Cd_0}')
