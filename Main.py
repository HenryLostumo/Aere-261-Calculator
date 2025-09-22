# Calculator used to aid design choices for Aere 261 Semester Project

import math
import time

from ambiance import Atmosphere


def reynolds(vel, char_len):
    return (alt.density * vel * char_len) / alt.dynamic_viscosity


def mach(vel):
    specific_gas = 287  # metric: J/kgK
    return vel / math.sqrt(1.4 * specific_gas * alt.temperature)


def ff_wing(xc_m, tc, sweep_angle):
    return ((1 + ((0.6 / xc_m) * tc) + (100 * (tc ** 4))) *
            ((1.34 * ((mach(velocity)) ** 0.18)) * (math.cos(math.radians(sweep_angle))) ** 0.28))


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
while True:
    try:
        print('\nCurrent Calculators:')
        for calc in calc_dict:
            print(f'{calc}: {calc_dict[calc]}')
        raw_calc_selection = input('\nSelect a calculator: ')
        int(raw_calc_selection)
    except ValueError:
        print('\nInvalid Selection')
    else:
        calc_selection = int(raw_calc_selection)
        print(f'You have selected the {calc_dict[calc_selection]} calculator\n')
        break

# Altitude Selection
while True:
    try:
        raw_user_altitude = input("Enter your altitude: ")
        float(raw_user_altitude)
    except ValueError:
        print('\nInvalid Selection')
    else:
        user_altitude = float(raw_user_altitude)
        alt = Atmosphere(user_altitude)
        break

# Calculators:

# Cd,0
if calc_selection == 1:

    Q_c = []  # initialize interference factors

    # C_f calculation
    print("\nEnter information for Reynolds Number: \n")
    velocity = float(input("Enter the free stream velocity: "))
    char_length = float(input("Enter the characteristic length, ex. chord length or fuselage length: "))

    print(alt.dynamic_viscosity)
    R = reynolds(velocity, char_length)
    M = mach(velocity)

    print(f"The Reynold's Number, R, is {R}")
    print(f"The Mach Number, M, is {M}")

    # Choose Flow Type
    if R > 100000:
        print("\nReynolds is greater than 100,000. Turbulent Flow has been selected")
        print(R)
        print(M)
        C_f = 0.455 / ((math.log10(15645290) ** 2.58) * ((1 + 0.144 * 0.282 ** 2) ** 0.65))
        #C_f = 0.455 / (((math.log10(R)) ** 2.58) * (1 + (0.114 * (M ** 2))) ** 0.65)
    elif R <= 100000:
        print("\nReynolds is less than 100,000. Laminar flow has been selected\n")
        C_f = 1.328 / math.sqrt(R)
    print(f'\nThe value for Skin Friction, C_f, is {C_f}')

    # FF_c calculation
    print(f'\nForm Factors:\n1: Wing, Tail, Strut, and Pylon\n2: Fuselage and Smooth Canopy\n3: '
          f'Nacelle and Smooth External Store')
    while True:
        try:
            raw_ff_selection = input('\nSelect your Form Factor: ')
            int(raw_ff_selection)
        except ValueError:
            print('\nInvalid Selection')
        else:
            ff_selection = int(raw_ff_selection)
            break

    # Wing, Tail, Strut, and Pylon
    if ff_selection == 1:
        print("\nYou have selected Wing, Tail, Strut, and Pylon")
        max_thick = float(input("\nEnter Pt of Max Thickness, (x/c)_m: "))
        thick_ratio = float(input("\nEnter the Thickness Ratio, (t/c): "))
        sweepback_angle = float(input("\nEnter the Sweepback Angle, Lambda_m: "))

        FF_c = ff_wing(max_thick, thick_ratio, sweepback_angle)
        print(f'\nFF_c is: {FF_c}')

    # Fuselage and Smooth Canopy
    elif ff_selection == 2:
        print("\nYou have selected Fuselage and Smooth Canopy")
        length = float(input("\nEnter the length of the component, l: "))
        cross_sectional_area_max = float(input("\nEnter the Maximum Cross-Sectional Area, A_max: "))

        FF_c = ff_fuse(length, cross_sectional_area_max)
        print(f'\nThe value for FF_c is: {FF_c}')

    # Nacelle and Smooth External Store
    elif ff_selection == 3:
        print("\nYou have selected Nacelle and Smooth External Store")
        length = float(input("\nEnter the length of the component, l: "))
        cross_sectional_area_max = float(input("\nEnter the Maximum Cross-Sectional Area, A_max: "))

        FF_c = ff_nacelle(length, cross_sectional_area_max)
        print(f'\nFF_c is: {FF_c}')
    else:
        print("Invalid Input")
        FF_c = 1

    # Interference Factor
    while True:
        while True:
            try:
                raw_int_factor = input('\nEnter your Interference Factor value, Q: ')
                float(raw_int_factor)
            except ValueError:
                print('\nInvalid Input')
            else:
                Q_c.append(float(raw_int_factor))
                break
        int_factor_continue = input("Do you have additional component interference factors y/n: ")
        if int_factor_continue == "y":
            continue
        elif int_factor_continue == "n":
            break
        else:
            print('Invalid Input')
            continue
    print(f'\nThe value(s) for the Component Interference Factor, Q, is {Q_c}')

    # Wetted Area Estimation
    while True:  # Allow the user to enter a custom wetted area gathered from somewhere like XFLR5
        try:
            print('\n1: Calculate Estimated Wetted Area\n2: Custom Wetted Area, ex. from XFLR5 ')
            raw_wetted_area_selection = int(input('Select 1 or 2: '))
        except ValueError:
            print('\nInvalid Selection')
        else:
            wetted_area_selection = int(raw_wetted_area_selection)
            break

    # Calculate an Estimated Wetted Area
    if wetted_area_selection == 1:
        print(f'\n1: Airfoil Shapes\n2: Fuselage')
        while True:
            try:
                raw_wetted_area_shape_selection = input('\nSelect your shape: ')
                int(raw_wetted_area_shape_selection)
            except ValueError:
                print('\nInvalid Selection')
            else:
                wetted_area_shape_selection = int(raw_wetted_area_shape_selection)
                break
        if wetted_area_shape_selection == 1:
            S_exposed = float(input('\nEnter your planform area, S_exposed: '))
            S_wet = s_wet_wing(S_exposed)
        elif wetted_area_shape_selection == 2:
            Area_top = float(input('\nEnter the Projected Area of the TOP of the fuselage, A_top: '))
            Area_side = float(input('\nEnter the Projected Area of the SIDE of the fuselage, A_side: '))
            S_wet = s_wet_fuse(Area_top, Area_side)
    elif wetted_area_selection == 2:
        S_wet = float(input('\nEnter your value for the Wetted Area, S_wet: '))
    else:
        print("Invalid input")
        S_wet = 0

    # Calculate S_wet / S_wing
    S_wing = float(input('\nEnter a value for the wing planform area: '))
    wet_area_ratio = S_wet / S_wing

    print(f'\nThe value of S_wet / S_wing is: {wet_area_ratio}')

    # Calculate Cd_0
    Cd_0 = C_f * FF_c * math.prod(Q_c) * wet_area_ratio

    print(f'\n\nThe value for Cd_0 is: {Cd_0}')
