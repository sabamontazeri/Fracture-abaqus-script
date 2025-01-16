import math
import csv
Thrusts=[]
cutting=[]
def calculate_cutting_force_and_thrust_force():
    print("Drill Cutting Force and Thrust Force Calculator")

    # Material properties dictionary
    materials = {
        "S235": {0.1: 5635, 0.2: 4990, 0.3: 4475, 0.5: 3930, 1.0: 3070},
        "E295": {0.1: 5850, 0.2: 5250, 0.3: 4700, 0.5: 4140, 1.0: 3230},
        "C45": {0.1: 4950, 0.2: 4350, 0.3: 3950, 0.5: 3450, 1.0: 2700},
        # Add more materials as needed
    }

    # Workpiece properties dictionary
    workpiece = {
        "Structured_steel": {3: 0.05, 6: 0.14, 10: 0.25, 16: 0.32, 25: 0.43},
        "Casehardenedsteel": {3: 0.04, 6: 0.1, 10: 0.15, 16: 0.22, 25: 0.3},
        "Toolssteel": {3: 0.024, 6: 0.05, 10: 0.1, 16: 0.15, 25: 0.2},
    }

    # Fixed inputs
    cutting_velocity = 28  # vc in m/min
    point_angle = 118  # sigma in degrees
    correction_factor = 1.3  # C

    # Iterate over all workpieces, drill sizes, materials, and chip thicknesses
    for workpiece_material, drill_sizes in workpiece.items():
        print(f"\nWorkpiece Material: {workpiece_material}")
        for drill_diameter, feed_rate in drill_sizes.items():
            for material, thicknesses in materials.items():
                print(f"  Material: {material}")
                for chip_thickness, kc in thicknesses.items():
                    # Calculations
                    spindle_speed = (1000 * cutting_velocity) / (math.pi * drill_diameter)  # N in rpm
                    point_angle_radians = math.radians(point_angle)
                    calculated_chip_thickness = (feed_rate / 2) * math.sin(point_angle_radians / 2)  # h = (f / 2) * sin(sigma / 2)
                    chip_section_area = (math.pi * drill_diameter * calculated_chip_thickness) / 2  # A = (pi * d * h) / 2
                    cutting_force = 1.2 * kc * chip_section_area * correction_factor  # Fc = 1.2 * Kc * A * C
                    crossarea=math.pi*drill_diameter*feed_rate
                    Ft=kc*crossarea
                    # Output results for each combination
                    print(f"    Thrust Force (Ft): {Ft:.2f} N")
                    Thrusts.append(Ft)
                    print(f"    Cutting Force (Fc): {cutting_force:.2f} N")
                    cutting.append(cutting_force)
   

if __name__ == "__main__":
    calculate_cutting_force_and_thrust_force()
print(cutting)
with open("cutting_force_results.csv", "w", newline="") as file:
        writer = csv.writer(file)
        for item in cutting:
            writer.writerow([item])

with open("thrust_force_results.csv", "w", newline="") as file:
        writer = csv.writer(file)
        for item in Thrusts:
            writer.writerow([item])


