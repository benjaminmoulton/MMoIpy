import sys
mmoipy_directory = '../mmoipy/'
sys.path.insert(1, mmoipy_directory)
from system import AircraftSystem as mmoi

if __name__ == "__main__":
    # test several inputs together
    ex_1 = mmoi("test_input.json")

    # report total mass properties
    ex_1.get_mass_properties(report=True)

    # report individual mass properties
    ex_1.get_mass_properties(report=True,individual=True)

    # report total mass properties with negative tensor definition
    ex_1.get_mass_properties(report=True,positive_tensor=False)

    # simple foam wings case study
    sfw = mmoi("simple_foam_wings.json")
    sfw.get_mass_properties(report=True,individual=True)

    # rotor case study
    rotor = mmoi("rotor.json")
    rotor.get_mass_properties(report=True)

    # CRM case study
    CRM = mmoi("CRM.json")
    CRM.get_mass_properties(report=True)

    # Horizon case study
    hzn = mmoi("horizon.json")
    hzn.get_mass_properties(report=True)
