import sys
mmoipy_directory = '../mmoipy/'
exampl_directory = '../examples/'
sys.path.insert(1, mmoipy_directory)
sys.path.insert(1, exampl_directory)
from system import AircraftSystem as mmoi

if __name__ == "__main__":
    # # test several inputs together
    # ex_1 = mmoi("dev_input.json")

    # # report total mass properties
    # ex_1.get_mass_properties(report=True,individual=True)

    # # visualize
    # ex_1.visualize()

    # # horizon
    # hzn = mmoi(exampl_directory + "horizon.json")
    # # hzn.get_mass_properties(report=True)
    # hzn.visualize()

    # # CRM
    # CRM = mmoi(exampl_directory + "CRM.json")
    # CRM.visualize()

    # test new cylinders
    ex_2 = mmoi("dev_cylinder_input.json")
    ex_2.get_mass_properties(report=True,individual=True)
    ex_2.visualize()
