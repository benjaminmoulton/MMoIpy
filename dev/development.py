import sys
mmoipy_directory = '../mmoipy/'
exampl_directory = '../examples/'
sys.path.insert(1, mmoipy_directory)
sys.path.insert(1, exampl_directory)
from system import AircraftSystem as mmoi

if __name__ == "__main__":
    # test several inputs together
    ex_1 = mmoi("dev_input.json")

    # # report total mass properties
    # ex_1.get_mass_properties(report=True)

    # visualize
    ex_1.visualize()

    # # horizon
    # hzn = mmoi("horizon.json")
    # hzn.visualize()
