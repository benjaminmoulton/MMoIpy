import sys
mmoipy_directory = '../mmoipy/'
docs_directory = '../docs/'
exampl_directory = '../examples/'
sys.path.insert(1, mmoipy_directory)
sys.path.insert(1, docs_directory)
sys.path.insert(1, exampl_directory)
from system import AircraftSystem as mmoi

if __name__ == "__main__":
    # RC f16
    f16 = mmoi("RC_F16.json")
    f16.get_mass_properties(report=True,individual=True)#mass_multiplier=32.1740485564305,length_multiplier=12.,individual=True)
    # f16.visualize()