import sys
mmoipy_directory = '../mmoipy/'
docs_directory = '../docs/'
exampl_directory = '../examples/'
sys.path.insert(1, mmoipy_directory)
sys.path.insert(1, docs_directory)
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

    # # test new cylinders
    # ex_2 = mmoi("dev_cylinder_input.json")
    # ex_2.get_mass_properties(report=True,individual=True)
    # # # ex_2.visualize()

    # # test ellipsoid
    # ex_3 = mmoi("dev_ellipsoid_input.json")
    # ex_3.get_mass_properties(report=True,individual=True)
    # ex_3.visualize(plot_ids=None)

    # save plots of objects
    ex_4 = mmoi("dev_input.json")
    keys = ex_4.components.keys()
    for i in keys:
        ex_4.visualize(no_color=True,plot_ids=[i],
            filename=docs_directory + ex_4.components[i].name + ".png")
