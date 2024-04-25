import sys
mmoipy_directory = '../mmoipy/'
docs_directory = '../docs/'
exampl_directory = '../examples/'
sys.path.insert(1, mmoipy_directory)
sys.path.insert(1, docs_directory)
sys.path.insert(1, exampl_directory)
from system import AircraftSystem as mmoi

if __name__ == "__main__":
    # input file
    # symmetric wing
    sym_wing = {
        "ID" : 1,
        "density": 0.3,
        "type" : "symmetric_airfoil",
        "side" : "right",
        "connect_to" : {
            "dx" :  2.0,
            "dy" :  3.0,
            "dz" :  -1.0
        },
        "semispan" : 4.0,
        "sweep" : 10.0,
        "chord" : [ [0.0, 1.0],
                    [1.0, 0.5]],
        "thickness" : [ [0.0, 0.08],
                        [1.0, 0.10]]
    }
    lan_wing = {**sym_wing}
    lan_wing["ID"] = 2
    lan_wing["type"] = "lanham_wing"
    inp_dict = {
        "components" : {
            "sym" : sym_wing,
            "lan" : lan_wing
        }
    }

    # run code
    test = mmoi(inp_dict)
    test.get_mass_properties(report=True,individual=True)
    test.visualize(plot_ids=None)

    # # save plots of objects
    # ex_4 = mmoi("dev_input.json")
    # keys = ex_4.components.keys()
    # for i in keys:
    #     ex_4.visualize(no_color=True,plot_ids=[i],
    #         filename=docs_directory + ex_4.components[i].name + ".png")
