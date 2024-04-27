import numpy as np
from matplotlib import pyplot as plt
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
        "mass": 0.3,
        "type" : "symmetric_airfoil",
        "side" : "right",
        "connect_to" : {
            "dx" :  2.0,
            "dy" :  3.0,
            "dz" :  -1.0
        },
        "semispan" : 4.0,
        "sweep" : 10.0,
        "dihedral" : -5.0,
        "chord" : [ [0.0, 1.0],
                    [1.0, 0.5]],
        "thickness" : [ [0.0, 0.08],
                        [1.0, 0.10]]
    }
    lan_wing = {**sym_wing}
    lan_wing["ID"] = 2
    lan_wing["type"] = "lanham_wing"
    prm_wing = {**sym_wing}
    prm_wing["ID"] = 3
    prm_wing["type"] = "prismoid"
    inp_dict = {
        "components" : {
            "sym" : sym_wing,
            "lan" : lan_wing,
            # "prsm" : prm_wing
        }
    }
    sym_dict = {
        "components" : {
            "sym" : sym_wing
        }
    }
    lan_dict = {
        "components" : {
            "lan" : lan_wing
        }
    }

    # run code
    test = mmoi(inp_dict)
    test.get_mass_properties(report=True,individual=True)
    # test.visualize(plot_ids=None)

    # run for various geometries
    num = 100
    ms  = np.zeros((2,num,))
    Vs  = np.zeros((2,num,))
    cgs = np.zeros((2,3,num))
    Is  = np.zeros((2,3,3,num))
    Ios = np.zeros((2,3,3,num))
    RT = np.linspace(0.0,1.0,num)
    for i in range(num):
        prop = "chord"
        ct = RT[i]*sym_dict["components"]["sym"][prop][0][1]
        sym_dict["components"]["sym"][prop][1][1] = \
            lan_dict["components"]["lan"][prop][1][1] = ct
        lan = mmoi(lan_dict,verbose=False)
        sym = mmoi(sym_dict,verbose=False)
        lan_props = lan.get_mass_properties()# individual=True)
        sym_props = sym.get_mass_properties()# individual=True)
        # print(lan_props,sym_props)
        ms[0,i] = sym_props["mass"]
        ms[1,i] = lan_props["mass"]
        Vs[0,i] = sym_props["volume"]
        Vs[1,i] = lan_props["volume"]
        cgs[0,:,i] = sym_props["cg_location"][:,0]
        cgs[1,:,i] = lan_props["cg_location"][:,0]
        Is [0,:,:,i] = sym_props["inertia_tensor"]
        Is [1,:,:,i] = lan_props["inertia_tensor"]
        Ios[0,:,:,i] = sym_props["origin_inertia_tensor"]
        Ios[1,:,:,i] = lan_props["origin_inertia_tensor"]
    
    
    # plot
    # change plot text parameters
    plt.rcParams["font.family"] = "Serif"
    plt.rcParams["font.size"] = 8.0
    plt.rcParams["axes.labelsize"] = 8.0
    plt.rcParams['lines.linewidth'] = 1.0
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams["xtick.direction"] = plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.bottom"] = plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.left"] = plt.rcParams["ytick.right"] = True
    plt.rcParams["xtick.major.width"] = plt.rcParams["ytick.major.width"] = 0.75
    plt.rcParams["xtick.minor.width"] = plt.rcParams["ytick.minor.width"] = 0.75
    plt.rcParams["xtick.major.size"] = plt.rcParams["ytick.major.size"] = 5.0
    plt.rcParams["xtick.minor.size"] = plt.rcParams["ytick.minor.size"] = 2.5
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    plt.rcParams['figure.dpi'] = 300.0
    # initialize plots
    plot_dict = dict(
        figsize = (3.25,3.5),
        constrained_layout = True,
        sharex=True
    )
    fig_mV ,axs_mV  = plt.subplots(2,1,**plot_dict)
    fig_cgs,axs_cgs = plt.subplots(3,1,**plot_dict)
    fig_Ims,axs_Ims = plt.subplots(3,1,**plot_dict)
    fig_Ips,axs_Ips = plt.subplots(3,1,**plot_dict)

    #   mass & volume
    axs_mV [0].grid(which="major",lw=0.6,ls="-",c="0.75")
    axs_mV [1].grid(which="major",lw=0.6,ls="-",c="0.75")
    axs_mV [0].plot(RT,ms[0],"k-",label="Analytic")
    axs_mV [0].plot(RT,ms[1],"k--",label="Lanham")
    axs_mV [1].plot(RT,Vs[0],"k-")
    axs_mV [1].plot(RT,Vs[1],"k--")
    axs_mV [0].set_ylabel(r"mass, $m$")
    axs_mV [1].set_ylabel(r"volume, $V$")
    axs_mV [1].set_xlabel("Taper Ratio, $c_t/c_r$")
    axs_mV [1].set_xlim((RT[0],RT[-1]))
    axs_mV [0].legend()
    #   cgs
    axs_cgs[0].grid(which="major",lw=0.6,ls="-",c="0.75")
    axs_cgs[1].grid(which="major",lw=0.6,ls="-",c="0.75")
    axs_cgs[2].grid(which="major",lw=0.6,ls="-",c="0.75")
    axs_cgs[0].plot(RT,cgs[0,0],"k-",label="Analytic")
    axs_cgs[0].plot(RT,cgs[1,0],"k--",label="Lanham")
    axs_cgs[1].plot(RT,cgs[0,1],"k-")
    axs_cgs[1].plot(RT,cgs[1,1],"k--")
    axs_cgs[2].plot(RT,cgs[0,2],"k-")
    axs_cgs[2].plot(RT,cgs[1,2],"k--")
    axs_cgs[0].set_ylabel(r"$\bar{x}$")
    axs_cgs[1].set_ylabel(r"$\bar{y}$")
    axs_cgs[2].set_ylabel(r"$\bar{z}$")
    axs_cgs[2].set_xlabel("Taper Ratio, $c_t/c_r$")
    axs_cgs[2].set_xlim((RT[0],RT[-1]))
    axs_cgs[0].legend()
    #   Ims
    axs_Ims[0].grid(which="major",lw=0.6,ls="-",c="0.75")
    axs_Ims[1].grid(which="major",lw=0.6,ls="-",c="0.75")
    axs_Ims[2].grid(which="major",lw=0.6,ls="-",c="0.75")
    axs_Ims[0].plot(RT,Is [0,0,0],"k-",label="Analytic")
    axs_Ims[0].plot(RT,Is [1,0,0],"k--",label="Lanham")
    axs_Ims[1].plot(RT,Is [0,1,1],"k-")
    axs_Ims[1].plot(RT,Is [1,1,1],"k--")
    axs_Ims[2].plot(RT,Is [0,2,2],"k-")
    axs_Ims[2].plot(RT,Is [1,2,2],"k--")
    # axs_Ims[0].plot(RT,Ios[0,0,0],"0.4",ls="-",label="origin")
    # axs_Ims[0].plot(RT,Ios[1,0,0],"0.4",ls="--")
    # axs_Ims[1].plot(RT,Ios[0,1,1],"0.4",ls="-")
    # axs_Ims[1].plot(RT,Ios[1,1,1],"0.4",ls="--")
    # axs_Ims[2].plot(RT,Ios[0,2,2],"0.4",ls="-")
    # axs_Ims[2].plot(RT,Ios[1,2,2],"0.4",ls="--")
    axs_Ims[0].set_ylabel(r"$I_{xx}$")
    axs_Ims[1].set_ylabel(r"$I_{yy}$")
    axs_Ims[2].set_ylabel(r"$I_{zz}$")
    axs_Ims[2].set_xlabel("Taper Ratio, $c_t/c_r$")
    axs_Ims[2].set_xlim((RT[0],RT[-1]))
    axs_Ims[0].legend()
    #   Ips
    axs_Ips[0].grid(which="major",lw=0.6,ls="-",c="0.75")
    axs_Ips[1].grid(which="major",lw=0.6,ls="-",c="0.75")
    axs_Ips[2].grid(which="major",lw=0.6,ls="-",c="0.75")
    axs_Ips[0].plot(RT,-Is [0,0,1],"k-",label="Analytic")
    axs_Ips[0].plot(RT,-Is [1,0,1],"k--",label="Lanham")
    axs_Ips[1].plot(RT,-Is [0,0,2],"k-")
    axs_Ips[1].plot(RT,-Is [1,0,2],"k--")
    axs_Ips[2].plot(RT,-Is [0,1,2],"k-")
    axs_Ips[2].plot(RT,-Is [1,1,2],"k--")
    # axs_Ips[0].plot(RT,-Ios[0,0,1],"0.4",ls="-",label="origin")
    # axs_Ips[0].plot(RT,-Ios[1,0,1],"0.4",ls="--")
    # axs_Ips[1].plot(RT,-Ios[0,0,2],"0.4",ls="-")
    # axs_Ips[1].plot(RT,-Ios[1,0,2],"0.4",ls="--")
    # axs_Ips[2].plot(RT,-Ios[0,1,2],"0.4",ls="-")
    # axs_Ips[2].plot(RT,-Ios[1,1,2],"0.4",ls="--")
    axs_Ips[0].set_ylabel(r"$I_{xy}$")
    axs_Ips[1].set_ylabel(r"$I_{xz}$")
    axs_Ips[2].set_ylabel(r"$I_{yz}$")
    axs_Ips[2].set_xlabel("Taper Ratio, $c_t/c_r$")
    axs_Ips[2].set_xlim((RT[0],RT[-1]))
    axs_Ips[0].legend()

    # set name
    save_figs = True
    if save_figs:
        file_end = "png" # "pdf" # 
        fig_mV .savefig("comp_mVs."+file_end,dpi=300.,transparent=True)
        fig_cgs.savefig("comp_cgs."+file_end,dpi=300.,transparent=True)
        fig_Ims.savefig("comp_Ims."+file_end,dpi=300.,transparent=True)
        fig_Ips.savefig("comp_Ips."+file_end,dpi=300.,transparent=True)
    show_plots = False
    if show_plots:
        plt.show()
    else:
        plt.close()