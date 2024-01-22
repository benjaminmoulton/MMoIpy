import sys
mmoipy_directory = '../mmoipy/'
docs_directory = '../docs/'
exampl_directory = '../examples/'
sys.path.insert(1, mmoipy_directory)
sys.path.insert(1, docs_directory)
sys.path.insert(1, exampl_directory)
from system import AircraftSystem as mmoi
from scipy.optimize import curve_fit
import numpy as np
import json
from matplotlib import pyplot as plt


def sinusoid(x,A,w,p,z):
    return A*np.sin(w*x + p) + z

def abs_sinusoid(x,A,w,p,z):
    return A*np.abs(np.sin(w*x + p)) + z


if __name__ == "__main__":
    # # RC f16
    # f16 = mmoi("RC_F16.json")
    # f16.get_mass_properties(report=True,mass_multiplier=32.1740485564305,length_multiplier=12.)#,individual=True)
    # f16.visualize()
    # RC f16 without ventral fins or arms mounts
    f16 = mmoi("RC_F16_wo_ventral_and_mounts.json")
    If16 = f16.get_mass_properties(report=True)["inertia_tensor"]#,mass_multiplier=32.1740485564305,length_multiplier=12.)#,individual=True)
    # f16.visualize()
    # RC BIRE
    bire = mmoi("RC_BIRE.json")
    bire.get_mass_properties(report=True)#,individual=True)#,mass_multiplier=32.1740485564305,length_multiplier=12.)#,individual=True)
    # bire.visualize()
    # # RC BIRE
    # bise = mmoi("RC_BIRE_stab_extended.json")
    # bise.get_mass_properties(report=True)#,individual=True)#,mass_multiplier=32.1740485564305,length_multiplier=12.)#,individual=True)
    # # bise.visualize()#show_legend=True)

    # run through BIRE angles by 5 degrees
    dBnum = 37
    dB = np.linspace(-90.,90.,dBnum)
    pieces_to_rotate = [
        "l_stab","r_stab",
        "L_fuselage_side_aft","R_fuselage_side_aft"
    ]
    filename = "RC_BIRE.json"
    g = 32.1740485564305
    input_dict = json.loads(open(filename).read())
    lf = abs(input_dict["components"]["R_fuselage_side_aft"]["connect_to"]["dy"])
    ls = abs(input_dict["components"]["r_stab"]["connect_to"]["dy"])
    Is = np.zeros((6,dBnum))
    for i in range(len(dB)):
        # rotate tail
        # print("running tail rotation {:> 5.1f} deg".format(dB[i]))
        for part in pieces_to_rotate:
            # # which side
            # if part[0] in ["l","L"]:
            #     dBi = dB[i]*1.
            # else:
            #     dBi = dB[i]*-1.
            if part[0] in ["l","L"]:
                si = -1.
            else:
                si = 1.
            dBi = dB[i]*1.
            if part[-4:] == "stab":
                input_dict["components"][part]["dihedral"] = -si*dBi*1.
                di = ls*1.
            else:
                input_dict["components"][part]["orientation"] = [dBi*1.,0.0,0.0]
                di = lf*1.
            dBi = np.deg2rad(dBi)
            input_dict["components"][part]["connect_to"]["dy"] = si*di*np.cos(dBi)
            input_dict["components"][part]["connect_to"]["dz"] = si*di*np.sin(dBi)
        
        # initialize aircraft
        craft = mmoi(input_dict,verbose=False)
        I = craft.get_mass_properties()["inertia_tensor"]
        # craft.visualize(view=(0.0,0.0))
        Is[0,i] =  I[0,0]
        Is[1,i] =  I[1,1]
        Is[2,i] =  I[2,2]
        Is[3,i] = -I[0,1]
        Is[4,i] = -I[0,2]
        Is[5,i] = -I[1,2]
    
    # determine fit
    dB_rad = np.deg2rad(dB)
    mid = int( (dBnum - 1)/2)
    p0 = [1.0,2.0,0.0,Is[1,mid]]
    Ixx_fit = [0.0,0.0,0.0,Is[0,0]]
    A = abs(np.max(Is[1,:]) - np.min(Is[1,:]))
    Iyy_fit = curve_fit(sinusoid,dB_rad,Is[1,:],p0=[A,2.,np.pi/2.,Is[1,mid]])[0]
    # print(Iyy_fit)
    A = abs(np.max(Is[2,:]) - np.min(Is[2,:]))
    Izz_fit = curve_fit(sinusoid,dB_rad,Is[2,:],p0=[-A,2.,np.pi/2.,Is[2,mid]])[0]
    Ixy_fit = [0.0,0.0,0.0,0.0]
    Ixz_fit = [0.0,0.0,0.0,Is[4,0]]
    A = abs(np.max(Is[5,:]) - np.min(Is[5,:]))
    Iyz_fit = curve_fit(sinusoid,dB_rad,Is[5,:],p0=[A,2.,0.,Is[5,mid]])[0]
    # print(Ixx_fit)
    # print(Iyy_fit)

    fits = [Ixx_fit,Iyy_fit,Izz_fit,Ixy_fit,Ixz_fit,Iyz_fit]
    funs = [sinusoid]*6

    # for i in range(6):
    #     plt.plot(dB,Is[i,:],"bo",lw=0.0)
    #     plt.plot(dB,funs[i](dB_rad,*fits[i]),"b")
    #     plt.show()

    # report
    names = ["$I_{xx_b}$","$I_{yy_b}$","$I_{zz_b}$",
    "$I_{xy_b}$","$I_{xz_b}$","$I_{yz_b}$"]
    print("{:<10s} & {:^11s} & {:^11s} & {:^11s} & {:^11s} \\\\".format(
        "","A","w","phi","z"
    ))
    print("BIRE Inertia ....")
    for j in range(6):
        print(("{:<10s}").format(names[j]),end="")
        for k in range(len(fits[j])):
            if k == 1:
                print(" & {:< 11.0f}".format(fits[j][k]),end="")
            else:
                if abs(fits[j][k]) <= 1e-15:
                    print(" & {:> 11.0f}".format(fits[j][k]),end="")
                else:
                    if k == 2:
                        print(" & {:< 11.4f}".format(fits[j][k]),end="")
                    else:
                        print(" & {:> 11.4e}".format(fits[j][k]),end="")
        print(" \\\\")
    print()

    print("F16  Inertia...")
    inds = [[0,0],[1,1],[2,2],[0,1],[0,2],[1,2]]
    # print(If16)
    for j in range(6):
        print(("{:<10s} [slugs-ft$^2$]").format(names[j]),end="")
        I = If16[inds[j][0],inds[j][1]]*1.
        if j > 2:
            I = -I
        if abs(I) <= 1e-15:
            print(" & {:> 11.0f}".format(I),end="")
        else:
            print(" & {:> 11.4e}".format(I),end="")
        print(" \\\\")



    # output
    #         _|      _|  _|      _|            _|_|_|                      
    #         _|_|  _|_|  _|_|  _|_|    _|_|      _|    _|_|_|    _|    _|  
    #         _|  _|  _|  _|  _|  _|  _|    _|    _|    _|    _|  _|    _|  
    #         _|      _|  _|      _|  _|    _|    _|    _|    _|  _|    _|  
    #         _|      _|  _|      _|    _|_|    _|_|_|  _|_|_|      _|_|_|  
    #                                                 _|              _|  
    #                                                 _|          _|_|    
            
    # TOTAL
    # =================================================================================================
    # Mass properties of RC_F16_wo_ventral_and_mounts.json

    # Mass =  0.22785334

    # Volume =  0.85382131

    # Center of mass:
    # 	X =    -2.36609626
    # 	Y =     0.00000000
    # 	Z =    -0.03192711

    # Angular momentum:
    # 	X =     0.00315636
    # 	Y =     0.00000000
    # 	Z =     0.00000000

    # Moment of inertia about the CG aligned with system coordinate frame:
    # 		Positive Tensor Formulation
    # 	Ixx =          0.02515302	Ixy =         -0.00000000	Ixz =         -0.00039819
    # 	Iyx =         -0.00000000	Iyy =          0.26480171	Iyz =          0.00000000
    # 	Izx =         -0.00039819	Izy =          0.00000000	Izz =          0.28230872

    # Moment of inertia about the origin aligned with system coordinate frame:
    # 		Positive Tensor Formulation
    # 	Ixx =          0.02538528	Ixy =         -0.00000000	Ixz =          0.01681445
    # 	Iyx =         -0.00000000	Iyy =          1.54065072	Iyz =         -0.00000000
    # 	Izx =          0.01681445	Izy =         -0.00000000	Izz =          1.55792547

    # =================================================================================================
                                                                
    #         _|      _|  _|      _|            _|_|_|                      
    #         _|_|  _|_|  _|_|  _|_|    _|_|      _|    _|_|_|    _|    _|  
    #         _|  _|  _|  _|  _|  _|  _|    _|    _|    _|    _|  _|    _|  
    #         _|      _|  _|      _|  _|    _|    _|    _|    _|  _|    _|  
    #         _|      _|  _|      _|    _|_|    _|_|_|  _|_|_|      _|_|_|  
    #                                                 _|              _|  
    #                                                 _|          _|_|    
            
    # TOTAL
    # =================================================================================================
    # Mass properties of RC_BIRE.json

    # Mass =  0.24783442

    # Volume =  0.84319231

    # Center of mass:
    # 	X =    -2.36020175
    # 	Y =     0.00000000
    # 	Z =    -0.03903927

    # Angular momentum:
    # 	X =     0.00315636
    # 	Y =     0.00000000
    # 	Z =     0.00000000

    # Moment of inertia about the CG aligned with system coordinate frame:
    # 		Positive Tensor Formulation
    # 	Ixx =          0.02460001	Ixy =         -0.00000000	Ixz =         -0.00052991
    # 	Iyx =         -0.00000000	Iyy =          0.28286493	Iyz =          0.00000000
    # 	Izx =         -0.00052991	Izy =          0.00000000	Izz =          0.30099468

    # Moment of inertia about the origin aligned with system coordinate frame:
    # 		Positive Tensor Formulation
    # 	Ixx =          0.02497772	Ixy =         -0.00000000	Ixz =          0.02230569
    # 	Iyx =         -0.00000000	Iyy =          1.66381726	Iyz =         -0.00000000
    # 	Izx =          0.02230569	Izy =         -0.00000000	Izz =          1.68156929

    # =================================================================================================
    #            &      A      &      w      &     phi     &      z      \\
    # BIRE Inertia ....
    # $I_{xx_b}$ &           0 &  0          &           0 &  2.4600e-02 \\
    # $I_{yy_b}$ & -1.0091e-03 &  2          &  1.5708     &  2.8387e-01 \\
    # $I_{zz_b}$ &  1.0091e-03 &  2          &  1.5708     &  2.9999e-01 \\
    # $I_{xy_b}$ &           0 &  0          &           0 &           0 \\
    # $I_{xz_b}$ &           0 &  0          &           0 & -5.2991e-04 \\
    # $I_{yz_b}$ &  1.0091e-03 &  2          &          -0 &           0 \\

    # F16  Inertia...
    # $I_{xx_b}$ [slugs-ft$^2$] &  2.5153e-02 \\
    # $I_{yy_b}$ [slugs-ft$^2$] &  2.6480e-01 \\
    # $I_{zz_b}$ [slugs-ft$^2$] &  2.8231e-01 \\
    # $I_{xy_b}$ [slugs-ft$^2$] &          -0 \\
    # $I_{xz_b}$ [slugs-ft$^2$] & -3.9819e-04 \\
    # $I_{yz_b}$ [slugs-ft$^2$] &           0 \\