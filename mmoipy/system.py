import numpy as np
import json
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as ax3

from wing import Wing
from component import Component, Cuboid, Cylinder, HalfCylinder, Sphere
from component import Ellipsoid, HalfEllipsoid, Rotor

class AircraftSystem:
    """A class calculating and containing the mass properties of an aircraft.

    Parameters
    ----------
    input_vars : string or dict , optional
        Must be a .json file or python dictionary.

    Raises
    ------
    TypeError
        If the input_vars type is not a dictionary or the file path to 
        a .json file
    """
    def __init__(self,input_vars={},verbose=True):

        # print cool logo
        print("""                                                               
        _|      _|  _|      _|            _|_|_|                      
        _|_|  _|_|  _|_|  _|_|    _|_|      _|    _|_|_|    _|    _|  
        _|  _|  _|  _|  _|  _|  _|    _|    _|    _|    _|  _|    _|  
        _|      _|  _|      _|  _|    _|    _|    _|    _|  _|    _|  
        _|      _|  _|      _|    _|_|    _|_|_|  _|_|_|      _|_|_|  
                                                _|              _|  
                                                _|          _|_|    
        """)

        # get info or raise error
        self._get_input_vars(input_vars)

        # retrieve info
        self._initialize_objects()


    def _get_input_vars(self,input_vars):
        # get info or raise error

        # determine if the input_vars is a file or a dictionary
        input_vars_type = type(input_vars)

        # dictionary
        if input_vars_type == dict:
            self.input_dict = input_vars
            self.file_name = "Aircraft"
        
        # json file
        elif input_vars_type == str and input_vars.split(".")[-1] == "json":
            self.file_name = input_vars
            self.input_dict = self._get_json(input_vars)

        # raise error
        else:
            raise TypeError("input_vars must be json file path, or " + \
                "dictionary, not {0}".format(input_vars_type))
    

    def _get_json(self,file_path):
        # import json file from file path
        json_string = open(file_path).read()

        # save to vals dictionary
        input_dict = json.loads(json_string)
        
        return input_dict


    def _initialize_objects(self):

        # store component input values
        components = self.input_dict.get("components",{})

        # check if mass in dictionary
        given_total_mass    = "mass"    in self.input_dict
        given_total_density = "density" in self.input_dict

        # check if we are given english units
        self.given_english_units = self.input_dict.get("units","English") == \
            "English"

        if given_total_mass:
            self.mass = self.input_dict.get("mass")
        elif given_total_density:
            self.density = self.input_dict.get("density")

        wing_types = ["prismoid","pseudo_prismoid","symmetric_airfoil",
        "diamond_airfoil"]
        comp_types = ["cuboid","cylinder","half_cylinder","sphere",
        "ellipsoid","half_ellipsoid","rotor"]

        # initialize get each id number
        ids = []; attach_ids = []; name = []
        for component in components:
            input_dict = components[component]
            ids.append(input_dict.get("ID"))
            input_dict["connect_to"] = input_dict.get("connect_to",{})
            attach_ids.append(input_dict["connect_to"].get("ID",0))
            name.append(component)
        
        # throw error if an id given is repeated
        unique_ids,counts = np.unique(ids,return_counts=True)
        repeated_ids = unique_ids[counts > 1]
        if len(unique_ids) != len(ids):
            raise TypeError("Repeated ID values: {}".format(repeated_ids))
               
        # reorganize attach_ids in order of run
        run_order = np.argsort(attach_ids).tolist()
        name_order = np.array(name)[run_order].tolist()
        attach_id_order = np.array(attach_ids)[run_order].tolist()

        # initialize components, save zero component
        self.components = {}
        self.components[0] = Component({"mass":0.0})
        self.components[0].type = "origin"


        for i in range(len(name_order)):

            component = name_order[i]

            # define component input dictionary
            input_dict = components[component]
            id_number = input_dict.get("ID")

            # overwrite mass and or densities if given
            if given_total_mass:
                input_dict["density"] = 1.0
                if "mass" in input_dict: input_dict.pop("mass")
            elif given_total_density:
                input_dict["density"] = self.density
                if "mass" in input_dict: input_dict.pop("mass")
            
            # shift to root or tip values if attached to a wing
            input_dict["connect_to"] = input_dict.get("connect_to",{})
            attach_loc = input_dict["connect_to"].get("location","tip")
            loc = self.components[attach_id_order[i]].locations[attach_loc]
            input_dict["connect_to"]["dx"] = \
                input_dict["connect_to"].get("dx",0.0) + loc[0,0]
            input_dict["connect_to"]["y_offset"] = \
                input_dict["connect_to"].get("y_offset",0.0) + loc[1,0]
            input_dict["connect_to"]["dz"] = \
                input_dict["connect_to"].get("dz",0.0) + loc[2,0]


            # initialize components
            if input_dict["type"] in wing_types:
                self.components[id_number] = Wing(input_dict)
            elif input_dict["type"] == "cuboid":
                self.components[id_number] = Cuboid(input_dict)
            elif input_dict["type"] == "cylinder":
                self.components[id_number] = Cylinder(input_dict)
            elif input_dict["type"] == "half_cylinder":
                self.components[id_number] = HalfCylinder(input_dict)
            elif input_dict["type"] == "sphere":
                self.components[id_number] = Sphere(input_dict)
            elif input_dict["type"] == "ellipsoid":
                self.components[id_number] = Ellipsoid(input_dict)
            elif input_dict["type"] == "half_ellipsoid":
                self.components[id_number] = HalfEllipsoid(input_dict)
            elif input_dict["type"] == "rotor":
                self.components[id_number] = Rotor(input_dict)
            
            # save "name"
            if input_dict["type"] in wing_types + comp_types:
                self.components[id_number].name = component
            self.components[id_number].type = input_dict["type"]
        
        # calculate total volume
        self.volume = 0.0
        for comp_id in self.components:
            self.volume += self.components[comp_id].volume
        
        # rewrite densities if given total mass
        if given_total_mass:
            # calculate total density, apply to each
            self.density = self.mass / self.volume
            for comp_id in self.components:
                self.components[comp_id].density = self.density * 1.
                self.components[comp_id].update_densities()


    def report_as_SolidWorks_report(self,info,positive_tensor=True,
        use_Lanham=False,name="",report_units=False):
        """Method which reports the mass and inertia properties as given in 
        SolidWorks.
        
        Parameters
        ----------
        info : dictionary
            The dictionary with 'mass', 'cg_location' and 'inertia_tensor'
            information.
        
        positive_tensor : boolean, optional
            Whether to report in positive tensor formulation. Defaults to true.
        """

        if name == "":
            partin = ""
        else:
            partin = name + " in "
        if use_Lanham:
            symbol = "-*"
            intro = "Lanham "
        else:
            symbol = "=="
            intro = ""
        print((symbol * 50)[:-3])
        print(intro+"Mass properties of",partin + self.file_name)
        print()

        # fix units if not given english units
        if not self.given_english_units:
            info["mass"] /= 14.59390
            info["cg_location"] /= 0.3048
            info["volume"] /= (0.3048**3.)
            info["inertia_tensor"] /= 14.59390 * (0.3048)**2.

        if report_units: units = " slugs"
        else: units = ""
        print(intro+"Mass = {:> 10.8f}{}".format(info["mass"],units))
        print()

        if report_units: units = " cubic feet"
        else: units = ""
        print(intro+"Volume = {:> 10.8f}{}".format(info["volume"],units))
        print()

        if report_units: units = " (feet)"
        else: units = ""
        print(intro+"Center of mass:{}".format(units))
        print("\tX = {:> 14.8f}".format(info["cg_location"][0,0]))
        print("\tY = {:> 14.8f}".format(info["cg_location"][1,0]))
        print("\tZ = {:> 14.8f}".format(info["cg_location"][2,0]))
        print()
        
        if report_units: units = " (slugs * square feet / seconds)"
        else: units = ""
        print(intro+"Angular momentum:{}".format(units))
        print("\tX = {:> 14.8f}".format(info["angular_momentum"][0,0]))
        print("\tY = {:> 14.8f}".format(info["angular_momentum"][1,0]))
        print("\tZ = {:> 14.8f}".format(info["angular_momentum"][2,0]))
        print()
        
        I = info["inertia_tensor"] * 1.0
        if positive_tensor:
            I[[0,0,1,1,2,2],[1,2,0,2,0,1]] *= -1.0
        [[Ixx,Ixy,Ixz],[Iyx,Iyy,Iyz],[Izx,Izy,Izz]] = I
        if report_units: units = "( slugs * square feet)"
        else: units = ""
        print(intro+"Moment of inertia about the CG" \
            + " aligned with system coordinate frame:{}".format(units))
        if positive_tensor:
            print("\t\tPositive Tensor Formulation")
        else:
            print("\t\tNegative Tensor Formulation")
        print("\tIxx = {:> 19.8f}\tIxy = {:> 19.8f}\tIxz = {:> 19.8f}".format(\
            Ixx,Ixy,Ixz))
        print("\tIyx = {:> 19.8f}\tIyy = {:> 19.8f}\tIyz = {:> 19.8f}".format(\
            Iyx,Iyy,Iyz))
        print("\tIzx = {:> 19.8f}\tIzy = {:> 19.8f}\tIzz = {:> 19.8f}".format(\
            Izx,Izy,Izz))
        print()
        
        I = info["origin_inertia_tensor"] * 1.0
        if positive_tensor:
            I[[0,0,1,1,2,2],[1,2,0,2,0,1]] *= -1.0
        [[Ixx,Ixy,Ixz],[Iyx,Iyy,Iyz],[Izx,Izy,Izz]] = I
        print(intro+"Moment of inertia about the origin" \
            + " aligned with system coordinate frame:{}".format(units))
        if positive_tensor:
            print("\t\tPositive Tensor Formulation")
        else:
            print("\t\tNegative Tensor Formulation")
        print("\tIxx = {:> 19.8f}\tIxy = {:> 19.8f}\tIxz = {:> 19.8f}".format(\
            Ixx,Ixy,Ixz))
        print("\tIyx = {:> 19.8f}\tIyy = {:> 19.8f}\tIyz = {:> 19.8f}".format(\
            Iyx,Iyy,Iyz))
        print("\tIzx = {:> 19.8f}\tIzy = {:> 19.8f}\tIzz = {:> 19.8f}".format(\
            Izx,Izy,Izz))
        print()
        print((symbol * 50)[:-3])


    def get_mass_properties(self,report=False,individual=False,use_Lanham=False,positive_tensor=True):

        # determine properties of each component
        for i in self.components:
            self.components[i].get_mass_properties()
        
        # determine total mass
        self.mass = 0.0
        for i in self.components:
            self.mass += self.components[i].get_mass()
        
        # determine lanham mass
        self.mass_lanham = 0.0
        for i in self.components:
            self.mass_lanham += self.components[i].get_mass(True)
        
        # determine lanham volume
        self.volume_lanham = 0.0
        for i in self.components:
            self.volume_lanham += self.components[i].get_volume(True)
        
        # determine total cg location
        self.cg_location = np.zeros((3,1))
        for i in self.components:
            self.cg_location += self.components[i].mass * \
                self.components[i].get_cg_location()
        self.cg_location /= self.mass
        
        # determine lanham cg location
        self.cg_location_lanham = np.zeros((3,1))
        for i in self.components:
            self.cg_location_lanham += self.components[i].get_mass(True) * \
                self.components[i].get_cg_location(True)
        self.cg_location_lanham /= self.mass_lanham
        
        # determine total angular momentum
        self.angular_momentum = 0.0
        for i in self.components:
            self.angular_momentum += self.components[i].angular_momentum
        
        # determine total inertia
        self.inertia_tensor = np.zeros((3,3))
        location = self.cg_location
        for i in self.components:
            self.inertia_tensor += \
                self.components[i].shift_properties_to_location(\
                location)["inertia_tensor"]
        
        # determine lanham inertia
        self.inertia_tensor_lanham = np.zeros((3,3))
        location = self.cg_location_lanham
        for i in self.components:
            self.inertia_tensor_lanham += \
                self.components[i].shift_properties_to_location(\
                location,True)["inertia_tensor"]
        
        # determine origin inertia
        self.origin_inertia_tensor = np.zeros((3,3))
        location = np.zeros((3,1))
        for i in self.components:
            self.origin_inertia_tensor += \
                self.components[i].shift_properties_to_location(\
                location)["inertia_tensor"]
        
        # determine lanham origin inertia
        self.origin_inertia_tensor_lanham = np.zeros((3,3))
        location = np.zeros((3,1))
        for i in self.components:
            self.origin_inertia_tensor_lanham += \
                self.components[i].shift_properties_to_location(\
                location,True)["inertia_tensor"]

        # return dictionary of values
        self.properties_dict = {
            "mass" : self.mass,
            "volume" : self.volume,
            "cg_location" : self.cg_location,
            "angular_momentum" : self.angular_momentum,
            "origin_inertia_tensor" : self.origin_inertia_tensor,
            "inertia_tensor" : self.inertia_tensor
        }

        # Lanham dictionary of values
        lanham = {
            "mass" : self.mass_lanham,
            "volume" : self.volume_lanham,
            "cg_location" : self.cg_location_lanham,
            "angular_momentum" : self.angular_momentum,
            "origin_inertia_tensor" : self.origin_inertia_tensor_lanham,
            "inertia_tensor" : self.inertia_tensor_lanham
        }

        # report
        if report:
            print("TOTAL")
            self.report_as_SolidWorks_report(self.properties_dict,\
                positive_tensor)
            if use_Lanham:
                self.report_as_SolidWorks_report(lanham,positive_tensor,\
                    use_Lanham)
            if individual:
                print("\nINDIVIDUAL")
                for i in self.components:
                    if i != 0:
                        info = self.components[i].properties_dict
                        shifted_dict = \
                            self.components[i].shift_properties_to_location(\
                            np.zeros((3,1)))
                        info["origin_inertia_tensor"] = \
                            shifted_dict["inertia_tensor"]
                        info["cg_location"] = \
                            shifted_dict["cg_location"]
                        shifted_dic2 = \
                            self.components[i].shift_properties_to_location(\
                            self.components[i].cg_location)
                        info["inertia_tensor"] = \
                            shifted_dic2["inertia_tensor"]
                        name = self.components[i].name
                        self.report_as_SolidWorks_report(info,positive_tensor,False,name)
                        if use_Lanham:
                            info = self.components[i].shift_properties_to_location(\
                                self.cg_location_lanham,True)
                            info["origin_inertia_tensor"] = \
                                self.components[i].shift_properties_to_location(\
                                np.zeros((3,1)),True)["inertia_tensor"]
                            self.report_as_SolidWorks_report(info,positive_tensor,use_Lanham,name)
        
        if use_Lanham:
            return self.properties_dict, lanham
        else:
            return self.properties_dict


    def get_mass_properties_about_point(self,point,report=False,individual=False,positive_tensor=True):

        # determine properties of each component
        for i in self.components:
            self.components[i].get_mass_properties()
        
        # determine total mass
        self.mass = 0.0
        for i in self.components:
            self.mass += self.components[i].mass
        
        # determine total cg location
        self.cg_location = np.zeros((3,1))
        for i in self.components:
            self.cg_location += self.components[i].mass * \
                self.components[i].get_cg_location()
        self.cg_location /= self.mass
        
        # determine total angular momentum
        self.angular_momentum = 0.0
        for i in self.components:
            self.angular_momentum += self.components[i].angular_momentum
        
        # determine total inertia
        self.inertia_tensor = np.zeros((3,3))
        location = point
        for i in self.components:
            self.inertia_tensor += \
                self.components[i].shift_properties_to_location(\
                location)["inertia_tensor"]

        # return dictionary of values
        self.properties_dict = {
            "mass" : self.mass,
            "volume" : self.volume,
            "cg_location" : self.cg_location,
            "angular_momentum" : self.angular_momentum,
            "inertia_tensor" : self.inertia_tensor
        }

        # report
        if report:
            print("TOTAL")
            self.report_as_SolidWorks_report(self.properties_dict,\
                positive_tensor)
            if individual:
                print("INDIVIDUAL")
                for i in self.components:
                    if i != 0:
                        info = self.components[i].properties_dict
                        info["origin_inertia_tensor"] = \
                            self.components[i].shift_properties_to_location(\
                            np.zeros((3,1)))["inertia_tensor"]
                        name = self.components[i].name
                        self.report_as_SolidWorks_report(info,positive_tensor,name)
        
        return self.properties_dict


    def _build_sphere(self,component):
        # create base arrays
        num = 100
        r = component.r2*np.ones((num,))
        t = np.linspace(0.,2.*np.pi,num=num)

        # create planar circles circle
        xy = np.block([ [r*np.cos(t)], [r*np.sin(t)], [0.*r] ])
        xz = np.block([ [r*np.cos(t)], [0.*r], [r*np.sin(t)] ])
        yz = np.block([ [0.*r], [r*np.cos(t)], [r*np.sin(t)] ])

        # create hollow lines
        hxy = xy*component.r1/component.r2
        hxz = xz*component.r1/component.r2
        hyz = yz*component.r1/component.r2

        return [xy,xz,yz],[hxy,hxz,hyz]


    def _build_ellipsoid(self,component):
        # create base arrays
        num = 100
        r = np.ones((num,))
        t = np.linspace(0.,2.*np.pi,num=num)

        # read in components
        ao = component._ao
        bo = component._bo
        co = component._co
        ai = component._ai
        bi = component._bi
        ci = component._ci

        # create planar circles circle
        xy = np.block([ [ao*r*np.cos(t)], [bo*r*np.sin(t)], [0.*r] ])
        xz = np.block([ [ao*r*np.cos(t)], [0.*r], [co*r*np.sin(t)] ])
        yz = np.block([ [0.*r], [bo*r*np.cos(t)], [co*r*np.sin(t)] ])

        # create hollow lines
        hxy = np.block([ [ai*r*np.cos(t)], [bi*r*np.sin(t)], [0.*r] ])
        hxz = np.block([ [ai*r*np.cos(t)], [0.*r], [ci*r*np.sin(t)] ])
        hyz = np.block([ [0.*r], [bi*r*np.cos(t)], [ci*r*np.sin(t)] ])

        return [xy,xz,yz],[hxy,hxz,hyz]


    def _build_half_ellipsoid(self,component):
        # create base arrays
        num = 100
        r = np.ones((num,))
        t = np.linspace(0.,2.*np.pi,num=num)
        h = np.linspace(-np.pi/2.,np.pi/2.,num=num)

        # read in components
        ao = component._ao
        bo = component._bo
        co = component._co
        ai = component._ai
        bi = component._bi
        ci = component._ci

        # create planar circles circle
        xy = np.block([ [ao*r*np.cos(h)], [bo*r*np.sin(h)], [0.*r] ])
        xz = np.block([ [ao*r*np.cos(h)], [0.*r], [co*r*np.sin(h)] ])
        yz = np.block([ [0.*r], [bo*r*np.cos(t)], [co*r*np.sin(t)] ])

        # create hollow lines
        hxy = np.block([ [ai*r*np.cos(h)], [bi*r*np.sin(h)], [0.*r] ])
        hxz = np.block([ [ai*r*np.cos(h)], [0.*r], [ci*r*np.sin(h)] ])
        hyz = np.block([ [0.*r], [bi*r*np.cos(t)], [ci*r*np.sin(t)] ])

        return [xy,xz,yz],[hxy,hxz,hyz]


    def _build_cylinder(self,component):
        # create base arrays
        num = 100
        r = np.ones((num,))
        t = np.linspace(0.,2.*np.pi,num=num)
        h = component._h
        bbo, bbi = component._bbo, component._bbi
        bto, bti = component._bto, component._bti
        cbo, cbi = component._cbo, component._cbi
        cto, cti = component._cto, component._cti

        # create planar circles circle
        fnt = np.block([ [0.*r  ], [bbo*r*np.cos(t)], [cbo*r*np.sin(t)] ])
        aft = np.block([ [0.*r+h], [bto*r*np.cos(t)], [cto*r*np.sin(t)] ])
        top = np.block([ 
            [0.,h],
            [0.,0.],
            [-cbo,-cto]
        ])
        bot = np.block([ 
            [0.,h],
            [0.,0.],
            [ cbo, cto]
        ])
        lef = np.block([ 
            [0.,h],
            [-bbo,-bto],
            [0.,0.]
        ])
        rig = np.block([ 
            [0.,h],
            [ bbo, bto],
            [0.,0.]
        ])

        # create hollow lines
        if bbi == bti == cbi == cti == 0.0:
            hfnt = np.zeros((3,2))
            haft = np.zeros((3,2))
            htop = np.zeros((3,2))
            hbot = np.zeros((3,2))
            hlef = np.zeros((3,2))
            hrig = np.zeros((3,2))
        else:
            hfnt = fnt*1.
            haft = aft*1.
            htop = top*1.
            hbot = bot*1.
            hlef = lef*1.
            hrig = rig*1.
            hfnt[1] = hfnt[1]*bbi/bbo
            hfnt[2] = hfnt[2]*cbi/cbo
            haft[1] = haft[1]*bbi/bbo
            haft[2] = haft[2]*cbi/cbo
            htop[1] = htop[1]*bbi/bbo
            htop[2] = htop[2]*cbi/cbo
            hbot[1] = hbot[1]*bbi/bbo
            hbot[2] = hbot[2]*cbi/cbo
            hlef[1] = hlef[1]*bbi/bbo
            hlef[2] = hlef[2]*cbi/cbo
            hrig[1] = hrig[1]*bbi/bbo
            hrig[2] = hrig[2]*cbi/cbo

        return [fnt,aft,top,bot,lef,rig],[hfnt,haft,htop,hbot,hlef,hrig]


    def _build_half_cylinder(self,component):
        # create base arrays
        num = 100
        r = np.ones((num,))
        t = np.linspace(-np.pi/2.,np.pi/2.,num=num)
        h = component._h
        bbo, bbi = component._bbo, component._bbi
        bto, bti = component._bto, component._bti
        cbo, cbi = component._cbo, component._cbi
        cto, cti = component._cto, component._cti

        # create planar circles circle
        fnt = np.block([ [0.*r  ], [bbo*r*np.cos(t)], [cbo*r*np.sin(t)] ])
        aft = np.block([ [0.*r+h], [bto*r*np.cos(t)], [cto*r*np.sin(t)] ])
        top = np.block([ 
            [0.,h],
            [0.,0.],
            [-cbo,-cto]
        ])
        bot = np.block([ 
            [0.,h],
            [0.,0.],
            [ cbo, cto]
        ])
        l0 = np.block([ 
            [0.,0.],
            [0.,0.],
            [ cbi, cbo]
        ])
        l1 = np.block([ 
            [0.,0.],
            [0.,0.],
            [-cbi,-cbo]
        ])
        l2 = np.block([ 
            [h,h],
            [0.,0.],
            [ cti, cto]
        ])
        l3 = np.block([ 
            [h,h],
            [0.,0.],
            [-cti,-cto]
        ])
        rig = np.block([ 
            [0.,h],
            [ bbo, bto],
            [0.,0.]
        ])

        # create hollow lines
        if bbi == bti == cbi == cti == 0.0:
            hfnt = np.zeros((3,2))
            haft = np.zeros((3,2))
            htop = np.zeros((3,2))
            hbot = np.zeros((3,2))
            hrig = np.zeros((3,2))
        else:
            hfnt = fnt*1.
            haft = aft*1.
            htop = top*1.
            hbot = bot*1.
            hrig = rig*1.
            hfnt[1] = hfnt[1]*bbi/bbo
            hfnt[2] = hfnt[2]*cbi/cbo
            haft[1] = haft[1]*bbi/bbo
            haft[2] = haft[2]*cbi/cbo
            htop[1] = htop[1]*bbi/bbo
            htop[2] = htop[2]*cbi/cbo
            hbot[1] = hbot[1]*bbi/bbo
            hbot[2] = hbot[2]*cbi/cbo
            hrig[1] = hrig[1]*bbi/bbo
            hrig[2] = hrig[2]*cbi/cbo

        return [fnt,aft,top,bot,l0,l1,l2,l3,rig],[hfnt,haft,htop,hbot,hrig]


    def _build_cuboid(self,component):
        # pull out components
        lx2 = component.lx2/2.
        ly2 = component.ly2/2.
        lz2 = component.lz2/2.
        lx1 = component.lx1/2.
        ly1 = component.ly1/2.
        lz1 = component.lz1/2.

        # create cuboid
        cub = [
            # x lines
            np.array([ [ lx2,-lx2],[ ly2, ly2],[-lz2,-lz2] ]),
            np.array([ [ lx2,-lx2],[-ly2,-ly2],[-lz2,-lz2] ]),
            np.array([ [ lx2,-lx2],[ ly2, ly2],[ lz2, lz2] ]),
            np.array([ [ lx2,-lx2],[-ly2,-ly2],[ lz2, lz2] ]),
            # front face
            np.array([ [ lx2, lx2],[-ly2, ly2],[-lz2,-lz2] ]),
            np.array([ [ lx2, lx2],[ ly2, ly2],[-lz2, lz2] ]),
            np.array([ [ lx2, lx2],[ ly2,-ly2],[ lz2, lz2] ]),
            np.array([ [ lx2, lx2],[-ly2,-ly2],[ lz2,-lz2] ]),
            # aft face
            np.array([ [-lx2,-lx2],[-ly2, ly2],[-lz2,-lz2] ]),
            np.array([ [-lx2,-lx2],[ ly2, ly2],[-lz2, lz2] ]),
            np.array([ [-lx2,-lx2],[ ly2,-ly2],[ lz2, lz2] ]),
            np.array([ [-lx2,-lx2],[-ly2,-ly2],[ lz2,-lz2] ])
        ]

        # hollow
        hcub = [
            # x lines
            np.array([ [ lx1,-lx1],[ ly1, ly1],[-lz1,-lz1] ]),
            np.array([ [ lx1,-lx1],[-ly1,-ly1],[-lz1,-lz1] ]),
            np.array([ [ lx1,-lx1],[ ly1, ly1],[ lz1, lz1] ]),
            np.array([ [ lx1,-lx1],[-ly1,-ly1],[ lz1, lz1] ]),
            # front face
            np.array([ [ lx1, lx1],[-ly1, ly1],[-lz1,-lz1] ]),
            np.array([ [ lx1, lx1],[ ly1, ly1],[-lz1, lz1] ]),
            np.array([ [ lx1, lx1],[ ly1,-ly1],[ lz1, lz1] ]),
            np.array([ [ lx1, lx1],[-ly1,-ly1],[ lz1,-lz1] ]),
            # aft face
            np.array([ [-lx1,-lx1],[-ly1, ly1],[-lz1,-lz1] ]),
            np.array([ [-lx1,-lx1],[ ly1, ly1],[-lz1, lz1] ]),
            np.array([ [-lx1,-lx1],[ ly1,-ly1],[ lz1, lz1] ]),
            np.array([ [-lx1,-lx1],[-ly1,-ly1],[ lz1,-lz1] ])
        ]

        return cub,hcub


    def _build_rotor(self,component):
        # pull in components
        Nb = component._Nb
        cr,ct = component._cr,component._ct
        tr,tt = component._tr,component._tt
        rr,rt = component._rr,component._rt
        u0 = component.u0_for_viz

        # create base arrays
        num = 100
        t = np.linspace(0.,2.*np.pi,num=num)

        # create height function
        r = np.linspace(rr,rt,num)
        c = (ct - cr)*(r - rr)/(rt - rr) + cr
        m = (tt - tr)*(r - rr)/(rt - rr) + tr
        h = Nb*m*c**2.*u0/2./np.pi/r

        # create planar circles circle
        frt = np.block([ [0.*t+h[ 0]/2.], [rr*np.cos(t)], [rr*np.sin(t)] ])
        art = np.block([ [0.*t-h[ 0]/2.], [rr*np.cos(t)], [rr*np.sin(t)] ])
        ftp = np.block([ [0.*t+h[-1]/2.], [rt*np.cos(t)], [rt*np.sin(t)] ])
        atp = np.block([ [0.*t-h[-1]/2.], [rt*np.cos(t)], [rt*np.sin(t)] ])

        draws = [frt,art,ftp,atp] #+ [fht,aht,fhb,ahb]

        # lines
        for i in range(Nb):
            p = np.deg2rad(i*360./Nb)
            draws.append( np.block([ [+h/2.], [r*np.sin(p)], [-r*np.cos(p)] ]))
            draws.append( np.block([ [-h/2.], [r*np.sin(p)], [-r*np.cos(p)] ]))


        return draws,None


    def _build_symmetric_airfoil(self,component):
        # pull in components
        cr,ct = component._cr,component._ct
        tr,tt = component._tr,component._tt
        a0 = component._a0
        a1 = component._a1
        a2 = component._a2
        a3 = component._a3
        a4 = component._a4
        b = component._b*component._delta

        # create base arrays
        num = 100
        xa = np.linspace(0.,1.,num=num)
        ta = a0*xa**0.5 + a1*xa + a2*xa**2. + a3*xa**3. + a4*xa**4.

        # root airfoil
        sw = component._b*np.tan(component._Lambda)
        rtu = np.block([ [(0.25-xa)*cr   ], [0.*xa  ], [-ta*tr*cr/2.] ])
        rtl = np.block([ [(0.25-xa)*cr   ], [0.*xa  ], [ ta*tr*cr/2.] ])
        ttu = np.block([ [(0.25-xa)*ct-sw], [0.*xa+b], [-ta*tt*ct/2.] ])
        ttl = np.block([ [(0.25-xa)*ct-sw], [0.*xa+b], [ ta*tt*ct/2.] ])
        lel = np.block([ [ 0.25*cr, 0.25*ct-sw], [0.,b], [0.0,0.0] ])
        tel = np.block([ [-0.75*cr,-0.75*ct-sw], [0.,b], [0.0,0.0] ])


        return [rtu,rtl,ttu,ttl,lel,tel],None


    def _build_diamond_airfoil(self,component):
        # pull in components
        cr,ct = component._cr,component._ct
        tr,tt = component._tr,component._tt
        xmt = component._xmt
        b = component._b*component._delta

        # create base arrays
        num = 101
        xa = np.linspace(0.,1.,num=num)
        ta = xa*0.
        ta[xa<=xmt] = xa[xa<=xmt]/xmt
        ta[xa>xmt] = (1. - xa[xa>xmt])/(1. - xmt)

        # root airfoil
        sw = component._b*np.tan(component._Lambda)
        rtu = np.block([ [(0.25-xa)*cr   ], [0.*xa  ], [-ta*tr*cr/2.] ])
        rtl = np.block([ [(0.25-xa)*cr   ], [0.*xa  ], [ ta*tr*cr/2.] ])
        ttu = np.block([ [(0.25-xa)*ct-sw], [0.*xa+b], [-ta*tt*ct/2.] ])
        ttl = np.block([ [(0.25-xa)*ct-sw], [0.*xa+b], [ ta*tt*ct/2.] ])
        lel = np.block([ [ 0.25*cr, 0.25*ct-sw], [0.,b], [0.0,0.0] ])
        tel = np.block([ [-0.75*cr,-0.75*ct-sw], [0.,b], [0.0,0.0] ])
        upl = np.block([ [(0.25-xmt)*cr,(0.25-xmt)*ct-sw], [0.,b], 
            [-tr*cr/2.,-tt*ct/2.] ])
        lol = np.block([ [(0.25-xmt)*cr,(0.25-xmt)*ct-sw], [0.,b], 
            [ tr*cr/2., tt*ct/2.] ])


        return [rtu,rtl,ttu,ttl,lel,tel,upl,lol],None


    def visualize(self,no_color=False,show_legend=False,plot_ids=None,
    filename=None):
        # initialize plot
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')

        # which ids to plot, or all
        if plot_ids == None:
            plot_ids = self.components.keys()
        

        # calculate component lines, plot
        ci = 0
        if no_color:
            colors = ["k"]
        else:
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        collen = len(colors)
        x_lims = [1.0e+100, -1.0e100]
        y_lims = [1.0e+100, -1.0e100]
        z_lims = [1.0e+100, -1.0e100]
        plottable_shapes = [
            "sphere","ellipsoid","half_ellipsoid",
            "cylinder","half_cylinder",
            "cuboid","rotor","symmetric_airfoil","diamond_airfoil"
        ]
        for i in plot_ids:
            # build object by type
            component = self.components[i]
            if component.type in plottable_shapes:
                if component.type in ["symmetric_airfoil","diamond_airfoil"]:
                    nums = len(component._components)
                else:
                    nums = 1
                    R_matrix = component.R
                    cg = component.get_cg_location()
                for j in range(nums):
                    if component.type == "sphere":
                        lines,hlines = self._build_sphere(component)
                    elif component.type == "ellipsoid":
                        lines,hlines = self._build_ellipsoid(component)
                    elif component.type == "half_ellipsoid":
                        lines,hlines = self._build_half_ellipsoid(component)
                    elif component.type == "cylinder":
                        lines,hlines = self._build_cylinder(component)
                    elif component.type == "half_cylinder":
                        lines,hlines = self._build_half_cylinder(component)
                    elif component.type == "cuboid":
                        lines,hlines = self._build_cuboid(component)
                    elif component.type == "rotor":
                        lines,hlines = self._build_rotor(component)
                    elif component.type == "symmetric_airfoil":
                        lines,hlines = self._build_symmetric_airfoil(
                            component._components[j])
                        R_matrix = component._components[j].R
                    elif component.type == "diamond_airfoil":
                        lines,hlines = self._build_diamond_airfoil(
                            component._components[j])
                        R_matrix = component._components[j].R
                    else:
                        lines = [[[],[],[]]]
                    
                    # plot lines
                    c = colors[ci % collen]
                    for linenum,line in enumerate(lines):
                        # rotate
                        line = np.matmul(R_matrix,line)

                        # shift by cg location
                        if component.type[-7:] == "airfoil":
                            cg = component._components[j].locations["root"]
                        elif component.type[-8:] == "cylinder" or \
                            component.type == "half_ellipsoid":
                            cg = component.origin
                        line[0] += cg[0]
                        line[1] += cg[1]
                        line[2] += cg[2]

                        # plot
                        if linenum == 0:
                            lbl = component.name
                        else:
                            lbl = ""
                        ax.plot(line[0],line[1],line[2],c=c,label=lbl)

                        # pull out max and min vals
                        x_lims[0] = min(np.min(line[0]),x_lims[0])
                        y_lims[0] = min(np.min(line[1]),y_lims[0])
                        z_lims[0] = min(np.min(line[2]),z_lims[0])
                        x_lims[1] = max(np.max(line[0]),x_lims[1])
                        y_lims[1] = max(np.max(line[1]),y_lims[1])
                        z_lims[1] = max(np.max(line[2]),z_lims[1])
                    
                    # plot hollow lines
                    if hlines is not None:
                        for line in hlines:
                            # rotate
                            line = np.matmul(R_matrix,line)

                            # shift by cg location

                            if component.type[-8:] == "cylinder" or \
                                component.type == "half_ellipsoid":
                                cg = component.origin
                            else:
                                cg = component.get_cg_location()
                            # cg = component.locations["root"]
                            line[0] += cg[0]
                            line[1] += cg[1]
                            line[2] += cg[2]

                            # plot
                            ax.plot(line[0],line[1],line[2],c=c,ls="--")

                            # pull out max and min vals
                            x_lims[0] = min(np.min(line[0]),x_lims[0])
                            y_lims[0] = min(np.min(line[1]),y_lims[0])
                            z_lims[0] = min(np.min(line[2]),z_lims[0])
                            x_lims[1] = max(np.max(line[0]),x_lims[1])
                            y_lims[1] = max(np.max(line[1]),y_lims[1])
                            z_lims[1] = max(np.max(line[2]),z_lims[1])

                # update colors index
                ci += 1


        # set axes lengths

        # Find out which axis has the widest limits
        x_diff = x_lims[1] - x_lims[0]
        y_diff = y_lims[1] - y_lims[0]
        z_diff = z_lims[1] - z_lims[0]
        max_diff = max([x_diff, y_diff, z_diff])

        # Determine the center of each set of axis limits
        x_cent = x_lims[0] + 0.5*x_diff
        y_cent = y_lims[0] + 0.5*y_diff
        z_cent = z_lims[0] + 0.5*z_diff

        # Scale the axis limits so they all have the same width as the widest set
        x_lims[0] = x_cent - 0.5*max_diff
        x_lims[1] = x_cent + 0.5*max_diff

        y_lims[0] = y_cent - 0.5*max_diff
        y_lims[1] = y_cent + 0.5*max_diff

        z_lims[0] = z_cent - 0.5*max_diff
        z_lims[1] = z_cent + 0.5*max_diff

        # Set limits so it is a right-handed coordinate system with z pointing down
        ax.set_xlim3d(x_lims[1], x_lims[0])
        ax.set_ylim3d(y_lims[0], y_lims[1])
        ax.set_zlim3d(z_lims[1], z_lims[0])

        # other plot parameters
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        white = ( 0.0,0.0,0.0,0.0)
        ax.w_xaxis.set_pane_color(white)
        ax.w_yaxis.set_pane_color(white)
        ax.w_zaxis.set_pane_color(white)
        plt.tight_layout()
        # ax.invert_xaxis()
        ax.view_init(30.,-140)
        if show_legend:
            ax.legend()

        # show plot
        if filename is not None:
            plt.savefig(filename)
            plt.close()
        else:
            plt.show()