import numpy as np
import json
import matplotlib.pyplot as plt

class Component:
    """A default class for calculating and containing the mass properties of an object.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # initialize cg location, inertia tensor, mass, density
        self.density = 0.0
        self.mass = 0.0
        self.volume = 0.0
        self.inertia_tensor = np.zeros((3,3))
        self.angular_momentum = np.zeros((3,1))
        self.cg_location = np.zeros((3,1))
        self.locations = {
            "root" : np.zeros((3,1)),
            "tip"  : np.zeros((3,1))
        }
        self.name = "empty_name"
        self.can_use_Lanham = False
        self._components = []

        # retrieve info
        self._read_in_info(input_dict)


    def _read_in_info(self,input_dict):
        """A function which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary
        # store mass / density variables
        if "density" in input_dict:
            self.density = input_dict.get("density")
            self._given_density_not_mass = True
        elif "mass" in input_dict:
            self.mass = input_dict.get("mass")
            self._given_density_not_mass = False
        else: # raise error if neither given
            self.mass = 0.0
            self._given_density_not_mass = False
            raise ValueError("No mass / density property given")

        # read in orientation
        self.orientation = input_dict.get("orientation",np.zeros((3,)))
        if len(self.orientation) == 3: # phi theta psi
            # define a rotation matrix
            rads = np.deg2rad(self.orientation)
            CF = np.cos(rads[0]); SF = np.sin(rads[0])
            CT = np.cos(rads[1]); ST = np.sin(rads[1])
            CS = np.cos(rads[2]); SS = np.sin(rads[2])
            self.R = np.array([
                [ CT*CS, SF*ST*CS - CF*SS, CF*ST*CS + SF*SS],
                [ CT*SS, SF*ST*SS + CF*CS, CF*ST*SS - SF*CS],
                [   -ST,            SF*CT,            CF*CT]
            ])
        else: # quaternion
            e0 , ex, ey, ez = self.orientation
            e02,ex2,ey2,ez2 = e0**2.,ex**2.,ey**2.,ez**2.
            exe0,exey,exez = ex*e0,ex*ey,ex*ez
            eye0     ,eyez = ey*e0      ,ey*ez
            eze0,ezey      = ez*e0,ez*ey
            # define a rotation matrix
            CG = np.cos(-self._delta*self._Gamma)
            SG = np.sin(-self._delta*self._Gamma)
            self.R = np.array([
                [ ex2+e02-ey2-ez2,  2.*(exey-eze0),  2.*(exez+eye0)],
                [  2.*(exey+eze0), ey2+e02-ex2-ez2,  2.*(eyez-exe0)],
                [  2.*(exez-eye0),  2.*(eyez+exe0), ez2+e02-ex2-ey2]
            ])
        self.R0 = self.R * 1.0 # initialize orientation for later use


    def update_densities(self):
        return 0


    def shift_properties_to_location(self,input_location,use_Lanham=False):
        """Method which determines the mass properties of the given component
        about a given location.
        
        Parameters
        ----------
        input_location : array
            The location about which to determine the mass properties. Must be
            formatted as input_location = [[x],[y],[z]].
        """

        # if input_location is a list, turn into numpy array
        if isinstance(input_location, list):
            input_location = np.array(input_location)
        
        # bring in mass, volume
        if use_Lanham and self.can_use_Lanham:
            mass = self.mass_lanham * 1.0
            volume = self.volume_lanham * 1.0
        else:
            mass = self.mass * 1.0
            volume = self.volume * 1.0
        
        # rotate inertia, cg location for rotation values
        if use_Lanham and self.can_use_Lanham:
            inertia_tensor = np.matmul(self.R0,np.matmul(self.inertia_tensor_lanham,self.R0.T))
        else:
            inertia_tensor = np.matmul(self.R,np.matmul(self.inertia_tensor,self.R.T))

        # shift cg by root location given to aircraft coordinate system
        new_cg_location = self.get_cg_location(use_Lanham)
        
        # determine new location
        s = input_location - new_cg_location

        # calculate mass shift from parallel axis theorem
        inner_product = np.matmul(s.T,s)[0,0]
        outer_product = np.matmul(s,s.T)
        I_shift = mass*( inner_product * np.eye(3) - outer_product )

        # calculate inertia tensor about new location
        Inew = inertia_tensor + I_shift

        # return dictionary of values
        output_dict = {
            "mass" : mass,
            "volume" : volume,
            "cg_location" : new_cg_location,
            "angular_momentum" : self.angular_momentum,
            "inertia_tensor" : Inew
        }
        return output_dict


    def get_mass_properties(self):
        """Method which returns mass, cg, I about cg rotated to total cframe"""

        self.properties_dict = {
            "mass" : self.mass,
            "volume" : self.volume,
            "cg_location" : self.cg_location,
            "angular_momentum" : self.angular_momentum,
            "inertia_tensor" : self.inertia_tensor
        }
        return self.properties_dict


    def get_mass(self,use_Lanham_approximations=False):
        """Method which returns the mass whether calculated using the method
        presented or the Lanham method."""
        if use_Lanham_approximations and self.can_use_Lanham:
            return self.mass_lanham
        else:
            return self.mass

    
    def get_volume(self,use_Lanham_approximations=False):
        """Method which returns the volume whether calculated using the
        method presented or the Lanham method."""
        if use_Lanham_approximations and self.can_use_Lanham:
            return self.volume_lanham
        else:
            return self.volume

    
    def get_cg_location(self,use_Lanham_approximations=False):
        """Method which returns the cg location whether calculated using the
        method presented or the Lanham method."""
        if use_Lanham_approximations and self.can_use_Lanham:
            return self.cg_location_lanham
        else:
            return self.cg_location

    
    def get_inertia_tensor(self,use_Lanham_approximations=False):
        """Method which returns the inertia tensor whether calculated using the
        method presented or the Lanham method."""
        if use_Lanham_approximations and self.can_use_Lanham:
            return self.inertia_tensor_lanham
        else:
            return self.inertia_tensor



class Cuboid(Component):
    """A default class for calculating and containing the mass properties of a
    Cuboid.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # invoke init of parent
        Component.__init__(self,input_dict)

        # retrieve additional info
        self._retrieve_info(input_dict)

        # calculate volume
        self.volume1 = self.lx1*self.ly1*self.lz1
        self.volume2 = self.lx2*self.ly2*self.lz2
        self.volume  = self.volume2 - self.volume1


    def _retrieve_info(self,input_dict):
        """A function which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary
        # store cg location
        connect_to = input_dict.get("connect_to",{})
        x_cg = connect_to.get("dx",0.0)
        y_cg = connect_to.get("dy",0.0)
        z_cg = connect_to.get("dz",0.0)
        self.cg_location = np.array([[x_cg],[y_cg],[z_cg]])

        # store lengths and hollow lengths
        self.lx2 = input_dict.get("x_length",1.0)
        self.ly2 = input_dict.get("y_length",1.0)
        self.lz2 = input_dict.get("z_length",1.0)
        self.lx1 = input_dict.get("x_hollow_length",0.0)
        self.ly1 = input_dict.get("y_hollow_length",0.0)
        self.lz1 = input_dict.get("z_hollow_length",0.0)


    def get_mass_properties(self):
        """Method which returns mass, cg, I about cg rotated to total cframe"""

        # calculate mass
        if self._given_density_not_mass:
            self.mass = self.density * self.volume
        
        # store values for ease of calculation
        lx1sq, ly1sq, lz1sq = self.lx1**2., self.ly1**2., self.lz1**2.
        lx2sq, ly2sq, lz2sq = self.lx2**2., self.ly2**2., self.lz2**2.
        v1, v2 = self.volume1, self.volume2

        # calculate inertia
        Ixxo = self.mass / 12. * (v2*(ly2sq+lz2sq) - v1*(ly1sq+lz1sq))/(v2-v1)
        Iyyo = self.mass / 12. * (v2*(lx2sq+lz2sq) - v1*(lx1sq+lz1sq))/(v2-v1)
        Izzo = self.mass / 12. * (v2*(lx2sq+ly2sq) - v1*(lx1sq+ly1sq))/(v2-v1)
        
        Ixyo = 0.0
        Ixzo = 0.0
        Iyzo = 0.0

        # create inertia tensor
        self.inertia_tensor = np.array([
            [ Ixxo,-Ixyo,-Ixzo],
            [-Ixyo, Iyyo,-Iyzo],
            [-Ixzo,-Iyzo, Izzo]
        ])

        self.properties_dict = {
            "mass" : self.mass,
            "volume" : self.volume,
            "cg_location" : self.cg_location,
            "angular_momentum" : self.angular_momentum,
            "inertia_tensor" : self.inertia_tensor
        }
        return self.properties_dict



class Cylinder(Component):
    """A default class for calculating and containing the mass properties of a
    Cylinder.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # invoke init of parent
        Component.__init__(self,input_dict)

        # retrieve additional info
        self._retrieve_info(input_dict)

        # calculate volume
        self.volume  = np.pi * self.l  * ( self.r2**2. - self.r1**2. )


    def _retrieve_info(self,input_dict):
        """A function which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary
        # store cg location
        connect_to = input_dict.get("connect_to",{})
        x_cg = connect_to.get("dx",0.0)
        y_cg = connect_to.get("dy",0.0)
        z_cg = connect_to.get("dz",0.0)
        self.cg_location = np.array([[x_cg],[y_cg],[z_cg]])

        # store lengths and hollow lengths
        self.l = input_dict.get("length",1.0)
        self.r2 = input_dict.get("radius",1.0)
        self.r1 = input_dict.get("hollow_radius",0.0)


    def get_mass_properties(self):
        """Method which returns mass, cg, I about cg rotated to total cframe"""

        # calculate mass
        if self._given_density_not_mass:
            self.mass = self.density * self.volume
        
        # store values for ease of calculation
        h, r1, r2 = self.l, self.r1, self.r2

        # calculate inertia
        Ixxo = self.mass /  2. *        (r2**2. + r1**2.)
        Iyyo = self.mass / 12. * ( 3. * (r2**2. + r1**2.) + h**2. )
        Izzo = self.mass / 12. * ( 3. * (r2**2. + r1**2.) + h**2. )
        
        Ixyo = 0.0
        Ixzo = 0.0
        Iyzo = 0.0

        # create inertia tensor
        self.inertia_tensor = np.array([
            [ Ixxo,-Ixyo,-Ixzo],
            [-Ixyo, Iyyo,-Ixzo],
            [-Ixzo,-Ixzo, Izzo]
        ])

        self.properties_dict = {
            "mass" : self.mass,
            "volume" : self.volume,
            "cg_location" : self.cg_location,
            "angular_momentum" : self.angular_momentum,
            "inertia_tensor" : self.inertia_tensor
        }
        return self.properties_dict



class Sphere(Component):
    """A default class for calculating and containing the mass properties of a
    Sphere.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # invoke init of parent
        Component.__init__(self,input_dict)

        # retrieve additional info
        self._retrieve_info(input_dict)

        # calculate volume
        self.volume  = 4. / 3. * np.pi * ( self.r2**3. - self.r1**3. )


    def _retrieve_info(self,input_dict):
        """A function which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary
        # store cg location
        connect_to = input_dict.get("connect_to",{})
        x_cg = connect_to.get("dx",0.0)
        y_cg = connect_to.get("dy",0.0)
        z_cg = connect_to.get("dz",0.0)
        self.cg_location = np.array([[x_cg],[y_cg],[z_cg]])

        # store lengths and hollow lengths
        self.r2 = input_dict.get("radius",1.0)
        self.r1 = input_dict.get("hollow_radius",0.0)


    def get_mass_properties(self):
        """Method which returns mass, cg, I about cg rotated to total cframe"""

        # calculate mass
        if self._given_density_not_mass:
            self.mass = self.density * self.volume
        
        # store values for ease of calculation
        r1, r2 = self.r1, self.r2

        # calculate inertia
        Ixxo = 2. / 5. * self.mass * (r2**5. - r1**5.) / (r2**3. - r1**3.)
        Iyyo = 2. / 5. * self.mass * (r2**5. - r1**5.) / (r2**3. - r1**3.)
        Izzo = 2. / 5. * self.mass * (r2**5. - r1**5.) / (r2**3. - r1**3.)
        
        Ixyo = 0.0
        Ixzo = 0.0
        Iyzo = 0.0

        # create inertia tensor
        self.inertia_tensor = np.array([
            [ Ixxo,-Ixyo,-Ixzo],
            [-Ixyo, Iyyo,-Ixzo],
            [-Ixzo,-Ixzo, Izzo]
        ])

        self.properties_dict = {
            "mass" : self.mass,
            "volume" : self.volume,
            "cg_location" : self.cg_location,
            "angular_momentum" : self.angular_momentum,
            "inertia_tensor" : self.inertia_tensor
        }
        return self.properties_dict



class Prismoid(Component):
    """A default class for calculating and containing the mass properties of a
    Psuedo Prismoid.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # invoke init of parent
        Component.__init__(self,input_dict)

        # retrieve additional info
        self.can_use_Lanham = True
        self._retrieve_info(input_dict)

        # initialize the volume
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        trcr, ttct = tr * cr, tt * ct
        ka = trcr * (2.*cr + ct) + ttct *(cr + 2.*ct)
        u0 = 1.0
        self.volume = self._b / 12. * ka * u0


    def _retrieve_info(self,input_dict):
        """A function which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary

        # store input values
        self._side = input_dict.get("side","right")
        if self._side == "right": self._delta = 1.0
        else: self._delta = -1.0
        
        # store other variables
        geometry = input_dict.get("geometry",{})
        self._b = geometry.get("span",1.0)
        self._cr,self._ct = geometry.get("chord",[1.0,1.0])
        self._tr,self._tt = geometry.get("thickness",[0.1,0.1])
        self._ur,self._ut = geometry.get("camber",[0.01,0.01])
        self._xmt = geometry.get("max_thickness_location",0.5)
        self._Lambda = np.deg2rad(input_dict.get("sweep",0.0))
        self._Gamma = np.deg2rad(input_dict.get("dihedral",0.0))
        root_location = input_dict.get("root_location",[0.0,0.0,0.0])
        tip_location  = input_dict.get( "tip_location",[0.0,0.0,0.0])
        self.locations = {
            "root" : np.array(root_location)[:,np.newaxis],
            "tip"  : np.array( tip_location)[:,np.newaxis]
        }

        # define a rotation matrix
        CG = np.cos(-self._delta*self._Gamma)
        SG = np.sin(-self._delta*self._Gamma)
        self.R = np.matmul( self.R0, np.array([
            [1.0, 0.0, 0.0],
            [0.0,  CG, -SG],
            [0.0,  SG,  CG]
        ]) )


    def _make_kappa_values(self):

        # initialize certian values for ease of calculation
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        # times values
        cr2, ct2 = cr**2., ct**2.
        crct = cr * ct
        cr3, ct3 = cr**3., ct**3.
        cr2ct, crct2 = cr2 * ct, cr * ct2
        trcr, ttct = tr * cr, tt * ct

        # calculate kappa values
        self._ka = trcr * (2.*cr + ct) + ttct *(cr + 2.*ct)

        one = 3.*cr2 + 2.*crct +    ct2
        two =    cr2 + 2.*crct + 3.*ct2
        self._kb = trcr * one + ttct * two
        
        self._kc = trcr * (   cr + ct) + ttct *(cr + 3.*ct)

        one = 3.*cr2 + 4.*crct + 3.*ct2
        two =    cr2 + 3.*crct + 6.*ct2
        self._kd = trcr * one + 2. * ttct * two

        one = 4.*cr3 + 3.*cr2ct + 2.*crct2 +    ct3
        two =    cr3 + 2.*cr2ct + 3.*crct2 + 4.*ct3
        self._ke = trcr * one + ttct * two

        self._kf = trcr * (2.*cr + 3.*ct) + ttct *(3.*cr + 12.*ct)

        one = 4.*cr +    ct
        two = 3.*cr + 2.*ct
        thr = 2.*cr + 3.*ct
        fou =    cr + 4.*ct
        oner = trcr**3.*one
        self._kg = oner + trcr**2.*ttct*two + trcr*ttct**2.*thr + ttct**3.*fou

        # fixes for constant multiples outside equation
        self._ka = self._ka *   1. /    2.
        self._kb = self._kb *  48. /   80.
        self._kc = self._kc *  12. /   60.
        self._kd = self._kd *   1. /    2.
        self._ke = self._ke * 960. / 1440.
        self._kf = self._kf *   1. /    2.
        self._kg = self._kg * 240. / 3360.


    def _make_upsilon_values(self):

        # calculate upsilon values
        self._u0 = 1.0
        self._u1 = 1.0
        self._u2 = 1.0
        self._u3 = 1.0


    def _get_lanham_mass_properties(self):

        # # determine properties for later use
        # b = self._b
        # c, ct = self._cr, self._ct
        # tr, tt = self._cr*self._tr, self._ct*self._tt
        # LL = np.arctan( np.tan(self._Lambda) + 0.25 / b * ( c - ct ) )
        # LT = np.arctan( np.tan(self._Lambda) - 0.75 / b * ( c - ct ) )
        # tLL, tLT = np.tan(LL), np.tan(LT)
        # ca = self._cr
        # cb = np.tan(self._Lambda) + 0.25 * ( c - ct )
        # cc = cb + ct

        # # save mass
        # m = self.mass
        # self.mass_lanham = m

        # # calculate volume
        # V = b * ( tr*(c + b/2.*(tLT-tLL)) - (tr-tt)*(c/2. + b/3.*(tLT-tLL)) )
        # self.volume_lanham = V

        # # calculate cg location
        # xcg = 2 * b * (-ca**2. + cb**2. + cc*cb + cc**2.) / (-ca + cb + cc)
        # ycg = b**2./V*(tr*(c/2.+b/3.*(tLT-tLL)) -(tr-tt)*(c/3.+b/4.*(tLT-tLL)))
        # self.cg_location_lanham = np.array([
        #     [xcg],
        #     [self._delta * ycg],
        #     [0.0]
        # ])
        # print(xcg,1/(-ca + cb + cc))

        # calculate inertia tensor
        # Ixx = m * b**3. / V * ( (tr-tt) * (c / 4. + b / 5. * (tLT-tLL)) \
        #     + tr * (c / 3. + b / 4. * (tLT-tLL)) )
        
        # Iyy = m * b / V * ( tr * (c**3. / 3. + b * c * tLT*(c/2. + b/3. * tLT)\
        #     + b**3. / 12. * (tLT**3. - tLL**3.)) \
        #         - (tr-tt) * (c**3. / 6. + b * c * tLT * (c /3. + b /4. * tLT) \
        #             + b**3. / 15 * (tLT**3. - tLL**3.)) )
        
        # Izz = Ixx + Iyy

        # one = c**2.*b**2./4. + c*b**3./3.*tLT + b**4./ 8.*(tLT**2. - tLL**2.)
        # two = c**2.*b**2./6. + c*b**3./4.*tLT + b**4./10.*(tLT**2. - tLL**2.)
        # thr = m / V * tr      * np.sin(self._Gamma) * one
        # fou = m / V * (tr-tt) * np.sin(self._Gamma) * two
        # Ixz = thr - fou

        self.mass_lanham = self.mass
        self.volume_lanham = 1.0
        self.cg_location_lanham = np.zeros((3,1))
        Ixx = Iyy = Izz = Ixz = 0.0
        Ixy = 0.0
        Iyz = 0.0

        self.inertia_tensor_lanham = np.array([
            [Ixx,Ixy,Ixz],
            [Ixy,Iyy,Iyz],
            [Ixz,Iyz,Izz]
        ])


    def get_mass_properties(self):
        """Method which returns mass, cg, I about cg rotated to total cframe"""

        # calculate kappa values
        self._make_kappa_values()

        # calculate upsilon values
        self._make_upsilon_values()

        # calculate mass
        if self._given_density_not_mass:
            self.mass = self.density * self.volume
        
        # calculate lanham properties
        self._get_lanham_mass_properties()
        
        # calculate center of gravity values
        num1 = 3. * self._kb * self._u1
        num = num1 + 4. * self._b * self._kc * self._u0 * np.tan(self._Lambda)
        xbar = - num / 20. / self._ka / self._u0
        ybar = self._delta * self._b * self._kc * self._u0/5./self._ka/self._u0
        self.cg_location = np.array([xbar,ybar,0.0])[:,np.newaxis]

        # calculate moments and products of inertia about the wing root c/4
        num = 56. * self._b**2. * self._kf * self._u0 + self._kg * self._u3
        Ixxo = self.mass * num / 280. / self._ka / self._u0
        
        one = 2. * self._b * self._kf * self._u0 * np.tan(self._Lambda)**2.
        two = self._kd * self._u1 * np.tan(self._Lambda)
        num1 = 84. * self._b * (one + two) + 49. * self._ke * self._u2
        num = num1 + 3. * self._kg * self._u3
        Iyyo = self.mass * num / 840. / self._ka / self._u0

        one = self._b * ( np.tan(self._Lambda)**2. + 1. ) * self._kf * self._u0
        two = self._kd * self._u1 * np.tan(self._Lambda)
        num = 12. * self._b * (2. * one + two) + 7. * self._ke * self._u2
        Izzo = self.mass * num / 120. / self._ka / self._u0

        num1 = 4. * self._b * self._kf * self._u0 * np.tan(self._Lambda)
        num = num1 + self._kd * self._u1
        Ixyo = -self._delta * self._b * self.mass * num / 20./self._ka/self._u0

        Ixzo = 0.0
        Iyzo = 0.0

        # create inertia tensor
        Io = np.array([
            [ Ixxo,-Ixyo,-Ixzo],
            [-Ixyo, Iyyo,-Iyzo],
            [-Ixzo,-Iyzo, Izzo]
        ])

        # calculate mass shift from parallel axis theorem
        s = self.cg_location
        inner_product = np.matmul(s.T,s)[0,0]
        outer_product = np.matmul(s,s.T)
        I_shift = self.mass*( inner_product * np.eye(3) - outer_product )

        # calculate inertia tensor about the cg
        self.inertia_tensor = Io - I_shift

        self.properties_dict = {
            "mass" : self.mass,
            "volume" : self.volume,
            "cg_location" : self.cg_location,
            "angular_momentum" : self.angular_momentum,
            "inertia_tensor" : self.inertia_tensor
        }
        return self.properties_dict


    def get_cg_location(self,use_Lanham_approximations=False):
        """Method which returns the cg location whether calculated using the
        method presented or the Lanham method."""
        if use_Lanham_approximations:
            cg_location = np.matmul(self.R,self.cg_location_lanham)
        else:
            cg_location = np.matmul(self.R,self.cg_location)

        # shift cg by root location given
        return cg_location + self.locations["root"]



class PseudoPrismoid(Prismoid):
    """A default class for calculating and containing the mass properties of a
    Psuedo Prismoid.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # invoke init of parent
        Prismoid.__init__(self,input_dict)

        # retrieve additional info
        self._retrieve_info(input_dict)

        # initialize the volume
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        cr2, ct2 = cr**2., ct**2.
        crct = cr * ct
        ka = tr * (3.*cr2 + 2.*crct + ct2) + tt *(cr2 + 2.*crct + 3.*ct2)
        u0 = 1.0
        self.volume = self._b / 12. * ka * u0


    def _make_kappa_values(self):

        # initialize certian values for ease of calculation
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        # times values
        cr2, ct2 = cr**2., ct**2.
        crct = cr * ct
        cr3, ct3 = cr**3., ct**3.
        cr2ct, crct2 = cr2 * ct, cr * ct2
        cr4, ct4 = cr**4., ct**4.
        cr3ct, cr2ct2, crct3 = cr3 * ct, cr2 * ct2, cr * ct3

        # calculate kappa values
        self._ka = tr * (3.*cr2 + 2.*crct + ct2) + tt *(cr2 + 2.*crct + 3.*ct2)
        
        one = 4.*cr3 + 3.*cr2ct + 2.*crct2 +    ct3
        two =    cr3 + 2.*cr2ct + 3.*crct2 + 4.*ct3
        self._kb = tr * one + tt * two

        one = 3.*cr2 + 4.*crct + 3.*ct2
        two =    cr2 + 3.*crct + 6.*ct2
        self._kc = tr * one + 2. * tt * two

        one = (cr + ct) * (2.*cr2 + crct + 2.*ct2)
        two = cr3 + 3.*cr2ct + 6.*crct2 + 10.*ct3
        self._kd = tr * one + tt * two

        one = 5.*cr4 + 4.*cr3ct + 3.*cr2ct2 + 2.*crct3 +    ct4
        two =    cr4 + 2.*cr3ct + 3.*cr2ct2 + 4.*crct3 + 5.*ct4
        self._ke = tr * one + tt * two

        one = cr2 + 2.*crct +  2.*ct2
        two = cr2 + 4.*crct + 10.*ct2
        self._kf = tr * one + tt * two

        one = 35.*cr4 + 20.*cr3ct + 10.*cr2ct2 +  4.*crct3 +     ct4
        two = 15.*cr4 + 20.*cr3ct + 18.*cr2ct2 + 12.*crct3 +  5.*ct4
        thr =  5.*cr4 + 12.*cr3ct + 18.*cr2ct2 + 20.*crct3 + 15.*ct4
        fou =     cr4 +  4.*cr3ct + 10.*cr2ct2 + 20.*crct3 + 35.*ct4
        self._kg = tr**3. *one + tr**2.*tt *two + tr*tt**2. *thr + tt**3. *fou


    def _make_upsilon_values(self):

        # calculate upsilon values
        self._u0 = 1.0
        self._u1 = 1.0
        self._u2 = 1.0
        self._u3 = 1.0



class SymmetricAirfoil(PseudoPrismoid):
    """A default class for calculating and containing the mass properties of a
    Symmetric Airfoil.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # invoke init of parent
        PseudoPrismoid.__init__(self,input_dict)

        # save a coefficients for NACA 4-digit thickness distribution
        avals = input_dict.get("thickness_distribution_coefficients","open_trailing_edge")
        if avals == "open_trailing_edge":
            avals = [2.969, -1.260, -3.516, 2.843, -1.015]
        elif avals == "closed_trailing_edge":
            avals = [2.969, -1.260, -3.516, 2.843, -1.036]
        elif avals == "closed_hunsaker":
            avals = [2.980, -1.320, -3.286, 2.441, -0.815]
        elif not isinstance(avals,(np.ndarray,list)) and not (len(avals)==5):
            raise ImportError("Incorrect thickness distribution coefficients")
        self._a0 = avals[0] * 1.0
        self._a1 = avals[1] * 1.0
        self._a2 = avals[2] * 1.0
        self._a3 = avals[3] * 1.0
        self._a4 = avals[4] * 1.0

        # initialize the volume
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        cr2, ct2 = cr**2., ct**2.
        crct = cr * ct
        ka = tr * (3.*cr2 + 2.*crct + ct2) + tt *(cr2 + 2.*crct + 3.*ct2)
        a0 = self._a0
        a1 = self._a1
        a2 = self._a2
        a3 = self._a3
        a4 = self._a4
        u0 = 1./ 60.*( 40.*a0+ 30.*a1+ 20.*a2+ 15.*a3+ 12.*a4)
        self.volume = self._b / 12. * ka * u0


    def _make_upsilon_values(self):

        # initialize certain values for ease of calculation
        a0 = self._a0
        a1 = self._a1
        a2 = self._a2
        a3 = self._a3
        a4 = self._a4
        a03, a13, a23, a33, a43 = a0**3., a1**3., a2**3., a3**3., a4**3.
        a02, a12, a22, a32, a42 = a0**2., a1**2., a2**2., a3**2., a4**2.

        # calculate upsilon values
        self._u0 = 1./ 60.*( 40.*a0+ 30.*a1+ 20.*a2+ 15.*a3+ 12.*a4)
        self._u1 = 1./ 60.*( 56.*a0+ 50.*a1+ 40.*a2+ 33.*a3+ 28.*a4)
        self._u2 = 1./980.*(856.*a0+770.*a1+644.*a2+553.*a3+484.*a4)
        one1 = 2./5.*a03 + 1*a02*a1 + 3./4.*a02*a2 + 3./5.*a02*a3
        one2 = 1./2.*a02*a4 + 6./7.*a0*a12 + 4./3.*a0*a1*a2
        one3 = 12./11.*a0*a1*a3 + 12./13.*a0*a1*a4
        two1 = 6./11*a0*a22 + 12./13.*a0*a2*a3 + 4./5.*a0*a2*a4 +2./5.*a0*a32
        two2 = 12./17.*a0*a3*a4 + 6./19.*a0*a42 + 1./4.*a13 + 3./5.*a12*a2
        thr1 = 1./2.*a12*a3 + 3./7.*a12*a4 + 1./2.*a1*a22 + 6./7.*a1*a2*a3
        thr2 = 3./4.*a1*a2*a4 + 3./8.*a1*a32 + 2./3.*a1*a3*a4 + 3./10.*a1*a42
        fou1 = 1./7.*a23 + 3./8.*a22*a3 + 1./3.*a22*a4 + 1./3.*a2*a32
        fou2 = 3./5.*a2*a3*a4 + 3./11.*a2*a42 + 1./10.*a33 + 3./11.*a32*a4
        fou3 = 1./4.*a3*a42 + 1./13.*a43
        one = one1 + one2 + one3
        self._u3 = one + two1 + two2 + thr1 + thr2 + fou1 + fou2 + fou3



class DiamondAirfoil(PseudoPrismoid):
    """A default class for calculating and containing the mass properties of a
    Diamond Airfoil.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # invoke init of parent
        PseudoPrismoid.__init__(self,input_dict)

        # initialize the volume
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        cr2, ct2 = cr**2., ct**2.
        crct = cr * ct
        ka = tr * (3.*cr2 + 2.*crct + ct2) + tt *(cr2 + 2.*crct + 3.*ct2)
        u0 = 1. / 2.
        self.volume = self._b / 12. * ka * u0


    def _make_upsilon_values(self):

        # calculate upsilon values
        self._u0 = 1. /  2.
        self._u1 = ( 4. * self._xmt     + 1. ) / 6.
        self._u2 = ( 8. * self._xmt**2. + 3. ) / 14.
        self._u3 = 1. / 4.



class Rotor(Component):
    """A default class for calculating and containing the mass properties of a
    Rotor.

    Parameters
    ----------
    input_vars : dict , optional
        Must be a python dictionary
    """
    def __init__(self,input_dict={}):

        # invoke init of parent
        Component.__init__(self,input_dict)

        # retrieve additional info
        self._retrieve_info(input_dict)

        # initialize the volume
        Nb = self._Nb
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        rr,rt = self._rr,self._rt
        cr2, ct2 = cr**2., ct**2.
        crct = cr * ct
        ka = tr * (3.*cr2 + 2.*crct + ct2) + tt *(cr2 + 2.*crct + 3.*ct2)
        a0 = self._a0
        a1 = self._a1
        a2 = self._a2
        a3 = self._a3
        a4 = self._a4
        u0 = 1./ 60.*( 40.*a0+ 30.*a1+ 20.*a2+ 15.*a3+ 12.*a4)
        self.volume = 1. / 12. * Nb * ka * u0 * (rt-rr)


    def _retrieve_info(self,input_dict):
        """A function which retrieves the information and stores it globally.
        """
        
        # store variables from file input dictionary
        connect_to = input_dict.get("connect_to",{})
        x_cg = connect_to.get("dx",0.0)
        y_cg = connect_to.get("dy",0.0)
        z_cg = connect_to.get("dz",0.0)
        self.cg_location = np.array([[x_cg],[y_cg],[z_cg]])

        # rotor geometry
        self._Nb = input_dict.get("number_blades",2)
        self._rr = input_dict.get("diameter_hub", 1.0) / 2.0
        self._rt = input_dict.get("diameter_rotor", 2.0) / 2.0
        self._wx = input_dict.get("rotation_rate",0.0)
        chord = input_dict.get("chord",[[0.0,1.0],[1.0,1.0]])
        self._cr,self._ct = chord[0][1], chord[1][1]
        thickness = input_dict.get("thickness",[[0.0,0.12],[1.0,0.12]])
        self._tr,self._tt = thickness[0][1], thickness[1][1]


        # save a coefficients for NACA 4-digit thickness distribution
        avals = input_dict.get("thickness_distribution_coefficients","open_trailing_edge")
        if avals == "open_trailing_edge":
            avals = [2.969, -1.260, -3.516, 2.843, -1.015]
        elif avals == "closed_trailing_edge":
            avals = [2.969, -1.260, -3.516, 2.843, -1.036]
        elif avals == "closed_hunsaker":
            avals = [2.980, -1.320, -3.286, 2.441, -0.815]
        elif not isinstance(avals,(np.ndarray,list)) and not (len(avals)==5):
            raise ImportError("Incorrect thickness distribution coefficients")
        self._a0 = avals[0] * 1.0
        self._a1 = avals[1] * 1.0
        self._a2 = avals[2] * 1.0
        self._a3 = avals[3] * 1.0
        self._a4 = avals[4] * 1.0


    def _make_gamma_values(self):

        # initialize certian values for ease of calculation
        cr,ct = self._cr,self._ct
        tr,tt = self._tr,self._tt
        rr,rt = self._rr,self._rt
        # times values
        cr2, ct2 = cr**2., ct**2.
        cr3, ct3 = cr2*cr, ct2*ct
        cr4, ct4 = cr3*cr, ct3*ct
        cr5, ct5 = cr4*cr, ct4*ct
        cr6, ct6 = cr5*cr, ct5*ct
        cr5ct, cr4ct2, cr3ct3 = cr5*ct, cr4*ct2, cr3*ct3
        cr2ct4, crct5, crct = cr2*ct4, cr*ct5, cr*ct
        tr2, tt2 = tr**2., tt**2.
        tr3, tt3 = tr2*tr, tt2*tt
        tr2tt, trtt2 = tr2*tt, tr*tt2
        ln = np.log(rr/rt)

        # initialize certain vectors for ease of calculation
        C = np.array([
            [cr6],
            [cr5ct],
            [cr4ct2],
            [cr3ct3],
            [cr2ct4],
            [crct5],
            [ct6]
        ])
        T = np.array([
            [tr3],
            [tr2tt],
            [trtt2],
            [tt3]
        ])
        

        self._ya = -280. * ct6 * tt3

        Db = np.diag([1., 6., 5., 20., 15., 2., 1.])
        Ab = np.array([
            [35., 15., 5., 1.],
            [5., 5., 3., 1.],
            [5., 9., 9., 5.],
            [1., 3., 5., 5.],
            [1., 5., 15., 35.],
            [5., 45., 315., -3123.],
            [5., 105., -3123., 4329.]
        ])
        Bb = np.array([
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 1.],
            [0., 0., 1., -3.]
        ])
        self._yb = - np.matmul(np.matmul(Db,C).T,np.matmul(Ab + 840.*ln*Bb,T))[0,0]

        Dc = np.diag([1., 2., 5., 20., 5., 2., 3.])
        Ac = np.array([
            [90., 40., 14., 3.],
            [40., 42., 27., 10.],
            [14., 27., 30., 20.],
            [3., 10., 20., 30.],
            [10., 60., 270., -1089],
            [20., 270., -3267., 996],
            [10., -363., 332., 840.]
        ])
        Bc = np.array([
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 1.],
            [0., 0., 3., -4.],
            [0., 1./3., -4./3., 0.]
        ])
        self._yc = 4.*np.matmul(np.matmul(Dc,C).T,np.matmul(Ac + 420.*ln*Bc,T))[0,0]

        Dd = np.diag([1., 2., 15., 20., 5., 6., 3.])
        Ad = np.array([
            [120., 56., 21., 5.],
            [56., 63., 45., 20.],
            [7., 15., 20., 20.],
            [5., 20., 60., -117.],
            [20., 180., -1053., -21.],
            [20., -351., -21., 280.],
            [-39., -7., 280., 280.]
        ])
        Bd = np.array([
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 1.],
            [0., 0., 9., -7.],
            [0., 3., -7., 0.],
            [1./3., -7./3., 0., 0.]
        ])
        self._yd = -14.*np.matmul(np.matmul(Dd,C).T,np.matmul(Ad + 60.*ln*Bd,T))[0,0]
        
        De = np.diag([1., 6., 5., 20., 15., 2., 1.])
        Ae = np.array([
            [168., 84., 35., 10.],
            [28., 35., 30., 20.],
            [35., 90., 180., -174.],
            [10., 60., -174., -33.],
            [20., -174., -99., 70.],
            [-174., -297., 630., 280.],
            [-33., 210., 280., 420.]
        ])
        Be = np.array([
            [0., 0., 0., 0.],
            [0., 0., 0., 0.],
            [0., 0., 0., 2.],
            [0., 0., 2., -1.],
            [0., 2., -3., 0.],
            [2., -9., 0., 0.],
            [-1., 0., 0., 0.]
        ])
        self._ye = 28.*np.matmul(np.matmul(De,C).T,np.matmul(Ae + 60.*ln*Be,T))[0,0]
        
        Df = np.diag([1., 10., 25., 300., 25., 10., 1.])
        Af = np.array([
            [126., 70., 35., 15.],
            [14., 21., 27., -12.],
            [7., 27., -36., -12.],
            [1., -4., -4., 1.],
            [-12., -36., 27., 7.],
            [-12., 27., 21., 14.],
            [15., 35., 70., 126.]
        ])
        Bf = np.array([
            [0., 0., 0., 0.],
            [0., 0., 0., 1./5.],
            [0., 0., 3./5., -1./5.],
            [0., 1./15., -1./15., 0.],
            [1./5., -3./5., 0., 0.],
            [-1./5., 0., 0., 0.],
            [0., 0., 0., 0.]
        ])
        self._yf = -70.*np.matmul(np.matmul(Df,C).T,np.matmul(Af + 60.*ln*Bf,T))[0,0]
        
        Dg = np.flip(De)
        Ag = np.flip(Ae)
        Bg = np.flip(Be)
        self._yg = 28.*np.matmul(np.matmul(Dg,C).T,np.matmul(Ag - 60.*ln*Bg,T))[0,0]

        Dh = np.flip(Dd)
        Ah = np.flip(Ad)
        Bh = np.flip(Bd)
        self._yh = -14.*np.matmul(np.matmul(Dh,C).T,np.matmul(Ah - 60.*ln*Bh,T))[0,0]

        Di = np.flip(Dc)
        Ai = np.flip(Ac)
        Bi = np.flip(Bc)
        self._yi = 4.*np.matmul(np.matmul(Di,C).T,np.matmul(Ai - 420.*ln*Bi,T))[0,0]

        Dj = np.flip(Db)
        Aj = np.flip(Ab)
        Bj = np.flip(Bb)
        self._yj = - np.matmul(np.matmul(Dj,C).T,np.matmul(Aj - 840.*ln*Bj,T))[0,0]

        self._yk = -280. * cr6 * tr3

        self._yl = -tr*(10.*cr2 + 4.*crct + ct2) - tt*(2.*cr2 + 2.*crct + ct2)

        self._ym = tr*(6.*cr2 - ct2) - tt*(2.*crct + 3.*ct2)

        self._yn = tr*(3.*cr2 + 2.*crct) + tt*(cr2 - 6.*ct2)

        self._yo = tr*(cr2 + 2.*crct + 2.*ct2) + tt*(cr2 + 4.*crct + 10.*ct2)


    def _make_upsilon_values(self):

        # initialize certain values for ease of calculation
        a0 = self._a0
        a1 = self._a1
        a2 = self._a2
        a3 = self._a3
        a4 = self._a4

        # calculate upsilon values
        self._u0 = 1./ 60.*( 40.*a0+ 30.*a1+ 20.*a2+ 15.*a3+ 12.*a4)


    def get_mass_properties(self):
        """Method which returns mass, cg, I about cg rotated to total cframe"""

        # calculate kappa values
        self._make_gamma_values()

        # calculate upsilon values
        self._make_upsilon_values()

        # calculate mass
        if self._given_density_not_mass:
            self.mass = self.density * self.volume

        # create values for later use
        rr  , rt   = self._rr, self._rt
        rr2 , rt2  = rr *rr, rt *rt
        rr3 , rt3  = rr2*rr, rt2*rt
        rr4 , rt4  = rr3*rr, rt3*rt
        rr5 , rt5  = rr4*rr, rt4*rt
        rr6 , rt6  = rr5*rr, rt5*rt
        rr7 , rt7  = rr6*rr, rt6*rt
        rr8 , rt8  = rr7*rr, rt7*rt
        rr9 , rt9  = rr8*rr, rt8*rt
        rr10, rt10 = rr9*rr, rt9*rt

        # calculate values for inertia use
        T000 = self.volume * 1.0

        fro = self._Nb**3. * self._u0**3.
        one = self._ya*rr10    + self._yb*rr9*rt  + self._yc*rr8*rt2
        two = self._yd*rr7*rt3 + self._ye*rr6*rt4 + self._yf*rr5*rt5
        thr = self._yg*rr4*rt6 + self._yh*rr3*rt7 + self._yi*rr2*rt8
        fou = self._yj*rr *rt9 + self._yk*rt10 
        den = 13440. * np.pi**2. * rr * rt * (rr-rt)**9.
        T200 = fro * (one + two + thr + fou) / den

        fro = 1. / 120. * self._Nb * self._u0
        T020 = fro*(rr3*self._yl+rr2*rt*self._ym+rr*rt2*self._yn+rt3*self._yo)
        # print("T200 =", T200)
        # print("T020 =", T020)

        # calculate moments and products of inertia about the wing root c/4
        Ixxr = self.mass * 2. * T020 / T000

        Iyyr = self.mass * (T200 + T020) / T000
        
        Izzr = Iyyr * 1.
        
        Ixyr = 0.0
        Ixzr = 0.0
        Iyzr = 0.0

        # create inertia tensor
        self.inertia_tensor = np.array([
            [ Ixxr,-Ixyr,-Ixzr],
            [-Ixyr, Iyyr,-Iyzr],
            [-Ixzr,-Iyzr, Izzr]
        ])

        # calculate angular momentum
        self.angular_momentum = np.array([
            [Ixxr * self._wx],
            [0.0],
            [0.0]
        ])

        self.properties_dict = {
            "mass" : self.mass,
            "volume" : self.volume,
            "cg_location" : self.cg_location,
            "angular_momentum" : self.angular_momentum,
            "inertia_tensor" : self.inertia_tensor
        }
        return self.properties_dict