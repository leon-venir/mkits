# Strcuture builder
import ase.io
import ase.build
import numpy as np
from mkits.globle import *
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon


"""
Class
-----
canvas:
    Ideal tiling canvas.

Functions
---------
refine_struct:
    Refine the structures
stack_struct:
    stack 
ngt_tiling:
    The NGT ideal tiling canvas.
carbontube:
    The builder of carbon nanotube.

"""


def stack_struct(substrate, support, distance=1.7, vacuum=15):
    """ 
    """
    zmin = support.positions[:, 2].min()
    zmax = substrate.positions[:, 2].max()
    support.positions += (0, 0, zmax - zmin + distance)
    c = support.positions[:, 2].max()
    interface = substrate + support
    interface.set_cell([interface.get_cell()[0].tolist(),
                        interface.get_cell()[1].tolist(), 
                        [0, 0, c]])
    ase.build.add_vacuum(interface, vacuum)
    return interface

class canvas(object):
    """

    """

    def __init__(self, center_pattern) -> None:
        """
        Parameters
        ----------
        pattern
            square_pos = np.loadtxt("./square.xyz", 
                                    skiprows=2, 
                                    usecols=[1,2,3])
            square_atom = np.loadtxt("./square.xyz", 
                                     skiprows=2, 
                                     usecols=[0], 
                                     dtype='str')
            pattern_square = {
                "corner": 4,
                "coord": square_pos,
                "atoms": square_atom
            }

        Attributes
        ----------
        self.position
            coordinates of the atom position
        self.atom
        self.edge_index
        self.edge_coordinates: array
            1x6 array  
        """
        self.coord = center_pattern["coord"]
        self.atom = center_pattern["atoms"]
        self.edge_list = np.array([0])
        self.edge_used = np.array([0])
        self.edge_max = 0
        self.flip = np.array([0])
        self.pattern_index = np.array([0])
        self.pattern_center = np.zeros(3)
        self.edge_coordinates = np.zeros(6)
        pos_head_tail = np.vstack((self.coord[:center_pattern["corner"]], 
                                   self.coord[0]))
        for i in range(int(center_pattern["corner"])):
            self.edge_max += 1
            self.edge_list = np.append(self.edge_list, self.edge_max)
            self.edge_coordinates = np.vstack((self.edge_coordinates, 
                                               np.hstack((pos_head_tail[i], pos_head_tail[i+1]))))
            self.flip = np.append(self.flip, 0)
        self.truncated = False
        self.truncated_edge = np.zeros(2)
        self.export_value = {}
        
    def add_pattern(self, pattern, edge_index):
        """
        Parameters
        ----------
        pattern: class
            pattern to add
        edge_index: int
            The position where adding the pattern
        reverse: bool
            Whether reverse the pattern or not
        clockwise: bool
            rotate in clockwise
        """
        edge_num = pattern["corner"]
        edge_init_position = self.edge_coordinates[edge_index]
        
        flip = self.flip[edge_index]
        if not flip:
            flip = 1
        else:
            flip = 0
        
        # determine derection
        if flip:
            reverse = True
            if edge_init_position[1] <= edge_init_position[4]:
                clockwise = True
            else:
                clockwise = False
        else:
            reverse = False
            if edge_init_position[1] <= edge_init_position[4]:
                clockwise = False
            else:
                clockwise = True

        if reverse:
            vr = np.array([[-1, 0, 0],
                           [0, 1, 0],
                           [0, 0, 1]])
            coord = (vr @ pattern["coord"].T).T
        else:
            coord = pattern["coord"]

        # rotation
        vector_init = (coord[1]-coord[0])[0:2]
        vector_fina = edge_init_position[3:5]-edge_init_position[0:2]
        angle = vector_angle(vector_init, vector_fina, "rad")
        if clockwise:
            vx = np.array([[np.cos(angle), np.sin(angle), 0],
                           [-np.sin(angle), np.cos(angle), 0],
                           [0,0,1]])
        else:
            vx = np.array([[np.cos(angle), -np.sin(angle), 0],
                           [np.sin(angle), np.cos(angle), 0],
                           [0,0,1]])
        coord_rotated = (vx @ coord.T).T

        # translation
        coord_translated = coord_rotated + np.append(edge_init_position[0:2], 0)
        self.coord = np.vstack((self.coord, coord_translated))
        self.atom = np.hstack((self.atom, pattern["atoms"]))

        pos_head_tail = np.vstack((coord_translated[0:edge_num], 
                                   coord_translated[0]))
        self.edge_used = np.append(self.edge_used, edge_index)

        for i in range(int(pattern["corner"])):
            self.edge_max += 1
            self.edge_list = np.append(self.edge_list, self.edge_max)
            self.edge_coordinates = np.vstack((self.edge_coordinates,
                                               np.hstack((pos_head_tail[i], 
                                                          pos_head_tail[i+1]))))
            self.flip = np.append(self.flip, flip)
        
        # pattern index
        self.pattern_index = np.append(self.pattern_index, 
                                       self.pattern_index[-1]+1)
        self.pattern_center = np.vstack((self.pattern_center,
                                         np.average(coord_translated[:pattern["corner"]], 
                                                    axis=0)))

    def plot_canvas(self, 
                    show_pattern_index=0,
                    show_edge_index=True,
                    save_fig="none"):
        """
        Parameters
        ----------
        show_pattern_index
            optional [-1, 0, 1]
            0 : don't show the indices
            1 : show the indices in positive sequence
            -1: show the indices in inverted sequence 
        save_fig:
            if provide, save the figure with specific size. eg: [20,20]
        """
        fig = plt.figure(figsize=[10,10])
        if save_fig != "none":
            fig = plt.figure(figsize=save_fig)

        for i in range(len(self.edge_coordinates) - 1):
            if i+1 in self.edge_used:
                continue
            else:
                plt.arrow(self.edge_coordinates[i+1, 0],
                          self.edge_coordinates[i+1, 1], 
                          self.edge_coordinates[i+1, 3] - self.edge_coordinates[i+1, 0],
                          self.edge_coordinates[i+1, 4] - self.edge_coordinates[i+1, 1],
                          width = 0.1)
                if show_edge_index:
                    plt.text((self.edge_coordinates[i+1, 3] + self.edge_coordinates[i+1, 0])/2,
                            (self.edge_coordinates[i+1, 4] + self.edge_coordinates[i+1, 1])/2,
                            str(self.edge_list[i+1]))
        plt.axis("equal")
        if self.truncated:
            plt.plot(self.truncated_edge[:, 0], 
                     self.truncated_edge[:, 1], 
                     "--", 
                     color="tab:red", 
                     linewidth=2.0)
            for i in range(len(self.truncated_edge)-1):
                plt.text(self.truncated_edge[i, 0],
                         self.truncated_edge[i, 1],
                         str(i),
                         color="tab:red",
                         fontsize=18)
                
        if show_pattern_index == 0:
            pass
        elif show_pattern_index == 1 or show_pattern_index == -1:
            for i in range(1, len(self.pattern_index)):
                plt.text(self.pattern_center[i, 0],
                         self.pattern_center[i, 1],
                         str(self.pattern_index[i * show_pattern_index]),
                         color="tab:red",
                         fontsize=16)
                
        if save_fig != "none":
            plt.savefig("./output.png")

    def savexyz(self, out):
        """
        """
        with open(out, "w") as f:
            f.write("%d\n" % len(self.coord))
            f.write("\n")
            for i in range(len(self.coord)):
                f.write("%s%23.15f%23.15f%23.15f\n" % ("{:<2}".format(self.atom[i]), 
                                                     self.coord[i, 0],
                                                     self.coord[i, 1],
                                                     self.coord[i, 2]))

    def points_in_polygon(self, polygon):
        """ 
        """
        atom_to_del = []
        polygon = Polygon(polygon)
        for i in range(1, len(self.atom)):
            points = Point(self.coord[i, 0:2])
            if not points.within(polygon):
                atom_to_del.append(i)
        self.atom = np.delete(self.atom, list(set(atom_to_del)))
        self.coord = np.delete(self.coord, list(set(atom_to_del)), axis=0)
    
    def truncate(self, method, inpara):
        """ 
        Truncate the tiling with different methods

        Parameters
        ----------
        method: str
            The way to truncate the tiling. Option [cube, plane]
            cube: Truncate with planes parallel to axis
                  specify xmin, xmax, ymin, ymax as a list
                  eg: inpara = [xmin, xmax, ymin, ymax]
            edge: Get the vertex from edge number and the head[0]/tail[1]
                  and offset
                  eg: inpara = ([1,2,3], 
                                [0,1,0], 
                                [-0.1, -0.1, 0.1])

        """
        self.truncated = True
        if method == "cube":
            xmin, xmax, ymin, ymax = inpara
            truncated_edge = np.array([[xmin, ymin],
                                       [xmax, ymin],
                                       [xmax, ymax],
                                       [xmin, ymax]])
            self.points_in_polygon(truncated_edge)
            # truncated edge
            truncated_edge = np.vstack((truncated_edge, truncated_edge[0]))
            self.truncated_edge = truncated_edge

        elif method == "edge":
            truncated_edge = np.zeros(2)
            offsetx = np.zeros(1)
            offsety = np.zeros(1)
            for i in range(len(inpara[0])):
                if inpara[1][i] == 0:
                    edge_points = self.edge_coordinates[inpara[0][i]][0:2]
                elif inpara[1][i] == 1:
                    edge_points = self.edge_coordinates[inpara[0][i]][3:5]
                truncated_edge = np.vstack((truncated_edge, edge_points))
                offsetx = np.append(offsetx, inpara[2][i])
                offsety = np.append(offsety, inpara[3][i])
            
            truncated_edge = truncated_edge + (np.vstack((offsetx, offsety))).T
            self.points_in_polygon(truncated_edge[1:])
            self.truncated_edge = np.vstack((truncated_edge[1:], 
                                             truncated_edge[1]))
    
    def get_edge_point(self, edge, initend):
        """
        """
        edge_point_coord = np.zeros(2)
        for i in range(len(edge)):
            if initend[i] == 0:
                edge_point_coord = np.vstack((edge_point_coord, 
                                              self.edge_coordinates[0:2]))
            elif initend[i] == 1:
                edge_point_coord = np.vstack((edge_point_coord, 
                                              self.edge_coordinates[3:5]))
            else:
                print("0 represents initial point, 1 represents end point")
                exit()
        return edge_point_coord
    
    def reflection_transform(self, edge, initend):
        """ 
        Parameters
        ----------
        """
        edge_points_coord = self.get_edge_point(edge, initend)
        edge_coordinates_reflected = trans_reflection_xy(self.edge_coordinates, 
                                                         "2ps", 
                                                         tuple(edge_points_coord))
        coord_reflected = trans_reflection_xy(self.coord, 
                                              "2ps", 
                                              tuple(edge_points_coord))
        self.atom = np.hstack((self.atom, self.atom))
        self.edge_coordinates = np.vstack((self.edge_coordinates, 
                                           edge_coordinates_reflected))
        self.coord = np.vstack((self.coord, coord_reflected))
    
    def build_crystal(self, out, vacuum, a_index, b_index, fmt="vasp", sort=True):
        """
        Parameters
        ----------
        out: str
            The name of output file
        vacuum: float
        a_index: int from 1
        b_index: int from 1
        fmt: str
            option [vasp]            
        """
        vector_a = self.truncated_edge[a_index[1], 0:2] - self.truncated_edge[a_index[0], 0:2]
        vector_b = self.truncated_edge[b_index[1], 0:2] - self.truncated_edge[b_index[0], 0:2]
        lattice_a = np.linalg.norm(vector_a)
        lattice_b = np.linalg.norm(vector_b)

        if a_index[0] != b_index[0]:
            print("Lattice a and b must start from the same point.")
            exit()
        angle = vector_angle(vector_a, np.array([1,0]), "rad")

        # translation
        self.coord = self.coord - np.append(self.truncated_edge[a_index[0]], 0)

        if vector_a[1] > 0:
            vx = np.array([[np.cos(angle), np.sin(angle), 0],
                           [-np.sin(angle), np.cos(angle), 0],
                           [0,0,1]])
        else:
            vx = np.array([[np.cos(angle), -np.sin(angle), 0],
                           [np.sin(angle), np.cos(angle), 0],
                           [0,0,1]])
        self.coord = (vx @ self.coord.T).T

        vector_a = np.append(vector_a, 0)
        vector_b = np.append(vector_b, 0)
        vector_a = (vx @ vector_a.T).T
        vector_b = (vx @ vector_b.T).T

        # delete too close pairs in periodic crystals
        toocloseatom = []
        # in a-direction
        for i in range(len(self.atom)):
            if abs(self.coord[i, 0]) < 1:
                for j in range(len(self.atom)):
                    if np.linalg.norm(self.coord[i] + np.array([vector_a[0],0,0] - self.coord[j])) < 1:
                        toocloseatom.append(j)
        # in b-direction
        for i in range(len(self.atom)):
            if abs(self.coord[i, 1]) < 1:
                for j in range(len(self.atom)):
                    if np.linalg.norm(self.coord[i] + np.array([0,vector_b[1],0] - self.coord[j])) < 1:
                        toocloseatom.append(j)
        self.atom = np.delete(self.atom, list(set(toocloseatom)))
        self.coord = np.delete(self.coord, list(set(toocloseatom)), axis=0)

        # sort the atoms
        if sort:
            self.sort_atom()

        # ase atoms
        cell = ase.Atoms(self.atom, 
                         positions=self.coord.tolist(), 
                         #cell=[lattice_a, lattice_b, vacuum],
                         cell=[vector_a.tolist(), vector_b.tolist(), [0, 0, vacuum]],
                         pbc=[1,1,0])
        
        self.export_value["crystal_lattice_a"] = vector_a
        self.export_value["crystal_lattice_b"] = vector_b
        
        ase.io.write(out, cell, format=fmt)
    
    def sort_atom(self):
        """ """
        atom_index = np.array([symbol_map[_] for _ in self.atom])
        atom_index = np.argsort(atom_index)
        self.atom = np.array(self.atom)[atom_index]
        self.coord = self.coord[atom_index]
        
    def export_val(self):
        """ 
        """
        self.export_value["atom"] = self.atom
        self.export_value["coord"] = self.coord
        self.export_value["edge"] = self.edge_coordinates
        return self.export_value
    
    def del_neighbor(self):
        """ 
        """
        toocloseatompairs = [[-1,-1]]
        for i in range(len(self.coord)):
            for j in range(len(self.coord)):
                if i != j:
                    dis = np.linalg.norm(self.coord[j] - self.coord[i])
                    if dis < 1 and self.atom[i] == self.atom[j]:
                        toocloseatompairs += [[i, j]]

        # delete reversed elements
        seen = set()
        toocloseatompairs = [x for x in toocloseatompairs \
                             if tuple(x[::-1]) not in seen and not seen.add(tuple(x))]

        new_atom = np.array(["x"])
        new_coord = np.zeros(3)
        old_atom_index = []
        for pair in toocloseatompairs:
            i, j = pair
            if self.atom[i] != self.atom[j]:
                print("The atoms in this pair %d is different" % i)
                exit()
            old_atom_index.append(i)
            old_atom_index.append(j)
            new_atom = np.append(new_atom, self.atom[i])
            new_coord = np.vstack((new_coord, (self.coord[i]+self.coord[j])/2))

        self.atom = np.delete(self.atom, list(set(old_atom_index)))
        self.coord = np.delete(self.coord, list(set(old_atom_index)), axis=0)

        for i in range(1, len(new_coord)):
            putatom = True
            for j in self.coord:
                if np.linalg.norm(new_coord[i]-j) < 1:
                    putatom = False
                else:
                    continue
            if putatom:
                self.atom = np.hstack((self.atom, new_atom[i]))
                self.coord = np.vstack((self.coord, new_coord[i]))


def carbontube(m:int=5, n:int=5, l:int=1, out:str="out", fmt:str="xyz", **kwargs):
    """
    A scripts to build carbon nano tube

    Properties -----
    fmt: str vasp xyz

    Output --------- 

    """

    # bond length
    cc_bond = 1.42165
    thick = cc_bond * np.sqrt(3)

    # parameters
    radium = 0
    vacuum = 15
    coor_tran = 0

    if fmt == "xyz":
        pass
    elif fmt == "vasp":
        print(kwargs)
        try: 
            vacuum  = float(kwargs["vacuum"])
        except:
            print("No vacuum specified, use 15 angstrom. If you want other value, specify with vaccum=15")
            exit()

    if m == n:
        """ """
        alpha = (np.pi/m)
        gamma = np.arctan((np.sin(alpha/2))/(1/2+np.cos(alpha/2)))
        beta = alpha-2*gamma
        radium = cc_bond / 2 / np.sin(gamma)
        coor_tran = vacuum/2 + radium

        
        cood_alpha = 0
        cood_xy = np.array([[coor_tran,coor_tran+radium,thick/2]])
        for i in range(m):
            # 1st point
            cood_alpha += 2*gamma
            x = np.sin(cood_alpha) * radium
            y = np.cos(cood_alpha) * radium
            cood_xy = np.vstack((cood_xy, np.array([[x,y,thick/2]])))
            # 2nd point
            cood_alpha += beta
            x = np.sin(cood_alpha) * radium
            y = np.cos(cood_alpha) * radium
            cood_xy = np.vstack((cood_xy, np.array([[x,y,0]])))
            # 3rd point
            cood_alpha += 2*gamma
            x = np.sin(cood_alpha) * radium
            y = np.cos(cood_alpha) * radium
            cood_xy = np.vstack((cood_xy, np.array([[x,y,0]])))
            # 4th point
            cood_alpha += beta
            x = np.sin(cood_alpha) * radium
            y = np.cos(cood_alpha) * radium
            cood_xy = np.vstack((cood_xy, np.array([[x,y,thick/2]])))


    # write to file
    with open("./tube_m%dn%dl%d.%s" % (m, n, l, fmt), "w") as f:
        if fmt == "vasp":
            lattice = np.array([[vacuum + 2*radium + 0.0001, 0, 0],
                                [0, vacuum + 2*radium + 0.0001, 0],
                                [0, 0, thick]])
            f.write("tube_m%dn%dl%d\n" % (m, n, l))
            f.write("1.000\n")
            np.savetxt(f, lattice, fmt="%20.15f%20.15f%20.15f")
            f.write("C\n")
            f.write("%d\n" % len(cood_xy[1:,:]))
            f.write("Cartesian\n")
            np.savetxt(f, cood_xy[1:,:] + np.array([coor_tran, coor_tran, 0]), fmt="%20.15f")
        elif fmt == "xyz":
            f.write("%d\n" % len(cood_xy[1:]))
            f.write("tube_m%dn%dl%d\n" % (m, n, l))
            for i in range(1, len(cood_xy)):
                f.write("%s%20.15f%20.15f%20.15f\n" % ("C", cood_xy[i][0], cood_xy[i][1], cood_xy[i][2]))
        else:
            pass


def ngt_til(pattern_square, pattern_triangle, pattern_rhombus):
    """
    Parameters
    ----------
    """
    ngt_tiling = canvas(pattern_square)
    ngt_tiling.add_pattern(pattern_triangle, 1)
    ngt_tiling.add_pattern(pattern_triangle, 2)
    ngt_tiling.add_pattern(pattern_triangle, 3)
    ngt_tiling.add_pattern(pattern_triangle, 4)
    ngt_tiling.add_pattern(pattern_triangle, 9)
    ngt_tiling.add_pattern(pattern_triangle, 16)
    ngt_tiling.add_pattern(pattern_square, 13)
    ngt_tiling.add_pattern(pattern_square, 12)
    ngt_tiling.add_pattern(pattern_triangle, 26)
    ngt_tiling.add_pattern(pattern_triangle, 28)
    ngt_tiling.add_pattern(pattern_square, 10)
    ngt_tiling.add_pattern(pattern_triangle, 6)
    ngt_tiling.add_pattern(pattern_square, 15)
    ngt_tiling.add_pattern(pattern_triangle, 7)
    ngt_tiling.add_pattern(pattern_triangle, 40)
    ngt_tiling.add_pattern(pattern_triangle, 52)
    ngt_tiling.add_pattern(pattern_rhombus, 18)
    ngt_tiling.add_pattern(pattern_triangle, 59)
    ngt_tiling.add_pattern(pattern_square, 62)
    ngt_tiling.add_pattern(pattern_triangle, 58)
    ngt_tiling.add_pattern(pattern_triangle, 69)
    ngt_tiling.add_pattern(pattern_square, 72)
    ngt_tiling.add_pattern(pattern_triangle, 65)
    ngt_tiling.add_pattern(pattern_triangle, 63)
    ngt_tiling.add_pattern(pattern_triangle, 82)
    ngt_tiling.add_pattern(pattern_triangle, 75)
    ngt_tiling.add_pattern(pattern_rhombus, 89)
    ngt_tiling.add_pattern(pattern_square, 80)
    ngt_tiling.add_pattern(pattern_triangle, 96)
    ngt_tiling.add_pattern(pattern_square, 100)
    ngt_tiling.add_pattern(pattern_triangle, 103)
    ngt_tiling.add_pattern(pattern_triangle, 66)
    ngt_tiling.add_pattern(pattern_triangle, 79)
    ngt_tiling.add_pattern(pattern_triangle, 97)
    ngt_tiling.add_pattern(pattern_triangle, 99)
    ngt_tiling.add_pattern(pattern_triangle, 104)
    ngt_tiling.add_pattern(pattern_triangle, 106)
    ngt_tiling.add_pattern(pattern_triangle, 88)
    ngt_tiling.add_pattern(pattern_triangle, 76)
    ngt_tiling.add_pattern(pattern_square, 109)
    ngt_tiling.add_pattern(pattern_square, 113)
    ngt_tiling.add_pattern(pattern_square, 115)
    ngt_tiling.add_pattern(pattern_square, 119)
    ngt_tiling.add_pattern(pattern_square, 121)
    ngt_tiling.add_pattern(pattern_square, 125)
    ngt_tiling.add_pattern(pattern_square, 127)
    ngt_tiling.add_pattern(pattern_square, 131)
    ngt_tiling.add_pattern(pattern_square, 35)
    ngt_tiling.add_pattern(pattern_triangle, 39)
    ngt_tiling.add_pattern(pattern_triangle, 134)
    ngt_tiling.add_pattern(pattern_triangle, 133)
    ngt_tiling.add_pattern(pattern_triangle, 138)
    ngt_tiling.add_pattern(pattern_triangle, 139)
    ngt_tiling.add_pattern(pattern_triangle, 142)
    ngt_tiling.add_pattern(pattern_triangle, 141)
    ngt_tiling.add_pattern(pattern_triangle, 146)
    ngt_tiling.add_pattern(pattern_triangle, 147)
    ngt_tiling.add_pattern(pattern_triangle, 150)
    ngt_tiling.add_pattern(pattern_triangle, 149)
    ngt_tiling.add_pattern(pattern_triangle, 154)
    ngt_tiling.add_pattern(pattern_triangle, 155)
    ngt_tiling.add_pattern(pattern_triangle, 158)
    ngt_tiling.add_pattern(pattern_triangle, 157)
    ngt_tiling.add_pattern(pattern_triangle, 162)
    ngt_tiling.add_pattern(pattern_triangle, 163)
    ngt_tiling.add_pattern(pattern_triangle, 166)
    ngt_tiling.add_pattern(pattern_triangle, 167)
    ngt_tiling.add_pattern(pattern_triangle, 29)
    ngt_tiling.add_pattern(pattern_rhombus, 224)
    ngt_tiling.add_pattern(pattern_triangle, 45)
    ngt_tiling.add_pattern(pattern_triangle, 234)
    ngt_tiling.add_pattern(pattern_rhombus, 236)
    ngt_tiling.add_pattern(pattern_triangle, 46)
    ngt_tiling.add_pattern(pattern_rhombus, 170)
    ngt_tiling.add_pattern(pattern_triangle, 247)
    ngt_tiling.add_pattern(pattern_triangle, 251)
    ngt_tiling.add_pattern(pattern_rhombus, 253)
    ngt_tiling.add_pattern(pattern_triangle, 256)
    ngt_tiling.add_pattern(pattern_triangle, 244)
    ngt_tiling.add_pattern(pattern_triangle, 248)
    ngt_tiling.add_pattern(pattern_triangle, 173)
    ngt_tiling.add_pattern(pattern_square, 237)
    ngt_tiling.add_pattern(pattern_square, 233)
    ngt_tiling.add_pattern(pattern_square, 254)
    ngt_tiling.add_pattern(pattern_square, 250)
    ngt_tiling.add_pattern(pattern_square, 175)
    ngt_tiling.add_pattern(pattern_triangle, 274)
    ngt_tiling.add_pattern(pattern_triangle, 293)
    ngt_tiling.add_pattern(pattern_triangle, 296)
    ngt_tiling.add_pattern(pattern_triangle, 272)
    ngt_tiling.add_pattern(pattern_square, 301)
    ngt_tiling.add_pattern(pattern_square, 295)
    ngt_tiling.add_pattern(pattern_square, 299)
    ngt_tiling.add_pattern(pattern_square, 218)
    ngt_tiling.add_pattern(pattern_rhombus, 209)
    ngt_tiling.add_pattern(pattern_rhombus, 197)
    ngt_tiling.add_pattern(pattern_rhombus, 194)
    ngt_tiling.add_pattern(pattern_rhombus, 185)
    ngt_tiling.add_pattern(pattern_rhombus, 182)
    ngt_tiling.add_pattern(pattern_triangle, 264)
    ngt_tiling.add_pattern(pattern_triangle, 277)
    ngt_tiling.add_pattern(pattern_triangle, 304)
    ngt_tiling.add_pattern(pattern_triangle, 305)
    ngt_tiling.add_pattern(pattern_triangle, 306)
    ngt_tiling.add_pattern(pattern_triangle, 310)
    ngt_tiling.add_pattern(pattern_triangle, 309)
    ngt_tiling.add_pattern(pattern_triangle, 308)
    ngt_tiling.add_pattern(pattern_triangle, 313)
    ngt_tiling.add_pattern(pattern_triangle, 314)
    ngt_tiling.add_pattern(pattern_triangle, 318)
    ngt_tiling.add_pattern(pattern_triangle, 317)
    ngt_tiling.add_pattern(pattern_triangle, 316)
    ngt_tiling.add_pattern(pattern_triangle, 321)
    ngt_tiling.add_pattern(pattern_triangle, 322)
    ngt_tiling.add_pattern(pattern_triangle, 206)
    ngt_tiling.add_pattern(pattern_triangle, 280)
    ngt_tiling.add_pattern(pattern_triangle, 266)
    ngt_tiling.add_pattern(pattern_square, 343)
    ngt_tiling.add_pattern(pattern_square, 353)
    ngt_tiling.add_pattern(pattern_square, 364)
    ngt_tiling.add_pattern(pattern_rhombus, 374)
    ngt_tiling.add_pattern(pattern_triangle, 373)
    ngt_tiling.add_pattern(pattern_triangle, 380)
    ngt_tiling.add_pattern(pattern_triangle, 382)
    ngt_tiling.add_pattern(pattern_triangle, 408)
    ngt_tiling.add_pattern(pattern_triangle, 410)
    ngt_tiling.add_pattern(pattern_triangle, 402)
    ngt_tiling.add_pattern(pattern_triangle, 398)
    ngt_tiling.add_pattern(pattern_triangle, 419)
    ngt_tiling.add_pattern(pattern_square, 385)
    ngt_tiling.add_pattern(pattern_triangle, 434)
    ngt_tiling.add_pattern(pattern_triangle, 435)
    ngt_tiling.add_pattern(pattern_triangle, 325)
    ngt_tiling.add_pattern(pattern_triangle, 326)
    ngt_tiling.add_pattern(pattern_triangle, 330)
    ngt_tiling.add_pattern(pattern_triangle, 329)
    ngt_tiling.add_pattern(pattern_triangle, 333)
    ngt_tiling.add_pattern(pattern_triangle, 334)
    ngt_tiling.add_pattern(pattern_triangle, 338)
    ngt_tiling.add_pattern(pattern_triangle, 337)
    ngt_tiling.add_pattern(pattern_triangle, 288)
    ngt_tiling.add_pattern(pattern_triangle, 289)
    ngt_tiling.add_pattern(pattern_square, 391)
    ngt_tiling.add_pattern(pattern_triangle, 476)
    ngt_tiling.add_pattern(pattern_triangle, 475)
    ngt_tiling.add_pattern(pattern_triangle, 474)
    ngt_tiling.add_pattern(pattern_rhombus, 485)
    ngt_tiling.add_pattern(pattern_square, 482)
    ngt_tiling.add_pattern(pattern_square, 444)
    ngt_tiling.add_pattern(pattern_square, 451)
    ngt_tiling.add_pattern(pattern_square, 460)
    ngt_tiling.add_pattern(pattern_square, 465)
    ngt_tiling.add_pattern(pattern_square, 442)
    ngt_tiling.add_pattern(pattern_triangle, 438)
    ngt_tiling.add_pattern(pattern_rhombus, 439)
    ngt_tiling.add_pattern(pattern_square, 420)
    ngt_tiling.add_pattern(pattern_square, 379)
    ngt_tiling.add_pattern(pattern_square, 416)
    ngt_tiling.add_pattern(pattern_triangle, 523)
    ngt_tiling.add_pattern(pattern_triangle, 527)
    ngt_tiling.add_pattern(pattern_triangle, 530)
    ngt_tiling.add_pattern(pattern_triangle, 515)
    ngt_tiling.add_pattern(pattern_triangle, 540)
    ngt_tiling.add_pattern(pattern_triangle, 534)
    ngt_tiling.add_pattern(pattern_triangle, 538)
    ngt_tiling.add_pattern(pattern_triangle, 531)
    ngt_tiling.add_pattern(pattern_triangle, 512)
    ngt_tiling.add_pattern(pattern_triangle, 513)
    ngt_tiling.add_pattern(pattern_triangle, 496)
    ngt_tiling.add_pattern(pattern_triangle, 447)
    ngt_tiling.add_pattern(pattern_triangle, 457)
    ngt_tiling.add_pattern(pattern_triangle, 462)
    ngt_tiling.add_pattern(pattern_square, 543)
    ngt_tiling.add_pattern(pattern_triangle, 558)
    ngt_tiling.add_pattern(pattern_rhombus, 546)
    ngt_tiling.add_pattern(pattern_square, 556)
    ngt_tiling.add_pattern(pattern_square, 581)
    ngt_tiling.add_pattern(pattern_square, 567)
    ngt_tiling.add_pattern(pattern_square, 571)
    ngt_tiling.add_pattern(pattern_square, 573)
    ngt_tiling.add_pattern(pattern_triangle, 508)
    ngt_tiling.add_pattern(pattern_triangle, 507)
    ngt_tiling.add_pattern(pattern_triangle, 604)
    ngt_tiling.add_pattern(pattern_triangle, 603)
    ngt_tiling.add_pattern(pattern_square, 611)
    ngt_tiling.add_pattern(pattern_triangle, 620)
    ngt_tiling.add_pattern(pattern_triangle, 623)
    ngt_tiling.add_pattern(pattern_triangle, 491)
    ngt_tiling.add_pattern(pattern_rhombus, 626)
    ngt_tiling.add_pattern(pattern_triangle, 492)
    ngt_tiling.add_pattern(pattern_square, 627)

    ngt_tiling.add_pattern(pattern_triangle, 641)
    ngt_tiling.add_pattern(pattern_triangle, 640)
    ngt_tiling.add_pattern(pattern_triangle, 639)
    ngt_tiling.add_pattern(pattern_square, 610)
    ngt_tiling.add_pattern(pattern_triangle, 652)
    ngt_tiling.add_pattern(pattern_triangle, 653)
    ngt_tiling.add_pattern(pattern_triangle, 654)
    ngt_tiling.add_pattern(pattern_square, 663)
    ngt_tiling.add_pattern(pattern_triangle, 665)
    ngt_tiling.add_pattern(pattern_triangle, 666)
    ngt_tiling.add_pattern(pattern_rhombus, 660)
    ngt_tiling.add_pattern(pattern_triangle, 647)
    ngt_tiling.add_pattern(pattern_square, 680)
    ngt_tiling.add_pattern(pattern_triangle, 682)
    ngt_tiling.add_pattern(pattern_triangle, 683)
    ngt_tiling.add_pattern(pattern_triangle, 684)
    ngt_tiling.add_pattern(pattern_square, 614)
    ngt_tiling.add_pattern(pattern_triangle, 599)
    ngt_tiling.add_pattern(pattern_triangle, 600)
    ngt_tiling.add_pattern(pattern_triangle, 601)
    ngt_tiling.add_pattern(pattern_triangle, 616)
    ngt_tiling.add_pattern(pattern_triangle, 699)
    ngt_tiling.add_pattern(pattern_rhombus, 708)
    ngt_tiling.add_pattern(pattern_square, 712)
    ngt_tiling.add_pattern(pattern_square, 365)
    ngt_tiling.add_pattern(pattern_triangle, 592)
    ngt_tiling.add_pattern(pattern_triangle, 596)
    ngt_tiling.add_pattern(pattern_triangle, 595)
    ngt_tiling.add_pattern(pattern_triangle, 706)
    ngt_tiling.add_pattern(pattern_triangle, 593)
    ngt_tiling.add_pattern(pattern_triangle, 669)
    ngt_tiling.add_pattern(pattern_triangle, 693)
    ngt_tiling.add_pattern(pattern_triangle, 690)
    ngt_tiling.add_pattern(pattern_square, 646)
    ngt_tiling.add_pattern(pattern_square, 726)
    ngt_tiling.add_pattern(pattern_square, 730)
    ngt_tiling.add_pattern(pattern_square, 550)
    ngt_tiling.add_pattern(pattern_square, 422)
    ngt_tiling.add_pattern(pattern_triangle, 756)
    ngt_tiling.add_pattern(pattern_triangle, 755)
    ngt_tiling.add_pattern(pattern_triangle, 760)
    ngt_tiling.add_pattern(pattern_triangle, 732)
    ngt_tiling.add_pattern(pattern_triangle, 762)
    ngt_tiling.add_pattern(pattern_square, 553)
    ngt_tiling.add_pattern(pattern_square, 431)
    ngt_tiling.add_pattern(pattern_square, 689)
    ngt_tiling.add_pattern(pattern_rhombus, 779)
    ngt_tiling.add_pattern(pattern_rhombus, 552)
    ngt_tiling.add_pattern(pattern_rhombus, 741)
    ngt_tiling.add_pattern(pattern_triangle, 759)
    ngt_tiling.add_pattern(pattern_triangle, 720)
    ngt_tiling.add_pattern(pattern_triangle, 696)
    ngt_tiling.add_pattern(pattern_triangle, 672)
    ngt_tiling.add_pattern(pattern_triangle, 794)
    ngt_tiling.add_pattern(pattern_square, 742)
    ngt_tiling.add_pattern(pattern_triangle, 787)
    ngt_tiling.add_pattern(pattern_triangle, 588)
    ngt_tiling.add_pattern(pattern_triangle, 576)
    ngt_tiling.add_pattern(pattern_triangle, 766)
    ngt_tiling.add_pattern(pattern_triangle, 767)
    ngt_tiling.add_pattern(pattern_triangle, 724)
    ngt_tiling.add_pattern(pattern_triangle, 361)
    ngt_tiling.add_pattern(pattern_rhombus, 838)
    ngt_tiling.add_pattern(pattern_rhombus, 340)
    ngt_tiling.add_pattern(pattern_rhombus, 580)
    ngt_tiling.add_pattern(pattern_triangle, 764)
    ngt_tiling.add_pattern(pattern_triangle, 789)
    ngt_tiling.add_pattern(pattern_triangle, 723)
    ngt_tiling.add_pattern(pattern_triangle, 260)
    ngt_tiling.add_pattern(pattern_triangle, 394)
    ngt_tiling.add_pattern(pattern_triangle, 395)
    ngt_tiling.add_pattern(pattern_triangle, 752)
    ngt_tiling.add_pattern(pattern_triangle, 751)
    ngt_tiling.add_pattern(pattern_triangle, 826)
    ngt_tiling.add_pattern(pattern_square, 841)
    ngt_tiling.add_pattern(pattern_square, 359)
    ngt_tiling.add_pattern(pattern_square, 870)
    ngt_tiling.add_pattern(pattern_triangle, 790)
    ngt_tiling.add_pattern(pattern_triangle, 763)
    ngt_tiling.add_pattern(pattern_triangle, 864)
    ngt_tiling.add_pattern(pattern_triangle, 888)
    ngt_tiling.add_pattern(pattern_triangle, 840)
    ngt_tiling.add_pattern(pattern_triangle, 893)
    ngt_tiling.add_pattern(pattern_triangle, 894)
    ngt_tiling.add_pattern(pattern_triangle, 399)
    ngt_tiling.add_pattern(pattern_triangle, 896)
    ngt_tiling.add_pattern(pattern_triangle, 897)
    ngt_tiling.add_pattern(pattern_triangle, 898)
    ngt_tiling.add_pattern(pattern_triangle, 388)
    ngt_tiling.add_pattern(pattern_triangle, 479)
    ngt_tiling.add_pattern(pattern_square, 907)
    ngt_tiling.add_pattern(pattern_square, 481)
    ngt_tiling.add_pattern(pattern_square, 931)
    ngt_tiling.add_pattern(pattern_square, 702)
    ngt_tiling.add_pattern(pattern_square, 747)
    ngt_tiling.add_pattern(pattern_square, 748)
    ngt_tiling.add_pattern(pattern_rhombus, 844)
    ngt_tiling.add_pattern(pattern_rhombus, 350)
    ngt_tiling.add_pattern(pattern_triangle, 874)
    ngt_tiling.add_pattern(pattern_square, 912)
    ngt_tiling.add_pattern(pattern_triangle, 964)
    ngt_tiling.add_pattern(pattern_triangle, 890)
    ngt_tiling.add_pattern(pattern_triangle, 889)
    ngt_tiling.add_pattern(pattern_triangle, 909)
    ngt_tiling.add_pattern(pattern_triangle, 956)
    ngt_tiling.add_pattern(pattern_triangle, 947)
    ngt_tiling.add_pattern(pattern_rhombus, 906)
    ngt_tiling.add_pattern(pattern_square, 985)
    ngt_tiling.add_pattern(pattern_triangle, 1002)
    ngt_tiling.add_pattern(pattern_triangle, 984)
    ngt_tiling.add_pattern(pattern_triangle, 975)
    ngt_tiling.add_pattern(pattern_triangle, 979)
    ngt_tiling.add_pattern(pattern_triangle, 915)
    ngt_tiling.add_pattern(pattern_square, 982)
    ngt_tiling.add_pattern(pattern_square, 1011)
    ngt_tiling.add_pattern(pattern_triangle, 955)
    ngt_tiling.add_pattern(pattern_triangle, 577)
    ngt_tiling.add_pattern(pattern_triangle, 859)
    ngt_tiling.add_pattern(pattern_triangle, 771)
    ngt_tiling.add_pattern(pattern_triangle, 952)
    ngt_tiling.add_pattern(pattern_triangle, 858)
    ngt_tiling.add_pattern(pattern_square, 1042)
    ngt_tiling.add_pattern(pattern_triangle, 1043)
    ngt_tiling.add_pattern(pattern_square, 1036)
    ngt_tiling.add_pattern(pattern_rhombus, 1049)
    ngt_tiling.add_pattern(pattern_triangle, 1052)
    ngt_tiling.add_pattern(pattern_square, 928)
    ngt_tiling.add_pattern(pattern_square, 972)
    ngt_tiling.add_pattern(pattern_square, 919)
    ngt_tiling.add_pattern(pattern_square, 922)
    ngt_tiling.add_pattern(pattern_square, 1007)
    ngt_tiling.add_pattern(pattern_square, 882)
    ngt_tiling.add_pattern(pattern_triangle, 1001)
    ngt_tiling.add_pattern(pattern_triangle, 941)
    ngt_tiling.add_pattern(pattern_triangle, 1079)
    ngt_tiling.add_pattern(pattern_triangle, 1020)
    ngt_tiling.add_pattern(pattern_triangle, 1021)
    ngt_tiling.add_pattern(pattern_triangle, 1073)
    ngt_tiling.add_pattern(pattern_triangle, 1069)
    ngt_tiling.add_pattern(pattern_triangle, 1068)
    ngt_tiling.add_pattern(pattern_triangle, 1067)
    ngt_tiling.add_pattern(pattern_triangle, 428)
    ngt_tiling.add_pattern(pattern_triangle, 877)
    ngt_tiling.add_pattern(pattern_triangle, 924)
    ngt_tiling.add_pattern(pattern_triangle, 948)
    ngt_tiling.add_pattern(pattern_triangle, 880)
    ngt_tiling.add_pattern(pattern_triangle, 636)
    ngt_tiling.add_pattern(pattern_triangle, 943)
    ngt_tiling.add_pattern(pattern_triangle, 493)
    ngt_tiling.add_pattern(pattern_triangle, 1025)
    ngt_tiling.add_pattern(pattern_triangle, 1024)
    ngt_tiling.add_pattern(pattern_square, 1100)
    ngt_tiling.add_pattern(pattern_rhombus, 1142)
    ngt_tiling.add_pattern(pattern_triangle, 944)
    ngt_tiling.add_pattern(pattern_rhombus, 933)
    ngt_tiling.add_pattern(pattern_triangle, 1096)
    ngt_tiling.add_pattern(pattern_triangle, 1072)
    ngt_tiling.add_pattern(pattern_triangle, 1076)
    ngt_tiling.add_pattern(pattern_triangle, 1106)
    ngt_tiling.add_pattern(pattern_rhombus, 1166)
    ngt_tiling.add_pattern(pattern_square, 1103)
    ngt_tiling.add_pattern(pattern_square, 1108)
    ngt_tiling.add_pattern(pattern_triangle, 1175)
    ngt_tiling.add_pattern(pattern_triangle, 1177)
    ngt_tiling.add_pattern(pattern_triangle, 1180)
    ngt_tiling.add_pattern(pattern_triangle, 1179)
    ngt_tiling.add_pattern(pattern_triangle, 1109)
    ngt_tiling.add_pattern(pattern_triangle, 959)
    ngt_tiling.add_pattern(pattern_square, 1111)
    ngt_tiling.add_pattern(pattern_square, 1112)
    ngt_tiling.add_pattern(pattern_square, 1196)
    ngt_tiling.add_pattern(pattern_square, 993)
    ngt_tiling.add_pattern(pattern_square, 1132)
    ngt_tiling.add_pattern(pattern_square, 1136)
    ngt_tiling.add_pattern(pattern_triangle, 1211)
    ngt_tiling.add_pattern(pattern_triangle, 1202)
    ngt_tiling.add_pattern(pattern_triangle, 1201)
    ngt_tiling.add_pattern(pattern_triangle, 1207)
    ngt_tiling.add_pattern(pattern_triangle, 1118)
    ngt_tiling.add_pattern(pattern_square, 1235)
    ngt_tiling.add_pattern(pattern_rhombus, 1229)
    ngt_tiling.add_pattern(pattern_triangle, 1246)
    ngt_tiling.add_pattern(pattern_triangle, 1245)
    ngt_tiling.add_pattern(pattern_triangle, 1242)
    ngt_tiling.add_pattern(pattern_triangle, 1241)
    ngt_tiling.add_pattern(pattern_square, 1251)
    ngt_tiling.add_pattern(pattern_square, 1255)
    ngt_tiling.add_pattern(pattern_rhombus, 1258)
    ngt_tiling.add_pattern(pattern_triangle, 1261)
    ngt_tiling.add_pattern(pattern_triangle, 1266)
    ngt_tiling.add_pattern(pattern_triangle, 1264)
    ngt_tiling.add_pattern(pattern_triangle, 1265)
    ngt_tiling.add_pattern(pattern_triangle, 1064)
    ngt_tiling.add_pattern(pattern_triangle, 1213)
    ngt_tiling.add_pattern(pattern_triangle, 1214)
    ngt_tiling.add_pattern(pattern_triangle, 1215)
    ngt_tiling.add_pattern(pattern_triangle, 1219)
    ngt_tiling.add_pattern(pattern_triangle, 1217)
    ngt_tiling.add_pattern(pattern_triangle, 1028)
    ngt_tiling.add_pattern(pattern_triangle, 1284)
    ngt_tiling.add_pattern(pattern_square, 1276)
    ngt_tiling.add_pattern(pattern_square, 1282)
    ngt_tiling.add_pattern(pattern_triangle, 1278)
    ngt_tiling.add_pattern(pattern_square, 1306)
    ngt_tiling.add_pattern(pattern_triangle, 1291)
    ngt_tiling.add_pattern(pattern_triangle, 1308)
    ngt_tiling.add_pattern(pattern_triangle, 1309)
    ngt_tiling.add_pattern(pattern_triangle, 1310)
    ngt_tiling.add_pattern(pattern_triangle, 1297)
    ngt_tiling.add_pattern(pattern_triangle, 1222)
    ngt_tiling.add_pattern(pattern_triangle, 1313)
    ngt_tiling.add_pattern(pattern_triangle, 1314)
    ngt_tiling.add_pattern(pattern_rhombus, 1330)
    ngt_tiling.add_pattern(pattern_rhombus, 1342)
    ngt_tiling.add_pattern(pattern_rhombus, 1294)
    ngt_tiling.add_pattern(pattern_rhombus, 1129)
    ngt_tiling.add_pattern(pattern_triangle, 1321)
    ngt_tiling.add_pattern(pattern_triangle, 1344)
    ngt_tiling.add_pattern(pattern_square, 1324)
    ngt_tiling.add_pattern(pattern_triangle, 1356)
    ngt_tiling.add_pattern(pattern_triangle, 1369)
    ngt_tiling.add_pattern(pattern_triangle, 1370)
    ngt_tiling.add_pattern(pattern_triangle, 1371)
    ngt_tiling.add_pattern(pattern_triangle, 1084)
    ngt_tiling.add_pattern(pattern_square, 1335)
    ngt_tiling.add_pattern(pattern_triangle, 1218)
    ngt_tiling.add_pattern(pattern_square, 1299)
    ngt_tiling.add_pattern(pattern_triangle, 1395)
    ngt_tiling.add_pattern(pattern_triangle, 1396)
    ngt_tiling.add_pattern(pattern_triangle, 1397)
    ngt_tiling.add_pattern(pattern_triangle, 1083)
    ngt_tiling.add_pattern(pattern_triangle, 1409)
    ngt_tiling.add_pattern(pattern_rhombus, 1411)
    ngt_tiling.add_pattern(pattern_triangle, 822)
    ngt_tiling.add_pattern(pattern_square, 1418)
    return ngt_tiling


def sigma_til(pattern_square, pattern_triangle, pattern_rhombus):
    """
    Parameters
    ----------
    """
    sigma = canvas(pattern_triangle)
    sigma.add_pattern(pattern_triangle, 3)
    sigma.add_pattern(pattern_square, 2)
    sigma.add_pattern(pattern_square, 1)
    sigma.add_pattern(pattern_triangle, 10)
    sigma.add_pattern(pattern_triangle, 12)
    sigma.add_pattern(pattern_square, 5)
    sigma.add_pattern(pattern_triangle, 14)
    sigma.add_pattern(pattern_square, 6)
    sigma.add_pattern(pattern_triangle, 8)
    return sigma


def approx1_til(pattern_square, pattern_triangle, pattern_rhombus):
    """
    Parameters
    ----------
    """
    approx1 = canvas(pattern_triangle)
    approx1.add_pattern(pattern_square, 2)
    approx1.add_pattern(pattern_square, 3)
    approx1.add_pattern(pattern_triangle, 5)
    approx1.add_pattern(pattern_triangle, 11)
    approx1.add_pattern(pattern_square, 13)
    approx1.add_pattern(pattern_square, 17)
    approx1.add_pattern(pattern_triangle, 25)
    approx1.add_pattern(pattern_square, 28)
    approx1.add_pattern(pattern_triangle, 20)
    approx1.add_pattern(pattern_triangle, 35)
    approx1.add_pattern(pattern_triangle, 24)
    approx1.add_pattern(pattern_triangle, 30)
    approx1.add_pattern(pattern_triangle, 23)
    approx1.add_pattern(pattern_triangle, 10)
    approx1.add_pattern(pattern_triangle, 6)
    approx1.add_pattern(pattern_triangle, 53)
    approx1.add_pattern(pattern_square, 1)
    approx1.add_pattern(pattern_triangle, 7)
    approx1.add_pattern(pattern_triangle, 62)
    approx1.add_pattern(pattern_triangle, 9)
    approx1.add_pattern(pattern_triangle, 60)
    return approx1


def approx2_til(pattern_square, pattern_triangle, pattern_rhombus):
    """
    Parameters
    ----------
    """
    approx2 = canvas(pattern_square)

    approx2.add_pattern(pattern_triangle, 1)
    approx2.add_pattern(pattern_triangle, 2)
    approx2.add_pattern(pattern_triangle, 3)
    approx2.add_pattern(pattern_triangle, 4)
    approx2.add_pattern(pattern_triangle, 13)
    approx2.add_pattern(pattern_triangle, 6)
    approx2.add_pattern(pattern_square, 12)
    approx2.add_pattern(pattern_square, 19)
    approx2.add_pattern(pattern_square, 9)
    approx2.add_pattern(pattern_square, 10)
    approx2.add_pattern(pattern_triangle, 30)
    approx2.add_pattern(pattern_triangle, 33)
    approx2.add_pattern(pattern_triangle, 32)
    approx2.add_pattern(pattern_triangle, 38)
    approx2.add_pattern(pattern_rhombus, 7)
    approx2.add_pattern(pattern_triangle, 24)
    approx2.add_pattern(pattern_triangle, 25)
    approx2.add_pattern(pattern_triangle, 26)
    approx2.add_pattern(pattern_rhombus, 63)
    approx2.add_pattern(pattern_triangle, 53)
    approx2.add_pattern(pattern_triangle, 67)
    approx2.add_pattern(pattern_triangle, 72)
    approx2.add_pattern(pattern_square, 69)
    approx2.add_pattern(pattern_square, 21)
    approx2.add_pattern(pattern_triangle, 54)
    approx2.add_pattern(pattern_triangle, 86)
    approx2.add_pattern(pattern_triangle, 37)
    approx2.add_pattern(pattern_triangle, 82)
    return approx2


def approx3_til(pattern_square, pattern_triangle, pattern_rhombus):
    """
    Parameters
    ----------
    """
    approx3 = canvas(pattern_square)
    approx3.add_pattern(pattern_triangle, 1)
    approx3.add_pattern(pattern_triangle, 2)
    approx3.add_pattern(pattern_triangle, 3)
    approx3.add_pattern(pattern_triangle, 4)
    approx3.add_pattern(pattern_triangle, 12)
    approx3.add_pattern(pattern_triangle, 13)
    approx3.add_pattern(pattern_square, 21)
    approx3.add_pattern(pattern_triangle, 24)
    approx3.add_pattern(pattern_rhombus, 29)
    approx3.add_pattern(pattern_rhombus, 18)
    approx3.add_pattern(pattern_triangle, 37)
    approx3.add_pattern(pattern_triangle, 32)
    approx3.add_pattern(pattern_square, 19)
    approx3.add_pattern(pattern_triangle, 35)
    approx3.add_pattern(pattern_triangle, 15)
    approx3.add_pattern(pattern_square, 7)
    approx3.add_pattern(pattern_square, 6)
    approx3.add_pattern(pattern_triangle, 46)
    approx3.add_pattern(pattern_triangle, 64)
    approx3.add_pattern(pattern_rhombus, 66)
    approx3.add_pattern(pattern_triangle, 56)
    approx3.add_pattern(pattern_triangle, 57)
    approx3.add_pattern(pattern_triangle, 59)
    approx3.add_pattern(pattern_triangle, 61)
    approx3.add_pattern(pattern_triangle, 25)
    approx3.add_pattern(pattern_triangle, 85)
    approx3.add_pattern(pattern_rhombus, 83)
    approx3.add_pattern(pattern_triangle, 91)
    return approx3


def bigapp_til(pattern_square, pattern_triangle, pattern_rhombus):
    """
    Parameters
    ----------
    """
    bigapprox = canvas(pattern_square)
    bigapprox.add_pattern(pattern_triangle, 1)
    bigapprox.add_pattern(pattern_triangle, 2)
    bigapprox.add_pattern(pattern_triangle, 3)
    bigapprox.add_pattern(pattern_triangle, 4)
    bigapprox.add_pattern(pattern_square, 6)
    bigapprox.add_pattern(pattern_square, 7)
    bigapprox.add_pattern(pattern_square, 9)
    bigapprox.add_pattern(pattern_triangle, 22)
    bigapprox.add_pattern(pattern_triangle, 23)
    bigapprox.add_pattern(pattern_triangle, 24)
    bigapprox.add_pattern(pattern_triangle, 18)
    bigapprox.add_pattern(pattern_triangle, 19)
    bigapprox.add_pattern(pattern_triangle, 20)
    bigapprox.add_pattern(pattern_triangle, 26)
    bigapprox.add_pattern(pattern_triangle, 27)
    bigapprox.add_pattern(pattern_triangle, 28)
    bigapprox.add_pattern(pattern_rhombus, 16)
    bigapprox.add_pattern(pattern_rhombus, 34)
    bigapprox.add_pattern(pattern_rhombus, 46)
    bigapprox.add_pattern(pattern_triangle, 55)
    bigapprox.add_pattern(pattern_triangle, 58)
    bigapprox.add_pattern(pattern_triangle, 59)
    bigapprox.add_pattern(pattern_triangle, 62)
    bigapprox.add_pattern(pattern_triangle, 63)
    bigapprox.add_pattern(pattern_triangle, 33)
    bigapprox.add_pattern(pattern_triangle, 43)
    bigapprox.add_pattern(pattern_triangle, 66)
    bigapprox.add_pattern(pattern_triangle, 67)
    bigapprox.add_pattern(pattern_triangle, 52)
    bigapprox.add_pattern(pattern_triangle, 75)
    bigapprox.add_pattern(pattern_triangle, 73)
    bigapprox.add_pattern(pattern_square, 69)
    bigapprox.add_pattern(pattern_square, 76)
    bigapprox.add_pattern(pattern_square, 78)
    bigapprox.add_pattern(pattern_rhombus, 102)
    bigapprox.add_pattern(pattern_square, 99)
    bigapprox.add_pattern(pattern_triangle, 88)
    bigapprox.add_pattern(pattern_triangle, 91)
    bigapprox.add_pattern(pattern_rhombus, 125)
    bigapprox.add_pattern(pattern_triangle, 118)
    bigapprox.add_pattern(pattern_triangle, 117)
    bigapprox.add_pattern(pattern_triangle, 110)
    bigapprox.add_pattern(pattern_triangle, 111)
    bigapprox.add_pattern(pattern_triangle, 121)
    bigapprox.add_pattern(pattern_triangle, 122)
    bigapprox.add_pattern(pattern_triangle, 114)
    bigapprox.add_pattern(pattern_square, 142)
    bigapprox.add_pattern(pattern_square, 141)
    bigapprox.add_pattern(pattern_square, 96)
    bigapprox.add_pattern(pattern_square, 97)
    bigapprox.add_pattern(pattern_triangle, 106)
    bigapprox.add_pattern(pattern_triangle, 169)
    bigapprox.add_pattern(pattern_triangle, 170)
    bigapprox.add_pattern(pattern_triangle, 164)
    bigapprox.add_pattern(pattern_triangle, 165)
    bigapprox.add_pattern(pattern_square, 37)
    bigapprox.add_pattern(pattern_square, 39)
    bigapprox.add_pattern(pattern_square, 90)
    bigapprox.add_pattern(pattern_square, 129)
    bigapprox.add_pattern(pattern_triangle, 162)
    bigapprox.add_pattern(pattern_triangle, 199)
    bigapprox.add_pattern(pattern_triangle, 200)
    bigapprox.add_pattern(pattern_triangle, 201)
    bigapprox.add_pattern(pattern_triangle, 189)
    bigapprox.add_pattern(pattern_triangle, 81)
    bigapprox.add_pattern(pattern_triangle, 150)
    bigapprox.add_pattern(pattern_triangle, 158)
    bigapprox.add_pattern(pattern_triangle, 160)
    bigapprox.add_pattern(pattern_rhombus, 204)
    bigapprox.add_pattern(pattern_triangle, 161)
    bigapprox.add_pattern(pattern_square, 84)
    bigapprox.add_pattern(pattern_square, 218)
    bigapprox.add_pattern(pattern_square, 173)
    bigapprox.add_pattern(pattern_triangle, 247)
    bigapprox.add_pattern(pattern_triangle, 245)
    bigapprox.add_pattern(pattern_triangle, 237)
    bigapprox.add_pattern(pattern_triangle, 238)
    bigapprox.add_pattern(pattern_triangle, 188)
    bigapprox.add_pattern(pattern_square, 126)
    bigapprox.add_pattern(pattern_triangle, 264)
    bigapprox.add_pattern(pattern_triangle, 266)
    bigapprox.add_pattern(pattern_triangle, 241)
    return bigapprox


def hexapp_til(pattern_square, pattern_triangle, pattern_rhombus):
    """
    Parameters
    ----------
    """
    hexapprox = canvas(pattern_rhombus)
    hexapprox.add_pattern(pattern_triangle, 1)
    hexapprox.add_pattern(pattern_triangle, 2)
    hexapprox.add_pattern(pattern_triangle, 3)
    hexapprox.add_pattern(pattern_triangle, 4)
    hexapprox.add_pattern(pattern_triangle, 7)
    hexapprox.add_pattern(pattern_triangle, 13)
    hexapprox.add_pattern(pattern_triangle, 18)
    hexapprox.add_pattern(pattern_triangle, 21)
    hexapprox.add_pattern(pattern_square, 15)
    hexapprox.add_pattern(pattern_square, 16)
    hexapprox.add_pattern(pattern_square, 9)
    hexapprox.add_pattern(pattern_square, 10)
    hexapprox.add_pattern(pattern_triangle, 22)
    hexapprox.add_pattern(pattern_triangle, 19)
    hexapprox.add_pattern(pattern_triangle, 35)
    hexapprox.add_pattern(pattern_triangle, 36)
    hexapprox.add_pattern(pattern_triangle, 30)
    hexapprox.add_pattern(pattern_triangle, 31)
    hexapprox.add_pattern(pattern_triangle, 43)
    hexapprox.add_pattern(pattern_triangle, 44)
    hexapprox.add_pattern(pattern_triangle, 38)
    hexapprox.add_pattern(pattern_triangle, 39)
    hexapprox.add_pattern(pattern_square, 47)
    hexapprox.add_pattern(pattern_square, 50)
    hexapprox.add_pattern(pattern_triangle, 76)
    hexapprox.add_pattern(pattern_triangle, 77)
    hexapprox.add_pattern(pattern_triangle, 78)
    hexapprox.add_pattern(pattern_triangle, 80)
    hexapprox.add_pattern(pattern_triangle, 81)
    hexapprox.add_pattern(pattern_triangle, 82)
    hexapprox.add_pattern(pattern_rhombus, 27)
    hexapprox.add_pattern(pattern_rhombus, 46)
    hexapprox.add_pattern(pattern_rhombus, 24)
    hexapprox.add_pattern(pattern_rhombus, 49)
    hexapprox.add_pattern(pattern_square, 93)
    hexapprox.add_pattern(pattern_square, 84)
    hexapprox.add_pattern(pattern_triangle, 122)
    hexapprox.add_pattern(pattern_triangle, 118)
    hexapprox.add_pattern(pattern_triangle, 120)
    hexapprox.add_pattern(pattern_triangle, 124)
    hexapprox.add_pattern(pattern_square, 52)
    hexapprox.add_pattern(pattern_square, 64)
    hexapprox.add_pattern(pattern_triangle, 56)
    hexapprox.add_pattern(pattern_triangle, 68)
    return hexapprox