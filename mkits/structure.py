# -*- coding: utf-8 -*

import numpy as np
import matplotlib.pyplot as plt
import spglib as spg
import os
import sys
import logging
from mkits.functions import *
from mkits.database import *
from mkits.qe import * 
from mkits import qe
from mkits import database
from mkits import functions


"""
Class
-----

Functions
---------


"""



class struct(object):
    """
    DESCRIPTION:
    ------------
    Structure Object.

    PARAMETERS:
    -----------
    inp, string
        The name of the input structure 

    Attributs:
    ----------
    inp
        The name of the input file.
    calculator: string
        The calculator of the class. Avaiable options:
        poscar - POSCAR Version 6.X.X structures
        qein  - 
        qeout - 
        wien2k - 
    total_atom: int
        The total number of atoms in structure.
    lattice9: 3x3 numpy array
        Cartesian lattice [xx xy xz yz yy yz zx zy zz].
    lattice6: 6x1 numpy array
        Direct lattice parameters: [a b c alpha beta gamma].
    position: 10*n numpy array
        [atomic_index, frac_x, _y, _z, cart_x, _y, _z, dyn_x, _y, _z]
    
    Functions:
    ----------
    __parse_calculator:
        Parse the calculator from the input file.
    __parse_structure:
        Parse the structure.


    return_dict
    """

    def __init__(self, inp): 
        self.inp = inp
        self.calculator = "none"
        self.title = "none"
        self.lattice9 = np.array([[]])
        self.lattice6 = np.array([])
        self.position = np.zeros([1, 10])
        self.total_atom = 0

        try:
            if type(inp) == dict:
                #self.struct_dict = inp              
                self.calculator = inp["calculator"]
            elif type(inp) == str and "\n" in inp:
                lexit("The string func hasn't been added yet.")
            elif inp == "none":
                pass
            else:
                with open(inp, "r") as f: 
                    _struct_lines = f.readlines()
                self.__parse_calculator(_struct_lines)
                self.__parse_structure(_struct_lines)
                self.sort_atoms()
        except:
            lexit("Cannot find the structure file: ", inp)

        #self.atoms_sequence = listcross(self.atom_type, self.atom_num)
        #self.atoms_index_seq = listcross(self.atom_index, self.atom_num)

    
    def __parse_dict(self):
        """

        """
        pass


    def __parse_calculator(self, lines):
        """ 
        Parse the calculator from the string list.

        Parameters:
        -----------

        Return: 
        -------
        A string contains calculator:
            poscar
            outcar
            vaspxml
            vasph5out
            wien
            qein
            qeout
        """

        # qeout
        if "open-source Quantum ESPRESSO" in lines[3]: 
            self.calculator = "qeout"
        # qein
        elif ("&SYSTEM" in "".join(lines[:]) or
              "&system" in "".join(lines[:])):
            self.calculator = "qein"
        # poscar
        elif len(lines[2].split()) == 3:
            self.calculator = "poscar"
        # cif 
        elif "_cell_length_a" in "".join(lines) and "_cell_length_b" in "".join(lines):
            self.calculator = "cif"


    def __parse_structure(self, lines):
        """
        DESCRIPTION:
        -----------
        Parse 

        PARAMETERS:
        -----------
        lines:string list
        
        UPDATE:
        -------
        """

        if self.calculator == "cif":
            pass

        elif self.calculator == "poscar":
            self.title = lines[0][:-1]
            _ratio = float(lines[1][:-1])
            _lattice_x = list(map(float, lines[2].split()))
            _lattice_y = list(map(float, lines[3].split()))
            _lattice_z = list(map(float, lines[4].split()))
            self.lattice9 = np.array([_lattice_x,
                                      _lattice_y,
                                      _lattice_z]) * _ratio
            self.lattice6 = lattice_conversion(self.lattice9)
            # read elements from 6th, 7th line, int -> vasp5, str -> vasp6
            try:
                _atom_type = lines[5].split()
                _atom_num = np.array(list(map(int, lines[6].split())))
                _atoms = listcross(_atom_type, _atom_num)
                self.total_atom = int(np.sum(_atom_num))
            except:
                lexit("POSCAR 5.X.X not supported yet.")
            _atom_index = np.array(
                [database.symbol_map[i] for i in _atoms]
            )
            _atom_index = _atom_index.reshape(len(_atom_index), 1)
            
            # check dyn on 8th line
            if (rmspace(lines[7]).upper()[0] == "S" and
                len(lines[9].split()) >= 6):
                
                # capitalize f t t and replace with 0/1
                for _ in range(self.total_atom):
                    lines[9+_] = (lines[9+_].upper()).replace("F", "0")
                    lines[9+_] = (lines[9+_].upper()).replace("T", "1")

                # read coordinates
                _position = np.array([list(map(float,(i[:-1].split())[:3])) \
                                for i in lines[9:9+self.total_atom]])
                _dyn = np.array([list(map(float,(i[:-1].split())[3:])) \
                                for i in lines[9:9+self.total_atom]])
            elif (rmspace(lines[7]).upper()[0] == "D" or 
                  rmspace(lines[7]).upper()[0] == "C"):
                _position = np.array([list(map(float, (i[:-1].split())[:3])) \
                                for i in lines[8:8+self.total_atom]])
                _dyn = np.ones([self.total_atom, 3])
            else:
                lexit("Error, cannot detemine dyn poscar.")
            
            # convert frac <-> cart
            if (rmspace(lines[7]).upper()[0] == "D" or
                rmspace(lines[8]).upper()[0] == "D"):
                _ = frac2cart(self.lattice9, _position)
                _ = np.hstack((_atom_index, _position, _, _dyn))
                self.position = np.vstack((self.position, _))
            elif (rmspace(lines[7]).upper()[0] == "C" or
                  rmspace(lines[8]).upper()[0] == "C"):
                _ = cart2frac(self.lattice6, _position)
                _ = np.hstack((_atom_index, _, _position, _dyn))
                self.position = np.vstack((self.position, _))
            else:
                lexit("Error, cannot detemine cartesian/fractional poscar.")

        elif self.calculator == "qein":
            _qeinblock = qe.parse_qein_block(lines)

            # cell parameters
            _lattice = _qeinblock["CELL_PARAMETERS"]
            self.lattice9 = np.array([list(map(float,(i[:-1].split())[:])) \
                                      for i in _lattice[1:]])
            if "bohr" in _qeinblock["CELL_PARAMETERS"]:
                self.lattice9 *= uc_bohr2ang
            self.lattice6 = lattice_conversion(self.lattice9)

            # coordinates
            _position = _qeinblock["ATOMIC_POSITIONS"]
            self.total_atom = int(len(_position)-1)

            # atomic index
            _atoms = []
            _dyn = np.array([])
            for i in _position[1:]:
                _atoms.append(i.split()[0])
            _atom_index = np.array([symbol_map[i] for i in _atoms])
            _atom_index = _atom_index.reshape(len(_atom_index), 1)

            # get coordinates
            if len(_position[1].split()) == 4:
                _position = np.array([list(map(float,(i[:-1].split())[1:4])) \
                                      for i in _position[1:]])
                _dyn = np.ones([len(_atom_index), 3])
            elif len(_position[1].split()) == 7:
                _position = np.array([list(map(float,(i[:-1].split())[1:4])) \
                                      for i in _position[1:]])
                _dyn = np.array([list(map(int,(i[:-1].split())[4:])) \
                                      for i in _position[1:]])
            else:
                lexit("Error, atomic postion block is not n*4 or n*7.")

            # check fractional or cartesian
            if "crystal" in _qeinblock["ATOMIC_POSITIONS"][0]:
                _ = frac2cart(self.lattice9, _position)
                _ = np.hstack((_atom_index, _position, _, _dyn))
                self.position = np.vstack((self.position, _))
            elif "angstrom" in _qeinblock["ATOMIC_POSITIONS"][0]:
                _ = cart2frac(self.lattice6, _position)
                _ = np.hstack((_atom_index, _, _position, _dyn))
            elif "bohr" in _qeinblock["ATOMIC_POSITIONS"][0]:
                _ = cart2frac(self.lattice6, _position * uc_bohr2ang)
                _ = np.hstack((_atom_index, _, _position, _dyn))
                self.position = np.vstack((self.position, _))
        
        elif self.calculator == "qeout":
            _qeout = qe.parse_qeout(lines)
            self.total_atom = _qeout["nat"]
            self.lattice9 = _qeout["cell_paramater"][-1]
            self.lattice6 = functions.lattice_conversion(self.lattice9)
            #print(_qeout["atomic_index"])
            

    def sort_atoms(self):
        """
        DESCRIPTION:
        -----------
        Sort the atomic coordinates by atomic index.     
        UPDATE:
        -------
        Update the self.position.
        """
        self.position = self.position[np.argsort(self.position[:, 0]), :]
    

    def add_dyn(self,
                xmin=-1e8, xmax=1e8,
                ymin=-1e8, ymax=1e8,
                zmin=-1e8, zmax=1e8,
                fix="none",
                move="none",
                frac=True):
        """
        DESCRIPTION:
        -----------
        Set enclosed atoms to free and the rest to fixed. 

        UPDATE:
        -------
        Update the dynrange of coordinates.
        """
        _coord = self.position[1:, 1:4]
        if not frac:
            _coord = self.position[1:, 1:4]

        for i in range(len(_coord)):
            if xmax > _coord[i, 0] > xmin and \
               ymax > _coord[i, 1] > ymin and \
               zmax > _coord[i, 2] > zmin:
                #print(zmax, _coord[i, 2], zmin)
                self.position[1+i, 7] = 0.0
                self.position[1+i, 8] = 0.0
                self.position[1+i, 9] = 0.0
            else:
                self.position[1+i, 7] = 1.0
                self.position[1+i, 8] = 1.0
                self.position[1+i, 9] = 1.0
        # fix by elements
        if fix == "none" and fix in database.symbol_map.keys():
            _fix_atomic_idx = database.symbol_map[fix]
            for _ in range(self.total_atom+1):
                if _fix_atomic_idx == self.position[_, 0]:
                    self.position[_, 7:] = 0
        elif move != "none" and move in database.symbol_map.keys():
            _fix_atomic_idx = database.symbol_map[move]
            for _ in range(self.total_atom+1):
                if _fix_atomic_idx == self.position[_, 0]:
                    self.position[_, 7:] = 1
    

    def scale_lattice9(
            self,
            scale = np.array([
                [1, 1, 1],
                [1, 1, 1],
                [1, 1, 1]
            ])
    ):
        """
        """
        self.lattice9 = self.lattice9 * scale


    def write_struct(self, 
                     fpath:str="./", 
                     fname:str="pwscf.in", 
                     calculator:str="qein",
                     dyn:bool=False,
                     frac:bool=True,
                     write2file:bool=True):
        """
        DESCRIPTION:
        -----------
        Write 

        PARAMETERS:
        -----------
        lines:string list
        
        UPDATE:
        -------
        """

        if calculator == "none":
            calculator == self.calculator
        
        _coord = self.position[1:, 1:4]
        if not frac:
            _coord = self.position[1:, 4:7]
        _coord = convert_array2strlist(_coord, 
                                       fmt="{:20.10f}")
        
        # atomic index to atomic type and atomic number
        _atomic_index = [int(_) for _ in self.position[1:,0]]
        _atomic_type = list(set(_atomic_index))
        _atomic_type = np.sort(np.array(_atomic_type))
        _atomic_types = [atom_data[_][1] for _ in _atomic_index]
        _atomic_num = [_atomic_index.count(_) for _ in _atomic_type]
        
        _lines = []
        _qesystemblock = {}

        if calculator == "poscar":
            _lines.append("%s\n" % self.title)
            _lines.append("1.0\n")
            _lines += convert_high2writeablelist(
                convert_array2strlist(self.lattice9))

            # merge the atomic type
            _atomic_type = [atom_data[_][1] for _ in _atomic_type]
            _lines.append("   ".join(_atomic_type)+"\n")
            _lines += convert_high2writeablelist(
                convert_array2strlist(_atomic_num, fmt="{:<6d}"))
            
            # write coordinates
            if dyn:
                _lines.append("Selective dynamics\n")
                _dyn = self.position[1:, 7:]
                _dyn = np.where(_dyn>0.5, "T", "F")
                _dyn = convert_array2strlist(_dyn, fmt="{:>5s}")
                if frac:
                    _lines.append("Direct\n")
                else:
                    _lines.append("Cartesian\n")
                #_lines += convert_high2writeablelist
                _lines += convert_high2writeablelist(
                    hstack_append_list(_coord, _dyn))
            else:
                if frac:
                    _lines.append("Direct\n")
                else:
                    _lines.append("Cartesian\n")
                _lines += convert_high2writeablelist(_coord)

        elif calculator == "qein":
            #
            _qesystemblock = {
                "qeblock": "system",
                "ibrav": 0,
                "nat": self.total_atom,
                "ntyp": len(_atomic_type)
            }

            if write2file:
                _lines = qe.qeblock2lines(_lines, _qesystemblock)
            else:
                pass

            _atomic_types = np.array(_atomic_types).reshape(self.total_atom, -1)
            _atomic_types = convert_array2strlist(_atomic_types, fmt="{:>6s}")

            if frac:
                _lines.append("ATOMIC_POSITIONS crystal\n")
            else:
                _lines.append("ATOMIC_POSITIONS angstrom\n")

            if dyn:
                _dyn = self.position[1:, 7:]
                _dyn = np.where(_dyn>0.5, "1.0", "0.0")
                _dyn = convert_array2strlist(_dyn, fmt="{:>5s}")
                _atomic_types = hstack_append_list(_atomic_types, _coord)
                _lines += convert_high2writeablelist(
                    hstack_append_list(_atomic_types, _dyn)
                    )               
            else:
                _lines += convert_high2writeablelist(
                    hstack_append_list(_atomic_types, _coord)
                    )
            
            # write CELL PARAMETERS
            _lines.appedope_sitend("CELL_PARAMETERS angstrom\n")
            _lines += convert_high2writeablelist(
                convert_array2strlist(self.lattice9)
            )
        
        if write2file:
            with open(fpath+"/"+fname, "w", newline="\n") as f:
                f.writelines(_lines)
        else:
            return _qesystemblock, _lines


    def get_rec_basis(self, with_2pi = False):
        """
        Return the reciprocal basis
        """
        if with_2pi:
            return np.linalg.inv(self.lattice9).T
        else:
            return np.linalg.inv(self.lattice9).T * 2 * np.pi


    def supercell(self, super_matrix=[2, 2, 1]):
        """
        DESCRIPTION:
        -----------
        Write 

        PARAMETERS:
        -----------
        super_matrix: list
            supercell matrix
        
        UPDATE:
        -------
        """
        _newlattice = self.lattice6[:3]*np.array(super_matrix)
        self.lattice6 = np.hstack((_newlattice, self.lattice6[3:]))


    def add_atom(self, atomic_symbo, position, is_frac=True):
        """

        """
        _frac_pos, _cart_pos = 0, 0
        if is_frac:
            _frac_pos = position
            _cart_pos = functions.frac2cart(self.lattice9, position)
        else:
            _frac_pos = functions.cart2frac_single(self.lattice9, position)
            _cart_pos = position

        self.total_atom += 1
        _new_line = np.array([
            functions.symbol_map[atomic_symbo],
             _frac_pos[0],
             _frac_pos[1],
             _frac_pos[2],
             _cart_pos[0],
             _cart_pos[1],
             _cart_pos[2],
             1,1,1
        ])
        self.position = np.vstack((
            self.position, _new_line
        ))
        self.sort_atoms()
    
    def replace_atom(
            self, 
            idx=0, 
            ranges="xmin=-1,xmax=-1,ymin=-1,ymax=-1,zmin=-1,zmax=-1,r=0.1", 
            atomic_symbo="Ti", 
            is_frac=True
        ):
        """
        Docstring for replace_atom
        
        :param self: Description
        :param idx: Description
        :param ranges: Description
        :param atomic_symbo: Description
        :param is_frac: Description
        """
        replace_range = {
            "xmin": -1e8, "ymin": -1e8, "zmin": -1e8, 
            "xmax": 1e8,  "ymax": 1e8,  "zmax": 1e8,
        }
            
        
        if idx != 0:
            if is_frac:
                pass
            else:
                pass

        else:
            _= functions.parser_inputpara(ranges)
            for key in _.keys():
                replace_range[key] = float(_[key])
        
            if is_frac:
                for _ in range(self.total_atom):
                    if self.position[_+1, 1] > replace_range["xmin"] and \
                       self.position[_+1, 1] < replace_range["xmax"] and \
                       self.position[_+1, 2] > replace_range["ymin"] and \
                       self.position[_+1, 2] < replace_range["ymax"] and \
                       self.position[_+1, 3] > replace_range["zmin"] and \
                       self.position[_+1, 3] < replace_range["zmax"]:
                        self.position[_+1, 0] = functions.symbol_map[atomic_symbo]
            else:
                pass
        
        self.sort_atoms()
                    
    

    def supercell(self, super_matrix=[2, 2, 1]):
        """
        DESCRIPTION:
        -----------
        Write 

        PARAMETERS:
        -----------
        super_matrix: list
            supercell matrix
        
        UPDATE:
        -------
        """

        if np.all(self.position[0] == 0):
            positions = self.position[1:]
        else:
            positions = self.position


        _expanded_x_positions = []
        for row in self.position[1:]:
            for i in range(super_matrix[0]):
                new_row = row.copy()
                new_row[4] += i * self.lattice6[0]
                new_row[1] = (new_row[1] + i) / super_matrix[0]
                _expanded_x_positions.append(new_row)

        # 在y方向扩展
        _expanded_positions = []
        for row in _expanded_x_positions:
            for j in range(super_matrix[1]):
                new_row = row.copy()
                new_row[5] += j * self.lattice6[1]
                new_row[2] = (new_row[2] + j) / super_matrix[1]
                _expanded_positions.append(new_row)

        # z 


        # 添加一行全 0 的行
        _expanded_positions.insert(0, np.zeros(10))

        self.position = np.array(_expanded_positions, dtype=object)

        _newlattice = self.lattice6[:3]*np.array(super_matrix)
        self.lattice6 = np.hstack((_newlattice, self.lattice6[3:]))
        self.lattice9 = functions.lattice_conversion(self.lattice6)

        self.sort_atoms()


class volumetric(struct):
    """
    DESCRIPTION:
    ------------
    Volumetric data object.

    PARAMETERS:
    -----------
    inp, string
        The name of the input structure
    
    """

    def __init__(self, inp):
        super().__init__(inp)
        self.grid = np.array([0, 0, 0])
        self.voldata = np.array([])
        self.__parse_volumetric()

    
    def __parse_volumetric(self):
        """
        """
        if self.calculator == "poscar":
            self.grid = np.loadtxt(
                self.inp,
                skiprows=self.total_atom+9,
                max_rows=1
            )
            _total_grid = int(self.grid[0]*self.grid[1]*self.grid[2])
            self.voldata = np.loadtxt(
                self.inp,
                skiprows=self.total_atom+10,
                max_rows=_total_grid//5
            ).flatten()
            if _total_grid % 5 != 0:
                self.voldata = np.hstack((
                    self.voldata,
                    np.loadtxt(
                        self.inp,
                        skiprows=self.total_atom+9+_total_grid//5,
                        max_rows=1
                    )
                ))
            self.voldata = self.voldata.reshape((
                int(self.grid[2]),
                int(self.grid[1]),
                int(self.grid[0])
            ))
        elif self.calculator == "qeout":
            pass
        else:
            pass
    
    def write_volumetric(
        self, 
        fpath = "./", 
        fname = "pwscf.in", 
        calculator = "qein", 
        dyn = False, 
        frac = True
        ):
        _, _lines = super().write_struct(
            fpath, fname, calculator, dyn, frac, write2file=False
        )
        _total_grid = int(self.grid[0]*self.grid[1]*self.grid[2])

        if calculator == "poscar":
            _lines.append("\n")
            with open(fpath+"/"+fname, "w", newline="\n") as f:
                f.writelines(_lines)
                np.savetxt(
                    f,
                    np.array([self.grid]),
                    fmt="%5d"
                )
                _voldata = self.voldata.flatten()[:_total_grid//5*5]
                np.savetxt(
                    f,
                    _voldata.reshape((-1, 5)),
                    fmt="%10.11e"
                )
                if _total_grid%5 != 0:
                    np.savetxt(
                        f,
                        self.voldata.flatten()[_total_grid//5*5:],
                        fmt="%10.11e"
                    )
        elif "qe" in  calculator:
            pass
    

    def set_value(
        self,
        value,
        xmin=0,
        xmax=1,
        ymin=0,
        ymax=1,
        zmin=0,
        zmax=1,
        frac=True
    ):
        """
        

        """
        if frac:
            xmin = int(xmin*self.grid[0])
            xmax = int(xmax*self.grid[0])
            ymin = int(ymin*self.grid[1])
            ymax = int(ymax*self.grid[1])
            zmin = int(zmin*self.grid[2])
            zmax = int(zmax*self.grid[2])
            self.voldata[zmin:zmax, ymin:ymax, xmin:xmax] = value
