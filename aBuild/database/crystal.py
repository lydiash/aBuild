from os import path
from aBuild import msg
import os
import ase
class Lattice:


    def __init__(self,specs,symops=None, clusters=None):

        # If the input is a dictionary, I'm probably passing all of the needed
        # information (lv, bv, etc) in explicitly
        if isinstance(specs, dict): 
            self._init_dict(specs)
        # If the input is a string, I'm just specifying a canonical lattice and I
        # expect the class to infer all the details
        elif isinstance(specs, str):
            self._init_string(specs)

        

    def _init_dict(self,dict):

        necessaryItems = ['lattice','basis','coordsys','name']

        if not all([x in dict for x in necessaryItems]):
            msg.fatal('Missing information when initializing Lattice object')

        if len(dict["lattice"]) == 3 and len(dict["lattice"][0]) == 3 and np.linalg.det(dict["lattice"]) != 0:
            self.lattice = dict["lattice"]
            self.lattice_name = "custom"
        else:
            msg.fatal("The lattice vectors must be a 3x3 matrix with the vectors as rows "
                        "and the vectors must be linearly independent.")
        self.basis = dict["basis"]
        self.coordsys = dict["coordsys"]
        self.lattice_name = dict["name"]
        self.nBasis = len(self.basis)
        
    def _init_string(self,string):
        lVLookupDict = {'sc':[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]] ,'fcc': [[0.5,0.5,0.0],[0.5,0.0,0.5],[0.0,0.5,0.5]],'bcc': [[-0.5,0.5,0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5]],'hcp': [[1,0,0],[.5, 0.866025403784439, 0],[0, 0, 1.6329931618554521]]}
        bVLookupDict = {'sc': [[0.0,0.0,0.0]],'fcc': [[0.0,0.0,0.0]],'bcc': [[0.0,0.0,0.0]],'hcp': [[0,0,0],[0.5,0.28867513459,0.81649658093]]}

        latDict = {}
        self.lattice = lVLookupDict[string]
        self.basis = bVLookupDict[string]
        self.coordsys = 'D'
        self.lattice_name = string
        self.nBasis = len(self.basis)
        
        

    @property
    def Bv_cartesian(self):
        from numpy import sum as nsum, array
        if self.coordsys[0].upper() != 'C':
            return [nsum([B[i]*self.Lv[i]*self.latpar for i in [0,1,2]], axis=0) for B in self.Bv]
        else:
            return self.Bv


        
    @property
    def Bv_direct(self):
        from numpy import sum as nsum, array
        if self.coordsys[0].upper() == 'D':
            return self.Bv

        from numpy.linalg import inv
        from numpy import transpose, array, equal
        inv_lattice = inv(self.Lv.transpose())        
        d_space_vector = [ list(dot(inv_lattice, array(b))) for b in self.Bv ]
#############################################################################################################
#        output = []
#        for i in d_space_vector:
#            if i not in output:
#                output.append(_chop_all(epsilon, i)) #no idea where this _chop_all function is. can't use it

#        return output
        return d_space_vector #just return the d_space_vector instead
#############################################################################################################  
    
    def Lv_strings(self,nolatpar = False):
        if nolatpar:
            return [' '.join(list(map(str,x))) for x in self.Lv * self.latpar]

        else:
            return [' '.join(list(map(str,x))) for x in self.Lv]
            
        
class Crystal(object):

    def __init__(self,crystalSpecs, systemSpecies,crystalSpecies = None,lFormat = 'mtpselect'):
        self.species = systemSpecies
        if isinstance(crystalSpecs,dict):  # When the crystal info is given in a dictionary
            self._init_dict(crystalSpecs)
        elif isinstance(crystalSpecs,str): # Read a file with (probably poscar)
            print('initializing from a file')
            self._init_file(crystalSpecs)
        elif isinstance(crystalSpecs,list): # Like when you read a "new_training.cfg" or a "structures.in"
            print('initializing from a list')
            self._init_lines(crystalSpecs, lFormat)
        self.results = None
        self.nAtoms = sum(self.atom_counts)
        self.nTypes = len(self.atom_counts)
        
        if self.nAtoms != len(self.basis):
            msg.fatal("We have a problem")

            #        print("Before, the atom counts was {}.".format(self.atom_counts)) 
        self._add_zeros(systemSpecies,crystalSpecies)
        #print("After, the atom counts was {}.".format(self.atom_counts)) 
            #''' Didn't read in a POTCAR that was associated with the crystal '''
            

        if len(self.species) != self.nTypes:
            msg.fatal('The number of atom types provided ({})is not consistent with the number of atom types found in the poscar ({})'.format(self.species,self.nTypes))
        if self.species == None:
            msg.fatal('I have to know what kind of atoms are in the crystal')
        if self.latpar is None and self.species is not None:
            self.set_latpar()
        #self.validateCrystal()        


    #  Sometimes a crystal object will be instantiated and the number of atom types is not consistent with
    # the order of the system being studied. For example, maybe you are studying a ternary system, but some 
    # of the configurations in your training set happened to be binary configurations. (Happens more often at higher order)
    #  In these cases, we need to add zeros to the atom_counts variable to clearly indicate which species are 
    # present and which are not.  In other words, our standard for the Crystal object is that the number of species
    # will *always* be equal to the order of the system. No exceptions!
    def _add_zeros(self,systemSpecies,crystalSpecies):
        self.species = systemSpecies
        if crystalSpecies is None:
            from numpy import array
            # If you don't tell me what atoms are in
            # the crystal, then I'll just riffle the system species
            # in, starting at the front.  For example, when I read in the prototype
            # files, zeros are not included in the list of 
            if len(systemSpecies) != self.nTypes:
                diff = len(systemSpecies) - self.nTypes
                self.atom_counts = array(list(self.atom_counts) + [0 for x in range(diff)])
                self.nTypes = len(self.atom_counts)
                #                self.species = systemSpecies
                #            else:
                #self.species = systemSpecies
        else:
            if len(crystalSpecies) != self.nTypes:
                msg.fatal("The number of species that was read in (POTCAR) does not agree with atom_counts (POSCAR)")
            #  This is the case when a crystal *with the atomic species* are read in
            # but it just so happens that the number of species in this crystal does not match
            # the number of species in the system being studied.  In this case, we need to specify which
            # atomic system species are missing from this particular crystal.  We do this by augmenting zeros
            # to atom_counts at the appropriate location.  
            elif len(crystalSpecies) != len(systemSpecies):
                from numpy import insert
                lacking = list(set(systemSpecies) - set(crystalSpecies))
                indices = [systemSpecies.index(x) for x in lacking]
                for idx,ele in enumerate(indices):
                    self.atom_counts = insert(self.atom_counts,ele + idx,0)
                if len(lacking) > 1:
                    print(" I haven't tested this case, can you verify that it's working the way it should")
                    print("The system species is {}, and the crystal species is {} and our new atom counts is {}.".format( systemSpecies,crystalSpecies,self.atom_counts))
                    import sys
                    sys.exit()
                    # self.species = systemSpecies
                self.nTypes = len(self.atom_counts)
                # else:
                #self.species = systemSpecies
        
    def _init_file(self,filepath):

        if 'poscar' in filepath.lower():
            self.from_poscar(filepath)
        elif 'input.in' in filepath.lower():
            self.from_lammpsin(filepath)
        else:
            msg.fatal("Not sure about the format of the file you want me to read")

            
    def _init_lines(self,lines,linesFormat):
        if linesFormat == 'mlpselect':
            self.fromMLPSelect(lines)

    def _init_dict(self,crystalDict):
        necessary = ['lattice','basis','atom_counts','coordsys','species']

        if not all([x in crystalDict for x in necessary]):
            msg.fatal("Some necessary information not set upon initialization of Crystal object")

        from numpy import array

        self.lattice = crystalDict["lattice"]
        self.basis = crystalDict["basis"]
        self.atom_counts = crystalDict["atom_counts"]
        self.coordsys = crystalDict["coordsys"]
        self.species = crystalDict["species"]
        if sorted(self.species,reverse = True) != self.species:
            msg.fatal("The order of your atomic species is not in reverse alphabetical order... OK?")
        
        if 'title' in crystalDict:
            self.title = crystalDict["title"]
        else:
            self.title = None

        if 'latpar' in crystalDict:
            self.latpar = crystalDict["latpar"]
        else:
            self.latpar = None

                

            #        self.directory = directory
        if len(self.species) != self.nTypes:
            msg.fatal('The number of atom types provided ({})is not consistent with the number of atom types found in the poscar ({})'.format(species,self.lattice.nTypes))
        try:
            self.strN = int(self.title.split()[-1])
        except:
            self.strN = 0
        self.calcResults = None



    def __str__(self):
        """Returns the string representation of the POSCAR lines to write
        to a file."""
        return '\n'.join(self.lines())

    # Checks to see if any basis vectors in the crystal
    # are outside of the first unit cell.  If they are, we
    # map them back inside the first cell.
    def validateCrystal(self):
        from numpy import array,any
        from math import floor
       # print(self.Bv_direct, 'D')
        if any(array(self.Bv_direct) > 1) or any(array(self.Bv_direct) < 0):
            basis_inside = []
            for j in self.Bv_direct:
                new_point = []
                for i in j:
                    if i < 0.0 or i > 1.0:
                        new_point.append(i - floor(i))
                    elif i == 1.0:
                        new_point.append(0.0)
                    else:
                        new_point.append(i)
                basis_inside.append(new_point)
            self.basis = basis_inside
            self.coordsys = 'D'
            msg.info('Crystal is fixed')
        else:
            msg.info("Crystal didn't need fixing")

    # Calculates the distances between all of the atoms and finds
    # the minimum value from all of them.
    @property
    def minDist(self):
        from numpy import array,dot,min,einsum,add,roll,column_stack
        from numpy.linalg import norm
        from itertools import product


        # Need to make sure that all of the atoms are inside the first
        # unit cell before we compile list of distances
        self.validateCrystal()
        
        #  Calculate all possible shifts.  These are vectors of integers
        # representing the amount of each lattice vector that we are going
        # to add to each basis atom.  We only do combinations of (-1,0,1) because
        # that should be enough to find all possible distances between atoms.
        #  We're not trying to get all atoms out to some cutoff radius, we just want to
        # make sure we get enough atoms in there to find the min separation.
        offsets = array([x for x in product(range(-1,2),repeat = 3)])

        # Now shift every basis atom by every shift previously calculated
        neighborsDirect  = array([self.Bv_direct + array(x) for x in offsets])
        #Convert list of atomic positions to cartesian coordinates
        neighborsCartesian = einsum('abc,cd',neighborsDirect,self.latpar * self.lattice)
        # Flatten the list down to a single list of position vectors
        neighborsCartesian.resize(len(offsets)* self.nAtoms,3)

        # Build a matrix where each row is a shifted version of atomic positions.
        rolledNeighbors = array([roll(neighborsCartesian,x,axis = 0) for x in range(len(neighborsCartesian))])
        #Now we can just subtract the first row (unshifted) from all of the other rows
        # and calculate the norm of each vector
        distances = norm(rolledNeighbors[0,:,:] - rolledNeighbors,axis = 2)
        print(self.latpar, 'latpar')
        # Return the min, excluding 0 distances.
        return min(distances[distances > 1e-5])

    @property
    def appMinDist(self):
        from aBuild.calculators import data
        return data.nnDistance(self.species,self.atom_counts)
        
    @property
    def recip_Lv(self):
        if len(self.lattice) == 0:
            raise ValueError("Lattice vectors are required for finding reciprocal"
                             " lattice vectors.")

        from numpy import array,cross,dot
        groupings = [[self.lattice[1], self.lattice[2]], [self.lattice[2], self.lattice[0]], 
                     [self.lattice[0], self.lattice[1]]]
        crosses = [cross(v[0], v[1]) for v in groupings]
        dots = [dot(self.lattice[i], crosses[i]) for i in [0,1,2]]
        return [crosses[i]/dots[i] for i in [0,1,2]]
################################################################################
##############
 
    @property
    def Lv_cartesian(self):
        '''returns the magnitude of each lattice vector... should maybe be called Lv_polar'''
        from numpy.linalg import norm
        #calculate magnitudes of each lattice vector
        params = [norm(self.lattice[i])*self.latpar for i in [0,1,2]]
        return params

    @property
    def Lv_angles(self):
        '''return the angles between lattice vectors'''
        from numpy import dot, array, arccos, pi
        from numpy.linalg import norm
        #calculate angles between lattice vectors
        groupings = [ [ self.lattice[1], self.lattice[2] ],
                      [ self.lattice[2], self.lattice[0] ], 
                      [ self.lattice[0], self.lattice[1] ] ]
        dots = [dot(v[0], v[1]) for v in groupings]
        mags = [norm(v[0])*norm(v[1]) for v in groupings]
        angles = [arccos(dots[i]/mags[i])*180/pi for i in [0,1,2]]
        return angles

    #copied from the Lattice object
    @property
    def Bv_direct(self):
        '''returns the basis in direct coordinates, whether it started as Direct or not'''
        from numpy import sum as nsum, array, dot
        if self.coordsys[0].upper() == 'D':
            return self.basis

        from numpy.linalg import inv
        from numpy import transpose, array, equal
        inv_lattice = inv(self.lattice.transpose())        
        d_space_vector = [ list(dot(inv_lattice, array(b))) for b in self.basis ]

#        output = []
#        for i in d_space_vector:
#            if i not in output:
#                output.append(_chop_all(epsilon, i))

#        return output
        return d_space_vector
#################################################################################################

    @property
    def volume(self):

        from numpy import cross, dot
        return dot(cross(self.lattice[0],self.lattice[1]), self.lattice[2])

    @property
    def concentrations(self):
        return [self.atom_counts[x]/sum(self.atom_counts) for x in range(self.nTypes)]
    
    @property
    def Bv_cartesian(self):
        from numpy import sum as nsum, array
        if self.coordsys[0].upper() != 'C':
            return [nsum([B[i]*self.lattice[i]*self.latpar for i in [0,1,2]], axis=0) for B in self.basis]
        else:
            return self.basis


    def findNeighbors(self):
        from numpy import array,dot,min,einsum,add,roll,column_stack,append,where,extract,argwhere,set_printoptions
        from numpy.linalg import norm
        from itertools import product

        # Get all of the offsets to translate the basis vectors to equivalent
        # positions
        offsets = array([x for x in product(range(-4,4),repeat = 3)])

        # Translate them (in direct coordinates)  Shape: nBasis x nOffset x nCoord (3)
        neighborsDirect  = array([x + offsets for x in self.Bv_direct])
        # Convert coordinates to cartesian
        neighborsCartesian = einsum('abc,cd',neighborsDirect,self.latpar * self.lattice)
        diffs = array([x  - neighborsCartesian for x in self.Bv_cartesian])
        distances = norm(diffs, axis = 3)
        keepIndices = argwhere(distances < self.rcut)

        self.neighbors = [[] for x in self.atom_types]
        for [centerAtom,neighborAtom,shift] in keepIndices:
            self.neighbors[centerAtom].append([neighborsCartesian[neighborAtom,shift], self.atom_types[neighborAtom] ])
        print(array(self.neighbors))


    @property
    def Bv_direct(self):
        from numpy import sum as nsum, array
        if self.coordsys[0].upper() == 'D':
            return self.basis
        print('getting direct vectors')
        from numpy.linalg import inv
        from numpy import transpose, array, equal
        inv_lattice = inv(self.lattice.transpose())        
        d_space_vector = [ list(dot(inv_lattice, array(b))) for b in self.basis ]

        output = []
        for i in d_space_vector:
            if i not in output:
                output.append(_chop_all(epsilon, i))

        return output


    @property
    def lattice_lines_LAMMPS(self):
        """Return \n joined lattice vector text lines."""
        lines = ' a1 '
        lines += ' '.join(map(str,self.lattice[0]))
        lines += ' a2 '
        lines += ' '.join(map(str,self.lattice[1]))
        lines += ' a3 '
        lines += ' '.join(map(str,self.lattice[2]))

        return lines
#        return zip(['a1','a2','a3'],[' '.join(map(str,x)) for x in self.lattice]) #'\n'.join(list(self.lattice))

    @property
    def lattice_lines(self):
        """Return \n joined lattice vector text lines."""
        
        return '\n'.join( [' '.join(map(str,x)) for x in self.lattice]) #'\n'.join(list(self.lattice))

    @property
    def lattice_lines_nolatpar(self):
        """Return \n joined lattice vector text lines."""
        
        return '\n'.join( [' '.join(map(str,x)) for x in self.lattice * self.latpar]) #'\n'.join(list(self.lattice))

    @property
    def basis_lines_LAMMPS(self):
        """Return \n joined lattice vector text lines."""
        lines = ''
        for i in self.basis:
            lines += ' basis '
            lines += ' '.join(map(str,i))

        return lines

    @property
    def basis_lines_ESPRESSO(self):
        """Return \n joined lattice vector text lines."""
        speciesList = []
        for i in range(self.nTypes):
            speciesList += [self.species[i] for x in range(self.atom_counts[i])]
            
        lines = ''
        for idx,i in enumerate(self.basis):
            lines += speciesList[idx] + ' '
            lines += ' '.join(map(str,i))
            lines += '\n'

        return lines[:-1]

######################################################################################
    @property
    def basis_lines_cif(self):
        ''' Returns basis vectors formatted for .cif files '''
        speciesList = []
        for i in range(self.nTypes):
            speciesList += [self.species[i] for x in range(self.atom_counts[i])]
        numList = []
        #this could probably be a list comprehension but I'm not sure how to do it
        for i in range(self.nTypes):
            for j in range(self.atom_counts[i]):
                numList.append(j+1)
                
        #write the lines that should look like "A1 A 0.0 0.0 0.0 1.0"
        lines = ''
        for idx,i in enumerate(self.Bv_direct): #CIF is written in direct coordinates
            lines += speciesList[idx] 
            lines += str(numList[idx])  + ' '
            lines += speciesList[idx] + ' '
            lines += ' '.join(map(str,i))
            lines += ' 1.0 \n'

        return lines
#######################################################################################

    @property
    def unknown(self):
        """Return \n joined lattice vector text lines."""
        atomLabels = [ x for sublist in [ [  i for k in range(self.atom_counts[i]) ]   for i in range(self.nTypes)] for x in sublist]
        lines = ""
        for i,v in enumerate(atomLabels):
            lines += ' basis '
            lines += str(i+1) + ' ' + str(v+1)

        return lines
    
    @property
    def basis_lines(self):
        """Return \n joined basis vector text lines."""
        return '\n'.join( [' '.join(map(str,x)) for x in self.basis])


    def mtpLines(self,relax = False):
        import numpy as np
        if not relax and self.results is None:
            msg.fatal("You want me to write result information but I don't have any.")
        result = []
        result.append('BEGIN_CFG')
        result.append('Size')
        result.append(str(self.nAtoms))
        result.append('SuperCell')
        for lv in self.latpar * self.lattice:
            result.append('{:12.6f} {:12.6f} {:12.6f}'.format(lv[0],lv[1],lv[2]  ))
        if not relax:
            result.append('   AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz')
        else:
            result.append('   AtomData:  id type       cartes_x      cartes_y      cartes_z')
        #counter = crystal.lattice.nTypes - 1
        
        # Took me a few minutes to figure this one out.  Very Pythonic
        #        atomLabels = [ x for sublist in [ [ counter - i for k in range(crystal.lattice.atom_counts[i]) ]   for i in range(counter + 1)] for x in sublist]
        atomLabels = [ x for sublist in [ [  i for k in range(self.atom_counts[i]) ]   for i in range(self.nTypes)] for x in sublist]
        for i in range(self.nAtoms):
            if not relax:
                forces = self.results["forces"][i]
            coords = self.Bv_cartesian[i]
            if not relax:
                result.append('{:16d} {:3d} {:16.6f} {:12.6f} {:12.6f} {:18.6f} {:10.6f} {:10.6f}'.format(i+ 1,atomLabels[i], coords[0],coords[1],coords[2], forces[0],forces[1],forces[2]  ))
            else:
                result.append('{:16d} {:3d} {:16.6f} {:12.6f} {:12.6f}'.format(i+ 1,atomLabels[i], coords[0],coords[1],coords[2] ))
            #line +=  ' '.join([map(str,crystal.lattice.Bv_cartesian[i]), crystal.forces[i]])

        if not relax:
            result.append('Energy')
            result.append(str(self.results["energyF"]) + '')
        
        
            result.append(' Stress:   xx          yy           zz            yz           xz           xy')
            s = self.results["stress"]
            stressesline = '{:16.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}'.format(s[0],s[1],s[2],s[3],s[4],s[5])
            result.append(stressesline)
        result.append(''.join([' Feature   conf_id ', '  '.join([self.symbol,self.title]),'']))
        result.append('END_CFG\n')
        return result
    
    def vasplines(self):
        """Returns a list of strings for each line in the POSCAR file.

        :arg vasp: when true, the atom_counts line is checked for zeros before
          it is created. Vasp can't handle zero for the number of atoms of a
          certain type; just remove it."""
        result = []
        result.append(self.title)
        result.append(str(self.latpar))
        result.append(self.lattice_lines)
        result.append(' '.join(map(str,self.atom_counts)))
        #        if vasp:
        #    result.append(' '.join([a for a in self.atom_counts if a != '0' and a != ' ']))
        #else:
        #    result.append(' '.join([a for a in self.atom_counts if a != ' ']))
        result.append(self.coordsys)
        result.append(self.basis_lines)
        return result

    def lines(self,fileformat):
        if fileformat.lower() == 'vasp':
            return self.vasplines()
        elif fileformat.lower() == 'mtptrain':
            return self.mtpLines()
        elif fileformat.lower() == 'mtprelax':
            return self.mtpLines(relax = True)

    def write(self, filename, fileformat = 'vasp'):
        """Writes the contents of this POSCAR to the specified file."""
        #fullpath = os.path.abspath(filepath)
        with open(filename, 'w') as f:
            f.write('\n'.join(self.lines(fileformat)))


    def scrambleAtoms(self,scrambleKey):
        from numpy import array
        Bvs = []
        start = 0
        end = self.atom_counts[0]
        Bvs.append(self.basis[start:end])
        for index,aC in enumerate(self.atom_counts[:-1]):
            start = start + aC
            end = end + self.atom_counts[index + 1]
            Bvs.append(self.basis[start:end])

        self.basis = [y for sublist in [Bvs[x] for x in scrambleKey] for y in sublist]
        self.atom_counts = [self.atom_counts[x] for x in scrambleKey]
        self.set_latpar()
        
    @property
    def symbol(self):
        symbol = ''
        for elem, count in zip(self.species,self.atom_counts):
                symbol += elem + '_' + str(count)
#        if self.calcResults is not None and "species" in self.calcResults:
#            for elem, count in zip(self.calcResults["species"],self.atom_counts):
#                symbol += elem + '_' + str(count) + '-'
#        else:
#            for elem, count in zip(self.species,self.atom_counts):
#                symbol += elem + '_' + str(count)
#        symbol += '     ' + self.title  + '    '
#        symbol += self.filepath
        return symbol
        
    def set_latpar(self):
        from aBuild.calculators import data
        #if self.latpar == 1.0 or self.latpar < 0:
            # We must first reverse sorte the species list so we get the right atom in the right place.
        if sorted(self.species,reverse = True) != self.species:
            msg.fatal("Your species are not in reverse alphabetical order... OK?")

#        print(self.volume,' volume')
#        print(self.lattice,' lattice vecs')
        self.latpar = data.vegardsVolume(self.species,self.atom_counts,self.volume)
#        previously = data.vegard(self.species,[float(x)/self.nAtoms for x in self.atom_counts])
#        print("setting latpar to {}. Previously it was set to {}".format(self.latpar,previously) )
#        print('-------------------------')
#        import sys
#        sys.exit()
    def from_poscar(self,filepath):
        """Returns an initialized Lattice object using the contents of the
        POSCAR file at the specified filepath.

        :arg strN: an optional structure number. If the label in the POSCAR doesn't
          already include the strN, it will be added to the title.
        """
        from aBuild.calculators.vasp import POSCAR
        lines = POSCAR(filepath)
        from numpy import array
        #First get hold of the compulsory lattice information for the class
        try:
            self.lattice = array([list(map(float, l.strip().split()[0:3])) for l in lines.Lv])
            self.basis = array([list(map(float, b.strip().split()[0:3])) for b in lines.Bv])
            self.atom_counts = array(list(map(int, lines.atom_counts.split() )  ))
            self.latpar = float(lines.latpar.split()[0])
            if self.latpar == 1.0 or self.latpar < 0:
                self.latpar = None
            #print(self.latpar, 'lat par read in')
            self.coordsys = lines.coordsys
#            self.title  = lines.label
            self.title = path.split(filepath)[0].split('/')[-1] + '_' + lines.label

        except:
            raise ValueError("Lv, Bv or atom_counts unparseable in {}".format(filepath))
            

    @property
    def reportline(self):
        line = [self.title , self.results["energyF"], self.results["fEnth"]]
        lineWrite = '{:35s}  {:9.4f}   {:9.4f}'
        
        for i in self.concentrations:
            line.append(i)
            lineWrite += '  {:4.2f}  '
        for i in self.atom_counts:
            line.append(i)
            lineWrite += '  {:4d}  '
        lineWrite += '\n'
            #        for i in range(self.knary):
#            line.append(self.pures[i].results["energy"])
#            lineWrite += '{:8.5f}'
#        for i in range(self.knary):
#            line.append(self.atom_counts[i])
#            lineWrite += '{:2d}'


        return lineWrite.format(*line)


    def from_lammpsin(self,filepath):
        from numpy import array
        with open(filepath,'r') as f:
            lines = f.readlines()

        for idx,line in enumerate(lines):
            if 'lattice' in line:
                self.lattice = array([map(float,line.split()[4:7]),map(float,line.split()[8:11]), map(float,line.strip('&').split()[12:15])])
                nBasis = lines[idx + 1].count('basis')
                self.basis = [map(float,lines[idx+1].strip().strip('basis').split()[i * 4:4 * (i+ 1) - 1]) for i in range(nBasis)]
                self.latpar = float(line.split()[2])
            if 'create_box' in line:
                self.nTypes = int(line.split()[1])

            if 'create_atoms' in line:
                self.atom_counts = [0 for i in range(self.nTypes) ]
                for idx,elem in enumerate(line.split()[5::3]):
                    self.atom_counts[int(elem)-1] += 1
        self.coordsys = 'D'
        self.title = path.split(filepath)[0].split('/')[-1] + lines[1].split('\n')[0]
        
    def fromMLPSelect(self,lines):
        from numpy import array
        nAtoms = int(lines[2].split()[0])
        latDict = {}
        self.lattice = array([list(map(float,x.split())) for x in lines[4:7]])
        self.basis = array([list(map(float,x.split()[2:5])) for x in lines[8:8 + nAtoms]])
        self.nAtoms = len(self.basis)
        self.coordsys = 'C'
        atoms = [int(x.split()[1]) for x in lines[8:8 + nAtoms]]
        self.atom_counts = array([ atoms.count(x) for x in range(max(atoms)+1)])
        titleindex = ['conf_id' in x for x in lines].index(True)
        self.title = ' '.join(lines[titleindex].split()[2:])
        
        self.latpar = 1.0
        if sum(self.atom_counts) != nAtoms:
            msg.fatal('atomCounts didn\'t match up with total number of atoms')
#        self.set_latpar()
        self.latpar = 1.0  # MLP files are formatted with no lattice parameter.  It's
                           # already built into the lattice vectors.
        print(self.latpar, 'latpar')
#        self.lattice = self.lattice / self.latpar  # I think I did this to ensure that the lattice
                                                    # vectors didn't change but I know that the lattice
                                                    # parameter is just 1.0 for MLP formatting.
        
    @staticmethod  # Needs fixed!!!
    def fromEnum(enumDict,structNum):


        enumLattice.generatePOSCAR(struct)
        result = Crystal.fromPOSCAR(enumLattice.root, self.species,
                                         filename = "poscar.{}.{}".format(lat,struct),
                                         title = ' '.join([lat," str #: {}"]).format(struct))

        return result


#    def setCalcResults(self):
#
#        from aBuild.calculators.vasp import VASP
#        from aBuild.utility import chdir
#
#        thiscalc = VASP(directory = self.directory)
#        #        with chdir(self.directory):
#        thiscalc.read_results(allIonic = False,allElectronic= False)
#        self.calcResults = thiscalc.results
#        if self.calcResults is not None and len(self.calcResults["forces"]) != self.lattice.nAtoms:
#            msg.fatal('The number of forces ({}) is not equal to the number of atoms ({})'.format(self.calcResults["forces"],self.lattice.nAtoms) )


    def mlpLines(self):
        pass
        
    def setAtypes(self):
        pass
        
    #    def writePOSCAR(self,filepath):
    #    from aBuild.calculators.vasp import POSCAR

        # Instantiate a POSCAR object from my crystal object so that it intializes to have
        # everything that's inside of crystal
        #    lines = POSCAR(self)
        #lines.write(filepath,vasp=True)

#################################################################
    #generate a cif file from a POSCAR
    def generate_cif(self):
        import shutil
        from jinja2 import Environment, PackageLoader  # Package for building files from a template
        from os import path
        import aBuild
                                                                 
        medpath = path.abspath(aBuild.__file__)
        reporoot = path.dirname(path.dirname(medpath))
        #build the library to fill the template with
        settings = {}
        settings["lpara"] = self.Lv_cartesian[0]
        settings["lparb"] = self.Lv_cartesian[1]
        settings["lparc"] = self.Lv_cartesian[2]
        settings["alpha"] = self.Lv_angles[0]
        settings["beta"] = self.Lv_angles[1]
        settings["gamma"] = self.Lv_angles[2]
        settings["bVs"] = self.basis_lines_cif
        settings["title"] = self.title
        #find the template (template.cif)
        env = Environment(loader=PackageLoader('aBuild', 'templates'))
        template = env.get_template("template.cif")
        # make the label for the file
        fileName = str(settings["title"].split("_")[0])+".cif"
        #write the .cif file
        with open(fileName,'w') as f:
            f.write(template.render(**settings))


            
    #If I want to make a crystal AFM, I need to assign spin to the specified atoms so they
    #line up in a specified plane
    def get_spin(self,plane,spinType):
        from numpy import dot, zeros, array
        
        #find the layers of the magnetic substance
        layer_values = self.AFM_layers( plane , spinType )
        #find the indices that correspond to atoms with spin
        indices = self.which_atoms_spin(spinType)

        #if it's an even number of layers, great, it preserved the periodicity of the crystal
        #if not, return an empty list of spins
        if layer_values == 0:
            spins = [] #the lattice vectors didn't work for AFM. return an empty list of spins
        elif len(layer_values) % 2 == 1:
            print("Not an even number of planes in the {} direction.".format(plane) )
            spins = [] #return an empty list of spins. This means we can't use this as an AFM structure
        else:           
            #convert to cartesian
            bases = self.Bv_cartesian
            #keep track of the spins. 0 if it's not the species specified
            spins = zeros(self.nAtoms)
            #for each layer value, assign appropriate spin
            for i in range(len(layer_values)):
                #loop over the bases with spin
                for j in range(indices[0],indices[1]):
                    basis = bases[j]
                    this_value = dot(basis,plane)
                    #if the atom is in the layer we are investigating
                    diff = abs(this_value - layer_values[i])
                    if diff < 1e-3:
                        #give it some spin
                        if i % 2 == 0:
                            spins[j] = 2.0
                        else:
                            spins[j] = -2.0
        self.spins = spins

        
    # cartesian values of the layers for AFM
    def AFM_layers(self,direction,spinType):
        from numpy import dot

        #first we want to check if the lattice vectors are friendly to AFM structures.
        test = [ dot(self.lattice[i],direction) for i in [0,1,2] ]
        tally = 0
        for i in range( len(test) ):
            #if the lattice points in the direction we are searching, add one to the tally
            if abs(test[i]) > 1e-3:
                tally += 1
        #if more than one lattice vector points in the direction we're searching, it won't work for AFM
        if tally != 1:
            print("These lattice vectors aren't nice for AFM.")
            return 0 #return 0. Must check for this case

        #if the lattice vectors are friendly, now calculate the layers
        #convert to cartesian
        bases = self.Bv_cartesian
        layer_values = []

        indices = self.which_atoms_spin(spinType)
        
        #loop over all the bases (of the atom type you want)
        for i in range(indices[0],indices[1]):
            basis = bases[i]
            this_value = dot(basis,direction) #the value of the basis in the direction we are investigating
            #add this_value to the list of layers if it's not already in the list
            if this_value not in layer_values:
                layer_values.append(this_value)

        #if there's values that are (very close to) the same
        for value in layer_values:
            for i in range(len(layer_values)):
                if abs(value - layer_values[i]) < 1e-3:
                    layer_values[i] = float(value) #make them equal eachother

        #now get rid of duplicates
        new_layer_values = []
        for value in layer_values:
            if value not in new_layer_values:
                new_layer_values.append(value)

        #sort the layer_values in ascending order
        new_layer_values.sort()

        return new_layer_values

    def which_atoms_spin(self,spinType):
        from numpy import array
        #figure out which bases need to have spin assigned
        if spinType == 0:
            basisNumStart = 0
            basisNumEnd = self.atom_counts[0]
        if spinType == 1:
            basisNumStart = self.atom_counts[0]-1
            basisNumEnd = basisNumStart + self.atom_counts[1]
        if spinType == 2:
            basisNumStart = self.atom_counts[0]+self.atom_counts[1]-1
            basisNumEnd = basisNumStart + self.atom_counts[2]
        return array( [basisNumStart,basisNumEnd] )
    
################################################################# 

        
        

       
