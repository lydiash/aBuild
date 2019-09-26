from os import path
from aBuild import msg
import os
import ase

def convert_direct(lattice,vec):
    from numpy import sum as nsum, array,dot
    from numpy.linalg import inv
    from numpy import transpose, array, equal
    inv_lattice = inv(lattice.transpose())
    d_space_vector = dot(inv_lattice, array(vec))

    return d_space_vector

def convert_cartesian(lattice,directVec):
    from numpy import sum as nsum, array
    return nsum([directVec[i]*lattice[i] for i in [0,1,2]], axis=0)


def vec_in_list(vec,veclist):
#    print("Considering if vector {} is in list {}".format(vec,veclist))
#    print([abs(vec - x) for x in veclist], 'diffs')
#    print([ all(abs(vec - x) < 1e-4) for x in veclist], 'check')
    if any([ all(abs(vec - x) < 1e-4) for x in veclist]):
        return True
    else:
        return False


def map_into_cell(vec):
    from math import floor
    from numpy import array
    new_point = []
    for i in vec:
        if i < 0.0 or i > 1.0:
         #   print(i,' i')
         #   print(floor(i),' floor')
         #   print(i - floor(i),' result')
            new_point.append(i - floor(i))
        elif i == 1.0:
            new_point.append(0.0)
        else:
            new_point.append(i)
    return array(new_point)

def orthogonality_defect(lattice):
    from numpy.linalg import norm,det
    from numpy import prod
    return prod([norm(x) for x in lattice])/abs(det(lattice))


def _chop(epsilon, const, i, j):
    """Sets the value of i[j] to exactly 'const' if its value already lies
    within 'epsilon' of 'const'."""
#    print('chopping')
#    print(i[j], 'to be chopped')
#    print(const, 'to go to')
#    print(epsilon, 'epsilon')
    if abs(const-i[j]) <= epsilon:
        i[j] = float(const)

def _chop_all(epsilon, i):
    """Performs chop() on the expected values of +- 0, 0.5, 1 for each
    element of i."""
    for j in range(len(i)):
        _chop(epsilon, 3, i, j)
        _chop(epsilon, -3, i, j)
        _chop(epsilon, 2, i, j)
        _chop(epsilon, -2, i, j)
        _chop(epsilon, 1, i, j)
        _chop(epsilon, -1, i, j)
        _chop(epsilon, 0, i, j)
        _chop(epsilon, 0.5, i, j)
        _chop(epsilon, -0.5, i, j)

    return i





class Lattice:


    def __init__(self,specs,symops=None, clusters=None):

        # If the input is a dictionary, I'm probably passing all of the needed
        # information (lv, bv, etc) in explicitly
<<<<<<< HEAD
        print(specs,' specs')
=======
>>>>>>> 03893ba98eddd6991de841e54e88613fa8b4165d
        if isinstance(specs["lattice"], list):
            self._init_dict(specs)
        # If the input is a string, I'm just specifying a canonical lattice and I
        # expect the class to infer all the details
        elif isinstance(specs["lattice"], str):
            self._init_string(specs["lattice"])

<<<<<<< HEAD
        

=======

            
>>>>>>> 03893ba98eddd6991de841e54e88613fa8b4165d
    def _init_dict(self,enumdict):
        import numpy as np
        necessaryItems = ['lattice','basis','coordsys','name']

        if not all([x in enumdict for x in necessaryItems]):
            msg.fatal('Missing information when initializing Lattice object')

        if len(enumdict["lattice"]) == 3 and len(enumdict["lattice"][0]) == 3 and np.linalg.det(enumdict["lattice"]) != 0:
            self.lattice = enumdict["lattice"]
            self.lattice_name = enumdict["name"]
        else:
            msg.fatal("The lattice vectors must be a 3x3 matrix with the vectors as rows "
                        "and the vectors must be linearly independent.")
        self.basis = enumdict["basis"]
        self.coordsys = enumdict["coordsys"]
        self.lattice_name = enumdict["name"]
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
        from aBuild.utility import _chop_all
        if self.coordsys[0].upper() == 'D':
            return self.Bv

        from numpy.linalg import inv
        from numpy import transpose, array, equal
        inv_lattice = inv(self.Lv.transpose()*self.latpar)        
        d_space_vector = [ list(dot(inv_lattice, array(b))) for b in self.Bv ]

        output = []
        for i in d_space_vector:
            if i not in output:
                output.append(_chop_all(epsilon, i)) 

        return output


    
    def Lv_strings(self,nolatpar = False):
        if nolatpar:
            return [' '.join(list(map(str,x))) for x in self.Lv * self.latpar]

        else:
            return [' '.join(list(map(str,x))) for x in self.Lv]
            
        
class Crystal(object):

    def __init__(self,crystalSpecs, systemSpecies,crystalSpecies = None,lFormat = 'mtpselect'):
        from numpy import array
        self.species = systemSpecies
        if isinstance(crystalSpecs,dict):  # When the crystal info is given in a dictionary
            self._init_dict(crystalSpecs)
        elif isinstance(crystalSpecs,str): # Read a file with (probably poscar)
            print('initializing from a file')
            self._init_file(crystalSpecs)
        elif isinstance(crystalSpecs,list): # Like when you read a "new_training.cfg" or a "structures.in"
            print('initializing from a list')
            self._init_lines(crystalSpecs, lFormat)
<<<<<<< HEAD
            
#        self.results = None
#        typesList = [[idx] * aCount for idx,aCount in enumerate(self.atom_counts)]
#        self.atom_types = []
#        for i in typesList:
#            self.atom_types += i
#        self.nAtoms = sum(self.atom_counts)
#        self.nTypes = len(self.atom_counts)
=======
        self.results = None
        typesList = [[idx] * aCount for idx,aCount in enumerate(self.atom_counts)]
        self.atom_types = []
        for i in typesList:
            self.atom_types += i
        self.nAtoms = sum(self.atom_counts)
        self.nTypes = len(self.atom_counts)
        self.basis += array([0,0,0])
        print(self.basis)
        print(self.nAtoms, "nAtoms", len(self.basis), "len of basis")
>>>>>>> 03893ba98eddd6991de841e54e88613fa8b4165d
        if self.nAtoms != len(self.basis):
            print("This is where we are.")
            msg.fatal("We have a problem")

            #        print("Before, the atom counts was {}.".format(self.atom_counts)) 
        self._add_zeros(systemSpecies,crystalSpecies)
        
        #print("After, the atom counts was {}.".format(self.atom_counts)) 
            #''' Didn't read in a POTCAR that was associated with the crystal '''
            
        self.AFMPlanes = None
        if len(self.species) != self.nTypes:
            msg.fatal('The number of atom types provided ({})is not consistent with the number of atom types found in the poscar ({})'.format(self.species,self.nTypes))
        if self.species == None:
            msg.fatal('I have to know what kind of atoms are in the crystal')
        if self.latpar is None and self.species is not None:
            self.set_latpar()
        #self.validateCrystal()
        #Just testing this function out.
#        print(self.atom_counts)
#        self.superPeriodics(3)

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
                msg.fatal("It appears that we need to insert zeros in atom_counts, but I don't know where to add them")
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


    def randomDisplace(self):
        from numpy.random import randn
        msg.warn("Warning:  You are about to move basis atoms around")

        
        for ibv,bv in enumerate(self.Bv_cartesian):
            disp = 0.05 * randn(3)
            print(disp, 'disp')
            self.basis[ibv] =  bv + disp
        self.coordsys = 'C'
        
    def _init_file(self,filepath):

        if 'poscar' in filepath.lower():
            self.from_poscar(filepath)
        elif 'input.in' in filepath.lower():
            self.from_lammpsin(filepath)
        else:
            msg.fatal("Not sure about the format of the file you want me to read")

            
    def _init_lines(self,lines,linesFormat):
        if linesFormat == 'mlp':
            self.fromMLP(lines)

    def _init_dict(self,crystalDict):
        necessary = ['lattice','basis','atom_counts','coordsys','species']

        if not all([x in crystalDict for x in necessary]):
            msg.fatal("Some necessary information not set upon initialization of Crystal object")

        from numpy import array

        self.lattice = crystalDict["lattice"]
        self.basis = crystalDict["basis"]
        self.atom_counts = crystalDict["atom_counts"]
        typesList = [[idx] * aCount for idx,aCount in enumerate(self.atom_counts)]
        self.atom_types = []
        for i in typesList:
            self.atom_types += i
        self.nAtoms = sum(self.atom_counts)
        self.nTypes = len(self.atom_counts)
 
        self.coordsys = crystalDict["coordsys"]
        self.species = crystalDict["species"]
        self.nTypes = crystalDict["nTypes"]
#########################################################
#        if not crystalDict["nAtoms"]:
#            self.nAtoms = sum(self.atom_counts)
#########################################################
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
#        if len(self.species) != self.nTypes:
#            msg.fatal('The number of atom types provided ({})is not consistent with the number of atom types found in the poscar ({})'.format(species,self.lattice.nTypes))
#        try:
#            self.strN = int(self.title.split()[-1])
#        except:
#            self.strN = 0
#        self.calcResults = None



    def __str__(self):
        """Returns the string representation of the POSCAR lines to write
        to a file."""
        return '\n'.join(self.lines())

    # Checks to see if any basis vectors in the crystal
    # are outside of the first unit cell.  If they are, we
    # map them back inside the first cell.
    def validateCrystal(self):
        from numpy import array,any,cross,dot,copy
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
<<<<<<< HEAD
           # msg.info('Crystal is fixed')
        else:
            pass
           # msg.info("Crystal didn't need fixing")

=======
####################################my changes################################
        if dot(cross(self.lattice[0],self.lattice[1]), self.lattice[2]) < 0:
            self.basis = self.Bv_cartesian
            self.coordsys = 'C'
            temp = copy(self.lattice[0])
            self.lattice[0] = copy(self.lattice[1])
            self.lattice[1] = temp
##############################################################################
            #msg.info('Crystal is fixed') ########I commented this out
        #else: ########I commented this out
            #msg.info("Crystal didn't need fixing") ########I commented this out
            
>>>>>>> 03893ba98eddd6991de841e54e88613fa8b4165d
    @property
    def orthogonality_defect(self):
        from numpy.linalg import norm,det
        from numpy import prod

        return prod([norm(x) for x in self.lattice])/abs(det(self.lattice))
<<<<<<< HEAD
=======
    
>>>>>>> 03893ba98eddd6991de841e54e88613fa8b4165d
    @property
    def cell_volume(self):
        from numpy import dot,cross
        return abs(dot(cross(self.lattice[0],self.lattice[1]),self.lattice[2]))

    # Calculates the distances between all of the atoms and finds
    # the minimum value from all of them.
    @property
    def minDist(self):
        from numpy import array,dot,min,einsum,add,roll,column_stack
        from numpy.linalg import norm
        from itertools import product
        import sys
        import numpy
        numpy.set_printoptions(threshold=sys.maxsize)


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
#        print(neighborsCartesian, 'nc')
#        print(rolledNeighbors,' rn')
#        print(rolledNeighbors.shape,' rn shape')
#        print(neighborsCartesian.shape,' nc shape')
#        import sys
#        sys.exit()
        #Now we can just subtract the first row (unshifted) from all of the other rows
        # and calculate the norm of each vector
        distances = norm(rolledNeighbors[0,:,:] - rolledNeighbors,axis = 2)
        # Return the min, excluding 0 distances.
        from numpy import count_nonzero,nonzero,transpose
        nZeroOccurrences = count_nonzero(distances==0)
        if nZeroOccurrences > self.nAtoms * len(offsets):
            ##import sys
            #import numpy
            #numpy.set_printoptions(threshold=sys.maxsize)
            #print(distances, 'distances')
            #print(transpose(nonzero(distances == 0)))
            msg.fatal("Atoms are on top of each other")
        return min(distances[distances > 1e-5])

#    @property
#    def appMinDist(self):
#        from aBuild.calculators import data
#        return data.nnDistance(self.species,self.atom_counts)
        
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

    def getAFMPlanes(self,direction):
        from numpy import sort,array,where,any
        from numpy.linalg import norm

        self.findNeighbors(0.75 * norm(self.latpar * self.lattice[0]+self.latpar * self.lattice[1]+self.latpar * self.lattice[2]) )
        # Find all of the A atoms
        aAtoms = where(array(self.atom_types) == 0)[0]
        # Get position vectors for all of the A atoms, 
        locationAatoms = array([self.neighbors[x][y] for x in aAtoms for y in range(len(self.neighbors[x])) ])
        # Build a dictionary, one list for every plane of atoms
        planesDict = {}
        for idx,atom in enumerate(locationAatoms):
            if atom[0][0] not in planesDict:
                planesDict[atom[0][0]] = [atom]
            else:
                planesDict[atom[0][0]].append(atom)

        # Check to see if every plane contains the same basis atom.
        foundPlanes = True
        for plane in list(planesDict.keys()):
            if any(abs(array([x[2] for x in planesDict[plane]]) - planesDict[plane][0][2]) > 1e-3):
                foundPlanes = False
                
        if foundPlanes:
            print([planesDict[x][0][2] for x in list(planesDict.keys())])
            return [planesDict[x][0][2] for x in list(planesDict.keys())]
        else:
            print("AFM doesn't work")
            return []
#        import sys
#        sys.exit()

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
        from aBuild.utility import _chop_all
        '''returns the basis in direct coordinates, whether it started as Direct or not'''
        from numpy import sum as nsum, array, dot
        if self.coordsys[0].upper() == 'D':
            return self.basis

        from numpy.linalg import inv
        from numpy import transpose, array, equal
        inv_lattice = inv(self.lattice.transpose()*self.latpar)        
        d_space_vector = [ list(dot(inv_lattice, array(b))) for b in self.basis ]

        output = []
        for i in d_space_vector:
            if i not in output:
                output.append(_chop_all(epsilon, i))

        return output

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


    def findNeighbors(self,rcut):
<<<<<<< HEAD
        from numpy import array,dot,min,einsum,add,roll,column_stack,append,where,extract,argwhere,set_printoptions,round
=======
        from numpy import array,dot,min,einsum,add,roll,column_stack,append,where,extract,argwhere,set_printoptions
>>>>>>> 03893ba98eddd6991de841e54e88613fa8b4165d
        from numpy.linalg import norm
        from itertools import product

        # Get all of the offsets to translate the basis vectors to equivalent
        # positions
        offsets = array([x for x in product(range(-1,2),repeat = 3)])

        # Translate them (in direct coordinates)  Shape: nBasis x nOffset x nCoord (3)
        neighborsDirect  = array([x + offsets for x in self.Bv_direct])
        # Convert coordinates to cartesian
        neighborsCartesian = einsum('abc,cd',neighborsDirect,self.latpar * self.lattice)
        diffs = array([x  - neighborsCartesian for x in self.Bv_cartesian])
        distances = norm(diffs, axis = 3)
        keepIndices = argwhere(distances < rcut)
<<<<<<< HEAD
        self.neighbors = [[] for x in self.atom_types]
        for [centerAtom,neighborAtom,shift] in keepIndices:
            self.neighbors[centerAtom].append([round(neighborsCartesian[neighborAtom,shift],3), self.atom_types[neighborAtom],neighborAtom])
=======

        self.neighbors = [[] for x in self.atom_types]
        for [centerAtom,neighborAtom,shift] in keepIndices:
            self.neighbors[centerAtom].append([neighborsCartesian[neighborAtom,shift], self.atom_types[neighborAtom] ]) 
#        print(array(self.neighbors))
>>>>>>> 03893ba98eddd6991de841e54e88613fa8b4165d



    # Warning: This routing is NOT general.  It is designed to
    # determine whether alternating (100) planes is possible for the given
    # crystal structure.  Need to generalize to any direction and any
    # alternating pattern.
    def getAFMPlanes(self,direction):
        from numpy import sort,array,where,any
        from numpy.linalg import norm

        if self.atom_counts[0] < 2:
            return []
        self.findNeighbors(0.75 * norm(self.latpar * self.lattice[0]+self.latpar * self.lattice[1]+self.latpar * self.lattice[2]) )
        # Find all of the A atoms
        aAtoms = where(array(self.atom_types) == 0)[0]
        # Get position vectors for all of the A atoms,
        locationAatoms = [basis[idy] for idx,basis in enumerate(self.neighbors) for idy,y in enumerate(basis) if basis[idy][2] in aAtoms]

        eps = 1e-1
        # Build a dictionary, one list for every plane of atoms
        planesDict = {}

        # Sometimes after a relaxation, atoms will have moved slightly off of their ideal
        # planes.  So we need to check to see if atoms are close to their planes.
        for idx,atom in enumerate(locationAatoms):
            if not any(abs(array([x for x in planesDict.keys()]) - atom[0][0]) < eps): #Found new plane
#                print('new plane found')
#                print(atom[0][0], 'x  coord')
#                print(atom, 'full atom to record')
#            if atom[0][0] not in planesDict:
                planesDict[atom[0][0]] = [atom]
            else:  #Adding an atom to a plane that I already found
 #               print('adding atom to previously found plane.')
                loc = abs(array([x for x in planesDict.keys()]) - atom[0][0]).argmin()
                planeAdd = array([x for x in planesDict.keys()])[loc]
 #               print(list(planesDict.keys()), "planes found thus far")
 #               print(array([x for x in planesDict.keys()]) - atom[0][0], 'diffs')
 ##               print(abs(array([x for x in planesDict.keys()]) - atom[0][0]).argmin(), 'location of min')
 #               print("Placing atom {} in plane {}.".format(atom,planeAdd))
                planesDict[planeAdd].append(atom)
        planes = sorted(list(planesDict.keys()))
  #      print(planes, 'here are the planes')
  #      print(planesDict,' planesDict')
        nPlanes = len(planes)
        planesAtoms = [set([x[2] for x in planesDict[y]]) for y in planes  ]
        evenPlanes = set([x for y in range(0,nPlanes,2) for x  in planesAtoms[y] ])
        oddPlanes = set([x for y in range(1,nPlanes,2) for x  in planesAtoms[y] ])
#        print(evenPlanes, 'even planes')
#        print(oddPlanes, 'odd planes')
        check = [x in evenPlanes for x in oddPlanes]
        if any(check):
            self.AFMPlanes = None
           # print(planes, 'here are the planes')
           # print(evenPlanes, 'even planes')
           # print(oddPlanes, 'odd planes')
           # print(planesDict,' planes dict')
            return False
        else:
            self.AFMPlanes = [0 for x in aAtoms]
            for i in evenPlanes:
                self.AFMPlanes[i] = 2
            for i in oddPlanes:
                self.AFMPlanes[i] = -2
            print(evenPlanes,oddPlanes, self.AFMPlanes,' check AFM results') 
            return True
    @property
    def Bv_direct(self):
<<<<<<< HEAD
        from numpy import sum as nsum, array,dot
=======
        from numpy import sum as nsum, array, dot
        from aBuild.utility import _chop_all
>>>>>>> 03893ba98eddd6991de841e54e88613fa8b4165d
        if self.coordsys[0].upper() == 'D':
            return self.basis
        from numpy.linalg import inv
        from numpy import transpose, array, equal
        inv_lattice = inv(self.lattice.transpose()*self.latpar)        
        d_space_vector = [ list(dot(inv_lattice, array(b))) for b in self.basis ]

        epsilon = 1e-4
        output = []
        for i in d_space_vector:
            if _chop_all(epsilon, i) not in output:
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
        #print("lines",self.lines(fileformat))
        with open(filename, 'w') as f:
            f.write('\n'.join(self.lines(fileformat)))


    def concsOK(self,concRestrictions=None):
        from numpy import all,array
        if concRestrictions == None:
            return True
        lowerlimit = array([x[0]/x[2] for x in concRestrictions])
        upperlimit = array([x[1]/x[2] for x in concRestrictions])
        if all(self.concentrations >= lowerlimit) and all(self.concentrations <= upperlimit):
            return True
        return False
        #  This routine switches A atoms for B atoms and vice versa
        #  For prototype structures, we only provide one crystal structure
        # excluding the concentration brothers.  
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
        typesList = [[idx] * aCount for idx,aCount in enumerate(self.atom_counts)]
        self.atom_types = []
        for i in typesList:
            self.atom_types += i
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

        self.latpar = data.vegardsVolume(self.species,self.atom_counts,self.volume)

    def from_poscar(self,filepath):
        """Returns an initialized Lattice object using the contents of the
        POSCAR file at the specified filepath.

        :arg strN: an optional structure number. If the label in the POSCAR doesn't
          already include the strN, it will be added to the title.
        """
        from aBuild.calculators.vasp import POSCAR
        print("filepath is ", filepath)
        lines = POSCAR(filepath)
        from numpy import array
        #First get hold of the compulsory lattice information for the class
        try:
            self.lattice = array([list(map(float, l.strip().split()[0:3])) for l in lines.Lv])
            self.basis = array([list(map(float, b.strip().split()[0:3])) for b in lines.Bv])
            self.atom_counts = array(list(map(int, lines.atom_counts.split() )  ))
            typesList = [[idx] * aCount for idx,aCount in enumerate(self.atom_counts)]
            self.atom_types = []
            for i in typesList:
                self.atom_types += i
        
            self.latpar = float(lines.latpar.split()[0])
            if self.latpar == 1.0 or self.latpar < 0:
                self.latpar = None
            #print(self.latpar, 'lat par read in')
            self.coordsys = lines.coordsys
#            self.title  = lines.label
            if len(lines.label) > 50:
                self.title = path.split(filepath)[0].split('/')[-1] + '_proto'
            else:
                self.title = path.split(filepath)[0].split('/')[-1] + '_' + lines.label
            self.nAtoms = sum(self.atom_counts)
            self.nTypes = len(self.atom_counts)
 
        except:
            raise ValueError("Lv, Bv or atom_counts unparseable in {}".format(filepath))
            

    @property
    def reportline(self):
        line = [self.title , self.results["energyF"],self.results["energyF"]/self.nAtoms, self.results["fEnth"]]
        lineWrite = '{:55s} {:9.4f} {:12.4f}{:12.4f}'
        
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
        
    def fromMLP(self,lines):
        from numpy import array
        import os
        from aBuild.calculators.vasp import VASP
        nAtoms = int(lines[2].split()[0])
        latDict = {}
        self.results = {}
        self.lattice = array([list(map(float,x.split())) for x in lines[4:7]])
        self.basis = array([list(map(float,x.split()[2:5])) for x in lines[8:8 + nAtoms]])
        self.nAtoms = len(self.basis)
        self.coordsys = 'C'
        self.atom_types = [int(x.split()[1]) for x in lines[8:8 + nAtoms]]
        self.results["forces"] = array([list(map(float,x.split()[5:8])) for x in lines[8:8 + nAtoms]])
        self.atom_counts = array([ self.atom_types.count(x) for x in set(self.atom_types)]) #range(max(atoms)+1)
        self.nTypes = len(self.atom_counts)
        # THe add_zeros function is here in case you run into a config with few atom types than the
        # system being studied. Like running into a pure config for a binary system.  Or a binary config
        # when studying a ternary system.
        self._add_zeros(self.species,[self.species[x] for x in list(set(self.atom_types))])
        titleindex = ['conf_id' in x for x in lines].index(True)
        self.title = ' '.join(lines[titleindex].split()[2:])
        if any(['Energy' in x for x in lines] ):
            energyindex = ['Energy' in x for x in lines].index(True)
            
            self.results["energyF"] = float(lines[energyindex + 1].strip().split()[0])
#            root = os.getcwd()
#            pures = [VASP(path.join(root,'training_set','pure' + x),systemSpecies = self.species)   for x in self.species]
#            puresDict = {}
#            for ispec,spec in enumerate(self.species):
#                #thispure = VASP(path.join(root,'training_set','pure' + spec,systemSpecies = self.species))
#                pures[ispec].read_results()
##                thispure.read_results()
#                puresDict[spec] = pures[ispec].crystal.results["energypatom"]
#            self.results["fEnth"] = self.results["energyF"]/self.nAtoms - sum(   [ pures[i].crystal.results["energyF"]/pures[i].crystal.nAtoms * self.concentrations[i] for i in range(self.nTypes)])
        else:
            self.results = None
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

<<<<<<< HEAD
    def superPeriodics(self,size):
        from numpy import dot,array,prod,einsum,any,all,copy
        from numpy.linalg import inv,det,norm
        from numpy import sum as nsum,transpose,matmul,cross
        from itertools import product,combinations

        crystals = []
        # If there are 0 occurrences of the first atom type. Just skip it.
        if self.atom_counts[0] == 0:
            return []
        # Generate all super-periodic lattice vectors
        dimension = len(self.lattice)
        multiplesOf = array([x for x in combinations([x for x in product(range(-2,3),repeat = dimension) if sum(abs(array(x))) !=0],dimension) ]) 
        lattices = einsum('abc,dc->abd',multiplesOf,transpose(self.lattice))
        #Keep only the lattices with i) positive triple product  ii) specified size increase or less  iii) orthogonality defect
        # less than 2, and sort them by their orthogonality defect.
        keeps = sorted([x  for x in lattices if det(x) > 1e-3 and abs(det(x))/abs(det(self.lattice)) <= size and orthogonality_defect(x)/orthogonality_defect(self.lattice) < 2], key = lambda k: orthogonality_defect(k))

#        sorted(keeps, key = lambda k: orthogonality_defect(k))
        # Done getting super-periodic lattice vectors.
        
        # Get all atoms by adding multiples of the parent lattice vectors
        combs =  [x for x in product(range(-4,4),repeat = 3)]
        latticeVecCombinations = einsum('ab,bd->ad', combs,self.lattice*self.latpar)
        basisAtoms = [x + latticeVecCombinations for x in self.Bv_cartesian]
#        print(self.Bv_cartesian, 'cartesian bv')
#        print(latticeVecCombinations, 'translates')
##        print(basisAtoms, 'result')
 #       import sys
 #       sys.exit()
        count = 1
        for lattice in keeps:

            count += 1
            basDirect = [map_into_cell(_chop_all(1e-4,map_into_cell(convert_direct(lattice*self.latpar, x)))) for x in self.Bv_cartesian]
            crystalDict = {"lattice":lattice, "basis":basDirect, "coordsys":'D', "atom_counts":[x for x in self.atom_counts],"species":self.species,"latpar":self.latpar,"title":self.title + '_super'}
            newCrystal = Crystal(crystalDict,self.species)
            if newCrystal.orthogonality_defect/self.orthogonality_defect > 2:
                print('Cell is too skew, not considering it')
                continue
            newCrystal.validateCrystal() #Maps all the atoms into the first unit cell
            sizeIncrease = round(newCrystal.cell_volume/self.cell_volume)
            if newCrystal.cell_volume < 1e-3:
                print("Found a zero volume cell!! Continuing without considering it")
                continue
            if abs(int(sizeIncrease) - sizeIncrease) > 1e-3:
                print(sizeIncrease,int(sizeIncrease), "Cell volume increase doesn't appear to be an integer multiple of the parent volume")
                import sys
                sys.exit()

            # We don't need to populate the cell if it's size 1 because we know there
            # are the same number of basis atoms as the original
#            if sizeIncrease == 1:
#                continue
            # Now populate the cell with all of the basis atoms.
            print("condisdering lattice {}".format(lattice))
            for idx, origbasis in enumerate(basisAtoms):
                #print(idx,' adding atoms of this type')
                track = 1
                for equiv in origbasis:
                    # First check to see if we have enough atoms of this type already.
                    if track == sizeIncrease:
                        break

                    candidate = map_into_cell(_chop_all(1e-4,map_into_cell(convert_direct(lattice*self.latpar, equiv))))
#                    candidate = _chop_all(1e-4,convert_direct(lattice*self.latpar, equiv))
#                    candidate = _chop_all(1e-4,map_into_cell(_chop_all(1e-4,convert_direct(lattice*self.latpar,equiv))))
                    if not vec_in_list(candidate,newCrystal.basis) and (not any(candidate >= 1) and not any(candidate < 0) ) :
                        newCrystal.basis.append(candidate)
                        atomType = self.atom_types[idx]
                        newCrystal.atom_types.append(atomType)
                        newCrystal.atom_counts[atomType] += 1
                        newCrystal.nAtoms += 1
                        track += 1
                
            if abs(newCrystal.nAtoms - sizeIncrease * self.nAtoms) > 1e-5:
                print("You don't have enough atoms")
                print(self.nAtoms, 'n primitive atoms')
                print(newCrystal.nAtoms, 'n current atoms')
                print(newCrystal.lattice, 'new lattice')
                print(self.cell_volume, 'size primitive')
                print(newCrystal.cell_volume, 'size current')
                print(det(transpose(lattice)), 'det')
                print(sizeIncrease, 'size increase')
                print(array(newCrystal.basis), 'new basis')
                print(newCrystal.atom_counts,'atom counts')
                print(newCrystal.atom_types,'atom types')
                print(len(crystals))
                import sys
                sys.exit()
            sortKey = sorted(range(newCrystal.nAtoms), key = lambda x: newCrystal.atom_types[x])
            newCrystal.basis = [newCrystal.basis[x] for x in sortKey]
            newCrystal.atom_types = sorted(newCrystal.atom_types)
            
            print(newCrystal.minDist, 'mindist)')
            with open('supers.out','a+') as f:
                f.writelines("Title\n")
                f.writelines("{}\n".format(sizeIncrease))
                f.writelines(newCrystal.lattice_lines)
                f.writelines("\n{} {}\n".format(newCrystal.atom_counts[0],newCrystal.atom_counts[1]))
                f.writelines(newCrystal.basis_lines)
                f.writelines('\n\n\n')

            if abs(newCrystal.minDist - self.minDist)    > 1e-3:
                msg.fatal("mindist for super-periodic is different than parent")
                
            print("Got Super-periodic.  Checking for AFM compatible")
            if newCrystal.getAFMPlanes([1,0,0]):
                return newCrystal
            print('not AFM compatible')
            #crystals.append( newCrystal )
        #import sys
        #sys.exit()
        return []
=======
#################################my changes#################################################3
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
    def get_spin(self,plane,spinType,eps=1e-3):
        from numpy import dot, zeros, zeros_like, array
        from numpy.linalg import norm
        from aBuild.utility import _chop_all

        bases = self.Bv_cartesian
        indices = self.which_atoms_spin(spinType)
        
        if self.atom_counts[spinType] <= 1:
            print("One or less atom with spin in the unit cell.")
            spin = []
        else: #if there's more than one, find layers...
            LVs = [ dot(self.lattice[i]*self.latpar,plane) for i in [0,1,2] ] #value of LVs in direction of planes
            maxLV = max( LVs )
            minLV = min( LVs )
            if maxLV < 0 and minLV < 0:
                maxLV = 0
            if maxLV > 0 and minLV > 0:
                minLV = 0

            layers = [ dot( bases[i],plane ) for i in range( indices[0],indices[1] ) ]
            self.findNeighbors(0.75 * norm(self.latpar * self.lattice[0]+self.latpar * self.lattice[1]+self.latpar * self.lattice[2]) ) #neighbor[i][j] holds the following info: [<cartesian vector>, <atom_type>]
            for i in range(indices[0],indices[1]): #loop over all the neighbors for the atom type we care about
                layers += [ dot(neighbor[0], plane) for neighbor in self.neighbors[i] ]

>>>>>>> 03893ba98eddd6991de841e54e88613fa8b4165d
        
            #if there's values that are (very close to) the same
            new_layers = _chop_all(eps,layers)
#        for i in range( len(layers) ):
#            for value2 in layers:
#                diff = abs(layers[i] - value2)
#                if diff < eps: #if they're within some epsilon of eachother
#                    layers[i] = float(value2) 

            #now get rid of duplicates
            new_layer_values = []
#        new_layer_values = set(new_layers)
#        new_layer_values = [value for value in layers_inside if value not in new_layer_values]
            for value in new_layers:
                if value not in new_layer_values:
                    new_layer_values.append(value)        
#        for value in layers:
#            if value not in new_layer_values:
#                new_layer_values.append(value)

            new_layer_values.sort()
        
            layers_inside = [ x for x in new_layer_values if x >= minLV and x <= maxLV ] #layers "inside" the unit cell
            if len(layers_inside) % 2 == 1: 
                print( "There's an odd number of layers in the {} direction.".format(plane) )
                spin = []
            else: #Cool, there's an even number of layers. Now check other conditions
                new_layer_values.sort() #paranoia check
                layers_inside.sort() #just to make sure...

                #check if lattice vectors keep periodicity
                tally = 0
                for LV in LVs:
                    for i in range( len(layers_inside) ):
                        diff = abs(LV - layers_inside[i])
                        if diff < eps and i % 2 == 0: 
                            tally += 1 #if the LV fits in one of the even layers, it's good.
                            break  #Stop checking this LV
                if tally < 3: #not all three LVs worked
                    print( "This lattice doesn't keep the periodicity for AFM." )
                    spin = []

                elif tally == 3: #Cool, the lattices work. Now check the bases..
                                 #If an atom is spin up in the unit cell,
                                 #it has to be spin up in all the neighboring cells.
                                 #assign spin to all the atoms in the unit cell
                    unit_cell_spin = zeros(self.nAtoms)
                    for i in range(indices[0],indices[1]):
                        for j in range( len(new_layer_values) ):
                            basis = bases[i]
                            this_value = dot(basis,plane)
                            diff = abs(this_value - new_layer_values[j])
                            if diff < eps: #if the atom is in the layer we are investigating
                                if j % 2 == 0: #if it's in an even layer, give it up spin
                                    unit_cell_spin[i] = 2.0
                                else: #if it's in an odd layer, give it down spin
                                    unit_cell_spin[i] = -2.0
                    #assign spin to all the neighboring atoms
                    neighbor_spin = [ zeros(len(self.neighbors[x])) for x in range( indices[0],indices[1] ) ] #list of lists
                    for i in range(indices[0],indices[1]):
                        for j in range( len(self.neighbors[i]) ):
                            for k in range( len(new_layer_values) ):
                                basis = bases[i]
                                this_value = dot(basis, plane)
                                diff = abs(this_value - new_layer_values[k])
                                if diff < eps:
                                    if k % 2 == 0:
                                        neighbor_spin[i][j] = 2.0
                                    else:
                                        neighbor_spin[i][j] = -2.0
                    #now check that for each atom, the same atom in the neighboring cell has the same spin
                    success = True
                    for i in range(indices[0],indices[1]):
                        this_success = True
                        for j in range( len(self.neighbors[i]) ):
                            if unit_cell_spin[i] != neighbor_spin[i][j]: #if even one neighbor fails this test
                                this_success = False #the whole structure fails
                        if not this_success: #if even one basis fails, the whole structure fails
                            success = False
                    if not success:
                        print( "This basis set doesn't keep the periodicity for AFM." )
                        spin = []
                    else:
                        print( "Success! Assigning spin..." )
                        spin = unit_cell_spin
        self.spins = spin
       
    #checks if a structure can be AFM
    def checkAFM(self,direction,spinType,eps=1e-3):
        from numpy import array
        self.get_spin(direction,spinType,eps)
        if self.spins != []:
            return True
        else:
            return False

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
    
################################################################################################# 
    def superPeriodics(self,size,special_settings=None): ############Added special_settings
        from numpy import dot,array,prod,einsum,any,all,copy
        from numpy.linalg import inv,det,norm
        from numpy import sum as nsum,transpose,matmul,cross
        from itertools import product,combinations
        from aBuild.utility import _chop_all,map_into_cell, convert_direct,vec_in_list
        crystals = []

        # Generate all super-periodic lattice vectors
        dimension = len(self.lattice)
        multiplesOf = array([x for x in combinations([x for x in product(range(-1,2),repeat = dimension) if sum(abs(array(x))) !=0],dimension) ]) 

        lattices = einsum('abc,dc->abd',multiplesOf,transpose(self.lattice))
        keeps = [x  for x in lattices if abs(det(transpose(x))) > 1e-3 and abs(det(transpose(x)))/abs(det(transpose(self.lattice))) <= size]
        # Done getting super-periodic lattice vectors.
        
        # Get all atoms by adding multiples of the parent lattice vectors
        combs =  [x for x in product(range(-2,2),repeat = 3)]
        latticeVecCombinations = einsum('ab,bd->ad', combs,self.lattice*self.latpar)
        basisAtoms = [x + latticeVecCombinations for x in self.Bv_cartesian]
        for lattice in keeps:
            basDirect = [_chop_all(1e-4,convert_direct(lattice*self.latpar, x)) for x in self.Bv_cartesian]
            crystalDict = {"lattice":lattice, "basis":basDirect, "coordsys":'D', "atom_counts":[x for x in self.atom_counts],"species":self.species,"latpar":self.latpar,"nTypes":self.nTypes}
            crystalDict["title"] = self.title + ' super'
#            print("initializing crystal")
#            print(self.atom_counts,' atom counts')
#            print(self.lattice, 'orig lattice')
#            print(crystalDict, 'dict')
            newCrystal = Crystal(crystalDict,self.species)
            if newCrystal.orthogonality_defect/self.orthogonality_defect > 2:
                print('Cell is too skew, not considering it')
                continue
#            print("intialized")
            newCrystal.validateCrystal() #Maps all the atoms into the first unit cell
            sizeIncrease = round(newCrystal.cell_volume/self.cell_volume)
#            print(sizeIncrease, "size Inc")
            if newCrystal.cell_volume < 1e-3:
                print("Found a zero volume cell!! Continuing without considering it")
                continue
            if abs(int(sizeIncrease) - sizeIncrease) > 1e-3:
                print(sizeIncrease,int(sizeIncrease), "Cell volume increase doesn't appear to be an integer multiple of the parent volume")
                import sys
                sys.exit()

            # We don't need to populate the cell if it's size 1 because we know there
            # are the same number of basis atoms as the original
#            if sizeIncrease == 1:
#                continue
            # Now populate the cell with all of the basis atoms.
            for idx, origbasis in enumerate(basisAtoms):
                track = 1
                for equiv in origbasis:
                    # First check to see if we have enough atoms of this type already.
                    if track == sizeIncrease:
                        break
                    candidate = _chop_all(1e-4,map_into_cell(_chop_all(1e-4,convert_direct(newCrystal.lattice*newCrystal.latpar,equiv))))
                    if not vec_in_list(candidate,newCrystal.basis) and (not any(candidate >= 1) and not any(candidate < 0) ) :
                        newCrystal.basis.append(candidate)
                        atomType = self.atom_types[idx]
                        newCrystal.atom_types.append(atomType)
                        newCrystal.atom_counts[atomType] += 1
                        newCrystal.nAtoms += 1
                        track += 1
                        
            if abs(newCrystal.nAtoms - sizeIncrease * self.nAtoms) > 1e-5:
                print("You don't have enough atoms")
                print(self.nAtoms, 'n primitive atoms')
                print(newCrystal.nAtoms, 'n current atoms')
                print(newCrystal.lattice, 'new lattice')
                print(self.cell_volume, 'size primitive')
                print(newCrystal.cell_volume, 'size current')
                print(det(transpose(lattice)), 'det')
                print(sizeIncrease, 'size increase')
                print(array(newCrystal.basis), 'new basis')
                print(newCrystal.atom_counts,'atom counts')
                print(newCrystal.atom_types,'atom types')
                print(len(crystals))
                import sys
                sys.exit()
            sortKey = sorted(range(newCrystal.nAtoms), key = lambda x: newCrystal.atom_types[x])
            newCrystal.basis = [newCrystal.basis[x] for x in sortKey]
            newCrystal.atom_types = sorted(newCrystal.atom_types)
#            with open('supers.out','a+') as f:
#                f.writelines("Title\n")
#                f.writelines(newCrystal.lattice_lines)
#                f.writelines("\n{} {}\n".format(newCrystal.atom_counts[0],newCrystal.atom_counts[1]))
#                f.writelines(newCrystal.basis_lines)
#                f.writelines('\n\n\n')
#            print('adding new crystal')
#            if newCrystal.getAFMPlanes([1,0,0]) != []:
#                print('Found one')
#                return newCrystal ##############I added this
#                print(newCrystal.lattice, 'lattice')
#                print(newCrystal.basis,'basis')
#                print(sizeIncrease, 'size increase')
#                import sys
#                sys.exit()
#            crystals.append( newCrystal )
#        print("Found all super-periodics")
#        print(len(crystals))
#        import sys
#        sys.exit()
#        return crystals
#####################my changes#############################################
            AFM = newCrystal.checkAFM( special_settings["AFM"]["plane"] , special_settings["AFM"]["spin_type"])# , special_settings["eps"] )
            if AFM:
                print("Found a superperiodic crystal that works for AFM")
                return newCrystal
        print("Did not find any superperiodic crystals that can be AFM.")
        return None
##############################################################################I commented the following out  
       

