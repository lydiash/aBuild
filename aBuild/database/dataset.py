
from aBuild import msg
from aBuild.utility import chdir, _get_reporoot

class dataset:

    def __init__(self,dset,systemSpecies,root=None,calculator = None,lFormat = 'mtpselect'):
        from os import path,makedirs
        from aBuild.database.crystal import Crystal
        from aBuild.calculators.vasp import VASP

        self.calculator = calculator
        
        if isinstance(dset,list):
            
            if isinstance(dset[0], dict):
                self.init_enum(dset,systemSpecies)
            elif isinstance(dset[0], Crystal):
                self.crystals = dset
                self.nCrystals = len(dset)
            elif isinstance(dset[0], VASP):
                self.calcs = dset
                self.nCalcs = len(dset)
            elif isinstance(dset[0],str):
                self.init_paths(dset,systemSpecies)
        elif isinstance(dset, str):
           self.init_file(dset,systemSpecies,lFormat)

        self.species = systemSpecies
        self.root = root

    # Used to be called 'buildFoldersFromEnum
    def init_enum(self,enumdicts,systemSpecies,runGetKpoints = True):
        from aBuild.enumeration import Enumerate
        from aBuild.calculators.vasp import VASP
        from aBuild.database.crystal import Crystal
        from aBuild.jobs import Job
        from random import randrange
        from aBuild.utility import chdir

        #        from crystal import Crystal
        from os import path
        import os

        #    if not path.isdir(self.root):
        #    os.mkdir(self.root)
        print("Building database from enumerations")
        self.crystals = []
        #        configIndex = startPoint = self._starting_point
        for eDict in enumdicts:
            enumController = Enumerate(eDict)
            if enumController.nEnumStructs == 0:
                msg.warn('There are no enumerated structures for lattice type {}.  Not building any VASP folders for them.'.format(self.enumDicts[index]["lattice"]))
                enumController.buildInputFile()

                enumController.enumerate()

            # Loop to generate random structures for a given lattice type
            for i in range(eDict["nconfigs"]):
                rStruct = randrange(1,enumController.nEnumStructs)
                print('Adding {} structure # {} to database'.format(eDict["lattice"],rStruct) )
                with open('structNums','a+') as f:
                    f.write(eDict["lattice"] + ' ' + str(rStruct) + '\n')
                    #print("Building VASP folder for {} structure #: {}".format(eDict["lattice"],rStruct))
                enumController.generatePOSCAR(rStruct)

                
                poscarpath = path.join(enumController.root,"poscar.{}.{}".format(eDict["lattice"],rStruct))
                thisCrystal = Crystal(poscarpath, systemSpecies = systemSpecies) #title = ' '.join([self.enumDicts[index]["lattice"]," str #: {}"]).format(rStruct)
                self.crystals.append(thisCrystal)

    # Sometimes an entire dataset is stored in one file.  I'd like to extract each crystal from the file to 
    # create a list of crystal objects
    def init_file(self,datafile,systemSpecies,linesformat):
        from aBuild.database.crystal import Crystal
        possibleFiles = {'new_training.cfg':'mlpadd','train.cfg': 'mlptrain','structures.in':'ce'}
        #selectedFile = path.join(self.root,'new_training.cfg')

        with open(datafile,'r') as f:
            lines = f.readlines()

        self.crystals = []
        nCrystals = 0
        for index,line in enumerate(lines):
            
            if line == 'BEGIN_CFG\n':
                #                nCrystals += 1
                #if numOfStructs is not 'all' and (nCrystals < start or nCrystals > start + numOfStructs):
                #    continue
                nAtoms = int(lines[index+2].split()[0])
                structlines = lines[index:index + 11 + nAtoms]
                thisCrystal = Crystal(structlines,systemSpecies,lFormat = linesformat)
                self.crystals.append(thisCrystal)



    def init_paths(self,paths,systemSpecies):
        from aBuild.database.crystal import Crystal
        from aBuild.calculators.vasp import VASP
        from aBuild.calculators.lammps import LAMMPS
        from os import path
        
        self.crystals = []
        for dirpath in paths:
            print(dirpath, 'path')
            if self.calculator == 'VASP':
                calc = VASP(dirpath,systemSpecies)
                calc.read_results()

            #Added for LAMMPS compatibility
            if self.calculator == 'LAMMPS':
                calc = LAMMPS(dirpath,systemSpecies)
                print('made it here!!')
                calc.read_results()
#            calc.add_to_results('energyLAMMPS',calcTwo.energy)
#            calc.add_to_results('energyLAMMPSPureA',calcTwo.pureA.energy)
#            calc.add_to_results('energyLAMMPSPureB',calcTwo.pureB.energy)
            if calc.crystal.results is not None:
                self.crystals.append(calc.crystal)

        print(path.split(paths[0]))
        for i in calc.crystal.species:
            purepath = path.join(path.split(paths[0])[0],'pure' + i)
            pure = LAMMPS(purepath,systemSpecies)
            pure.read_results()
            self.crystals.append(pure.crystal)
            

    def starting_point(self,folderpath):
        from os import path
        from glob import glob
        from aBuild.utility import chdir

        with chdir(folderpath):
            dirs = glob('E.*')
        prevCalcs = [int(x.split('.')[1])  for x in dirs]
        prevCalcs.sort()
        if prevCalcs != []:
            return prevCalcs[-1] + 1
        else:
            return 1

        #    def write(self):

        
    def buildFolders(self,buildpath,calculator,runGetKpoints = True,foldername = 'E'):
        from os import path
        from aBuild.calculators.vasp import VASP
        from aBuild.calculators.lammps import LAMMPS
        from aBuild.calculators.espresso import ESPRESSO
        from aBuild.jobs import Job

        import os
        print("Building folders in {}".format(buildpath))
        if not path.isdir(buildpath):
                os.mkdir(buildpath)
                print('Made path:',buildpath)
        configIndex = startPoint = self.starting_point(buildpath)
        for crystal in self.crystals:
            if 'vasp' in calculator["name"]:
                vaspspecs = {"incar":calculator["incar"],"kpoints":calculator["kpoints"], 'potcar':calculator["potcars"],"crystal":crystal}
                thisVASP = VASP(vaspspecs,self.species)
            if 'lammps' in calculator["name"]:
                print('lammps triggered')
                specsDict = {"crystal":crystal, "potential":calculator["potential"]}
                thisLAMMPS = LAMMPS(specsDict,self.species)
                
            if 'qe' in calculator["name"]:
                print('espresso triggered')
                specsDict = {"crystal":crystal, "pseudopotentials":calculator["pseudopotentials_qe"]}
                thisESPRESSO = ESPRESSO(specsDict,self.species)
            

                
            runpath = path.join(buildpath,foldername + ".{}".format(configIndex) )
            if not path.isdir(runpath):
                os.mkdir(runpath)
            else:
                msg.fatal("I'm gonna write over top of a current directory. ({})  I think I'll stop instead.".format(runpath)) 
            print("Building folder for structure: {}".format(crystal.title) )
            with chdir(runpath):
                print(calculator["name"])
                if 'vasp' in calculator["name"]:
                    print('building for vasp')
                    thisVASP.buildFolder(runGetKPoints = runGetKpoints)
                if 'lammps' in calculator["name"]:
                    print('building for lammps')
                    thisLAMMPS.buildFolder()
                if 'qe' in calculator["name"]:
                    print('building for espresso')
                    thisESPRESSO.buildFolder()
            configIndex += 1

        exdir = path.join(buildpath,'E.')
        mljob = Job(calculator["execution"],exdir,calculator["execution"]["exec_path"], arrayStart = startPoint,arrayEnd = configIndex - 1)
        with chdir(buildpath):
            print('Building job file')
            mljob.write_jobfile()






        
    


    def build_relax_select_input(self):
        from os import remove,path
        from aBuild.enumeration import Enumerate
        from aBuild.database.crystal import Crystal
        from aBuild.fitting.mtp import MTP
        from aBuild.utility import unpackProtos,getAllPerms
        from glob import glob
        fittingRoot = path.join(self.root,'fitting','mtp')
        
        for ilat  in range(self.nEnums):
            lat = self.enumDicts[ilat]["lattice"]
            enumLattice = Enumerate(self.enumDicts[ilat])

            if lat == 'protos':
                structures = getProtoPaths()
                for struct in structures:
                    scrambleOrder = getAllPerms(self.knary,justCyclic = 'uniqueUnaries' in struct)
                    for scramble in scrambleOrder:
                        thisCrystal = Crystal(struct,species = self.species)
                        thisCrystal.scrambleAtoms(scramble)
                        thisMTP = MTP(fittingRoot,dataSet = [thisCrystal],forRelax=True)
                        with open(path.join(fittingRoot,'to-relax.cfg'),'a+') as f:
                            f.writelines(thisMTP.lines)
                
            else:
                for struct in range(1,enumLattice.nEnumStructs+1):
                    enumLattice.generatePOSCAR(struct)
                    thisCrystal = Crystal.fromPOSCAR(enumLattice.root, self.species,
                                                     filename = "poscar.{}.{}".format(lat,struct),
                                                     title = ' '.join([lat," str #: {}"]).format(struct))
                    thisMTP = MTP(fittingRoot,dataSet = [thisCrystal],forRelax=True)
                    with open(path.join(fittingRoot,'to-relax.cfg'),'a+') as f:
                        f.writelines(thisMTP.lines)

                    delpath = path.join(enumLattice.root,"poscar.{}.{}".format(lat,struct))
                    remove(delpath)
                
        thisMTP.write_relaxin()



    def writeReport(self):
        import datetime
        nAtoms = len(self.crystals[0].species)
        with open('dataReport_' + self.calculator + '.txt', 'w') as f:
            f.write(self.calculator + ' REPORT\n')
            f.write(str(datetime.datetime.now()) + '\n')
            f.write("# S. Number               F. Enthalpy                Conc.     S. Energy.   E(A)    E(B)   A atoms      B atoms\n")
            f.write('------------------------------------------------------------------------------------------------------------------\n')
            for crystal in self.crystals[:-nAtoms]:
                f.write(crystal.reportline)
                for i in range(-nAtoms,0,1):
                    print('writing pure energies')
                    f.write(" {:8.5f}".format(self.crystals[i].results["energy"]))
                #f.write(" {:8.5f}".format(self.crystals[-nAtoms].results["energy"]))
                for i in range(nAtoms):
                    f.write(" {:2d}".format(crystal.atom_counts[i]))
                    #f.write(" {:2d}".format(crystal.atom_counts[1]))
                f.write('\n')










                
