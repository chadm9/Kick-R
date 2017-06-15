# Kick-R v. 1.0
# Author: W. C. McKee
#This algorithm performs a random search for isomers of a given composition inside of a box  set by the user 
#in the 'Kick-R.in' file. All atoms in the search are randomly displaced inside of the box away from their
#user determined original positions.  Atomic starting coordinates are gennerally the origin or those of a molecule.
#A usage example is 'python Kick-R.py'.

import math, random

class Atom:

    def __init__(self, symbol, x, y, z):
        self.symbol = symbol
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.CovR1 = 0.0
        self.CovR3 = 0.0

    def getSymbol(self):
        return self.symbol

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def getZ(self):
        return self.z

    def getCovR1(self):
        return self.CovR1

    def getCovR3(self):
        return self.CovR3

    def setCovR1(self):
        covr1_dictionary = {"H":0.32, "He":0.46, "Li":1.33, "Be":1.02, "B":0.85, "C":0.75, "N":0.71,
        "O":0.63, "F":0.64, "Ne":0.96, "Na":1.60, "Mg":1.39, "Al":1.26, "Si":1.16, "P":1.11, "S":1.03,
        "Cl":0.99, "Ar":1.07, "K":1.96, "Ca":1.71, "Sc":1.48, "Ti":1.36, "V":1.34, "Cr":1.22, "Mn":1.19,
        "Fe":1.16, "Co":1.11, "Ni":1.10, "Cu":1.20, "Zn":1.20, "Ga":1.24, "Ge":1.21, "As":1.21,
        "Se":1.16, "Br":1.14, "Kr":1.21, "Rb":2.1, "Sr":1.85, "Y":1.63, "Zr":1.54, "Nb":1.47,
        "Mo":1.38, "Tc":1.28, "Ru":1.25, "Rh":1.25, "Pd":1.20, "Ag":1.39, "Cd":1.44, "In":1.46,
        "Sn":1.40, "Sb":1.40, "Te":1.36, "I":1.33, "Xe":1.35, "Cs":2.32, "Ba":1.96, "Hf":1.52,
        "Ta":1.46, "W":1.37, "Re":1.31, "Os":1.29, "Ir":1.22, "Pt":1.23, "Au":1.24, "Hg":1.42,
        "Tl":1.50, "Pb":1.44, "Bi":1.51, "Po":1.45, "At":1.47, "Rn":1.45, "Fr":2.23, "Ra":2.01,
        "Rf":1.57, "Db":1.49, "Sg":1.43, "Bh":1.41, "Hs":1.34, "Mt":1.29, "Ds":1.28, "Rg":1.21,
        "La":1.8, "Ce":1.63, "Pr":1.76, "Nd":1.74, "Pm":1.73, "Sm":1.72, "Eu":1.68, "Gd":1.69,
        "Tb":1.68, "Dy":1.67, "Ho":1.66, "Er":1.65, "Tm":1.64, "Yb":1.7, "Lu":1.62, "Ac":1.86,
        "Th":1.75, "Pa":1.69, "U":1.7, "Np":1.71, "Pu":1.72, "Am":1.66, "Cm":1.66, "Bk":1.68,
        "Cf":1.68, "Es":1.65, "Fm":1.67, "Md":1.73, "No":1.76, "Lr":1.61}

        if covr1_dictionary.has_key(self.symbol):
            self.CovR1 = covr1_dictionary[self.symbol]
        else:
            print "No single bond covalent radius found for atom " + self.symbol + ". Terminating Kick-R."
            exit()

    def setCovR3(self):
        covr3_dictionary = {"H":0.32, "He":0.46, "Li":1.24, "Be":0.85, "B":0.73, "C":0.60, "N":0.54,
        "O":0.53, "F":0.53, "Ne":0.67, "Na":1.55, "Mg":1.27, "Al":1.11, "Si":1.02, "P":0.94, "S":0.95,
        "Cl":0.93, "Ar":0.96, "K":1.93, "Ca":1.33, "Sc":1.14, "Ti":1.08, "V":1.06, "Cr":1.03, "Mn":1.03,
        "Fe":1.02, "Co":0.96, "Ni":1.01, "Cu":1.12, "Zn":1.18, "Ga":1.17, "Ge":1.11, "As":1.06,
        "Se":1.07, "Br":1.09, "Kr":1.08, "Rb":2.02, "Sr":1.39, "Y":1.24, "Zr":1.21, "Nb":1.16,
        "Mo":1.13, "Tc":1.1, "Ru":1.03, "Rh":1.06, "Pd":1.12, "Ag":1.28, "Cd":1.36, "In":1.36,
        "Sn":1.3, "Sb":1.27, "Te":1.21, "I":1.25, "Xe":1.22, "Cs":2.09, "Ba":1.49, "Hf":1.22,
        "Ta":1.19, "W":1.15, "Re":1.1, "Os":1.09, "Ir":1.07, "Pt":1.1, "Au":1.21, "Hg":1.33,
        "Tl":1.42, "Pb":1.35, "Bi":1.35, "Po":1.29, "At":1.38, "Rn":1.33, "Fr":2.18, "Ra":1.59,
        "Rf":1.31, "Db":1.26, "Sg":1.21, "Bh":1.19, "Hs":1.18, "Mt":1.13, "Ds":1.12, "Rg":1.18,
        "La":1.39, "Ce":1.31, "Pr":1.28, "Nd":1.37, "Pm":1.35, "Sm":1.34, "Eu":1.34, "Gd":1.32,
        "Tb":1.35, "Dy":1.33, "Ho":1.33, "Er":1.33, "Tm":1.31, "Yb":1.29, "Lu":1.31, "Ac":1.40,
        "Th":1.36, "Pa":1.29, "U":1.18, "Np":1.16, "Pu":1.35, "Am":1.35, "Cm":1.36, "Bk":1.39,
        "Cf":1.4, "Es":1.4, "Fm":1.67, "Md":1.39, "No":1.76, "Lr":1.41}

        if covr3_dictionary.has_key(self.symbol):
            self.CovR3 = covr3_dictionary[self.symbol]
        else:
            print "No triple bond covalent radius found for atom " + self.symbol + ". Terminating Kick-R."
            exit()

    def setX(self, value):
        self.x = float(value)

    def setY(self, value):
        self.y = float(value)

    def setZ(self, value):
        self.z = float(value)


def getKeyWords(input):
    multiplicity, charge, kick_number, basis_set, dft_type, spin_method = 1, 0, 250, "LANL2DZ", "B3LYP", "R"
    for i in input:
        if "multiplicity" in str(i):
            multiplicity = int(input[input.index(i)].split()[1])

        if "charge" in str(i):
            charge = int(input[input.index(i)].split()[1])

        if "kick_number" in str(i):
            kick_number = int(input[input.index(i)].split()[1])

        if "basis_set" in str(i):
            basis_set = str(input[input.index(i)].split()[1].rstrip())

        if "dft_type" in str(i):
            dft_type = str(input[input.index(i)].split()[1].rstrip())

        if multiplicity > 1:
            spin_method = "U"

    return multiplicity, charge, kick_number, basis_set, dft_type, spin_method

def structure_parameters(input):
    for i in input:
        if "structure" in str(i):
            no_of_atoms = int(input[input.index(i)].split()[1])
            structure_index = input.index(i) + 1

    try:
        return no_of_atoms, structure_index
    except UnboundLocalError:
        print "The 'structure' keyword is required, but was not found.  Terminating Kick-R."
        exit()

def set_molecule(structure_index, no_of_atoms, input):
     molecule = []

     for i in range(structure_index, structure_index + no_of_atoms):
        molecule.append(Atom(input[i].split()[0], input[i].split()[1], input[i].split()[2], input[i].split()[3]))

     for atoms in molecule:
        atoms.setCovR1()
        atoms.setCovR3()
     return molecule

def set_BoxSize(input, molecule):
     for i in input:
        if "box_size" in str(i):
            box_x = float(input[input.index(i)].split()[1].split("x")[0])
            box_y = float(input[input.index(i)].split()[1].split("x")[1])
            box_z = float(input[input.index(i)].split()[1].split("x")[2])
     try:
        return box_x, box_y, box_z
     except UnboundLocalError:
        box_x = box_y = box_z = 2*CovRad1Sum(molecule)
        return box_x, box_y, box_z

def CovRad1Sum(molecule):
    sum = 0.0
    for atoms in molecule:
        sum = sum + atoms.getCovR1()
    return sum

def PerformKick(box_x, box_y, box_z, structure_index, no_of_atoms, input):
    trial_structure = set_molecule(structure_index, no_of_atoms, input)

    for atoms in trial_structure:
        atoms.setX(atoms.getX()+ 0.5*box_x*random.triangular(-1, 1))
        atoms.setY(atoms.getY() + 0.5*box_y*random.triangular(-1, 1))
        atoms.setZ(atoms.getZ() + 0.5*box_z*random.triangular(-1, 1))
    return trial_structure

def InteratomicDistance(atom1, atom2):
    distance = math.sqrt((atom1.getX()-atom2.getX())**2+(atom1.getY()-atom2.getY())**2+(atom1.getZ()-atom2.getZ())**2)
    return distance

def CovR3Sum(atom1, atom2):
    sum = atom1.getCovR3() + atom2.getCovR3()
    return sum

def CovR1Sum(atom1, atom2):
    sum = atom1.getCovR1() + atom2.getCovR1()
    return sum

def DistanceCriteria(molecule):
    passed, i, j = True, 1, 1
    while i < len(molecule):
        while j < len(molecule):
            if InteratomicDistance(molecule[i-1], molecule[j]) < 0.9*CovR3Sum(molecule[i-1], molecule[j]):
                passed = False
                break
            j = j + 1
        if passed == False:
            break
        i = i + 1
        j = i
    return passed

def IsBonded(atom1, atom2):
    if InteratomicDistance(atom1, atom2) < 1.1*CovR1Sum(atom1, atom2):
        return True
    else:
        return False

def getGraph(molecule):
    graph = {}
    for atoms in range(len(molecule)):
        graph[atoms+1] = []
    for i in range(len(molecule)):
        for j in range(len(molecule)):
            if i!= j and IsBonded(molecule[i], molecule[j]):
                graph[i+1].append(j+1)
    return graph

def isConnected(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return True
    for node in graph[start]:
        if node not in path:
            newpath = isConnected(graph, node, end, path)
            if newpath: return True
    return False

def ConnectedCriteria(graph):
    passed, i, j = True, 1, 2
    while i < len(graph.keys()):
        while j <= len(graph.keys()):
            if not isConnected(graph, i, j):
                passed = False
                break
            j = j + 1
        if passed == False:
            break
        i = i + 1
        j = i + 1
    return passed

def writeKick(multiplicity, charge, basis_set, dft_type, spin_method, trial_structure, kick):
    output = open(str(kick) + '.kick', 'w')
    output.write("%mem=4gb\n%lindaworkers=LINDA\n")
    output.write("# " + spin_method + dft_type + "/" + basis_set + " Opt=maxcycles=50\n")
    output.write("\nTitle\n\n")
    output.write(str(charge)+ " " + str(multiplicity) + "\n")
    for atoms in trial_structure:
        output.write(str(atoms.getSymbol()) + "  " + str(atoms.getX()) + "\t" + str(atoms.getY()) + "\t" + str(atoms.getZ())+"\n")
    output.write("\nEND")
    output.close()

def writeSummary(multiplicity, charge, kick_number, basis_set, dft_type,ShortDistanceFailures, ConnectedFailures,
                 TotalFailures, molecule, box_x, box_y, box_z):
    output = open('Kick-R.out', 'w')
    output.write("Input Parameters\n")
    output.write("DFT Method: " + str(dft_type) + "\n")
    output.write("Basis Set: " + str(basis_set) + "\n")
    output.write("Charge: " + str(charge) + "\n")
    output.write("Multiplicity: " + str(multiplicity) + "\n")
    output.write("Requested Number of Kicks: " + str(kick_number) + "\n")
    output.write("Number of Atoms: " + str(len(molecule)) + "\n")
    output.write("Box Dimensions: " + str(box_x) + " x " + str(box_y) + " x " + str(box_z) + "\n")
    output.write("2x(Sum of Covalent Radii): " + str(2*CovRad1Sum(molecule)) + "\n")
    output.write("\nStarting Structure\n")
    for atoms in molecule:
        output.write(str(atoms.getSymbol()) + "  " + str(atoms.getX()) + "\t" + str(atoms.getY()) + "\t" + str(atoms.getZ())+"\n")
    output.write("\nKick-R Statistics\n")
    output.write("Number of kicks rejected due to short interatomic distances: " + str(ShortDistanceFailures) + "\n")
    output.write("Number of kicks rejected because not all atoms were connected: " + str(ConnectedFailures)+ "\n")
    output.write("Total number of kicks rejected: " + str(TotalFailures) + "\n")
    output.close()





def main():

    read_in = open('Kick-R.in', 'r')
    input = read_in.readlines()

    multiplicity, charge, kick_number, basis_set, dft_type, spin_method = getKeyWords(input)
    no_of_atoms, structure_index = structure_parameters(input)
    molecule = set_molecule(structure_index, no_of_atoms, input)
    box_x, box_y, box_z = set_BoxSize(input, molecule)
    ShortDistanceFailures, ConnectedFailures, TotalFailures, kick = 0, 0, 0, 1

    while kick <= kick_number:
        while True:
            trial_structure = PerformKick(box_x, box_y, box_z, structure_index, no_of_atoms, input)
            graph = getGraph(trial_structure)
            distance_passed, connected_passed = False, False

            if DistanceCriteria(trial_structure):
                distance_passed = True
            else:
                ShortDistanceFailures = ShortDistanceFailures + 1
            if ConnectedCriteria(graph):
                connected_passed = True
            else:
                ConnectedFailures = ConnectedFailures + 1
            if distance_passed and connected_passed:
                writeKick(multiplicity, charge, basis_set, dft_type, spin_method, trial_structure, kick)
                kick = kick + 1
                break
            else:
                TotalFailures = TotalFailures + 1

    writeSummary(multiplicity, charge, kick_number, basis_set, dft_type,ShortDistanceFailures, ConnectedFailures,
                 TotalFailures, molecule, box_x, box_y, box_z)





main()
