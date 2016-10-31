#!/usr/bin/env python3

def read_quartz(fname='file07.dat'):

    with open(fname,'r') as fin:
        lines = fin.readlines()

    xyzfrac = []
    names = []
    import numpy as np
    for iline,line in enumerate(lines):
        sline = line.split()
        if len(sline[0]) == 30: # Coordinate line
            xyzfrac.append([float(sline[i]) for i in range(2,5)])
            names.append(splist[int(sline[-1])-1])
            if names[-1] == 'SI':
                names[-1] = 'Si'
        elif iline == 7:
            box  = np.array([float(sline[i]) for i in range(1,4)])
            angle = np.array([float(sline[i]) for i in range(4,7)])
        elif iline == 2:
            splist = [sline[i] for i in range(len(sline))]

    xyzfrac = np.asarray(xyzfrac)
    # Convert fractional positions to Angstroms, taking into account possible non-orthoganal axes.
    xyz = fractional_to_cartesian(box,angle,xyzfrac)

    zeros = np.zeros(len(names))
    ones = np.ones(len(names),dtype=int)

    dictlist = {'box': box, 'pos': xyz, 'names': names, 'file': fname, 'chg': zeros, 'mol': ones, 'angle': angle}

    # Calculate density - TODO : Update for non-diagonal boxes and for new definition of "angle" (ie, it's now cos(angle))
    #mass = {'Si': 28.085, 'O': 15.999} # g/mol
    #avogadro = 6.022141e23
    #totmass = np.sum([mass[aname] for aname in names]) / avogadro / 1000. # mass in kg
    #ang = 2.*np.pi*angle[2]/360.
    #volume = box[0] * np.sin(ang)*box[2] * box[1] * 1.0e-30
    #density = totmass / volume
    #print 'Density: {0:7.2f} kg/m3'.format(density)

    return Struct(**dictlist)


def fractional_to_cartesian(box,angle,xyzfrac):
    aa, bb, cc = box[0], box[1], box[2]
    cosalpha, cosbeta, cosgamma = angle[0], angle[1], angle[2]
    import numpy as np
    xyz = np.zeros([len(xyzfrac),3])
    from numpy import sin, arccos, sqrt
    singamma = sin(arccos(cosgamma))
    for iatom in range(len(xyzfrac)):
        ix, iy, iz = xyzfrac[iatom,0], xyzfrac[iatom,1], xyzfrac[iatom,2]
        xyz[iatom,0] = ix * aa + iy * bb * cosgamma + iz * cc * cosbeta
        xyz[iatom,1] = iy * bb * singamma + iz * (cc * (cosalpha - cosbeta * cosgamma)) / singamma
        ww = sqrt(1. - cosalpha**2 - cosbeta**2 - cosgamma**2 + 2. * cosalpha * cosbeta * cosgamma)
        xyz[iatom,2] = iz * cc * ww / singamma

    return xyz


class Struct(object): # http://stackoverflow.com/questions/1305532/convert-python-dict-to-object
    def __init__(self, **entries): 
        self.__dict__.update(entries)


if __name__ == '__main__':
    coords = read_quartz()
    print(coords.box)
    print(coords.pos)
    print(coords.names)
    print(coords.file)
    print(coords.chg)
    print(coords.mol)
    print(coords.angle)

