# pyQbox - Input/Output-Tools for Qbox
# Copyright (c) 2014, Andreas W. Hauser
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# The views and conclusions contained in the software and documentation are those
# of the authors and should not be interpreted as representing official policies, 
# either expressed or implied, of the FreeBSD Project.

#####################################################################
#                                                                   #
#                 pyQbox - Tools and Subroutines                   #
#                                                                   #
#####################################################################

# This file contains a few hidden methods for easy file reading with
# one reading method for each type of array.

# Additionally, tools for comparing geometries are included for
# cartesian objects


######################### STANDARD MODULES ##########################

from copy import deepcopy
import numpy as _np

############################# MODULES ###############################

# Let's start with importing the inputfile classes. 
# Hidden classes have to be imported explicitely.

from input_classes import *
from input_classes import _unsupported_array
from input_classes import _array
import pybel
from constants import dict_of_atomic_abbr

def _readpybel(extension,filename):
    infile = pybel.readfile(extension,filename).next()
    periodic=False
    try:
        if (infile.unitcell.GetValue()=='UnitCell'):
            periodic=True
    except:
        periodic=False
    cvec=_np.zeros(9).reshape(3,3)
    if periodic:
        cmat=infile.unitcell.GetCellMatrix()
        for i in xrange(3):
            cvec[i][0]=cmat.GetRow(i).GetX()
            cvec[i][1]=cmat.GetRow(i).GetY()
            cvec[i][2]=cmat.GetRow(i).GetZ()

    Natoms = len(infile.atoms)
    #atoms=_np.array([[dict_of_atomic_abbr[i.atomicnum]] for i in infile.atoms])    
    atoms=_np.array([[i.atomicnum] for i in infile.atoms])
    xyzs=_np.array([list(i.coords) for i in infile.atoms]).reshape(Natoms,3)
    atom_list=_np.concatenate((atoms,xyzs),axis=1)
    atom_list=atom_list.tolist()
    for num,i in enumerate(atom_list):
        atom_list[num][0]=dict_of_atomic_abbr[i[0]]

    re_file = inputfile()
    re_file.add(cartesian(atom_list=atom_list))    

    cell=cell_array()
    cell.set_cell(cvec)
    re_file.add(cell)

    return deepcopy(re_file)
    del re_file
    
def _readinput(file_input,silent=False):
    import re
    regex = re.compile('[^a-zA-Z]')
    #Check input type
    if type(file_input)==list:
        full_content = file_input
    else:
        infile = open(file_input,"r")
        full_content = infile.readlines()      

    blocks = []
    for num, line in enumerate(full_content):
        if line.strip().lower().startswith("run"):
            blocks.append(num+1)

    istart=0
    ifinish=len(full_content)
    if len(blocks)>1:
        re_file=multifile()        
        for i in blocks:
            ifinish=i
            ifile=_readinput(full_content[istart:ifinish])
            istart=i
            re_file.append(ifile)
            return deepcopy(re_file)
    else:
        re_file = inputfile() 
        carts=cartesian()
        pseudos =pseudo_array()
        commands=command_array()
        comments=comment_array()
        variables=variable_array()
        cell=cell_array()

        for i in full_content:
            i=i.strip()
            if i.lower().startswith("set cell"):
                cell.set_cell(_np.array(i.split()[2:],dtype='f'))
            elif i.lower().startswith("species"):
                pseudos.add(i.split()[1],i.split()[2])
            elif i.lower().startswith("atom"):
                j=i.split()
                k=regex.sub("",j[1])
                loc=_np.array(j[3:],dtype='f')
                carts.add_atom(k,loc=loc)
            elif i.startswith("#"):
                continue #keeping comments is for the "future"
            #elif i.lower().startswith("set"):
                #command_array
        re_file.add(cell)
        re_file.add(carts)
        re_file.add(pseudos)
        re_file.add(commands)
        re_file.add(variables)
        # Return the final inputfile object
        return deepcopy(re_file)
        del re_file
    

#Utilities for comparing xyz files taken from https://github.com/charnley/rmsd

#Copyright (c) 2013, Jimmy Charnley Kromann <jimmy@charnley.dk> & Lars Bratholm
#All rights reserved.

#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:

#1. Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

def fit(P, Q):
    """ Varies the distance between P and Q, and optimizes rotation for each step
    until a minimum is found.
    """
    step_size = P.max(0)
    threshold = step_size*1e-9
    rmsd_best = kabsch(P, Q)
    while True:
        for i in range(3):
            temp = numpy.zeros(3)
            temp[i] = step_size[i]
            rmsd_new = kabsch(P+temp, Q)
            if rmsd_new < rmsd_best:
                rmsd_best = rmsd_new
                P[:,i] += step_size[i]
            else:
                rmsd_new = kabsch(P-temp, Q)
                if rmsd_new < rmsd_best:
                    rmsd_best = rmsd_new
                    P[:,i] -= step_size[i]
                else:
                    step_size[i] /= 2
        if (step_size<threshold).all():
            break
    return rmsd_best
    
def rmsd(V, W):
    """ Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return numpy.sqrt(rmsd/N)
    
def kabsch(P, Q):
    """ The Kabsch algorithm

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    The algorithm starts with two sets of paired points P and Q.
    P and Q should already be centered on top of each other.

    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    The optimal rotation matrix U is then used to
    rotate P unto Q so the RMSD can be caculated
    from a straight forward fashion.

    """

    # Computation of the covariance matrix
    C = numpy.dot(numpy.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = numpy.linalg.svd(C)
    d = (numpy.linalg.det(V) * numpy.linalg.det(W)) < 0.0

    if(d):
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    # Create Rotation matrix U
    U = numpy.dot(V, W)

    # Rotate P
    P = numpy.dot(P, U)

    return rmsd(P,Q)

def kabsch2(Q, P):
    """ The Kabsch algorithm

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    The algorithm starts with two sets of paired points P and Q.
    P and Q should already be centered on top of each other.

    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    The optimal rotation matrix U is then used to
    rotate P unto Q so the RMSD can be caculated
    from a straight forward fashion.

    """

    # Computation of the covariance matrix
    C = numpy.dot(numpy.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = numpy.linalg.svd(C)
    d = (numpy.linalg.det(V) * numpy.linalg.det(W)) < 0.0

    if(d):
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    # Create Rotation matrix U
    U = numpy.dot(V, W)

    # Rotate P
    P = numpy.dot(P, U)

    return rmsd(P,Q),P
        
#End of utilities for comparing xyz files taken from https://github.com/charnley/rmsd    


#Let's add some utilities for reorganizing xyz files
def swap(cart,i,j):
    cart2=deepcopy(cart)
    temp=deepcopy(cart2.xyzs[i])
    cart2.xyzs[i]=deepcopy(cart2.xyzs[j])
    cart2.xyzs[j]=temp
    temp=deepcopy(cart2.list_of_atoms[i])
    cart2.list_of_atoms[i][0]=cart2.list_of_atoms[j][0]
    cart2.list_of_atoms[i][1]=cart2.list_of_atoms[j][1]
    cart2.list_of_atoms[i][2]=cart2.list_of_atoms[j][2]
    cart2.list_of_atoms[i][3]=cart2.list_of_atoms[j][3]
    cart2.list_of_atoms[j][0]=temp[0]
    cart2.list_of_atoms[j][1]=temp[1]
    cart2.list_of_atoms[j][2]=temp[2]
    cart2.list_of_atoms[j][3]=temp[3]
    return cart2

def order(P,Q):
    """Reordering Q to match structure P as best as possible"""
    natoms=len(P.xyzs)
    dist,Q.xyzs=kabsch2(P.xyzs,Q.xyzs)
    Q.move([0,0,0])
    for i in xrange(natoms):
        for j in xrange(natoms):
            if i==j:
                continue
            tmp=kabsch(P.xyzs,swap(Q,i,j).xyzs)
            if tmp<dist:
                Q=swap(Q,i,j)
                dist,Q.xyzs=kabsch2(P.xyzs,Q.xyzs)
                Q.move([0,0,0])
                print "Swapping atoms ",i," and ",j
    return Q

def orderset(cart_list):
    cnum=len(cart_list)
    for i in xrange(cnum):
        cart_list[i].move(-cart_list[i].com)
    for i in xrange(cnum-1):
        print "Comparing ",i," and ",i+1
        cart_list[i+1]=order(cart_list[i],cart_list[i+1])
        print "RMSD=",rmsd(cart_list[i].xyzs,cart_list[i+1].xyzs)
        print "Kabsch RMSD=",kabsch(cart_list[i].xyzs,cart_list[i+1].xyzs)
    return





