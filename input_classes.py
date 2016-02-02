# pyQbox - Input/Output-Tools for Qbox
# Copyright (c) 2014, Andreas W. Hauser, Matthew Goldey
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
#                      pyQbox - Input Classes                      #
#                                                                   #
#####################################################################

# This file contains all input file classes and their methods.

import numpy as _np
import pybel
import constants
import running_scripts

############################# RUNDATA  ###############################

class _rundata(object):
    def __init__(self):
        self.name=''
        self.qbox=''
        self.nt=1
        self.np=1
        self.timestamp=False

    def __str__(self):
        ret_str = "Submission status summary:\n" + 26*"-" + "\n\n"
        if self.name!='':
            ret_str +=  "Filename is " + self.name + "\n"
        else:
            ret_str += "No filename provided, will use timestamp instead\n"
        if self.qbox!='':
            ret_str += "Qbox version is " + self.qbox + "\n"
	if self.nt!=1:
		ret_str += "Using " + str(self.nt) + " threads\n"
	if self.np!=1:
		ret_str += "Using " + str(self.np) + " processors\n"
        return ret_str

    def info(self):
        print self


########################### MULTIFILE  ##############################

class multifile(object):

    def __init__(self, jobs=[]):
        self.list_of_jobs=[]
        self.list_of_content=[]
        for k in jobs:
            self.add(k)

    def add(self,new_job):
        ''' Adds an inputfile to your batch object.'''
        if type(new_job) == type(inputfile()):
            self.list_of_jobs.append(new_job)
        else:
            print "Can only add inputfiles."

    def remove(self,position=0): #if not specified delete last
        ''' Removes an inputfile from your batch object. If no other specified the last is removed.'''
        del self.list_of_content[position]
        del self.list_of_jobs[position]

    def __str__(self):
        if self.list_of_jobs==[]:
            ret_str =  "empty"
        else:
            ret_str = self.list_of_jobs[0].__str__()
        if len(self.list_of_jobs)>1:
            for k in self.list_of_jobs[1:]:
                ret_str += "# job spacer \n\n" + k.__str__()
        return ret_str

    def write(self,filename):
        '''Writes the batch jobfile to disk.'''
        f = open(filename,'w')
        str_ret = self.__str__()
        print >>f, str_ret
        f.close()

    def run(self,name='',qbox='',nt=1,np=1,timestamp=False):
        '''Makes qbox process the given batch inputfile object. Optional parameters are

        name  ...... filename (without file extension, will be \".in\" and \".out\" by default)
        nt ......... number of threads
        np ......... number of processors.
        timestamp... adds a timestamp to input and output if set True.

        If nothing specified, pyQbox will fall back on information in the corresponding runinfo object.'''

        running_scripts._run(self,name,qbox,nt,np,timestamp)


########################### INPUTFILE  ##############################

class inputfile(object):

    def __init__(self, arrays=[]):
        self.list_of_arrays=[]
        self.list_of_content=[]
        self.runinfo=_rundata()
        self.__jtype="undef"
        for k in arrays:
            self.add(k)

    def add(self,new_array):
        ''' Adds an array to your inputfile object.'''
        if type(new_array) == type(variable_array()):
            self.variables = new_array
            if "variables" in self.list_of_content:
                index = self.list_of_content.index("variables")
                self.list_of_arrays[index]=new_array
            else: #change this to allow sets of commands and variables
                self.list_of_content.append("variables")
                self.list_of_arrays.append(new_array)

        elif type(new_array) == type(command_array()):
            self.commands = new_array
            if "commands" in self.list_of_content:
                index = self.list_of_content.index("commands")
                self.list_of_arrays[index]=new_array
            else:
                self.list_of_content.append("commands")
                self.list_of_arrays.append(new_array)

        elif type(new_array) == type(cell_array()):
            self.cell = new_array
            if "cell" in self.list_of_content:
                index = self.list_of_content.index("cell")
                self.list_of_arrays[index]=new_array
            else:
                self.list_of_content.append("cell")
                self.list_of_arrays.append(new_array)

        elif type(new_array) == type(ref_cell_array()):
            self.ref_cell = new_array
            if "ref_cell" in self.list_of_content:
                index = self.list_of_content.index("ref_cell")
                self.list_of_arrays[index]=new_array
            else:
                self.list_of_content.append("ref_cell")
                self.list_of_arrays.append(new_array)

        elif type(new_array) == type(cartesian()):
            self.cartesian = new_array
            if "cartesian" in self.list_of_content:
                index = self.list_of_content.index("cartesian")
                self.list_of_arrays[index]=new_array
            else:
                self.list_of_content.append("cartesian")
                self.list_of_arrays.append(new_array)
        elif type(new_array) == type(pseudo_array()):
            self.pseudo = new_array
            if "pseudo" in self.list_of_content:
                index = self.list_of_content.index("pseudo")
                self.list_of_arrays[index]=new_array
            else:
                self.list_of_content.append("pseudo")
                self.list_of_arrays.append(new_array)
        elif type(new_array) == type(comment_array()):
            self.comments = new_array
            if "comments" in self.list_of_content:
                index = self.list_of_content.index("comments")
                self.list_of_arrays[index]=new_array
            else:
                self.list_of_content.append("comments")
                self.list_of_arrays.append(new_array)
        elif type(new_array) == type(_unsupported_array()):
            self.list_of_content.append(str(new_array.type))
            self.list_of_arrays.append(new_array)
        else:
            print "Array type unknown."

    def remove(self,position=0): #if not specified delete last
        ''' Removes an array from your inputfile object. If no other specified the last is removed.'''
        del self.list_of_content[position]
        del self.list_of_arrays[position]

    def __str__(self):
        ret_str = ""
        try:
            ret_str +=self.pseudo.__str__()+"\n"
        except:
            ret_str +="# No pseudo definitions"+"\n"
        try:
            ret_str +=self.comments.__str__()+"\n"
        except:
            ret_str +="# No comment" +"\n"
        try:
            ret_str +=self.cell.__str__()+"\n"
        except:
            ret_str +="# No cell definition"+"\n"

        try:
            ret_str +=self.ref_cell.__str__()+"\n"
        except:
            ret_str +="# No ref_cell definition"+"\n"

        try:
            ret_str +=self.cartesian.__str__()+"\n"
        except:
            ret_str +="# No cartesian definition"+"\n"
        try:
            ret_str +=self.variables.__str__()+"\n"
        except:
            ret_str +="# No variable definitions"+"\n"
        try:
            ret_str +=self.commands.__str__()+"\n"
        except:
            ret_str +="# No commands"+"\n"

        #        for k in self.list_of_arrays:
        #            ret_str += k.__str__() + "\n"
        return ret_str

    def write(self,filename):
        f = open(filename,'w')
        str_ret = self.__str__()
        print >>f, str_ret
        f.close()

    def info(self):
        '''A quick overview of your inputfile.''' # Health check could be put here
        if "command" in self.list_of_content:
            status = "valid enough"
        else:
            status = "invalid"

        print "Type: inputfile"
        print "Status: " + status

    def run(self,name='',qbox='',nt=1,np=1,timestamp=False):
        '''Makes Qbox process the given batch inputfile object. Optional parameters are

        name  ...... filename (without file extension, will be \".in\" and \".out\" by default)
        nt ......... number of threads
        np ......... number of processors.
        timestamp... adds a timestamp to input and output if set True.

        If nothing specified, pyQbox will fall back on information in the corresponding runinfo object.'''

        running_scripts._run(self,name,qbox,nt,np,timestamp)

    def __add__(self,other):
        #autoadd subarrays - works
        import copy
        new=copy.deepcopy(self)
        new.add(other)
        return new

#    def __radd__(self,other):
#	#not working right, but general form should be this
#	import copy
#	new=copy.deepcopy(self)
#	new.add(other)
#	return new

######################## INPUT FRAGMENTS ############################

class _array(object):

    def __init__(self,content="undef"):
        self.content = content

    def __str__(self):
        ret_str += str(self.content) + "\n"
        return ret_str

    def write(self,filename):
        f = open(filename,'w')
        str_ret = self.__str__()
        print >>f, str_ret
        f.close()

    def __add__(self,other):
        if isinstance(other,_array):
            a=inputfile()
            a.add(self)
            a.add(other)
            return a
        if isinstance(other,inputfile):
            return other+self
    def __radd__(self,other):
	if isinstance(other,_array):
	    a=inputfile()
	    a.add(self)
	    a.add(other)
	    return a

##################### UNSUPPORTED FRAGMENT ##########################

class _unsupported_array(_array):

    def __init__(self,arraytype="undef"):
        self.content = []
        self.type = arraytype

    def add_line(self,line):
        self.content.append(line.strip())

    def __str__(self):
        for k in self.content:
            ret_str += k + "\n"
        return ret_str

####################### CARTESIAN FRAGMENT ##########################

class cartesian(_array):
    def __init__(self,title="",atom_list=[]):
        import copy
        self.__title = title
        self.__Natoms = 0
        self.xyzs=[]
        self.com=_np.array([0.0,0.0,0.0])
        self.centroid=_np.array([0.0,0.0,0.0])
        self.list_of_atoms = copy.deepcopy(atom_list)
        if atom_list!=[]:
            self.__Natoms=len(atom_list)
            for i in xrange(self.__Natoms):
                x=self.list_of_atoms[i][1]
                y=self.list_of_atoms[i][2]
                z=self.list_of_atoms[i][3]
                self.xyzs.append(_np.array([float(x),float(y),float(z)]))
            self.xyzs=_np.array(self.xyzs)
            self.__center_of_mass()

    def fix(self):
        """This fixes any odd errors resulting from modifying the number of atoms"""
        self.__Natoms=len(self.list_of_atoms)
        if self.__Natoms==0:
            return
        self.list_of_atoms.sort()
        self.xyzs=[]
        for i in xrange(self.__Natoms):
            x=self.list_of_atoms[i][1]
            y=self.list_of_atoms[i][2]
            z=self.list_of_atoms[i][3]
            self.xyzs.append(_np.array([float(x),float(y),float(z)]))
        self.xyzs=_np.array(self.xyzs)
        self.__center_of_mass()

    def __center_of_mass(self):
        """This computes the centroid and center of mass using standard atomic masses"""
        #print self.xyzs, self.__Natoms
        self.com=_np.array([0.0,0.0,0.0])
        self.centroid=_np.array([0.0,0.0,0.0])
        if len(self.xyzs)==0:
            return
        total_mass=0.0
        self.centroid=sum(self.xyzs)/len(self.xyzs)
        wts=[constants.dict_of_atomic_masses[self.list_of_atoms[i][0].replace("@","")]  for i in xrange(self.__Natoms)]
        for i,atom in enumerate(self.xyzs):
            wt=wts[i]
            total_mass=total_mass+wt
            self.com=self.com+atom*wt
        self.centroid=_np.array([i/self.__Natoms for i in self.centroid])
        self.com=_np.array([i/total_mass for i in self.com])

    def title(self,title="show"):
        if title == "show":
            return self.__title
        else:
            self.__title=title

    def add_atom(self,name="H",x="0",y="0",z="0",loc=[]):
        if len(loc)!=0:
            loc=[ str(i) for i in loc]
            self.list_of_atoms.append([name,loc[0],loc[1],loc[2]])
        else:
            self.list_of_atoms.append([name,x,y,z])
        self.fix()
        self.__center_of_mass()

    def remove_atom(self,position=0):
        del self.list_of_atoms[position]  # First atom is atom 1
        self.fix()
        self.__center_of_mass()

    def atoms(self):
        for i, k in enumerate(self.list_of_atoms):
            print  k[0] + "\t" + str(k[1]) + "\t" + str(k[2]) + "\t" + str(k[3])

    def atomic_distance(self,a,b):
	"""Gives the pair-wise distance between two atoms (counting from 0)"""
	from math import sqrt
	d=self.xyzs[b]-self.xyzs[a]
	return sqrt(d.dot(d))

    def print_centroid(self):
        print str(self.centroid[0])+'\t'+str(self.centroid[1])+'\t'+str(self.centroid[2])

    def print_center_of_mass(self):
        print str(self.com[0])+'\t'+str(self.com[1])+'\t'+str(self.com[2])

    def move(self,dir,amt=1.0):
        dir=_np.array(dir)
        for i in xrange(self.__Natoms):
            self.xyzs[i]=self.xyzs[i]+dir*amt
            self.list_of_atoms[i][1]=str(self.xyzs[i][0])
            self.list_of_atoms[i][2]=str(self.xyzs[i][1])
            self.list_of_atoms[i][3]=str(self.xyzs[i][2])
        self.__center_of_mass()

    def __str__(self,angstrom=False):
        #str_ret = str(self.__Natoms) + "\n" + self.__title + "\n"
        str_ret=""
        atoms=_np.array(self.list_of_atoms)
        if angstrom==True:
            for k in atoms:
                str_ret += k[0] + " " +str(k[1])+ " " +str(k[2])+ " " +str(k[3])+ "\n"
            return str_ret
        names=atoms.T[0].tolist()
        kinds=set()
        for i in names:
            kinds.add(i)
        kinds=list(kinds)
        kinds.sort()
        counts=[names.count(i) for i in kinds]
        counters=_np.zeros(len(kinds),dtype=_np.int)

        for k in atoms:
            myctr=counters[kinds.index(k[0])]
            str_ret += "atom " + k[0] + str(myctr+1) + " " + constants.dict_abbr_to_name[k[0]].lower() + "   " + str(float(k[1])*1.8897543761) + "    " + str(float(k[2])*1.8897543761) + "    " + str(float(k[3])*1.8897543761) + "\n"
            counters[kinds.index(k[0])]+=1
        return str_ret

    def __sub__(self,other):
        if type(other)==type([]):  #let's move the atoms
            self.move(other,-1.0)
        if type(other)==type(self.com):  #let's move the atoms using a numpy array
            self.move(other,-1.0)

    def __add__(self,other):
        if type(other)==type([]):  #let's move the atoms
            self.move(other,1.0)
        if type(other)==type(self.com):  #let's move the atoms using a numpy array
            self.move(other,1.0)
        if type(other)==type(self):                      #merge two cartesians
            atoms=self.list_of_atoms+other.list_of_atoms
            return cartesian(atom_list=atoms)
        if isinstance(other,_array):
            return other+mol_array(self)

    def __radd__(self,other):  #reverse of above
        if type(other)==type([]):
            self.move(other)
        if type(other)==type(self):
            atoms=self.list_of_atoms+other.list_of_atoms
            return cartesian(atom_list=atoms)
        if isinstance(other,_array):
            return other+mol_array(self)

######################### CELL FRAGMENT ############################

class cell_array(_array):

    def __init__(self):
        self.cell=_np.array(_np.zeros(9).reshape((3,3)))

    def setABC(self,A,B,C):
        self.cell[0][0]=A
        self.cell[1][1]=B
        self.cell[2][2]=C

    def set_cell(self,cell):
        self.cell=_np.array(cell).reshape(3,3)

    #add more fancy things with slants, maybe bravais lattices if time

    def __str__(self):
        ret_str="set cell "
        for i in self.cell:
            for j in i:
                ret_str+=str(j*1.8897543761)+' '
        ret_str+="\n"
        return ret_str


class ref_cell_array(cell_array):

    def __str__(self):
        ret_str="set ref_cell "
        for i in self.cell:
            for j in i:
                ret_str+=str(j*1.8897543761)+' '
        ret_str+="\n"
        return ret_str

######################### PSEUDO FRAGMENT ############################

class pseudo_array(_array):

    def __init__(self):
        self.dict_of_atoms = {}

    def add(self,atom,line):
        if constants.dict_abbr_to_name.has_key(atom):
            atom=constants.dict_abbr_to_name[atom]
        if self.dict_of_atoms.has_key(atom):
            self.dict_of_atoms[atom]=line
        else:
            self.dict_of_atoms[atom]=line

    def remove(self,atom):
        if self.dict_of_atoms.has_key(atom):
            del self.dict_of_atoms[atom]

    def __str__(self):
        ret_str=""
        for key in self.dict_of_atoms:
            ret_str += "species " + str(key).lower() +" "+self.dict_of_atoms[key]+"\n"
        return ret_str

####################### COMMENT FRAGMENT ############################

class comment_array(_array):

    def __init__(self,content=""):
        self.content = content+'\n'

    def add(self,value):
        '''\nFor values without documenation herein, please add keyword and value manually'''
        self.content+=value+'\n'

    def clear(self):
        '''Removes all keywords from array.'''
        self.content=""

    def __str__(self):
        ret_str=""
        lines = filter(bool,self.content.split("\n")) #remove empty list entries
        for line in lines:
            ret_str += "#   " + line.strip() + "\n"
        return ret_str

######################### VARIABLE FRAGMENT ##############################

class variable_array(_array):

    __tabstop = 30
    def __init__(self,variable_init=""):
        self.dict_of_keywords = {}
        variable_init=variable_init.splitlines()
        if len(variable_init)!=0:
            for i in variable_init:
                i=i.split(" ")
                if len(i)==0:
                    i=i.split("=")
                if i[0].startswith("$"):
                    continue
                self.add(i[0],i[1])
    def add(self,keyword,value):
        '''\nFor values without documenation herein, please add keyword and value manually'''
        self.dict_of_keywords[keyword.lower()]=value.upper()

    def remove(self,keyword):
        del self.dict_of_keywords[keyword.upper()]

    def clear(self):
        '''Removes all keywords from array.'''
        self.dict_of_keywords.clear()

    def __str__(self):
        ret_str=""
        for key,value in self.dict_of_keywords.iteritems():
            ret_str += "set " +key.lower() + " " + value + "\n"
        return ret_str

    def info(self):
        print "Type: variable array"
        print "Keywords: " + str(len(self.dict_of_keywords))



######################### COMMAND FRAGMENT ##############################

class command_array(_array):

    __tabstop = 30

    def __init__(self,command_init=""):
        self.list_of_keywords = []
        command_init=command_init.splitlines()
    	if len(command_init)!=0:
    		for i in command_init:
    			j=i.split(" ")
    			if len(j)==0:
    				i=i.split("=")
    			if j[0].startswith("$"):
    				continue
    			self.add(j[0],i.replace(j[0],""))

    def add(self,keyword,value):
        '''\nFor values without documenation herein, please add keyword and value manually'''
        self.list_of_keywords.append([keyword.lower(),value])

    def clear(self):
        '''Removes all keywords from array.'''
        self.list_of_keywords=[]

    def __str__(self):
        ret_str=""
        for key,value in self.list_of_keywords:
            ret_str += key.lower() + " " + value + "\n"
        return ret_str

    def info(self):
        print "Type: command array"
        print "Keywords: " + str(len(self.list_of_keywords))
