# IPython log file

import pyqbox as pq
import os
from copy import deepcopy

cell=pq.cell_array()
cell.setABC(8,8,8)
pseudo=pq.pseudo_array()

pseudo.add("Ar","http://fpmd.ucdavis.edu/potentials/Ar/Ar_HSCV_PBE-1.0.xml")
pseudo.add("C", "http://fpmd.ucdavis.edu/potentials/C/C_HSCV_PBE-1.0.xml")
pseudo.add("H", "http://fpmd.ucdavis.edu/potentials/H/H_HSCV_PBE-1.0.xml")
pseudo.add("N", "http://fpmd.ucdavis.edu/potentials/N/N_HSCV_PBE-1.0.xml")
pseudo.add("O", "http://fpmd.ucdavis.edu/potentials/O/O_HSCV_PBE-1.0.xml")
pseudo.add("B", "http://fpmd.ucdavis.edu/potentials/B/B_HSCV_PBE-1.0.xml")
pseudo.add("F", "http://fpmd.ucdavis.edu/potentials/F/F_HSCV_PBE-1.0.xml")

var=pq.variable_array()
var.add("xc","PBE")
var.add("wf_dyn","PSDA")
var.add("ecut","40")
var.add("ecutprec","5")
var.add("scf_tol","1e-8")

cmd=pq.command_array()
job_list=[]
for i in os.popen("ls geoms/*.xyz").read().splitlines():
    inp=pq.read(i,silent=True)
    name=i.split("/")[-1].replace(".xyz","")
    cmd.clear()
    cmd.add("run","0 500")
    cmd.add("save",name+".xml")
    inp.runinfo.name=name
    inp.runinfo.np=1
    inp.runinfo.nt=1
    inp.runinfo.qbox='~/Programs/qbox-vdw/src/qb'
    inp=inp+pseudo+cell+var+deepcopy(cmd)
    job_list.append(inp)
pq.queue(job_list,num_workers=8)

