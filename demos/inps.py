import pyqbox as pq
import os

cell=pq.cell_array()
cell.setABC(20,20,20)
pseudo=pq.pseudo_array()
pseudo.add("Ar","pseudos/Ar_HSCV_PBE-1.0.xml")
pseudo.add("C", "pseudos/C_HSCV_PBE-1.0.xml")
pseudo.add("H", "pseudos/H_HSCV_PBE-1.0.xml")
pseudo.add("N", "pseudos/N_HSCV_PBE-1.0.xml")
pseudo.add("O", "pseudos/O_HSCV_PBE-1.0.xml")

#pseudo.add("Ar","http://fpmd.ucdavis.edu/potentials/Ar/Ar_HSCV_PBE-1.0.xml")
#pseudo.add("C", "http://fpmd.ucdavis.edu/potentials/C/C_HSCV_PBE-1.0.xml")
#pseudo.add("H", "http://fpmd.ucdavis.edu/potentials/H/H_HSCV_PBE-1.0.xml")
#pseudo.add("N", "http://fpmd.ucdavis.edu/potentials/N/N_HSCV_PBE-1.0.xml")
#pseudo.add("O", "http://fpmd.ucdavis.edu/potentials/O/O_HSCV_PBE-1.0.xml")

var=pq.variable_array()
var.add("xc","PBE")
var.add("wf_dyn","PSDA")
var.add("ecut","85")
var.add("ecutprec","5")
var.add("scf_tol","1e-8")

cmd=pq.command_array()
for i in os.popen("ls geoms/*.xyz").read().splitlines():
    inp=pq.read(i,silent=True)
    cmd.clear()
    cmd.add("run","0 500")
    cmd.add("save",i.split("/")[-1].replace(".xyz",".xml"))
    inp=inp+pseudo+cell+var+cmd
    inp.write(i.split("/")[-1].replace(".xyz",".in"))

