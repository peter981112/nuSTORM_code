# load pyg4ometry
import pyg4ometry
import numpy as np

def deg_to_rad(deg_ang):
    return deg_ang*np.pi/180

# registry to store gdml data
reg  = pyg4ometry.geant4.Registry()

# world solid and logical
r0 = 9.80493e4
rmin = 9.7049e4
rmax = 9.9049e4
ang_min = 2.352-1.700
ang_max = 2.352+1.700
pdz = r0*deg_to_rad(1.7*2)
ws   = pyg4ometry.geant4.solid.Tubs("ws", rmin, rmax, pdz, deg_to_rad(ang_min), deg_to_rad(ang_max), reg)#pyg4ometry.geant4.solid.Box("ws",400,400,400,reg), pLowNorm, pHighNorm
wl   = pyg4ometry.geant4.LogicalVolume(ws,"G4_Galactic","wl",reg)
reg.setWorld(wl.name)

# box placed at origin
b1   = pyg4ometry.geant4.solid.Tubs("b1", rmin, rmax, pdz, deg_to_rad(ang_min), deg_to_rad(ang_max), reg)#pyg4ometry.geant4.solid.Box("b1",300,300,300,reg)
b1_l = pyg4ometry.geant4.LogicalVolume(b1,"G4_Galactic","b1_l",reg)
b1_p = pyg4ometry.geant4.PhysicalVolume([0,0,0],[0,-r0,0],b1_l,"b1_p",wl,reg)#PhysicalVolume(rotation, position, logicalVolume, name, motherVolume, registry=None, copyNumber=0, addRegistry=True, scale=None)
print("geomatry created")

# Write geometry to a gdml file
out_file_name = "file_test.gdml"
w = pyg4ometry.gdml.Writer()
w.addDetector(reg)
w.write(out_file_name)
print("geomatry written to : "+out_file_name)

# visualise geometry
v = pyg4ometry.visualisation.VtkViewer()
v.addLogicalVolume(wl)
v.addAxes(20)
v.view()