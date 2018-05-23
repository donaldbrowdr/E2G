from dolfin import *
from fenics import *
import numpy as np
import mshr


t_sleeve=0.188
t_wall=0.188
t_gap=0.02
L_sleeve=10*t_wall
L_wall=1.5*L_sleeve

domain_vertices = [Point(0,0), Point(L_wall,0), Point(L_wall,t_wall), Point(0,t_wall)]
domain          = mshr.Polygon(domain_vertices)
domain_vertices = [Point(L_wall-L_sleeve,t_wall),Point(L_wall-L_sleeve,t_wall+t_gap) ,Point(L_wall,t_wall+t_gap), Point(L_wall,t_wall+t_gap+t_sleeve),
                       Point(L_wall-L_sleeve,t_wall+t_gap+t_sleeve),Point(L_wall-L_sleeve-t_gap-t_sleeve,t_wall)]
domain          =  domain+mshr.Polygon(domain_vertices)
#domain_vertices = [Point(L_wall-L_sleeve,t_wall+t_gap), Point(L_wall,t_wall+t_gap),Point(L_wall,t_wall+t_gap+t_sleeve) , Point(L_wall-L_sleeve,t_wall+t_gap+t_sleeve)]
#domain          =  domain+mshr.Polygon(domain_vertices)
# mesh resolution: higher numbers give more elements
mesh_resolution = 10
mesh            = mshr.generate_mesh(domain,mesh_resolution)
bmesh = BoundaryMesh(mesh, "exterior", True)
bmesh_coord=bmesh.coordinates()
print(bmesh_coord[1,1])
# Create VTK file for saving solution
vtkfile = File('heat_weld/mesh.pvd')

vtkfile<< mesh

#File("my_mesh.xml.gz") << mesh






#domain_vertices = [Point(0,0), Point(L_wall,0), Point(L_wall,t_wall), Point(0,t_wall)]
#domain          = mshr.Polygon(domain_vertices)
#domain_vertices = [Point(L_wall-L_sleeve,t_wall), Point(L_wall-L_sleeve,t_wall+t_gap+t_sleeve),Point(L_wall-L_sleeve-t_gap-t_sleeve,t_wall)]
#domain          =  domain+mshr.Polygon(domain_vertices)
#domain_vertices = [Point(L_wall-L_sleeve,t_wall+t_gap), Point(L_wall,t_wall+t_gap),Point(L_wall,t_wall+t_gap+t_sleeve) , Point(L_wall-L_sleeve,t_wall+t_gap+t_sleeve)]
#domain          =  domain+mshr.Polygon(domain_vertices)
