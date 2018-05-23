from __future__ import print_function
from fenics import *
import numpy as np
import mshr

T = 10.0            # final time
num_steps = 1000     # number of time steps
dt = T / num_steps # time step size
tol = 1E-14
# Physical Constants
c_p=.119;rho=.284;k=31.95/(3600.0*12.0)
T_Proc=325.0;h_Proc=48.0/(144.0*3600.0)
T_Amb=70.0;h_Amb=9.0/(144.0*3600.0)

# Boundary Size
t_wall=8
L_wall=8

# Create mesh and define function space
domain_vertices = [Point(0,0), Point(L_wall,0), Point(L_wall,t_wall), Point(0,t_wall)]
domain          = mshr.Polygon(domain_vertices)
# mesh resolution: higher numbers give more elements
mesh_resolution = 40
mesh            = mshr.generate_mesh(domain,mesh_resolution)

#Define function space (FEM Space)
V = FunctionSpace(mesh, 'P', 1)


# Construct facet markers
bndry = FacetFunction("size_t", mesh)
for f in facets(mesh):
    mp = f.midpoint()
    if near(mp[1], t_wall,tol): # Ambient
        bndry[f] = 1
    elif near(mp[1], 0.0,tol): # Process
        bndry[f] = 2
#Boundary Measure, ds(1) is Ambient, ds(2) is Process
ds = Measure("ds", subdomain_data=bndry)

# Define initial  condition and forcing
T_init = Constant(0.5*(T_Amb+T_Proc) )
f = Expression('(x[0] >= 0.45*L_wall  + tol && x[0] <= 0.55*L_wall  + tol) && x[1] >= 0.5*t_wall+tol ? 2700*0.5*(1 + cos((t+5)*pi/5) ) : 0.0', degree=2,
               tol=tol,L_wall=L_wall,t_wall=t_wall, t=0)


# Previous and current solution
T0 = interpolate (T_init , V )
T1 = Function ( V )

# Define variational problem
T = TrialFunction(V)
v = TestFunction(V)

#Define the variational form and right hand side. 
a = c_p*rho*T*v*dx +dt*( k*inner ( grad ( T ) , grad ( v ) )*dx+h_Amb*T*v*ds(1)+h_Proc*T*v*ds(2))
L = c_p*rho*T0*v*dx + dt*(f*v*dx+h_Amb*T_Amb*v*ds(1)+h_Proc*T_Proc*v*ds(2))

# Create VTK file for saving solution
vtkfile = File('heat_weld/solution.pvd')

t = 0
for n in range(num_steps):
    # Update forcing and solve
    f.t = t
    solve ( a == L , T1 )
    # Update
    T0.assign ( T1 )
    t += float ( dt )
    # Save to file
    vtkfile << (T1, t)
    

print('All is Done Boss')
