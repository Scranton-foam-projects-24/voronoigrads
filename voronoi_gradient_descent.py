# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 13:43:38 2024

@author: Joseph Klobusicky
"""



#This is an attempt at Voronoi coarsening when the boundary edges (including edges that only have one vertex on the bounary)
#are ignored

    



import numpy as np
import matplotlib.pyplot as plt
#import random as z
import scipy.interpolate as interpolate
import scipy as sp
import scipy.special as spec
import scipy.stats as stats
import csv
import pandas as pd
import scipy.integrate as integrate
import pylab
import random



import pyvoro

import matplotlib
import matplotlib.pyplot as plt


np.random.seed([1938430])



# dots_num = 75
# colors = np.random.rand(dots_num, 3) # or predefined 
# points = np.random.rand(dots_num, 2)



# points = [ [.4,.4], [.4, .6],[.6,.4], [.7, .6]]
# points = np.array(points)
# colors = np.random.rand(2, 3) # or predefined 




# #For hexagonal ICs
# points = []
# spacing = 2
# vec1 = [np.sqrt(3)/2, 1/2]
# vec2 = [np.sqrt(3)/2, -1/2]
# pt = [0,0]
# noiz = .1
# for j in range(-spacing,spacing):
#     for k in range(-spacing,spacing):
#         pt = np.array([.5+((j+k)*np.sqrt(3)/2)/spacing+ noiz*np.random.uniform(), .5+((j-k)/2)/spacing+ noiz*np.random.uniform()])
#         if pt[0] >= 0 and pt[0] <= 1 and pt[1] >= 0 and pt[1] <= 1 :
#             points.append(pt)
            
# points = np.array(points)
            
        
# print(len(points))
# dots_num = len(points)
# colors = np.random.rand(dots_num, 3) 




# #For random point ics

# points = [[.5, .48], [.5, .62], [.61, .53], [.39, .5], [.52, .4]]
points = [.51, .52], [.57, .51], [.42, .56], [.41, .39]
points = np.array(points)
dots_num = len(points)
colors = np.random.rand(dots_num, 3) 





# #two cloud initial points:
   
# noiz = .03
# numpoints = 30
# pt1ct = [1/3, 1/3]
# pt2ct = [2/3,2/3]
# points1 = [[pt1ct[0] + noiz*np.random.normal(), pt1ct[1] + noiz*np.random.normal()] for i in range(numpoints)]
# points2 = [[pt2ct[0] + noiz*np.random.normal(), pt2ct[1] + noiz*np.random.normal()] for i in range(numpoints)]
# points = points1 + points2
# points = np.array(points)
# dots_num = 2*numpoints
# colors = np.random.rand(dots_num, 3) 


 

# make color map (unnecessary for just random colorization)
color_map = {tuple(coords):color for coords, color in zip(points, colors)}

cells = pyvoro.compute_2d_voronoi(
points, # point positions, 2D vectors this time.
[[0.0, 1.0], [0.0, 1.0]], # box size
2.0, # block size
periodic = [False, False]
)

# colorize
for i, cell in enumerate(cells):    
    polygon = cell['vertices']
    plt.fill(*zip(*polygon),  color = 'black', alpha=0.1)
    # plt.fill(*zip(*polygon), color = color_map[tuple(cell['original'])], alpha=0.5)

plt.plot(points[:,0], points[:,1], 'ko')
plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 1.1)

plt.show()











#This just gives some definitions of distance adjusted for periodic cell


def pdisttwo(x,y):
    return(  np.sqrt( (x[0]-y[0])**2 + (x[1]-y[1])**2))



def pnormtwo(x):
    return( x[0]**2+ x[1]**2)


# def pdiff(v,w):
#     z = v-w
    
#     if (abs(z[0])>1/2):
#         if z[0]>0:
#             z[0] = z[0]-1
#         else:
#             z[0] = z[0]+1

#     if (abs(z[1])>1/2):
#         if z[1]>0:
#             z[1] = z[1]-1
#         else:
#             z[1] = z[1]+1
#     return v-w






# def matdiff(v,w):
#     z = v-w
    
#     if (abs(z[0,0])>1/2):
#         z[0,0] = z[0,0]-1
#     if (abs(z[1])>1/2):
#         z[1] = z[1]-1
#     return z


def indy(i,j):
    return np.ix_(  [2*i, 2*i+1] , [2*j, 2*j+1])   



def perp(v):
    return np.array([-v[1],v[0]])



def vco(i):
    return vertices[edges[i]]





#Exponent on potential function

potexp = 3
def gradufun(v,w):
    return -(v-w)/  (pdisttwo(v,w)**potexp)





stepz = 3000


ptmatrix = np.zeros( shape = (stepz+1, len(points), 2)  )

ptmatrix[0,:,:] = points

Energy = 0

envec = []




typlength = .01
alphaparm = (typlength)**2

time = np.zeros(stepz)

dt =   (1/10)*(typlength/100) 


for kk in range(stepz):


    cells = pyvoro.compute_2d_voronoi(
    points, # point positions, 2D vectors this time.
    [[0.0, 1.0], [0.0, 1.0]], # box size
    2.0, # block size
    periodic = [False, False]
    )

    
    print(kk)
    print("energy is", Energy)
    
    
    Energy = 0
    
    
    
    edges = set()
    
    
    vertices = [C['vertices'] for C in cells]
    
    
    
    #flatten and pluck out boundary vertices
    vertices = [val for sublist in vertices for val in sublist]
    vertices = [   [float('%.8f'%(item[0]% 1)), float('%.8f'%(item[1]% 1))] for item in vertices]
    
    vertices = [p for p in vertices if p[0] != 0 and p[0] != 1 and p[1] != 0 and p[1] != 1]
    
    

         
    # #unique, sorted, and arrayed
    vertices = sorted(list(set(tuple(p) for p in vertices)))
    
    vertices = [np.array(q) for q in vertices]
    
    tripoint = [set() for i in range(len(vertices))]
    
    

    #this just gets x-coords for vertices
    E = [j[0] for j in vertices]
    
    if min(np.diff(E)) == 0:
        
        
        points = np.array([ [p[0]+ (10**(-7))*np.random.normal(), p[1]+ + (10**(-7))*np.random.normal()*np.random.normal()] for p in points])
        print('4vertex!')
        continue 

        
    
    
    
    gradu = np.zeros(2*len(points))
    
    

    

    
    
    #get vertex labels, as well as Delaunay lengths.
    
    
    
        
    
    for i, C in enumerate(cells):   
        
        
        
        #this is very wrong, need to look over each vertex one by one
        vc = C['vertices']
        vc =  np.array([   [float('%.8f'%(p[0]% 1)), float('%.8f'%(p[1]% 1))] for p in vc ])
        isint = np.array([(p[0] != 0 and p[0] != 1 and p[1] != 0 and p[1] != 1) for p in vc])
        vc = [np.searchsorted(E, q[0]) for q in vc]
        newedge = { frozenset([vc[j], vc[j+1]])  for j in range(len(vc)-1) if (isint[j] and isint[j+1])}     
        if isint[0] and isint[-1] and sum(isint)> 1:
            newedge = newedge.union({ frozenset( [vc[-1], vc[0]])})

        
    
        edges = edges.union(newedge)
        
        #add to grad u, here we are using a potential of 1/r, with a gradient of (xi-xj)/r^2
        
        vz = [C['faces'][q]['adjacent_cell'] for q in range(len(C['faces']))]
        vz = [p for p in vz if p >= 0]
        gradu[(2*i):(2*i+2)] = sum(   [    gradufun(points[i],points[q])   for q in vz]      )
        Energy += np.sum(   [    1/pdisttwo(points[i],points[q])   for q in vz]      )
        
 
    
        #add to tripoint
        
        #contract vc to only contain good points
        vc = np.array(vc)[isint]
        for ind in vc:
            tripoint[ind] = tripoint[ind].union({i})
    
    #Adjust to alphaparm
    Energy = Energy*alphaparm
    
    edges = list(edges)
    tripoint = [list(q) for q in tripoint]
    

    #great, can we now import over the old matrix functions and see if they work?
    
    
    #computing the total perimeter
    
    
    Energy += sum( [pdisttwo(vertices[list(a)[0]],vertices[list(a)[-1]]) for a in edges ]  )
    
    envec += [Energy]
    
    gradv = np.zeros(shape = (2*len(vertices)))
    
    
    #need to adjust for periodicity in subtraction
    for e in edges:
        ee = list(e)
        a = ee[0]
        b = ee[-1]
        
        if a != b:
            gradv[(2*a):(2*a+2)] += (vertices[a]-vertices[b])/ pdisttwo(vertices[a],vertices[b])
            gradv[(2*b):(2*b+2)] += (vertices[b]-vertices[a])/ pdisttwo(vertices[a],vertices[b])
    
    
    QT = np.zeros((2*len(vertices), 2*len(vertices)))
    
    
    #cycle through edges and add contributions
    
    
    for e in edges:
        ee = list(e)
        v = vertices[ee[0]]
        w = vertices[ee[-1]]
        
        if ee[0] != ee[-1]:
            QT[indy(ee[0],ee[1])] += (1/3)*(1/pdisttwo( v,w  ))*np.outer(perp(v-w),perp(v-w))
            QT[indy(ee[1],ee[0])] += (1/3)*(1/pdisttwo( v,w  ))*np.outer(perp(v-w),perp(v-w))
            QT[indy(ee[0],ee[0])] += (1/3)*(1/pdisttwo( v,w  ))*np.outer(perp(v-w),perp(v-w))
            QT[indy(ee[1],ee[1])] += (1/3)*(1/pdisttwo( v,w  ))*np.outer(perp(v-w),perp(v-w))
        
        
    
    
    
    #What about the Jacobian?

    
    P = np.zeros((2*len(vertices), 2*len(points)))
    
    
    #loop over vertices
    
    
    def trishift(a,b,c):
        q = [a[0], b[0], c[0], a[1], b[1], c[1]]
        # if max([  abs(a[0]-b[0]) ,  abs(a[0]-c[0]), abs(b[0]-c[0])])      > .5 :
        #     q[0] += (a[0]<.3)
        #     q[1] += (b[0]<.3)
        #     q[2] += (c[0]<.3)
            
            
        # if max([  abs(a[1]-b[1]) ,  abs(a[1]-c[1]), abs(b[1]-c[1])])      > .5 :
        #     q[3] += (a[1]<.3)
        #     q[4] += (b[1]<.3)
        #     q[5] += (c[1]<.3)
            
        return([[q[0], q[3]], [q[1], q[4]], [q[2], q[5]] ])
            
    
    for l in range(len(vertices)):
        
        #Don't loop over values in tripoint, just do it one at a time, three times    
        zz = trishift(points[tripoint[l][0]], points[tripoint[l][1]], points[tripoint[l][2]])
        pi = zz[0]
        pj = zz[1]
        pk = zz[2]
        
        
        #Now plug in values into Jacobian
            
        #easiest case
        
        P[indy(l,tripoint[l][0])]   = ( ( np.dot(pk, pk)- np.dot(pj, pj)) *np.array([[0,-1],[1,0]]) + 2*  np.outer(perp((np.array(pj)-np.array(pk))), pi)   ) /  (  2* (np.dot(pi, perp(pj)) + np.dot(pj, perp(pk)) + np.dot(pk, perp(pi)) )  )  -    np.outer(   ( np.dot(pk, pk)- np.dot(pj, pj)) *perp(pi)  +  (np.dot(pi, pi)- np.dot(pk, pk)) *perp(pj) +  (np.dot(pj, pj)- np.dot(pi, pi)) *perp(pk)   , perp( np.array(pj)-np.array(pk) )    )      /   (  2* ((np.dot(pi, perp(pj)) + np.dot(pj, perp(pk)) + np.dot(pk, perp(pi)) ) )**2 )  
        
        
        
        P[indy(l,tripoint[l][1])]   = ( ( np.dot(pi, pi)- np.dot(pk, pk)) *np.array([[0,-1],[1,0]]) + 2*  np.outer(perp((np.array(pk)-np.array(pi))), pj)   ) /  (  2* (np.dot(pi, perp(pj)) + np.dot(pj, perp(pk)) + np.dot(pk, perp(pi)) )  )  -    np.outer(   ( np.dot(pk, pk)- np.dot(pj, pj)) *perp(pi)  +  (np.dot(pi, pi)- np.dot(pk, pk)) *perp(pj) +  (np.dot(pj, pj)- np.dot(pi, pi)) *perp(pk)   , perp( np.array(pk)-np.array(pi) )    )      /   (  2* ((np.dot(pi, perp(pj)) + np.dot(pj, perp(pk)) + np.dot(pk, perp(pi)) ) )**2 )      
        
        P[indy(l,tripoint[l][2])]   = ( ( np.dot(pj, pj)- np.dot(pi, pi)) *np.array([[0,-1],[1,0]]) + 2*  np.outer(perp((np.array(pi)-np.array(pj))), pk)   ) /  (  2* (np.dot(pi, perp(pj)) + np.dot(pj, perp(pk)) + np.dot(pk, perp(pi)) )  )  -    np.outer(   ( np.dot(pk, pk)- np.dot(pj, pj)) *perp(pi)  +  (np.dot(pi, pi)- np.dot(pk, pk)) *perp(pj) +  (np.dot(pj, pj)- np.dot(pi, pi)) *perp(pk)   , perp( np.array(pi)-np.array(pj) )    )      /   (  2* ((np.dot(pi, perp(pj)) + np.dot(pj, perp(pk)) + np.dot(pk, perp(pi)) ) )**2 )  
                              
    
        
        
    
            
    #On to matrix multiplication: ought to be (Pt QT P)^-1   PT (gradv)t
    
    
    #something's off with gradE
    
    # gradE =      np.linalg.inv(np.transpose(P).dot(QT).dot(P)).dot(np.transpose(P)).dot(np.transpose(gradv))
    
    # gradE = np.reshape(gradE, (dots_num, 2))
    
    gradnog = np.transpose(P).dot(np.transpose(gradv))
        
    
    gradnog = np.reshape(gradnog, (dots_num, 2))
    
    gradu = np.reshape(gradu, (dots_num, 2))
    
    
    fullgrad = gradnog+ alphaparm*gradu
    
    ptmatrix[kk+1,:,:] = points
    
    
    #depending on whether to include perimeter, potentials, or both
    
    # points = (points - dt* gradnog) % 1
    
    # points = (points - dt* gradu) % 1
    
    points = (points - dt* fullgrad) % 1

    









# print(gradnog)
plt.plot(envec)




#plot points again
# make color map (unnecessary for just random colorization)
color_map = {tuple(coords):color for coords, color in zip(ptmatrix[0,:,:], colors)}

cells = pyvoro.compute_2d_voronoi(
ptmatrix[0,:,:], # point positions, 2D vectors this time.
[[0.0, 1.0], [0.0, 1.0]], # box size
2.0, # block size
periodic = [False, False]
)

# colorize
for i, cell in enumerate(cells):    
    polygon = cell['vertices']
    plt.fill(*zip(*polygon), color = color_map[tuple(cell['original'])], alpha=0.5)

for kk in range(stepz):
    plt.plot(ptmatrix[kk,:,0], ptmatrix[kk,:,1], 'ko', color = 'red', markersize = .1)

plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 1.1)

plt.show()



# r =6
# for kk in range(stepz):
#     plt.plot(ptmatrix[kk,r,0], ptmatrix[kk,r,1], 'ko', color = 'red', markersize = .1)

# plt.xlim(-0.1, 1.1)
# plt.ylim(-0.1, 1.1)

# plt.show()





#Now to plot site point trajectories




# make color map (unnecessary for just random colorization)
color_map = {tuple(coords):color for coords, color in zip(points, colors)}

cells = pyvoro.compute_2d_voronoi(
points, # point positions, 2D vectors this time.
[[0.0, 1.0], [0.0, 1.0]], # box size
2.0, # block size
periodic = [False, False]
)

# colorize
for i, cell in enumerate(cells):    
    polygon = cell['vertices']
    plt.fill(*zip(*polygon),  color = 'black', alpha=0.1)
    # plt.fill(*zip(*polygon), color = color_map[tuple(cell['original'])], alpha=0.5)

plt.plot(points[:,0], points[:,1], 'ko')
plt.xlim(-0.1, 1.1)
plt.ylim(-0.1, 1.1)

plt.show()
