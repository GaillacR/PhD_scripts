#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np

# Atomic masses from Blue Obelisk's elements.xml
atMasses = { "X": 0.0000, "H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182, "B": 10.811, "C": 12.0107, "N": 14.0067,
  "O": 15.9994, "F": 18.9984032, "Ne": 20.1797, "Na": 22.989770, "Mg": 24.3050, "Al": 26.981538, "Si": 28.0855, "P": 30.973761,
  "S": 32.065, "Cl": 35.453, "Ar": 39.948, "K": 39.0983, "Ca": 40.078, "Sc": 44.955910, "Ti": 47.867, "V": 50.9415, "Cr":
  51.9961, "Mn": 54.938049, "Fe": 55.845, "Co": 58.933200, "Ni": 58.6934, "Cu": 63.546, "Zn": 65.409, "Ga": 69.723, "Ge": 72.64,
  "As": 74.92160, "Se": 78.96, "Br": 79.904, "Kr": 83.798, "Rb": 85.4678, "Sr": 87.62, "Y": 88.90585, "Zr": 91.224, "Nb":
  92.90638, "Mo": 95.94, "Tc": 98, "Ru": 101.07, "Rh": 102.90550, "Pd": 106.42, "Ag": 107.8682, "Cd": 112.411, "In": 114.818,
  "Sn": 118.710, "Sb": 121.760, "Te": 127.60, "I": 126.90447, "Xe": 131.293, "Cs": 132.90545, "Ba": 137.327, "La": 138.9055,
  "Ce": 140.116, "Pr": 140.90765, "Nd": 144.24, "Pm": 145, "Sm": 150.36, "Eu": 151.964, "Gd": 157.25, "Tb": 158.92534, "Dy":
  162.500, "Ho": 164.93032, "Er": 167.259, "Tm": 168.93421, "Yb": 173.04, "Lu": 174.967, "Hf": 178.49, "Ta": 180.9479, "W":
  183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217, "Pt": 195.078, "Au": 196.96655, "Hg": 200.59, "Tl": 204.3833, "Pb": 207.2,
  "Bi": 208.98038, "Po": 209, "At": 210, "Rn": 222, "Fr": 223, "Ra": 226, "Ac": 227, "Th": 232.0381, "Pa": 231.03588, "U":
  238.02891, "Np": 237, "Pu": 244, "Am": 243, "Cm": 247, "Bk": 247, "Cf": 251, "Es": 252, "Fm": 257, "Md": 258, "No": 259, "Lr":
  262, "Rf": 261, "Db": 262, "Sg": 266, "Bh": 264, "Hs": 277, "Mt": 268, "Ds": 281, "Rg": 272 }



###########################################
# Romain Atoms and distance code
###########################################
# Define class Distance (between two atoms)
class Distance(object):
  r = float(0)
  elt1 = 'H'
  elt2 = 'H'

  def __init__(self, r, elt1, elt2):
    self.r = r
    self.elt1 = elt1
    self.elt2 = elt2

def make_distance(r,elt1,elt2):
  d = Distance(r,elt1,elt2)
  return d

def print_distance(d):
  print "Distance value:", d.r, "First element:", d.elt1, "Second element:", d.elt2


# Define class Point

class Point(object):
  x = float(0)
  y = float(0)
  z = float(0)
  elt = 'H'

  def __init__(self, x, y, z, elt):
    self.x = x
    self.y = y
    self.z = z
    self.elt = elt

def make_point(x,y,z,elt):
  point = Point(x,y,z,elt)
  return point

def print_point(point):
  print point.elt, point.x, point.y, point.z
  return 0

def get_coord(point):
  coord=[]
  coord.append(point.x)
  coord.append(point.y)
  coord.append(point.z)
  return coord

# Define class Cell and some functions
class Cell(object):
  a = float(10)
  b = float(10)
  c = float(10)
  alpha = float(np.pi/2)
  beta = float(np.pi/2)
  gamma = float(np.pi/2)
  C=[0 for  x in range(6)]

  # The class "constructor" - It's actually an initializer
  def __init__(self, a, b, c, alpha, beta, gamma):
    self.a = a
    self.b = b
    self.c = c
    self.alpha = alpha
    self.beta = beta
    self.gamma = gamma


def make_cell(a, b, c, alpha, beta, gamma):
    cell = Cell(a, b, c, alpha, beta, gamma)
    return cell

def make_cell_C(cell):
  cell.C[0]=cell.a
  cell.C[1]=cell.b * np.cos(cell.gamma)
  cell.C[2]=cell.b * np.sin(cell.gamma)
  cell.C[3]=cell.c * np.cos(cell.beta)
  cell.C[4]=cell.c * (np.cos(cell.alpha) - np.cos(cell.beta)*np.cos(cell.gamma)) / np.sin(cell.gamma)
  cell.C[5]=cell.c * math.sqrt(1 - (cell.C[3]/cell.c)**2 - (cell.C[4]/cell.c)**2)

  return cell.C

def red2cart(cart, red, cell):
  cart[0] = red[0]*cell.C[0] + red[1]*cell.C[1] + red[2]*cell.C[3]
  cart[1] = red[1]*cell.C[2] + red[2]*cell.C[4]
  cart[2] = red[2]*cell.C[5]

def cart2red(cart, red, cell):
  red[0] = (cart[0]/cell.C[0]) - ((cell.C[1])/(cell.C[0]*cell.C[2]))*cart[1] + ((cell.C[1]*cell.C[4]-cell.C[2]*cell.C[3])/(cell.C[0]*cell.C[2]*cell.C[5]))*cart[2]
  red[1] = -(cart[2]*cell.C[4])/(cell.C[2]*cell.C[5]) + cart[1]/cell.C[2]
  red[2] = cart[2]/cell.C[5]

def distance_cart(r1,r2,cell):
  r3 = ClosestImage(r1,r2,cell)
  dr = [r3[i]-r1[i] for i in range(3)]
  return np.sqrt(np.dot(dr,dr))

def distance_red(r1,r2,cell):
  c=[]
  append = c.append
  for i in range(3):
    append(r2[i]-r1[i])
  for i in range(3):
    c[i]=c[i]-round(c[i])
  red2cart(dr,c,cell)
  d_red=math.sqrt(np.dot(dr,dr))
  return d_red

def ClosestImage(r1,r2,cell):
  r2_tmp = r2[:]
  v1 = [0.0 for p in range(3)]
  v2 = [0.0 for p in range(3)]
  dr = [r2[i]-r1[i] for i in range(3)]
  c = [0.0 for i in range(3)]
  cart2red(dr,c,cell)
  cart2red(r1,v1,cell)
  cart2red(r2,v2,cell)
  d_cart = [0.0 for k in range(8)]
  op = c[:]
  for x in range(len(c)):
    if (abs(c[x]) <= 1):
      if (c[x] >= 0):
        op[x] = -1
      else:
        op[x] = 1
    else:
      op[x] = -round(c[x])

  sol = -1
  if (cell.alpha==cell.beta and cell.beta==cell.gamma and abs(cell.alpha-(np.pi/2.0)) < 1e-5):
    v2[0] = v2[0] - round(c[0])
    v2[1] = v2[1] - round(c[1])
    v2[2] = v2[2] - round(c[2])
    sol = 0
  else:
    c2 = [c[0],c[1],c[2]]
    red2cart(dr,c2,cell)
    d_cart[0]=math.sqrt(np.dot(dr,dr))

    c2 = [c[0]+op[0],c[1],c[2]]
    red2cart(dr,c2,cell)
    d_cart[1]=math.sqrt(np.dot(dr,dr))

    c2 = [c[0]+op[0],c[1]+op[1],c[2]]
    red2cart(dr,c2,cell)
    d_cart[2]=math.sqrt(np.dot(dr,dr))

    c2 = [c[0]+op[0],c[1]+op[1],c[2]+op[2]]
    red2cart(dr,c2,cell)
    d_cart[3]=math.sqrt(np.dot(dr,dr))

    c2 = [c[0],c[1]+op[1],c[2]]
    red2cart(dr,c2,cell)
    d_cart[4]=math.sqrt(np.dot(dr,dr))

    c2 = [c[0],c[1]+op[1],c[2]+op[2]]
    red2cart(dr,c2,cell)
    d_cart[5]=math.sqrt(np.dot(dr,dr))

    c2 = [c[0],c[1],c[2]+op[2]]
    red2cart(dr,c2,cell)
    d_cart[6]=math.sqrt(np.dot(dr,dr))

    c2 = [c[0]+op[0],c[1],c[2]+op[2]]
    red2cart(dr,c2,cell)
    d_cart[7]=math.sqrt(np.dot(dr,dr))

    sol = d_cart.index(min(d_cart))
  if (sol==1):
    v2[0] = v2[0]+op[0]
  elif(sol==2):
    v2[0] = v2[0]+op[0]
    v2[1] = v2[1]+op[1]
  elif(sol==3):
    v2[0] = v2[0]+op[0]
    v2[1] = v2[1]+op[1]
    v2[2] = v2[2]+op[2]
  elif(sol==4):
    v2[1] = v2[1]+op[1]
  elif(sol==5):
    v2[1] = v2[1]+op[1]
    v2[2] = v2[2]+op[2]
  elif(sol==6):
    v2[2] = v2[2]+op[2]
  elif(sol==7):
    v2[0] = v2[0]+op[0]
    v2[2] = v2[2]+op[2]
  elif(sol<0):
    print sol
    print "Impossible value for sol, should be between 0 and 7"
    exit(1)

  red2cart(r2,v2,cell)
  dr = [r2[i]-r1[i] for i in range(3)]

  if (abs(np.sqrt(np.dot(dr,dr))) > min([cell.a,cell.b,cell.c])):
    print np.sqrt(np.dot(dr,dr))
    print min(d_cart)
    print sol
    print cell.alpha, cell.beta, cell.gamma
    exit(1)

  return r2

def NeiAtoms(A, atoms, cell, n=2, dmax=20):
  Nei = [make_point(0.0,0.0,0.0,'H') for i in range(n+1)]
  Nei[0]=A
  for k in range(1,n+1):
    if (len(atoms)>0):
      b = [get_coord(atoms[i]) for i in range(len(atoms))]
      d = [0.0 for i in range(len(atoms))]
      a = get_coord(A)
      for i in range(len(b)):
        d[i] = distance_cart(a,ClosestImage(a,b[i],cell),cell)
      sortdist = sorted(d)
      if (sortdist[0] < dmax):
        neiat = ClosestImage(a,b[d.index(sortdist[0])],cell)
        Nei[k] = make_point(neiat[0],neiat[1],neiat[2], atoms[d.index(sortdist[0])].elt)
        A = atoms[d.index(sortdist[0])]
        atoms.remove(A)
      else:
        if (k==0): print "No", k+1, "st neighbour found within", dmax, "Angström"
        elif (k==1): print "No", k+1, "nd neighbour found within", dmax, "Angström"
        elif (k==2): print "No", k+1, "rd neighbour found within", dmax, "Angström"
        else: print "No", k+1, "th neighbour found within", dmax, "Angström"
        print get_coord(A)
        exit(1)
  if (make_point(0.0,0.0,0.0,'H') in Nei): Nei.remove(make_point(0.0,0.0,0.0,'H'))
  return Nei

def Barycenter(list):
  coords = [get_coord(list[i]) for i in range(len(list))]
  elts = [list[i].elt for i in range(len(list))]
  coord = [0.0 for x in range(3)]
  Mass = 0
  for k in range(len(list)):
    coord = [coord[i]+coords[k][i]*atMasses[elts[k]] for i in range(3)]
    Mass = Mass + atMasses[elts[k]]
  coord = [coord[i]/Mass for i in range(3)]
  return coord

def Angle(A,B,C,cell):
  pos = [get_coord(A), get_coord(B), get_coord(C)]
  v1 = pos[0]
  v2 = ClosestImage(v1,pos[1],cell)
  v3 = ClosestImage(v2,pos[2],cell)
  V1 = [v2[h]-v1[h] for h in range(3)]
  V2 = [v2[h]-v3[h] for h in range(3)]
  Theta = math.acos(np.dot(V1,V2)/(math.sqrt(np.dot(V1,V1))*math.sqrt(np.dot(V2,V2))))*180.0/np.pi
  return Theta

def Norm(v):
  return math.sqrt(v[0]**2+v[1]**2+v[2]**2)


def FormFact(q,elt):
  felt = 0.0
  if (elt=="H"):
    a = [0.489918,0.262003,0.196767,0.049879]
    b = [20.6593,7.74039,49.5519,2.20159]
    c = 0.001305
  elif (elt=="C"):
    a = [2.31,1.02,1.5886,0.865]
    b = [20.8439,10.2075,0.5687,51.6512]
    c = 0.2156
  elif (elt=="N"):
    a = [12.2126,3.1322,2.0125,1.1663]
    b = [0.0057,9.8933,28.9975,0.5826]
    c = -11.529
  elif (elt=="Zn"):
    a = [14.0743,7.0318,5.1652,2.41]
    b = [3.2655,0.2333,10.3163,58.7097]
    c = 1.3041
  for i in range(len(a)):
    felt = felt + a[i]*np.exp(-b[i]*(q/(4*np.pi))**2)
  felt = felt + c
  if (abs(felt) > abs(sum(a)+c)):
    print felt
    print a
    print b
    print q
    print abs(sum(a))+c
    print elt
    exit(1)

  return felt



# Define function to get used elements in a set of atoms
def get_elements(list_points):
  used_elements=[]
  append = used_elements.append
  for i in range(len(list_points)):
    if (i==0):
      append(list_points[i].elt)
    if (i>0):
      if (list_points[i].elt not in used_elements):
        append(list_points[i].elt)
  return used_elements
