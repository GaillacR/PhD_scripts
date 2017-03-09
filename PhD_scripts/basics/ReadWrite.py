#!/usr/bin/env python
# -*- coding: utf-8 -*-

import PBC_RG as pbc
import numpy as np
import gzip

# Define function to get lines of a file
def getlines(Cij_file):
  lines=[]
  append = lines.append
  tmp=""
  for line in Cij_file:
    if(line!=tmp) : append(line)
    tmp=line
  return(lines)

def readFromXYZFile(filename,step=1,first=1,last=int(1e15)):
  traj = list()
  try:
    f = open (filename, "r")
  except IOError as e:
    f = gzip.open(filename+".gz", "r")
  i = 0
  #nat = int(f.readline())
  while (i+1<=last):
    line = f.readline()
    if (len(line) == 0): break
    if (i==0): nat = int(line)
    f.readline()
    snap = list()
    if (i+1 >= first and (i+1-first)%step==0):
      traj.append(snap)
      for at in range(nat):
        l2 = f.readline()
        [elt, x, y, z] = l2.split()
        atom = pbc.make_point(float(x),float(y),float(z),elt)
        snap.append(atom)
    else:
      for at in range(nat):
        l2 = f.readline()
    i = i + 1
  return traj

def writeTrajAsXYZFile(filename, traj):
  f = open (filename, "w")
  count = -1
  for conf in traj:
    count = count + 1
    f.write(str(len(conf)) + "\n")
    f.write("i = "+ str(count) + "\n")
    for at in conf:
      f.write(at.elt + " " + str(at.x)
                       + " " + str(at.y)
                       + " " + str(at.z) + "\n")
  f.close()

def appendTrajAsXYZFile(filename, traj, first):
  try:
    f = open (filename, "a")
  except IOError as e:
    f = gzip.open(filename+".gz", "a")
  count = first-1
  for conf in traj:
    count = count + 1
    f.write(str(len(conf)) + "\n")
    f.write("i = "+ str(count) + "\n")
    for at in conf:
      f.write(at.elt + " " + str(at.x)
                       + " " + str(at.y)
                       + " " + str(at.z) + "\n")
  f.close()

def readFromDistanceFile(filename, nconf):
  traj = list()
  f = open (filename, "r")
  while True:
    last_pos = f.tell()
    line = f.readline()
    if (len(line) == 0): break
    if (len(line)>2):
      snap = list()
      traj.append(snap)
      f.seek(last_pos)
      for at in range(nconf):
        [t, d] = f.readline().split()
        snap.append(float(d))
  return np.array(traj)


def nconfFromDist(filename):
  foo = 0
  nconf = -1
  with open(filename,"r") as dist:
    for line in dist:
      sp = line.split()
      if (len(sp)>0 and int(sp[0])+1>nconf): nconf=int(sp[0])+1
      if (len(sp)==0): foo = foo+1
      if (foo==2): break
  return nconf

def readCell(filename):
  with open(filename,"r") as cell_file:
    for line in cell_file:
      sp=line.split()
      if ("a" in sp):
        a=float(sp[len(sp)-1])
      if ("b" in sp):
        b=float(sp[len(sp)-1])
      if ("c" in sp):
        c=float(sp[len(sp)-1])
      if ("alpha" in sp):
        alpha=float(sp[len(sp)-1])*np.pi/180
      if ("beta" in sp):
        beta=float(sp[len(sp)-1])*np.pi/180
      if ("gamma" in sp):
        gamma=float(sp[len(sp)-1])*np.pi/180

  Cell = pbc.make_cell(a, b, c, alpha, beta, gamma)
  Cell.C = pbc.make_cell_C(Cell)
  return [Cell, Cell.C]

def readFromColfile(filename,header=0,d=None):
  cl = 0
  v = list()
  with open(filename,"r") as fin:
    for k in range(header):
      next(fin)
    for line in fin:
      sp = line.split(d)
      if (cl==0):
        v = [[] for x in range(len(sp))]
      for p in range(len(sp)):
        v[p].append(float(sp[p]))
      cl = cl + 1
  return v


