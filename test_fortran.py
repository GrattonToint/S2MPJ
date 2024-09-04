#  Compare problems in Fortran (as output from Pycutest) and Python
#  (to be run in the docker)
#
#  Programming: S. Gratton and Ph. Toint
#  This version: 29 V 2024
#

import numpy as np
import pycutest
import time
from numpy import linalg as LA
import sys
import subprocess

###############################################################################################

start = 249
stop  = 10000

use_frt_x0  = True
printres    = True

if printres:
    fidres = open("test_fortran.data", "a")
    
#
#    ===============  Get the list of problems ==================
#

start    = start-1  # to cope with Python indeces starting at 0

file1    = open('fullproblist', 'r')  # the list of problems

Lines    = file1.readlines()
problems = np.array([])
invalue  = np.array([])
for line in Lines:
    fields   = line.split()
    problems = np.append( problems, fields[0] )
    invalue  = np.append( invalue, fields[1] )
fid   = ' '

#problems = [ "LUKVLI4" ]
#invalue  = [ "10" ]

#
#    ===============  Loop over all  problems ==================
#

for i in range(start,min(stop,len(problems))):

#
#   Get the problem name, copy the relevant SIF fir to the docker mastsif and remove
#   related files from the caches
#
    
    name  = problems[i]
    oname = name
    if name[0] == "#":
        continue
    if name[0] == "n":
        name = name[ 1:len(name) ]
    inval = invalue[i]
    print( "=== Checking problem ", i+1, ": ", oname, "  (", inval, ")", fid )

    sifname = name+".SIF"
    subprocess.run( [ "cp", sifname, "/cutest/mastsif" ] )
    try:
        pycutest.clear_cache( name )
    except:
        pass

    start_time = time.time()
    p = pycutest.import_problem( name )
    end_time = time.time()
    duration = end_time - start_time
    print(f"   Fortran setup        : {duration:.2e} secondes.")
    m_frt = p.m
    n_frt = p.n
    v     = np.ones((n_frt,))
    if m_frt == 0:
        start_time = time.time()
        f_frt, g_frt = p.obj(p.x0, gradient=True)  # Function and gradient
        end_time = time.time()
        duration = end_time - start_time
        print(f"   Fortran evaluation   : {duration:.2e} secondes.")
        g_frt = np.array(g_frt).reshape(-1,1)
        L_frt = f_frt
        
        start_time = time.time()
        Hv    = p.hprod(v,x=p.x0.reshape((n_frt,)))
        end_time = time.time()
        duration = end_time - start_time
        print(f"   Fortran prod H*v     : {duration:.2e} secondes.")
    else:
        y     = np.ones((m_frt,1))
        start_time = time.time()
        f_frt = p.obj(p.x0)   # Function
        c_frt = p.cons(p.x0)  # Constraints
        L_frt = (f_frt + y.T.dot(c_frt))[0]
        C     = p.lagjac(p.x0)
        g_frt = C[1].T.dot(y)+np.array(C[0]).reshape(-1,1)
        end_time = time.time()
        duration = end_time - start_time
        print(f"   Fortran evaluation   : {duration:.2e} secondes.")
#        g_new = p.lag(x,y)   # does not work ("no lag method") ???

        y  = np.ones((m_frt,))
        start_time = time.time()
        Hv = p.hprod(v,x=p.x0.reshape((n_frt,)), v=y)
        end_time = time.time()
        duration = end_time - start_time
        print(f"   Fortran prod H*v     : {duration:.2e} secondes.")

    if printres:
        print( "=== Checking problem  %s :  %s ( n = %d  m = %d )"%
             (str(i), oname, p.n, p.m ), file = fidres )
        for j in range(0,p.n,5):
            print( "x0 = %+18.15e" % ( p.x0[j] ),  end = ' ', file = fidres )
            if  j+1 < n_frt:
                print( "%+18.15e" % ( p.x0[j+1] ), end = ' ', file = fidres )
            if  j+2 < n_frt:
                print( "%+18.15e" % ( p.x0[j+2] ), end = ' ', file = fidres )
            if  j+3 < n_frt:
                print( "%+18.15e" % ( p.x0[j+3] ), end = ' ', file = fidres )
            if  j+4 < n_frt:
                print( "%+18.15e" % ( p.x0[j+4] ), end = ' ', file = fidres )
            print( "", file = fidres )
        print( "L0 = %+18.15e" % ( L_frt ), file = fidres )
        for j in range(0,p.n,5):
            print( "g0 = %+18.15e" % ( g_frt[j,0] ),  end = ' ', file = fidres )
            if j+1 < n_frt:
                print( "%+18.15e" % ( g_frt[j+1,0] ), end = ' ', file = fidres )
            if j+2 < n_frt:
                print( "%+18.15e" % ( g_frt[j+2,0] ), end = ' ', file = fidres )
            if j+3 < n_frt:
                print( "%+18.15e" % ( g_frt[j+3,0] ), end = ' ', file = fidres )
            if j+4 < n_frt:
                print( "%+18.15e" % ( g_frt[j+4,0] ), end = ' ', file = fidres )
            print( "", file = fidres );
        for j in range(0,p.n,5):
             print( "Hv = %+18.15e" % ( Hv[j] ),  end = ' ', file = fidres )
             if  j+1 < n_frt:
                 print( "%+18.15e" % ( Hv[j+1] ), end = ' ', file = fidres )
             if  j+2 < n_frt:
                 print( "%+18.15e" % ( Hv[j+2] ), end = ' ', file = fidres )
             if  j+3 < n_frt:
                 print( "%+18.15e" % ( Hv[j+3] ), end = ' ', file = fidres )
             if  j+4 < n_frt:
                 print( "%+18.15e" % ( Hv[j+4] ), end = ' ', file = fidres )
             print( "", file = fidres );               
            
if printres:
    fidres.close()
file1.close()


