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

start = 1
stop  = 10000

use_fortran = True
use_python  = True
use_frt_x0  = True

#printfg = True
printfg = False

#
#    ===============  Get the list of problems ==================
#

start   = start-1  # to cope with Python indeces starting at 0

file1   = open('fullproblist', 'r')  # the list of problems

Lines   = file1.readlines()
problems = np.array([])
invalue = np.array([])
for line in Lines:
    fields  = line.split()
    problems = np.append( problems, fields[0] )
    invalue = np.append( invalue, fields[1] )
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
    if name[0] == "#":
        continue
    inval = invalue[i]
    print( "=== Checking problem ", i+1, ": ", name, "  (", inval, ")", fid )

    sifname = name+".SIF"
    subprocess.run( [ "cp", sifname, "/cutest/mastsif" ] )
    try:
        pycutest.clear_cache( name )
    except:
        pass
    try:
       cachefile = "__pycache__/"+name+".cpython-311.pyc"
       subprocess.run( [ "rm", cachefile ] )
    except:
        pass
#
#       ======================  FORTRAN =======================
#
    if use_fortran:
        frt_ok = True
        try:
           p = pycutest.import_problem( name )
           m_frt = p.m
           n_frt = p.n
#           print("n_frt = %i , m_frt = %i " % (n_frt, m_frt ) )
           start_time = time.time()
           if m_frt == 0:
               f_frt, g_frt = p.obj(p.x0, gradient=True)  # Function and gradient
               g_frt = np.array(g_frt).reshape(-1,1)
           else:
               y     = np.ones((m_frt,1))
               f_frt = p.obj(p.x0)   # Function
               c_frt = p.cons(p.x0)  # Constraints
               L_frt = (f_frt + y.T.dot(c_frt))[0]
               C     = p.lagjac(p.x0)
               g_frt = C[1].T.dot(y)+np.array(C[0]).reshape(-1,1)
           end_time = time.time()
           duration = end_time - start_time
           print(f"   Fortran evaluation   : {duration:.2e} secondes.")
           if printfg:
              print('x0_frt = ', p.x0)
              print("f_frt  = ",f_frt)
              if p.m > 0:
                  print("c_frt  = ",c_frt.T )
                  print("c_frt[0] = ", c_frt[0] )
                  print("L_frt = ", L_frt )
              print("g_frt  = ",g_frt.T)
        except:
            frt_ok = False
            print('pb in fortran', fid )
#
#      ====================== PYTHON =======================
#

    if use_python:
        python_ok = True
        try:
            import_command = f"from {name} import *"
            eval( compile(import_command, '<string>', 'exec') )
            exec('PB = '+name+'('+inval+')')
            m = PB.pb.m
            n = PB.pb.n
            active = np.array(~(PB.pb.xlower == PB.pb.xupper) )
            fixed  = np.array(  PB.pb.xlower == PB.pb.xupper)
            new_n  = active.sum()
            dimOK = True
            if use_fortran:
               if new_n < PB.pb.n:
                   dimOK = False
                   print("n_frt = ", n_frt, ", n_py = ", PB.pb.n, "--> new_n = ", new_n )
               if not((m == m_frt) and (new_n == n_frt)):
                    print( "===> Nonmatching sizes ", fid )
                    raise Exception("Size pb")
            start_time = time.time()
            if use_frt_x0 and use_fortran and dimOK:
                myx0 = p.x0.reshape(-1,1)
            else:
                myx0 = PB.pb.x0
            if new_n < PB.pb.n:
               myx0[fixed] = PB.pb.xlower[fixed]
            if m == 0:
                f, g = PB.fgx(myx0)
            else:
                y = np.ones((m,1))
                f = PB.fx(myx0)
                c = PB.cx(myx0)
                L, g = PB.Lgxy(myx0, y)
            end_time = time.time()
            duration = end_time - start_time
            print(f"   Python evaluation    : {duration:.2e} secondes.")
            gactive =  np.array(g[active]).reshape(-1,1)
            if printfg:
                print("x0 = ", PB.pb.x0.T )
                print("f  = ",f)
                if PB.pb.m > 0:
                    print("c  = ",c.T )
                    print("L = ", L )
                print("g  = ",gactive.T )
        except:
            python_ok = False
            print('pb in python',  fid )
#
#           ================================================================
#

    if use_fortran and use_python:
        if python_ok and frt_ok :
            diffx0 = LA.norm(p.x0-myx0[active])
            print("   Error in x0          : absolute =  %1.6e" %diffx0,   " relative   = %1.6e " %(LA.norm(p.x0-myx0[active]) /(1e-15+LA.norm(p.x0.reshape(-1, 1)))), fid )
            if printfg:
                print("   Error in function    = %1.6e " %abs((f_frt-float(f))/(1e-12+abs(f_frt))), fid )
                if p.m:
                    print("   Error in Lagrangian = %1.6e " %abs((L_frt-float(L))/(1e-12+abs(L_frt))), fid )
                print("difference in grad  = ", (g_frt-gactive).T )
            diffnorm = LA.norm((g_frt-gactive))
            print("   Error in Lag gradient: absolute =  %1.6e" %diffnorm, " relative   = %1.6e " %( diffnorm /(1e-15+LA.norm(g_frt))), fid)
        else:
            print("Comparison not possible ", fid )
            print("n_frt = %i , n = %i " % (n_frt,n ) )
            print("m_frt = %i , m = %i " % (m_frt,m ) )


sys.stdout = sys.__stdout__

exit()
