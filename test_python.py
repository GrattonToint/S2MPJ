#  Test the setup and evaluation of Matlab versions of problems
#
#  Programming : S. Gratton and Ph. Toint
#  This version: 15 V 2024
#

import numpy as np
import time
import sys

sys.path.insert(1, './python_problems')

start = 0
stop  = 10000

second_evaluation = True
eval_matvec       = True

file1 = open('fullproblist', 'r')

Lines = file1.readlines()

problem = np.array([])
invalue = np.array([])
for line in Lines:
    fields = line.split()
    problem = np.append( problem,  fields[0] )
    invalue = np.append( invalue,  fields[1] )

problem = [ "TAX1" ]
invalue = [ "1,3,3,2,2" ]
    
for i in range(start,min(stop,len(problem))):
    theproblem = problem[i]
    thevalue   = invalue[i]
    if ( theproblem[0] == "#" ):
       continue
    
    print("=== Checking problem "+str(i+1)+": "+theproblem)
    exec("from "+theproblem+" import *")
    tic = time.perf_counter()
    exec("PB = "+theproblem+"("+thevalue+")")
    toc = time.perf_counter()
    print(f"   Python setup       :  {toc - tic:0.4f} seconds")

    print( "=========")
    print(PB.pb)
    print( "=========")
    print(PB.pbm)

#    exit()#D

#    if hasattr( PB.pbm, "A" ):
#        print("A = ", PB.pbm.A )
#    if hasattr( PB.pbm, "H" ):
#        print("H = ", PB.pbm.H )
#    print( "=========")

    v = 1.5*np.ones((PB.pb.n,1))
    if PB.pb.m:
        y = np.ones((PB.pb.m,1))
    else:
        y = np.array([])
        PB.pbm.congrps = np.array([])
    
    tic = time.perf_counter()
    fx,gx,Hx = PB.fgHx( PB.pb.x0 )
    print( "Hessian = ", Hx)
    exit()#D
    Lxy, Lgxy, LHxy = PB.LgHxy( PB.pb.x0, y )
    toc = time.perf_counter()
    print(f"   Python evaluation  :  {toc - tic:0.4f} seconds")

    print("Lxy = ", Lxy)
#    print("Lgxy = ", Lgxy )
#    print("LgHxy =", LgHxy.toarray() )

    if ( second_evaluation ):
        tic = time.perf_counter()
        Lxy, Lgxy, LHxy = PB.LgHxy( PB.pb.x0, y )
        toc = time.perf_counter()
        print(f"   Python evaluation  :  {toc - tic:0.4f} seconds")


    if ( eval_matvec ):
        v     = np.ones((PB.pb.n,1))
        tic   = time.perf_counter()
        LHxyv = PB.LHxyv( PB.pb.x0, y, v )
        toc   = time.perf_counter()
        print(f"   Python prod H*v    :  {toc - tic:0.4f} seconds")

#    print("LHxyv =", LHxyv )
