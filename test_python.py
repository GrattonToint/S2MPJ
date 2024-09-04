#  Test the setup and evaluation of Python versions of problems
#
#  Programming : S. Gratton and Ph. Toint
#  This version: 2 IX 2024
#

import numpy as np
import time
import sys

sys.path.insert(1, './python_problems')

start = 1
stop  = 10000

second_evaluation = True
eval_matvec       = True
printtimes        = True # print timings in test_python.times
printres          = True  # print values of x0, L0, g0, and Hv in test_python.data

start = start - 1

if printres:
    fidres   = open("test_python.data", "a")
if printtimes:
    fidtimes = open("test_python.times", "a")
    
file1 = open("fullproblist", 'r')

Lines = file1.readlines()

problem = np.array([])
invalue = np.array([])
for line in Lines:
    fields  = line.split()
    problem = np.append( problem,  fields[0] )
    invalue = np.append( invalue,  fields[1] )

#problem = [ "TAX1" ]
#invalue = [ "1,3,3,2,2" ]
    
for i in range(start,min(stop,len(problem))):
    theproblem = problem[i]
    thevalue   = invalue[i]
    if ( theproblem[0] == "#" ):
        continue
    
    print("=== Checking problem "+str(i+1)+": "+theproblem)
    if printtimes:
       print("=== Checking problem "+str(i+1)+": "+theproblem, file = fidtimes )
    exec("from "+theproblem+" import *")
    tic = time.perf_counter()
    exec("PB = "+theproblem+"("+thevalue+")")
    toc = time.perf_counter()
    print(f"   Python setup       :  {toc - tic:0.4f} seconds")
    if printtimes:
        print(f"   Python setup       :  {toc - tic:0.4f} seconds", file=fidtimes)
    nm1 = PB.n - 1
    if printres:
        print( "=== Checking problem  %s :  %s ( n = %d  m = %d )"%
               (str(i+1), theproblem, PB.n, PB.m ), file = fidres )
        for j in range(0,PB.n,5):
             print( "x0 = %+18.15e" % ( PB.x0[j] ),  end = ' ', file = fidres )
             if  j+1 <= nm1:
                 print( "%+18.15e" % ( PB.x0[j+1] ), end = ' ', file = fidres )
             if  j+2 <= nm1:
                 print( "%+18.15e" % ( PB.x0[j+2] ), end = ' ', file = fidres )
             if  j+3 <= nm1:
                 print( "%+18.15e" % ( PB.x0[j+3] ), end = ' ', file = fidres )
             if  j+4 <= nm1:
                 print( "%+18.15e" % ( PB.x0[j+4] ), end = ' ', file = fidres )
             print( "", file = fidres )

    ifixed  = find(PB.xupper-PB.xlower, lambda x:x==0 )
    if len( ifixed ) > 0:
        PB.x0[ifixed] = PB.xlower[ifixed]
        
    v = np.ones((PB.n,1))
    if PB.m:
        y = np.ones((PB.m,1))
    else:
        y = np.array([])
        PB.congrps = np.array([])
    
    
#    fx,gx,Hx = PB.fgHx( PB.x0 )
#    print( "Hessian = ", Hx)
#    exit()#D

    tic = time.perf_counter()
    Lxy, Lgxy, LHxy = PB.LgHxy( PB.x0, y )
    
    toc = time.perf_counter()
    print(f"   Python evaluation  :  {toc - tic:0.4f} seconds")
    if printtimes:
        print(f"   Python evaluation  :  {toc - tic:0.4f} seconds", file = fidtimes )
    if printres:
        print( "L0 = %+18.15e" % ( Lxy ), file = fidres )
        for j in range(0,PB.n,5):
             print( "g0 = %+18.15e" % ( Lgxy[j] ),  end = ' ', file = fidres )
             if  j+1 <= nm1:
                 print( "%+18.15e" % ( Lgxy[j+1] ), end = ' ', file = fidres );
             if  j+2 <= nm1:
                 print( "%+18.15e" % ( Lgxy[j+2] ), end = ' ', file = fidres );
             if  j+3 <= nm1:
                 print( "%+18.15e" % ( Lgxy[j+3] ), end = ' ', file = fidres );
             if  j+4 <= nm1:
                 print( "%+18.15e" % ( Lgxy[j+4] ), end = ' ', file = fidres );
             print( "", file = fidres );

#    print("Lxy = ", Lxy)
#    print("Lgxy = ", Lgxy )
#    print("LgHxy =", LgHxy.toarray() )

    if ( second_evaluation ):
        tic = time.perf_counter()
        Lxy, Lgxy, LHxy = PB.LgHxy( PB.x0, y )
        toc = time.perf_counter()
        print(f"   Python evaluation  :  {toc - tic:0.4f} seconds")
        if printtimes:
           print(f"   Python evaluation  :  {toc - tic:0.4f} seconds", file = fidtimes )


    if ( eval_matvec ):
        v     = np.ones((PB.n,1))
        if len(ifixed) > 0:
            v[ifixed] = np.zeros((len(ifixed),1))
        tic   = time.perf_counter()
        LHxyv = PB.LHxyv( PB.x0, y, v )
        toc   = time.perf_counter()
        print(f"   Python prod H*v    :  {toc - tic:0.4f} seconds")
        if printtimes:
            print(f"   Python prod H*v    :  {toc - tic:0.4f} seconds", file = fidtimes )
        if printres:
            for j in range(0,PB.n,5):
                 print( "Hv = %+18.15e" % ( LHxyv[j] ),  end = ' ', file = fidres )
                 if  j+1 <= nm1:
                     print( "%+18.15e" % ( LHxyv[j+1] ), end = ' ', file = fidres );
                 if  j+2 <= nm1:
                     print( "%+18.15e" % ( LHxyv[j+2] ), end = ' ', file = fidres );
                 if  j+3 <= nm1:
                     print( "%+18.15e" % ( LHxyv[j+3] ), end = ' ', file = fidres );
                 if  j+4 <= nm1:
                     print( "%+18.15e" % ( LHxyv[j+4] ), end = ' ', file = fidres );
                 print( "", file = fidres );

#    print("LHxyv =", LHxyv )

if printres:
    fidres.close()
file1.close()
