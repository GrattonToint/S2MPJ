#  Test the setup and evaluation of Matlab versions of problems
#
#  Programming : S. Gratton and Ph. Toint
#  This version: 21 VI 2024
#

import numpy as np
import time
import sys

sys.path.insert(1, './python_problems')

start = 1
stop  = 10000

second_evaluation = True
eval_matvec       = True
printres          = True

start = start - 1

if printres:
    fidres = open("test_python.res", "a")

file1 = open("fullproblist", 'r')

Lines = file1.readlines()

problem = np.array([])
invalue = np.array([])
for line in Lines:
    fields = line.split()
    problem = np.append( problem,  fields[0] )
    invalue = np.append( invalue,  fields[1] )

#problem = [ "TAX1" ]
#invalue = [ "1,3,3,2,2" ]
    
for i in range(start,min(stop,len(problem))):
    print(" i= ",i)
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
    nm1 = PB.pb.n - 1
    if printres:
        print( "=== Checking problem  %s :  %s ( n = %d  m = %d )"%
               (str(i+1), theproblem, PB.pb.n, PB.pb.m ), file = fidres )
        for j in range(0,PB.pb.n,5):
             print( "x0 = %+18.15e" % ( PB.pb.x0[j] ),  end = ' ', file = fidres )
             if  j+1 <= nm1:
                 print( "%+18.15e" % ( PB.pb.x0[j+1] ), end = ' ', file = fidres )
             if  j+2 <= nm1:
                 print( "%+18.15e" % ( PB.pb.x0[j+2] ), end = ' ', file = fidres )
             if  j+3 <= nm1:
                 print( "%+18.15e" % ( PB.pb.x0[j+3] ), end = ' ', file = fidres )
             if  j+4 <= nm1:
                 print( "%+18.15e" % ( PB.pb.x0[j+4] ), end = ' ', file = fidres )
             print( "", file = fidres )

#    print( "=========")
#    print(PB.pb)
#    print( "=========")
#    print(PB.pbm)

#    exit()#D

#    if hasattr( PB.pbm, "A" ):
#        print("A = ", PB.pbm.A )
#    if hasattr( PB.pbm, "H" ):
#        print("H = ", PB.pbm.H )
#    print( "=========")

    ifixed  = find(PB.pb.xupper-PB.pb.xlower, lambda x:x==0 )
    if len( ifixed ) > 0:
        PB.pb.x0[ifixed] = PB.pb.xlower[ifixed]
        
    v = np.ones((PB.pb.n,1))
    if PB.pb.m:
        y = np.ones((PB.pb.m,1))
    else:
        y = np.array([])
        PB.pbm.congrps = np.array([])
    
    
#    fx,gx,Hx = PB.fgHx( PB.pb.x0 )
#    print( "Hessian = ", Hx)
#    exit()#D

    tic = time.perf_counter()
    Lxy, Lgxy, LHxy = PB.LgHxy( PB.pb.x0, y )
    
    toc = time.perf_counter()
    print(f"   Python evaluation  :  {toc - tic:0.4f} seconds")
    if printres:
        print( "L0 = %+18.15e" % ( Lxy ), file = fidres )
        for j in range(0,PB.pb.n,5):
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
        Lxy, Lgxy, LHxy = PB.LgHxy( PB.pb.x0, y )
        toc = time.perf_counter()
        print(f"   Python evaluation  :  {toc - tic:0.4f} seconds")


    if ( eval_matvec ):
        v     = np.ones((PB.pb.n,1))
        if len(ifixed) > 0:
            v[ifixed] = np.zeros((len(ifixed),1))
        tic   = time.perf_counter()
        LHxyv = PB.LHxyv( PB.pb.x0, y, v )
        toc   = time.perf_counter()
        print(f"   Python prod H*v    :  {toc - tic:0.4f} seconds")
        if printres:
            for j in range(0,PB.pb.n,5):
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
