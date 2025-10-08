#
#  This version handles missing derivatives
#
####################################################################################################
#####################################################################################################
#
#                                S2MPJ library for Python
#
#   Performs the runtime actions specific to S2MPJ, irrespective of the problem at hand.
#   Also contains the problem selection tool.
#
#   Programming: S. Gratton and Ph. Toint (this version 14 VII 2025)
#
#####################################################################################################
#####################################################################################################

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse import csr_matrix
from pprint import pprint
import re
import os

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class structtype():
    
    def __str__(self):
        pprint(vars(self))
        return ''
    pass
    def __repr__(self):
        pprint(vars(self))
        return ''

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   The generic class for CUTEst problems, including evaluation methods
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class CUTEst_problem:

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    #   Extract the values of the global element's and group's parameters, if any.
    #
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def getglobs( self ):
        try:
            eval( 'self.'+'e_globs(self)' )
        except:
            pass
        try:
            eval( 'self.'+'g_globs(self)' )
        except:
            pass

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    #   Define the main evaluations actions for the problem.
    # 
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    def fx( self, x ):              # input = ( x )
        if ( hasattr( self, "objgrps" ) and len( self.objgrps) ) or hasattr( self, "H" ):
            x = x.reshape(-1,1)
            self.getglobs()
            return self.evalgrsum( True, self.objgrps, x, 1 )
        else:
            print( " ")
            print( "ERROR: problem "+self.name+" has no objective function!" )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " ")
        
    def fgx( self, x ):             # input = ( x )
        if ( hasattr( self, "objgrps" ) and len( self.objgrps) ) or hasattr( self, "H" ):
            x = x.reshape(-1,1)
            self.getglobs()
            return self.evalgrsum( True, self.objgrps, x, 2 )
        else:
            print( " " )
            print( "ERROR: problem "+self.name+" has no objective function!" )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " " )
        
    def fgHx( self, x ):            # input = ( x )
        if ( hasattr( self, "objgrps" ) and len( self.objgrps) ) or hasattr( self, "H" ):
            x = x.reshape(-1,1)
            self.getglobs()
            return self.evalgrsum( True, self.objgrps, x, 3 )
        else:
            print( " " )
            print( "ERROR: problem "+self.name+" has no objective function!" )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " " )
        
    def fHxv( self, x, v ):          # input = ( x, v )
        if ( hasattr( self, "objgrps" ) and len( self.objgrps) ) or hasattr( self, "H" ):
            x = x.reshape(-1,1)
            v = v.reshape(-1,1)
            self.getglobs()
            return self.evalHJv( "Hv", self.objgrps, x, v, [] )
        else:
            print( " " )
            print( "ERROR: problem "+self.name+" has no objective function!" )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " " )
        
    def cx( self, x ):               # input = ( x )
        if hasattr( self, "congrps" ) and len( self.congrps):
            x = x.reshape(-1,1)
            self.getglobs()
            return self.evalgrsum( False, self.congrps, x, 1 )
        else:
            print( " " )
            print( 'ERROR: problem '+self.name+' has no constraint!' )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " " )
        
    def cJx( self, x ):               # input = ( x )
        if hasattr( self, "congrps" ) and len( self.congrps):
            x = x.reshape(-1,1)
            self.getglobs()
            return self.evalgrsum( False, self.congrps, x, 2 )
        else:
            print( " " )
            print( 'ERROR: problem '+self.name+' has no constraint!' )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " " )
        
    def cJHx( self, x ):             # input = ( x )
        if hasattr( self, "congrps" ) and len( self.congrps):
            x = x.reshape(-1,1)
            self.getglobs()
            return self.evalgrsum( False, self.congrps, x, 3 )
        else:
            print( " ")
            print( 'ERROR: problem '+self.name+' has no constraint!' )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " " )

    def cJxv( self, x, v ):          # input = ( x, v )
        if hasattr( self, "congrps" ) and len( self.congrps):
            x = x.reshape(-1,1)
            v = v.reshape(-1,1)
            self.getglobs()
            return self.evalHJv( "Jv", self.congrps, x, v, [] )
        else:
            print( " " )
            print( 'ERROR: problem '+self.name+' has no constraint!' )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " ")
        
    def cJtxv( self, x, v ):          # input = ( x, v )
        if hasattr( self, "congrps" ) and len( self.congrps):
            x = x.reshape(-1,1)
            v = v.reshape(-1,1)
            self.getglobs()
            return self.evalHJv( "Jtv", self.congrps, x, v, [] )
        else:
            print( " " )
            print( 'ERROR: problem '+self.name+' has no constraint!' )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " ")
        
    def cIx( self, x, clist ):       # input = ( x, clist )
        if hasattr( self, "congrps" ) and len( self.congrps) and len( clist ):
            x      = x.reshape(-1,1)
            iclist = [ self.congrps[i] for i in clist ]
            self.getglobs()
            return self.evalgrsum( False, iclist , x, 1 )
        else:
            print( " " )
            print( 'ERROR: problem '+self.name+' has no constraint!' )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " ")
        
    def cIJx( self, x, clist ):          # input = ( x , clist )
        if hasattr( self, "congrps" ) and len( self.congrps) > 0 and len( clist ):
            x      = x.reshape(-1,1)
            iclist = [ self.congrps[i] for i in clist ]
            self.getglobs()
            return self.evalgrsum( False, iclist, x, 2 )
        else:
            print( " " )
            print( 'ERROR: problem '+self.name+' has no constraint!' )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " ")
        
    def cIJHx( self, x , clist ):        # input = ( x, clist )
        if hasattr( self, "congrps" ) and len( self.congrps) > 0 and len( clist ):
            x      = x.reshape(-1,1)
            iclist = [ self.congrps[i] for i in clist ]
            self.getglobs()
            return self.evalgrsum( False, iclist, x, 3 )
        else:
            print( " " )
            print( 'ERROR: problem '+self.name+' has no constraint!' )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " ")
        
    def cIJxv( self, x, v, clist ):      # input = ( x, v, clist )
        if hasattr( self, "congrps" ) and len( self.congrps) > 0 and len( clist ):
            x      = x.reshape(-1,1)
            v      = v.reshape(-1,1)
            iclist = [ self.congrps[i] for i in clist ]
            self.getglobs()
            return self.evalHJv( "Jv", iclist, x, v, [] )
        else:
            print( " " )
            print( 'ERROR: problem '+self.name+' has no constraint!' )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " ")
        
    def cIJtxv( self, x, v, clist ):      # input = ( x, v, clist )
        if hasattr( self, "congrps" ) and len( self.congrps) > 0 and len( clist ):
            x      = x.reshape(-1,1)
            v      = v.reshape(-1,1)
            iclist = [ self.congrps[i] for i in clist ]
            self.getglobs()
            return self.evalHJv( "Jvt", iclist, x, v, [] )
        else:
            print( " " )
            print( 'ERROR: problem '+self.name+' has no constraint!' )
            print( "       Please refer to the problem classification for checking a problem's type." )
            print( " ")
        
    def Lxy( self, x, y ):           # input = ( x, y )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self, "congrps" ):
            return self.evalLx( self.objgrps, self.congrps, x, y, 1 )
        else:
            return self.evalgrsum( True, self.objgrps, x, 1 )
        
    def Lgxy( self, x, y ):          # input = ( x, y )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self, "congrps" ):
            return self.evalLx( self.objgrps, self.congrps, x, y, 2 )
        else:
            return self.evalgrsum( True, self.objgrps, x, 2 )
        
    def LgHxy( self, x, y ):         # input = ( x, y )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self, "congrps" ):
            return self.evalLx( self.objgrps, self.congrps, x, y, 3 )
        else:
            return self.evalgrsum( True, self.objgrps, x, 3 )
        
    def LHxyv( self, x, y, v ):      # input = ( x, y, v )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        v = v.reshape(-1,1)
        self.getglobs()
        if hasattr( self, "congrps" ):
            return self.evalLHxyv( self.objgrps, self.congrps, x, y, v )
        else:
            print( "Hv" )
            return self.evalHJv( "Hv", self.objgrps, x, v, [] )
        
    def LIxy( self, x, y, clist ):       # input = ( x, y, clist )
        x  = x.reshape(-1,1)
        y  = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self, "congrps" ):
            return self.evalLx( self.objgrps, [ self.congrps[i] for i in clist ], x, y, 1 )
        else:
            return self.evalgrsum( True, self.objgrps, x, 1 )
        
    def LIgxy( self, x, y, clist ):      # input = ( x, y, clist )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self, "congrps" ):
            return self.evalLx( self.objgrps, [ self.congrps[i] for i in clist ], x, y, 2 )
        else:
            return self.evalgrsum( True, self.objgrps, x, 2 )
        
    def LIgHxy( self, x, y, clist  ):    # input = ( x, y, clist )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self, "congrps" ):
            return self.evalLx( self.objgrps, [ self.congrps[i] for i in clist ], x, y, 3 )
        else:
            return self.evalgrsum( True, self.objgrps, x, 3 )
        
    def LIHxyv( self, x, y, v, I ):  # input = ( x, y, v, clist )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        v = v.reshape(-1,1)
        self.getglobs()
        if hasattr( self, "congrps" ):
            return self.evalLHxyv( self.objgrps, [ self.congrps[i] for i in clist ], x, y, v )
        else:
            return self.evalHJv( "Hv", self.objgrps, x, v, [] )


    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    #   Evaluate the value of a sum of groups (and, if requested, that of of its gradient and Hessian)
    #   at x, given the problem data available in the self struct.
    #
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def evalgrsum( self, isobj, glist, x, nargout ):

        debug = False#D
        
        #   Initializations
        
        n   = len( x )
        m   = len( glist )
        if isobj:
            fx = 0.0
        else:
            cx   = np.zeros(( m, 1 ))
            ic   = -1
            if hasattr( self, "conderlvl" ):
               lder = len( self.conderlvl )

        if nargout > 1:
            if isobj:
                if hasattr( self, "objderlvl" ):
                    if self.objderlvl >= 1:
                        gx = np.zeros(( n, 1 ))
                    else:
                        gx = np.zeros(( n, 1 ))
                        gx[0]  = np.nan
                else:
                    gx = np.zeros(( n, 1 ))
                    
            else:
                if hasattr( self, "conderlvl" ):
                   if any( x >=1 for x in  self.conderlvl ):
                       Jx = lil_matrix(( m, n ))
                   else:
                       Jx = lil_matrix(( m, n ))
                       Jx[0,0] = np.nan
                else:
                    Jx = lil_matrix(( m, n ))
                    
            if nargout > 2:
                if ( isobj ):
                    if hasattr( self, "objderlvl" ):
                       if self.objderlvl >= 2:
                          Hx = lil_matrix(( n, n ))
                       else:
                          Hx = lil_matrix(( n, n ))
                          Hx[0,0] = np.nan
                    else:
                        Hx = lil_matrix(( n, n ))
                else:
                    Hx = []

        #  Check for the presence and size of a linear term

        if hasattr( self, "A" ):
            sA1, sA2 = self.A.shape
            has_A = True
        else:
            has_A = False

        #  Evaluate the quadratic term, if any.

        if isobj and hasattr( self, "H" ):
            
            Htimesx = self.H .dot(x)
            if nargout ==  1:
                fx += 0.5 * x.T .dot(Htimesx)
            elif nargout == 2:
                gx += Htimesx
                fx += 0.5 * x.T .dot(Htimesx)
            elif nargout == 3:
                Htimesx = self.H .dot(x);
                gx += Htimesx
                fx += 0.5 * x.T .dot(Htimesx)
                Hx += self.H
                
        if debug: #D
            if isobj:
                print( "fx(quadratic) = ", fx )
                #print( "gx = ", gx )
            else:
                print( "cx(quadratic) = ", cx )
        
        #    Loop on the groups list 

        for iig in range( len( glist )):
            ig = int( glist[ iig ] )
            

            #  Find the level of available derivatives for the group.

            if  isobj:
                if hasattr( self, "objderlvl" ):
                   derlvl = self.objderlvl
                else:
                   derlvl = 2
            else:
                if hasattr( self, "conderlvl" ):
                    if lder == 1:
                        derlvl = self.conderlvl[ 0 ]
                    else:
                        derlvl = self.conderlvl[ np.where( self.congrps == ig )[0][0] ]
                else:
                    derlvl = 2
            nout = min( nargout, derlvl + 1 );

            #  Find the group's scaling.
            
            if hasattr(self,"gscale"):
                if ig < len(self.gscale) and not self.gscale[ig] is None and abs( self.gscale[ig] ) > 1.0e-15:
                    gsc = self.gscale[ig]
                else:
                    gsc = 1.0
            else:
                gsc = 1.0

            #  Evaluate the linear term, if any.

            if hasattr(self,"gconst") and ig < len(self.gconst) and not self.gconst[ig] is None:
                fin = float(-self.gconst[ig])
            else:
                fin = 0
            if has_A and ig < sA1:
                gin           = np.zeros( (n, 1) )
                gin[:sA2, :1] = self.A[ ig, :sA2 ].T.toarray()
                fin           = float( fin + gin.T .dot(x) )
            elif nargout >= 2:
                gin =  np.zeros(( n, 1 ))

            if nargout > 2:
                Hin = lil_matrix(( n, n ))

            if debug:
                print( "ig = ", ig, "  fin(linear)", fin )

            if hasattr( self, "grelt" ) and ig < len( self.grelt ) and not self.grelt[ ig ] is None:
                for iiel in range(len( self.grelt[ ig ] ) ):  #  loop on elements
                    iel    = self.grelt[ ig ][ iiel ]         #  the element's index
                    efname = self.elftype[ iel ]              #  the element's ftype
                    irange = [iv for iv in self.elvar[ iel ]] #  the elemental variable's indeces 
                    xiel   = x[ np.array(irange) ]            #  the elemental variable's values

                    if  hasattr( self, 'grelw' ) and ig <= len( self.grelw ) and not self.grelw[ig] is None :
                        has_weights = True;
                        wiel        = self.grelw[ ig ][ iiel ]
                    else:
                        has_weights = False

                    # Only the value is requested.
                    
                    if nout == 1:
                        fiel = eval('self.'+efname +'( self, 1, xiel, iel )')
                        if ( has_weights ):
                            fin += wiel * fiel
                        else:
                            fin += fiel
                        
                    #  The value and its gradient are requested.
                    
                    elif nout == 2:
                        fiel, giel = eval('self.'+efname +'( self, 2, xiel, iel)')
                        if  has_weights:
                            fin += wiel * fiel
                            for ir in range(len(irange)):
                                ii = irange[ ir ]
                                gin[ ii ] += wiel * giel[ ir ]
                        else:
                            fin = fin + fiel;
                            for ir in range(len(irange)):
                                ii = irange[ ir ]
                                gin[ ii ] += giel[ ir ]

                    elif nout == 3:
                        fiel, giel, Hiel = eval('self.'+efname +'( self, 3, xiel, iel )')
                        if has_weights:
                            fin += wiel * fiel
                            for ir in range(len(irange)):
                                ii = irange[ ir ]
                                gin[ ii ] += wiel * giel[ ir ]
                                for jr in range(len( irange )):
                                    jj  = irange[ jr ]
                                    Hin[ ii, jj ] += wiel * Hiel[ ir, jr ]
                        else:
                            fin = fin + fiel;
                            for ir in range(len(irange)):
                                ii = irange[ ir ]
                                gin[ ii ] += giel[ ir ]
                                for jr in range(len( irange )):
                                    jj  = irange[ jr ]
                                    Hin[ ii, jj ] += Hiel[ ir, jr ]
            if debug: #D
                if isobj:
                    print( "ig = ", ig, "  fin(nonlinear) = ", fin, "  fx = ", fx  )
                    #print( "gx = ", gx )
                else:
                    print( "ig = ", ig, "  fin(nonlinear) = ", fin, "  cx = ", cx  )
                    
            #  Evaluate the group function.
            
            #  1) the non-TRIVIAL case

            if hasattr(self,"grftype") and ig < len( self.grftype ) and not self.grftype[ ig ] is None:
                egname = self.grftype[ ig ]
            else:
                egname = "TRIVIAL"
            if egname!='TRIVIAL' and egname is not None:
                if isobj:
                    if nargout == 1:
                        fx += eval('self.'+egname+'( self, 1, fin, ig )') / gsc
                    elif nargout == 2:
                        [ fa, grada ] = eval('self.'+ egname+'( self, 2, fin, ig )')
                        fx += fa / gsc
                        if derlvl >= 1:
                            gx += grada * gin / gsc
                        else:
                            gx = np.nan * np.ones(( n, 1 ))
                    elif nargout == 3:
                        [ fa, grada, Hessa ] = eval('self.'+egname+'( self, 3, fin, ig )')
                        fx   += fa / gsc
                        if derlvl >= 1:
                            gx += grada * gin / gsc
                        else:
                            gx = np.nan * np.ones(( n, 1 ))
                        if derlvl >= 2:
                            sgin  = lil_matrix(gin)
                            Hx   += (Hessa * sgin.dot(sgin.transpose())+ grada * Hin) / gsc
                        else:
                            Hx      = lil_matrix(( n, n ))
                            Hx[0,0] = np.nan
                else:
                    ic = ic + 1
                    if nargout == 1:
                        fa = eval('self.'+egname+'( self, 1, fin, ig )')
                        cx[ ic ] = fa / gsc
                    elif nargout == 2:
                        fa, grada = eval('self.'+egname+'( self, 2, fin, ig )')
                        cx[ ic ]  = fa / gsc
                        if derlvl >= 1:
                            sgin      = lil_matrix( gin )
                            Jx[ic,:]  = grada * sgin.T / gsc
                        else:
                            Jx[ic,:]  = np.nan*np.ones(( 1, n ))
                    elif nargout == 3:
                        fa, grada, Hessa = eval('self.'+egname+'( self, 3, fin, ig )') 
                        cx[ ic ] = fa / gsc
                        if derlvl >= 1:
                            sgin      = lil_matrix( gin )
                            Jx[ic,:]  = grada * sgin.T / gsc
                        else:
                            Jx[ic,:]  = np.nan*np.ones(( 1, n ))
                        if derlvl >= 2:
                            Hx.append( ( Hessa * sgin.dot( sgin.transpose() )+ grada * Hin ) / gsc )
                        else:
                            Hxi = lil_matrix(( n, n ))
                            Hxi[0,0] = np.nan
                            Hx.append( Hxi )
                            
            #  2) the TRIVIAL case: the group function is the identity
            
            else:
                if isobj:
                    if nargout == 1:
                        fx += fin / gsc 
                    if nargout == 2:
                        fx += fin / gsc
                        if derlvl >= 1:
                            gx += gin / gsc
                        else:
                            gx = np.nan * np.ones(( n, 1 ))
                    if nargout == 3:
                        fx += fin / gsc
                        if derlvl >= 1:
                            gx += gin / gsc
                        else:
                            gx = np.zeros(( n, 1 ))
                            gx[0] = np.nan
                        if derlvl >= 2:
                            Hx += Hin / gsc
                        else:
                            Hx = lil_matrix(( n, n ))
                            Hx[0,0] = np.nan
                else:
                    ic += + 1
                    if nargout == 1:
                        cx[ic] = fin / gsc
                    elif nargout == 2:
                        cx[ic] = fin / gsc
                        if derlvl >= 1:
                            Jx[ic,:] = gin.transpose() / gsc
                        else:
                            Jx[ic,:] = np.nan * np.ones(( 1, n ))
                    elif nargout == 3:
                        cx[ic] = fin / gsc
                        if derlvl >= 1:
                            Jx[ic,:] = gin.transpose() / gsc
                        else:
                            Jx[ic,:] = np.nan * np.ones(( 1, n ))
                        if derlvl >= 2:
                           Hx.append( Hin / gsc )
                        else:
                           Hin = lil_matrix(( n, n ))
                           Hin[0,0] = np.nan
                           Hx.append( Hin )

            if debug:
                if isobj:
                    print( "ig = ", ig, "  fx(final) = ", fx )#D
                    #print( "gx(final) = ", gx )
                    #print( "Hx(final) = ", Hx )
                else:
                    print( "ig = ", ig, "  cx(final) = ", cx )#D
        if isobj:
            if nargout == 1:
                return float( fx )
            elif nargout == 2:
                return float( fx ), gx.reshape(-1,1)
            elif nargout == 3:
                return float( fx ), gx.reshape(-1,1), Hx
        else:
            if nargout == 1:
                return cx
            elif nargout == 2:
                return cx, Jx
            elif nargout == 3:
                return cx, Jx, Hx

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    #     Depending on mode:
    #     mode = "Hv"  : evaluate the product of the objective's Hessian times v (glist = obj groups)
    #     mode = "HIv" : evaluate the product of the constraints' Hessian  times v time the multiplier y
    #                 (glist = cons groups)
    #     mode = "Jv"  : evaluate the product of the constraints' Jacobian times v (glist = cons groups)
    #     The vector y is unused (and unreferenced) for modes "Hv" and "Jv".
    #
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def evalHJv( self, mode, glist, x, v, y ):

        debug = False#D
        
        #   Initializations
        
        n   = len( x );
        if mode == "Hv":
            if hasattr( self, "objderlvl" ):
               if self.objderlvl < 2:
                  HJv      = np.zeros(( n, 1 ))
                  HJv[0,0] = np.nan
                  return HJv.reshape(-1,1)
               else:
                  HJv    = np.zeros((n,1))
                  derlvl = 2
            else:
                derlvl = 2
        elif mode == "HIv":
            if hasattr( self, "conderlvl" ):
               m     = len( glist )
               if any( x < 2 for x in  self.conderlvl ) :
                  HJv    = np.zeros(( n, 1 ))
                  HJv[0] = np.nan
                  return HJv.reshape(-1,1)
               else:
                  HJv    = np.zeros(( n, 1 ))
                  ic     = -1
                  derlvl = 2
            else:
                derlvl = 2
        else:
            m     = len( glist )
            if hasattr( self, "conderlvl" ):
                if any( x < 1 for x in  self.conderlvl ) :
                    if mode == "Jv":
                        HJv    = np.zeros(( m, 1 ))
                        HJv[0] = np.nan
                        return HJv.reshape(-1,1)
                    else:
                        HJv    = np.zeros(( n, 1 ))
                        HJv[0] = np.nan
                        return HJv.reshape(-1,1)
                else:
                    if mode == "Jv":
                        HJv  = np.zeros(( m, 1 ))
                        ic   = -1
                        lder = len( self.conderlvl )
                    else:
                        HJv  = np.zeros(( n, 1 ))
                        ic   = -1
                        lder = len( self.conderlvl )
            else:
                derlvl = 2

        #  Check for the presence and size of a linear term

        if hasattr( self, "A" ):
            sA1, sA2 = self.A.shape
            has_A = True
        else:
            has_A = False

        #  Evaluate the quadratic term, if any.

        if mode == "Hv" and hasattr( self, "H" ):
            HJv += self.H.dot(v);
        
        if debug: #D
            print( "HJv(quadratic) = ", HJv )

        for iig in range(len( glist )):
            ig = int( glist[ iig ] );

            #  Find the level of available derivatives for the group.

            if mode == "Jv" or mode == "Jtv":
                if hasattr( self, "conderlvl" ):
                    if lder == 1:
                        derlvl = self.conderlvl[0]
                    else:
                        derlvl = self.conderlvl[ np.where( self.congrps == ig )[0][0] ]
                else:
                    derlvl = 2

                #  Avoid computation for group ig if its first derivative is missing

                if derlvl < 1:
                    ic += 1
                    HJv[ic] = np.nan
                    continue
                    
            #  Find the group's scaling.
            
            if hasattr(self,"gscale"):
                if ig < len( self.gscale ) and not self.gscale[ ig ] is None and abs( self.gscale[ ig ] ) > 1.0e-15:
                    gsc = self.gscale[ ig ];
                else:
                    gsc = 1.0;
            else:
                gsc = 1.0;

            #  Evaluate the linear term, if any.

            if hasattr( self, "gconst" ) and ig < len( self.gconst ) and not self.gconst[ ig ] is None:
                fin = float(-self.gconst[ ig ] )
            else:
                fin = 0.0
            gin = np.zeros((n,1))
            if has_A and ig < sA1:
                gin[:sA2, :1] = self.A[ ig, :sA2 ].T.toarray()
                fin           = float(fin + gin.T .dot(x))

            if debug:
                print( "ig = ", ig, "  fin(linear)", fin ) #D

            Hinv = np.zeros((n,1))

            if hasattr( self, "grelt" ) and ig < len( self.grelt ) and not self.grelt[ ig ] is None:
                for iiel in range( len( self.grelt[ ig ] ) ):  #  loop on elements
                    iel    = self.grelt[ ig ][ iiel ]          #  the element's index
                    efname = self.elftype[ iel ];              #  the element's ftype
                    irange = [iv for iv in self.elvar[ iel ]]  #  the elemental variable's indeces
                    xiel   = x[ np.array(irange) ]             #  the elemental variable's values

                    if  hasattr( self, 'grelw' ) and ig <= len( self.grelw ) and not self.grelw[ig] is None :
                        has_weights = 1;
                        wiel        = self.grelw[ ig ][ iiel ]
                    else:
                        has_weights = 0;

                    #  The group is an objective group
                    
                    if mode == "Hv" or mode == "HIv":
                        fiel, giel, Hiel = eval('self.'+efname +'( self, 3, xiel, iel )')
                        if has_weights:
                            fin += wiel * fiel;
                            for ir in range( len( irange ) ):
                                ii = irange[ ir ]
                                gin[ ii ] = gin[ ii ] + wiel * giel[ ir ]
                                for jr in range(len( irange )):
                                    jj  = irange[ jr ]
                                    Hinv[ ii ] +=  wiel * Hiel[ ir, jr ] * v[ jj ]
                        else:
                            fin = fin + fiel;
                            for ir in range( len( irange ) ):
                                ii = irange[ ir ]
                                gin[ ii ] = gin[ ii ] + giel[ ir ]
                                for jr in range(len( irange )):
                                    jj  = irange[ jr ];
                                    Hinv[ ii ] +=  Hiel[ ir, jr ].dot( v[ jj ] )

                     #   The group is a constraint group.

                    elif derlvl >= 1:

                        fiel, giel = eval('self.'+efname +'( self, 2, xiel, iel)')
                        if has_weights:
                            fin = fin + wiel * fiel;
                            for ir in range( len( irange ) ):
                                ii = irange[ ir ]
                                gin[ ii ] += wiel * giel[ ir ]
                        else:
                            fin = fin + fiel;
                            for ir in range( len( irange ) ):
                                ii = irange[ ir ]
                                gin[ ii ] += giel[ ir ]
            
            if debug: #D
                print( "ig = ", ig, "  Hinv(nonlinear) = ", Hinv, "  HJv = ", HJv  )#D

            #  Include contribution from the group function.
            
            #  1) the non-TRIVIAL case
            
            if hasattr(self,"grftype") and ig < len( self.grftype ) and not self.grftype[ ig ] is None:
                egname = self.grftype[ ig ]; #  the group's ftype
            else:
                egname = "TRIVIAL"
                
            if mode == "Hv":
                if egname == "TRIVIAL":
                    HJv += Hinv / gsc
                else:
                    fa, grada, Hessa = eval('self.'+egname+'( self, 3, fin, ig )' )
                    sgin = lil_matrix(gin);
                    HJv  += ( ( Hessa * sgin) * (sgin.transpose().dot(v)) + grada * Hinv) / gsc
            elif mode == "HIv":
                if egname == "TRIVIAL":
                    ic  += 1
                    HJv += y[ic] * Hinv / gsc
                else:
                    fa, grada, Hessa = eval('self.'+egname+'( self, 3, fin, ig )' )
                    sgin = lil_matrix(gin);
                    ic  += 1
                    HJv += y[ic] * ( ( Hessa * sgin) * (sgin.transpose().dot(v)) + grada * Hinv) / gsc
            elif  mode == "Jv":
                ic += 1
                if derlvl >= 1:
                    if egname == "TRIVIAL":
                        HJv[ic] = gin.transpose().dot( v ) / gsc
                    else:
                        fa, grada = feval( self.name, egname, fin, ig )
                        HJv[ic] = grada * gin.transpose().dot( v ) / gsc
                else:
                    HJv[ic] = np.nan
            elif mode == "Jtv":
                if derlvl >= 1:
                    if egname == "TRIVIAL":
                        HJv += gin * v[ic] / gsc
                    else:
                        fa, grada = feval( self.name, egname, fin, ig )
                        HJv += grada * gin * v[ic] / gsc
                else:
                    HJv[ic] = np.nan
                    
            if debug:
                print( "ig = ", ig, "  HJv(final) = ", HJv )#D
        #end for iig

        return HJv.reshape(-1,1)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    #   Evaluate the value and derivatives of the Lagrangian function.
    #
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def evalLx( self, gobjlist, gconlist, x, y, nargout ):
        
        if nargout == 1:
            if len( gobjlist ) or hasattr( self, "H" ):
               Lxy = self.evalgrsum( True, gobjlist, x, 1 )
            else:
               Lxy = 0.
            if len( gconlist ):
                c   = self.evalgrsum( False, gconlist, x, 1 )
                Lxy = Lxy + y.T.dot(c)
            return float(Lxy)
        elif nargout == 2:
            if len( gobjlist ) or hasattr( self, "H" ):
                Lxy, Lgxy = self.evalgrsum( True, gobjlist, x, 2 )
            else:
                Lxy  = 0.
                Lgxy = np.zeros((len(x),1))
            if len( gconlist ):
                c, J = self.evalgrsum( False, gconlist, x, 2 )
                Lxy  = Lxy + y.T.dot(c)
                Lgxy = Lgxy  + J.T.dot(y)
            return float(Lxy), Lgxy
        elif nargout == 3:
            if len( gobjlist ) or hasattr( self, "H" ):
                Lxy, Lgxy, LgHxy = self.evalgrsum( True, gobjlist, x, 3 )
            else:
                n = len( x )
                Lxy   = 0.
                Lgxy  = np.zeros((n,1))
                LgHxy = lil_matrix((n,n))
            if len( gconlist ):
                c, J, cHi = self.evalgrsum( False, gconlist, x, 3 )
                Lxy  = Lxy + y.T.dot(c)
                Lgxy = Lgxy  + J.T.dot(y)
                for ig in range( len( gconlist ) ):
                    LgHxy = LgHxy + y[ig,0] * cHi[ig]
            return float(Lxy), Lgxy.reshape(-1,1), LgHxy
        
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    #   Evaluate the product of the Lagrangian's Hessian times a vector v.
    #
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def evalLHxyv( self, gobjlist, gconlist, x, y, v ):

        if len( gobjlist ) or hasattr( self, "H" ):
            LHxyv = self.evalHJv( "Hv", gobjlist, x, v, [] )
        else:
            LHxyv = np.zeros((len(x),1))
        if len( gconlist ):
#            for ig in range( len( gconlist ) ):
#                LHxyv += self.evalHJv( "HIv", gconlist[ig:ig+1], x, v, y )
            LHxyv += self.evalHJv( "HIv", gconlist, x, v, y )
           
        return LHxyv

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Computes the effective index name in List, or add name to List if not in there already. Return
#   the index of name in List and new = 1 if the List has been enlarged or 0 if name was already
#   present in List at the call.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def s2mpj_ii( name, List ):
    
    if name in List:
       idx = List[ name ]
       new = 0
    else:
       idx = len( List)
       List[ name ] = idx
       new = 1
    return idx, List, new
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Get the index of a nonlinear variable.  This implies adding it to the variables' dictionary ix_
#   if it is a new one, and adjusting the bounds, start point and types according to their default
#   setting.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def s2mpj_nlx( self, name, List, getxnames=None, xlowdef=None, xuppdef=None, x0def=None ):

    iv, List, newvar = s2mpj_ii( name, List );
    if( newvar ):
        self.n = self.n + 1;
        if getxnames:
            self.xnames = arrset( self.xnames, iv, name )
        if hasattr( self, "xlower" ):
            thelen = len( self.xlower )
            if ( iv <= thelen ):
                self.xlower = np.append( self.xlower, np.full( (iv-thelen+1,1), float(0.0) ), 0 )
            if not xlowdef is None:
                self.xlower[iv] = xlowdef
        if hasattr( self, "xupper" ):
            thelen = len( self.xupper )
            if ( iv <= thelen ):
                self.xupper = np.append( self.xupper, np.full( (iv-thelen+1,1), float('Inf') ), 0 )
            if not xuppdef is None:
                self.xupper[iv] = xuppdef
        try:
            self.xtype  = arrset( self.xtype, iv, 'r' )
        except:
            pass
        thelen = len( self.x0 )
        if ( iv <= thelen ):
            self.x0 = np.append( self.x0, np.full( (iv-thelen+1,1), 0.0 ), 0 )
        if not x0def is None:
            self.x0[iv] =  x0def
    return iv, List

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   An emulation of the Matlab find() function for everything that can ne enumerated
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def find( lst, condition ):
    
    return np.array([i for i, elem in enumerate(lst) if condition(elem)])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Set the elements indexed by index of an np.array (arr) to value.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def arrset( arr, index, value ):
    
    if isinstance( index, np.ndarray):
        maxind = np.max( index )
    else:
        maxind = index
    if len(arr) <= maxind:
        arr= np.append( arr, np.full( maxind - len( arr ) + 1, None ) )
    arr[index] = value
    return arr

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Set the (i,j)-th element of a list of arrays (loa) to value.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def loaset( loa, i, j, value ):
    
    if len(loa) <= i:
       loa.extend( [None] * ( i - len( loa ) + 1 ) )
    if loa[i] is None:
       loa[i]= np.full( j + 1, None )
    if len(loa[i]) <= j:
       loa[i]= np.append(loa[i],np.full( j - len( loa[i] ) + 1, None ) )
    loa[i][j] = value
    return loa

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   This tool consider all problems in list_of_python_problems (whose files are in
#   the ./python_problems directory) and selects those whose SIF classification matches
#   that given by the input string classif. Matching is in the sense of  regular expressions
#   (regexp).
#   If varargin is empty (i.e. only classif is used as input argument), the function prints
#   the list of matching problems on standard output. Otherwise, the list is output in the
#   file whose name is a string passed as varargin{1} (Further components of varargin are
#   ignored).
#
#   If the input string is 'help'  or 'h', a message is printed on the standard output
#   describing the SIF classification scheme and an brief explanation of how to use the tool.
#
#   Thanks to Greta Malaspina (Firenze) for an inital implementation in Matlab.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
def s2mpjlib_select( classif, *args ):

    if classif in [ "help", "h" ]:
    
        print( "  " )
        print( " === The classification scheme ===" )
        print( "  " )
        print( " A problem is classified by a string of the form" )
        print( "    X-XXXXr-XX-n-m" )
        print( " The first character in the string identifies the problem collection" )
        print( " from which the problem is extracted. Possible values are" )
        print( "    C the CUTEst collection;" )
        print( "    S the SPARCO collection; and" )
        print( "    N none of the above." )
        print( " The character immediately following the first hyphen defines the type" )
        print( " of variables occurring in the problem. Its possible values are" )
        print( "    C the problem has continuous variables only;" )
        print( "    I the problem has integer variables only;" )
        print( "    B the problem has binary variables only; and" )
        print( "    M the problem has variables of different types." )
        print( " The second character after the first hyphen defines the type" )
        print( " of the problem''s objective function. Its possible values are" )
        print( "    N no objective function is defined;" )
        print( "    C the objective function is constant;" )
        print( "    L the objective function is linear;" )
        print( "    Q the objective function is quadratic;" )
        print( "    S the objective function is a sum of squares; and" )
        print( "    O the objective function is none of the above." )
        print( " The third character after the first hyphen defines the type of" )
        print( " constraints of the problem. Its possible values are" )
        print( "    U the problem is unconstrained;" )
        print( "    X the problem’s only constraints are fixed variables;" )
        print( "    B the problem’s only constraints are bounds on the variables;" )
        print( "    N the problem’s constraints represent the adjacency matrix of a (linear)" )
        print( "      network;" )
        print( "    L the problem’s constraints are linear;" )
        print( "    Q the problem’s constraints are quadratic; and" )
        print( "    O the problem’s constraints are more general than any of the above alone." )
        print( " The fourth character after the first hyphen indicates the smoothness of" )
        print( " the problem. There are two possible choices" )
        print( "    R the problem is regular, that is, its first and second derivatives " )
        print( "      exist and are continuous everywhere; or" )
        print( "    I the problem is irregular." )
        print( " The integer (r) which corresponds to the fourth character of the string is" )
        print( " the degree of the highest derivatives provided analytically within the problem" )
        print( " description. It is restricted to being one of the single characters O, 1, or 2." )
        print( " The character immediately following the second hyphen indicates the primary" )
        print( " origin and/or interest of the problem. Its possible values are" )
        print( "    A the problem is academic, that is, has been constructed specifically by" )
        print( "      researchers to test one or more algorithms;" )
        print( "    M the problem is part of a modeling exercise where the actual value of the" )
        print( "      solution is not used in a genuine practical application; and" )
        print( "    R the problem’s solution is (or has been) actually used in a real")
        print( "      application for purposes other than testing algorithms." )
        print( " The next character in the string indicates whether or not the problem" )
        print( " description contains explicit internal variables. There are two possible" )
        print( " values, namely," )
        print( "    Y the problem description contains explicit internal variables; or" )
        print( "    N the problem description does not contain any explicit internal variables." )
        print( " The symbol(s) between the third and fourth hyphen indicate the number of" )
        print( " variables in the problem. Possible values are" )
        print( "    V the number of variables in the problem can be chosen by the user; or" )
        print( "    n a positive integer giving the actual (fixed) number of problem variables." )
        print( " The symbol(s) after the fourth hyphen indicate the number of constraints" )
        print( " (other than fixed variables and bounds) in the problem. Note that fixed" )
        print( " variables are not considered as general constraints here. The two possible" )
        print( " values are" )
        print( "    V the number of constraints in the problem can be chosen by the user; or" )
        print( "    m a nonnegative integer giving the actual (fixed) number of constraints." )
        print( "  " )
        print( " === Using the problem selection tool ===" )
        print( "  " )
        print( " In order to use the selection too, you should first open Python in the parent" )
        print( " of the directory containing the Python problem files, then import the library" )
        print( " by issuing the command" )
        print( "    from s2mpjlib import *" )
        print( " or, more specifiaclly," )
        print( "    from s2mpjlib import s2mpjlib_select")
        print( " The select tool may then be called with its first argument being a string ")
        print( " which specifies the class of problems of interest.  This string is constructed" )
        print( " by replacing by a dot each character in the classification string for which" )
        print( " all possible values are acceptable (the dot is a wildcard character)." )
        print( " For instance" )
        print( "    s2mpjlib_select( \"C-CSU..-..-2-0\" ) ")
        print( " lists all CUTEst unconstrained ""sum-of-squares"" problems in two continuous" )
        print( " variables, while " )
        print( "    s2mpjlib_select( ""C-C....-..-V-V"" ) " )
        print( " lists all CUTEst problems with variable number of continuous variables and" )
        print( " variable number of constraints." )
        print( " The classification strings \"unconstrained\", \"bound-constrained\", " )
        print( " \"fixed-variables\", \"general-constraints\", \"variable-n\" and " )
        print( " \"variable-m\" are also allowed." )
        print( " NOTE: any regular expression may be used as the first argument of select " )
        print( "       to specify the problem class, so that, for instance, the previous " )
        print( "       selection can also be achieved by s2mpjlib_select( \"C-C.*V-V\" ) ")
        print( " Writing the list of selected problems to a file is obtained by specifying" )
        print( " the name of the file as a second argument of select, as in ")
        print( "    s2mpjlib_select( \"C-C....-..-V-V\", filename )" )

    else:
    
        #  Modify the filter to cope with fixed numbers of variables/constraints with more
        #  than one digit.

        if classif == "unconstrained":
            classif = ".-..U.*"
        elif classif == "bound-constrained":
            classif = ".-..B.*"
        elif classif == "fixed-variables":
            classif = ".-..X.*"
        elif classif == "general-constraints":
            classif = ".-..[LNQO].*"
        elif classif == "variable-n":
            classif = ".-..B..-..-V-[V0-9]*"
        elif classif == "variable-m":
            classif = ".-..B..-..-[V0-9]*-V"
        else:
            lencl = len( classif )
            if lencl > 11 and classif[11] == ".":
                oclassif = classif
                classif  = classif[0:11] + "[V0-9]*"
                if lencl> 12:
                    classif = classif + oclassif[12:lencl]
            lenclm1 = len( classif ) - 1
            if classif[ lenclm1 ] == ".":
                classif = classif[0:lenclm1] + "[V0-9]*"
        filter_pattern = f'classification = .*{classif}'
      
        list_of_problems = "./list_of_python_problems"
        python_problems  = "./python_problems/"

        if len(args) > 0:
            fid = open( args[0], "w" )
        else:
            fid = None

        filter_pattern = f'classification = .*{classif}'
        with open( list_of_problems, 'r' ) as f:
            allprobs = f.readlines()

        for theprob in allprobs:
            theprob = theprob.strip()
            problem_path = os.path.join( python_problems, theprob )
            if os.path.isfile( problem_path ):
                with open( problem_path, 'r' ) as prob_file:
                    content = prob_file.read()
                if re.search( filter_pattern, content ):
                    if fid:
                        fid.write(f'{theprob}\n')
                    else:
                        print( theprob )

        if fid:
            fid.close()

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

