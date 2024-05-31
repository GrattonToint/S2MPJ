####################################################################################################
#####################################################################################################
#
#                                S2X library for Python
#
#   Performs the runtime actions specific to S2X, irrespective of the problem at hand.
#
#   Programming: S. Gratton and Ph. Toint (this version 17 IV 2024)
#
#####################################################################################################
#####################################################################################################

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse import csr_matrix
from pprint import pprint

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
            eval( 'self.'+'e_globs(self.pbm)' )
        except:
            pass
        try:
            eval( 'self.'+'g_globs(self.pbm)' )
        except:
            pass

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    #   Define the main evaluations actions for the problem.
    # 
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    def fx( self, x ):              # input = ( x )
        x = x.reshape(-1,1)
        if hasattr( self.pbm, "objgrps" ) or hasattr( self.pbm, "H" ):
            self.getglobs()
            return self.evalgrsum( True, self.pbm.objgrps, x, 1 )
        else:
            print( "ERROR: no objective groups in "+self.pbm.name+"!" )
        
    def fgx( self, x ):             # input = ( x )
        x = x.reshape(-1,1)
        if hasattr( self.pbm, "objgrps" ) or hasattr( self.pbm, "H" ):
            self.getglobs()
            return self.evalgrsum( True, self.pbm.objgrps, x, 2 )
        else:
            print( "ERROR: no objective groups in "+self.pbm.name+"!" )
        
    def fgHx( self, x ):            # input = ( x )
        x = x.reshape(-1,1)
        if hasattr( self.pbm, "objgrps" ) or hasattr( self.pbm, "H" ):
            self.getglobs()
            return self.evalgrsum( True, self.pbm.objgrps, x, 3 )
        else:
            print( "ERROR: no objective groups in "+self.pbm.name+"!" )
        
    def fHxv( self, x, v ):          # input = ( x, v )
        x = x.reshape(-1,1)
        v = v.reshape(-1,1)
        if hasattr( self.pbm, "objgrps" ) or hasattr( self.pbm, "H" ):
            self.getglobs()
            return self.evalHJv( "Hv", self.pbm.objgrps, x, v, [] )
        else:
            print( "ERROR: no objective groups in "+self.pbm.name+"!" )
        
    def cx( self, x ):               # input = ( x )
        x = x.reshape(-1,1)
        if hasattr( self.pbm, "congrps" ):
            self.getglobs()
            return self.evalgrsum( False, self.pbm.congrps, x, 1 )
        else:
            print( 'ERROR: no constraint groups in '+self.pbm.name+'!' )
        
    def cJx( self, x ):              # input = ( x )
        x = x.reshape(-1,1)
        if hasattr( self.pbm, "congrps" ):
            self.getglobs()
            return self.evalgrsum( False, self.pbm.congrps, x, 2 )
        else:
            print( 'ERROR: no constraint groups in '+self.pbm.name+'!' )
        
    def cJHx( self, x ):             # input = ( x )
        x = x.reshape(-1,1)
        if hasattr( self.pbm, "congrps" ):
            self.getglobs()
            return self.evalgrsum( False, self.pbm.congrps, x, 3 )
        else:
            print( 'ERROR: no constraint groups in '+self.pbm.name+'!' )

    def cJxv( self, x, v ):          # input = ( x, v )
        x = x.reshape(-1,1)
        v = v.reshape(-1,1)
        if hasattr( self.pbm, "congrps" ):
            self.getglobs()
            return evalHJv( "Jv", self.pbm.congrps, x, v, [] )
        else:
            print( 'ERROR: no constraint groups in '+self.pbm.name+'!' )
        
    def cIx( self, x, clist ):       # input = ( x, clist )
        x      = x.reshape(-1,1)
        iclist = [ self.pbm.congrps[i] for i in clist ]
        if hasattr( self.pbm, "congrps" ) and len( iclist ):
            self.getglobs()
            return self.evalgrsum( False, iclist , x, 1 )
        else:
            print( 'ERROR: empty list of constraints for '+self.pbm.name+'!' )
        
    def cIJx( self, x, clist ):          # input = ( x , clist )
        x      = x.reshape(-1,1)
        iclist = [ self.pbm.congrps[i] for i in clist ]
        if hasattr( self.pbm, "congrps" ) and len( iclist  ):
            self.getglobs()
            return self.evalgrsum( False, iclist, x, 2 )
        else:
            print( 'ERROR: empty list of constraints for '+self.pbm.name+'!' )
        
    def cIJHx( self, x , clist ):        # input = ( x, clist )
        x      = x.reshape(-1,1)
        iclist = [ self.pbm.congrps[i] for i in clist ]
        if hasattr( self.pbm, "congrps" ) and len(iclist ):
            self.getglobs()
            return self.evalgrsum( False, iclist, x, 3 )
        else:
            print( 'ERROR: empty list of constraints for '+self.pbm.name+'!' )
        
    def cIJxv( self, x, v, clist ):      # input = ( x, v, clist )
        x      = x.reshape(-1,1)
        v      = v.reshape(-1,1)
        iclist = [ self.pbm.congrps[i] for i in clist ]
        if hasattr( self.pbm, "congrps" ) and len( iclist ):
            self.getglobs()
            return self.evalHJv( "Jv", iclist, x, v, [] )
        else:
            print( 'ERROR: empty list of constraints for '+self.pbm.name+'!' )
        
    def Lxy( self, x, y ):           # input = ( x, y )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self.pbm, "congrps" ):
            return self.evalLx( self.pbm.objgrps, self.pbm.congrps, x, y, 1 )
        else:
            return self.evalgrsum( True, self.pbm.objgrps, x, 1 )
        
    def Lgxy( self, x, y ):          # input = ( x, y )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self.pbm, "congrps" ):
            return self.evalLx( self.pbm.objgrps, self.pbm.congrps, x, y, 2 )
        else:
            return self.evalgrsum( True, self.pbm.objgrps, x, 2 )
        
    def LgHxy( self, x, y ):         # input = ( x, y )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self.pbm, "congrps" ):
            return self.evalLx( self.pbm.objgrps, self.pbm.congrps, x, y, 3 )
        else:
            return self.evalgrsum( True, self.pbm.objgrps, x, 3 )
        
    def LHxyv( self, x, y, v ):      # input = ( x, y, v )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        v = v.reshape(-1,1)
        self.getglobs()
        if hasattr( self.pbm, "congrps" ):
            return self.evalLHxyv( self.pbm.objgrps, self.pbm.congrps, x, y, v )
        else:
            return self.evalHJv( "Hv", self.pbm.objgrps, x, v, [] )
        
    def LIxy( self, x, y, clist ):       # input = ( x, y, clist )
        x  = x.reshape(-1,1)
        y  = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self.pbm, "congrps" ):
            return self.evalLx( self.pbm.objgrps, [ self.pbm.congrps[i] for i in clist ], x, y, 1 )
        else:
            return self.evalgrsum( True, self.pbm.objgrps, x, 1 )
        
    def LIgxy( self, x, y, clist ):      # input = ( x, y, clist )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self.pbm, "congrps" ):
            return self.evalLx( self.pbm.objgrps, [ self.pbm.congrps[i] for i in clist ], x, y, 2 )
        else:
            return self.evalgrsum( True, self.pbm.objgrps, x, 2 )
        
    def LIgHxy( self, x, y, clist  ):    # input = ( x, y, clist )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        self.getglobs()
        if hasattr( self.pbm, "congrps" ):
            return self.evalLx( self.pbm.objgrps, [ self.pbm.congrps[i] for i in clist ], x, y, 3 )
        else:
            return self.evalgrsum( True, self.pbm.objgrps, x, 3 )
        
    def LIHxyv( self, x, y, v, I ):  # input = ( x, y, v, clist )
        x = x.reshape(-1,1)
        y = y.reshape(-1,1)
        v = v.reshape(-1,1)
        self.getglobs()
        if hasattr( self.pbm, "congrps" ):
            return self.evalLHxyv( self.pbm.objgrps, [ self.pbm.congrps[i] for i in clist ], x, y, v )
        else:
            return self.evalHJv( "Hv", self.pbm.objgrps, x, v, [] )
            
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    #   Evaluate the value of a sum of groups (and, if requested, that of of its gradient and Hessian)
    #   at x, given the problem data available in the self.pbm struct.
    #
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def evalgrsum( self, isobj, glist, x, nargout ):

        debug = False#D
        
        #   Initializations
        
        pbm = self.pbm
        n   = len( x )
        m   = len( glist )
        if isobj:
            fx = 0.0
        else:
            cx = np.zeros(( m, 1 ))
            ic = -1

        if nargout > 1:
            if isobj:
                gx = np.zeros(( n, 1))
            else:
                Jx = lil_matrix(( m, n ))

            if nargout > 2:
                if ( isobj ):
                    Hx = lil_matrix(( n, n ))
                else:
                    Hx = []

        #  Check for the presence and size of a linear term

        if hasattr( pbm, "A" ):
            sA1 = pbm.Ashape[0]
            sA2 = pbm.Ashape[1]
            has_A = True
        else:
            has_A = False

        #  Evaluate the quadratic term, if any.

        if isobj and hasattr( pbm, "H" ):
            
            Htimesx = pbm.H .dot(x)
            if nargout ==  1:
                fx += 0.5 * x.T .dot(Htimesx)
            elif nargout == 2:
                gx += Htimesx
                fx += 0.5 * x.T .dot(Htimesx)
            elif nargout == 3:
                Htimesx = pbm.H .dot(x);
                gx += Htimesx
                fx += 0.5 * x.T .dot(Htimesx)
                Hx += pbm.H
                
        if debug: #D
            if isobj:
                print( "fx(quadratic) = ", fx )
                #print( "gx = ", gx )
            else:
                print( "cx(quadratic) = ", cx )
        
        #    Loop on the groups list 

        for iig in range( len( glist )):
            ig = glist[ iig ]
            
            #  Find the group's scaling.
            
            if hasattr(pbm,"gscale"):
                if ig < len(pbm.gscale) and not pbm.gscale[ig] is None and abs( pbm.gscale[ig] ) > 1.0e-15:
                    gsc = pbm.gscale[ig]
                else:
                    gsc = 1.0
            else:
                gsc = 1.0

            #  Evaluate the linear term, if any.

            if has_A and ig < sA1:# and pbm.A[ ig, :sA2 ].nnz:
                gin           = np.zeros( (n, 1) )
                gin[:sA2, :1] = pbm.A[ ig, :sA2 ].T.toarray()
                fin           = float(-pbm.gconst[ ig ] + gin.T .dot(x))
            else:
                fin = float(-pbm.gconst[ ig ])
                if nargout == 2 or nargout == 3:
                    gin =  np.zeros(( n, 1 ))

            if nargout > 2:
                Hin = lil_matrix(( n, n ))

            if debug:
                print( "ig = ", ig, "  fin(linear)", fin ) #D

            if hasattr(pbm,"grelt") and ig < len(pbm.grelt) and not pbm.grelt[ig] is None:
                for iiel in range(len( pbm.grelt[ ig ] ) ):  #  loop on elements
                    iel    = pbm.grelt[ ig ][ iiel ]         #  the element's index
                    efname = pbm.elftype[ iel ]              #  the element's ftype
                    irange = [iv for iv in pbm.elvar[ iel ]] #  the elemental variable's indeces
                    xiel   = x[ np.array(irange) ]           #  the elemental variable's values

#                    print("iel lib = ", iel, " efname = ", efname )#D
#                    print( "irange = ", irange )#D
#                    print( "xiel = = ", xiel )#D
                    
                    if  hasattr( pbm, 'grelw' ) and ig <= len( pbm.grelw ) and not pbm.grelw[ig] is None :
                        has_weights = True;
                        wiel        = pbm.grelw[ ig ][ iiel ]
                    else:
                        has_weights = False

                    # Only the value is requested.
                    
                    if nargout == 1:
                        fiel = eval('self.'+efname +'( self.pbm, 1, xiel, iel )')
                        if ( has_weights ):
                            fin += wiel * fiel
                        else:
                            fin += fiel
                        
                    #  The value and its gradient are requested.
                    
                    elif nargout == 2:
                        fiel, giel = eval('self.'+efname +'( self.pbm, 2, xiel, iel)')
#                        print( "irange = ", irange, "  x[irange] = ", x[irange], "  efname = ", efname, "  iel = ", iel,"  fiel = ", fiel )
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

                    elif nargout == 3:
                        fiel, giel, Hiel = eval('self.'+efname +'( self.pbm, 3, xiel, iel )')
#                        print("fiel = ", fiel )#D
#                        print("giel = ", giel )#D
#                        print("Hiel = ", Hiel )#D
#                        print("efname =", efname," irange = ", irange, ", ig = ", ig," iel = ", iel )#D
                        if has_weights:
                            fin += wiel * fiel
                            for ir in range(len(irange)):
                                ii = irange[ ir ]
                                gin[ ii ] += wiel * giel[ ir ]
                                for jr in range(len( irange )):
                                    jj  = irange[ jr ]
#                                    print( "Hiel", Hiel, " efname = ", efname )#D
                                    Hin[ ii, jj ] += wiel * Hiel[ ir, jr ]
                        else:
                            fin = fin + fiel;
                            for ir in range(len(irange)):
                                ii = irange[ ir ]
                                gin[ ii ] += giel[ ir ]
                                for jr in range(len( irange )):
                                    jj  = irange[ jr ]
                                    Hin[ ii, jj ] += Hiel[ ir, jr ]
                                #end for jr
                            #end for ir
                        #end has_weights
                    #end  if nargout
                #end for iiel
            #end if hasattr(pbm,"grelt")

            if debug: #D
                if isobj:
                    print( "ig = ", ig, "  fin(nonlinear) = ", fin, "  fx = ", fx  )
                    #print( "gx = ", gx )
                else:
                    print( "ig = ", ig, "  fin(nonlinear) = ", fin, "  cx = ", cx  )
                    
            #  Evaluate the group function.
            
            #  1) the non-TRIVIAL case

            if hasattr(pbm,"grftype") and ig < len( pbm.grftype ) and not pbm.grftype[ ig ] is None:
                egname = pbm.grftype[ ig ]
            else:
                egname = "TRIVIAL"
            if egname!='TRIVIAL' and egname is not None:
                if isobj:
                    if nargout == 1:
                        fx += eval('self.'+egname+'( self.pbm, 1, fin, ig )') / gsc
                    elif nargout == 2:
                        [ fa, grada ] = eval('self.'+ egname+'( self.pbm, 2, fin, ig )')
                        fx += fa / gsc
                        gx += grada * gin / gsc
                    elif nargout == 3:
                        [ fa, grada, Hessa ] = eval('self.'+egname+'( self.pbm, 3, fin, ig )')
                        fx   += fa / gsc
                        gx   += grada * gin / gsc
                        sgin  = lil_matrix(gin)
                        Hx   += (Hessa * sgin.dot(sgin.transpose())+ grada * Hin) / gsc
                else:
                    ic = ic + 1
                    if nargout == 1:
                        fa = eval('self.'+egname+'( self.pbm, 1, fin, ig )')
                        cx[ ic ] = fa / gsc
                    elif nargout == 2:
                        fa, grada = eval('self.'+egname+'( self.pbm, 2, fin, ig )')
                        cx[ ic ]  = fa / gsc
                        sgin      = lil_matrix( gin )
                        Jx[ic,:]  = grada * sgin.T / gsc
                    elif nargout == 3:
                        fa, grada, Hessa = eval('self.'+egname+'( self.pbm, 3, fin, ig )') 
                        cx[ ic ] = fa / gsc
                        sgin     = lil_matrix( gin )
                        Jx[ic,:] = grada * sgin.T / gsc
                        Hx.append( ( Hessa * sgin.dot( sgin.transpose() )+ grada * Hin ) / gsc )
                        
            #  2) the TRIVIAL case: the group function is the identity
            
            else:
                if isobj:
                    if nargout == 1:
                        fx += fin / gsc 
                    if nargout == 2:
                        fx += fin / gsc
                        gx += gin / gsc
                    if nargout == 3:
                        fx += fin / gsc
                        gx += gin / gsc
                        Hx += Hin / gsc
                else:
                    ic = ic + 1
                    if nargout == 1:
                        cx[ic] = fin / gsc
                    elif nargout == 2:
                        cx[ic] = fin / gsc
                        Jx[ic,:] = gin.transpose() / gsc
                    elif nargout == 3:
                        cx[ic] = fin / gsc
                        Jx[ic,:] = gin.transpose() / gsc
                        Hx.append( Hin / gsc )

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
    #                    (glist = cons groups)
    #     mode = "Jv"  : evaluate the product of the constraints' Jacobian times v (glist = cons groups)
    #     The vector y is unused (and unreferenced) for modes "Hv" and "Jv".
    #
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    def evalHJv( self, mode, glist, x, v, y ):

        debug = False#D
        
        #   Initializations
        
        pbm = self.pbm
        n   = len( x );
        m   = len( glist );
        if mode == "Hv":
            HJv = np.zeros((n,1))
        elif mode == "HIv":
            HJv = np.zeros((n,1))
            ic  = -1
        else:
            HJv = np.zeros((m,1))
            ic = -1;

        #  Check for the presence and size of a linear term

        if hasattr( pbm, "A" ):
            sA1 = pbm.Ashape[0]
            sA2 = pbm.Ashape[1]
            has_A = True
        else:
            has_A = False

        #  Evaluate the quadratic term, if any.

        if mode == "Hv" and hasattr( pbm, "H" ):
            HJv += pbm.H.dot(v);
        
        if debug: #D
            print( "HJv(quadratic) = ", HJv )

        for iig in range(len( glist )):
            ig = glist[ iig ];
            
            #  Find the group's scaling.
            
            if hasattr(pbm,"gscale"):
                if ig < len(pbm.gscale) and not pbm.gscale[ig] is None and abs( pbm.gscale[ig] ) > 1.0e-15:
                    gsc = pbm.gscale[ig];
                else:
                    gsc = 1.0;
            else:
                gsc = 1.0;

            #  Evaluate the linear term, if any.

            fin = 0.0
            if has_A and ig < sA1:
                gin =  np.zeros((n,1))
                gin[:sA2, :1] = pbm.A[ ig, :sA2 ].T.toarray()
                fin           = float(-pbm.gconst[ ig ] + gin.T .dot(x))
            elif mode == "Hv" or mode == "HIv":
                fin = float( -pbm.gconst[ ig ] )
                gin = np.zeros((n,1))
            else:
                gin = np.zeros((n,1))

            if debug:
                print( "ig = ", ig, "  fin(linear)", fin ) #D

            Hinv = np.zeros((n,1))

            if hasattr(pbm,"grelt") and ig < len(pbm.grelt) and not pbm.grelt[ig] is None:
                for iiel in range(len( pbm.grelt[ ig ] ) ):  #  loop on elements
                    iel    = pbm.grelt[ ig ][ iiel ]         #  the element's index
                    efname = pbm.elftype[ iel ];             #  the element's ftype
                    irange = [iv for iv in pbm.elvar[ iel ]] #  the elemental variable's indeces
                    xiel   = x[ np.array(irange) ]           #  the elemental variable's values

                    if  hasattr( pbm, 'grelw' ) and ig <= len( pbm.grelw ) and not pbm.grelw[ig] is None :
                        has_weights = 1;
                        wiel        = pbm.grelw[ ig ][ iiel ]
                    else:
                        has_weights = 0;

                    #  The group is an objective group
                    
                    if mode == "Hv" or mode == "HIv":
                        fiel, giel, Hiel = eval('self.'+efname +'( self.pbm, 3, xiel, iel )')
                        if has_weights:
                            fin += wiel * fiel;
                            for ir in range(len(irange)):
                                ii = irange[ ir ]
                                gin[ ii ] = gin[ ii ] + wiel * giel[ ir ]
                                for jr in range(len( irange )):
                                    jj  = irange[ jr ]
                                    Hinv[ ii ] +=  wiel * Hiel[ ir, jr ] * v[ jj ]
                        else:
                            fin = fin + fiel;
                            for ir in range(len(irange)):
                                ii = irange[ ir ]
                                gin[ ii ] = gin[ ii ] + giel[ ir ]
                                for jr in range(len( irange )):
                                    jj  = irange[ jr ];
                                    Hinv[ ii ] +=  Hiel[ ir, jr ].dot( v[ jj ] )

                     #   The group is a constraint group.

                    else:

                        fiel, giel = eval('self.'+efname +'( self.pbm, 2, xiel, iel)')
                        if  has_weights:
                            fin = fin + wiel * fiel;
                            for ir in range(len(irange)):
                                ii = irange[ ir ]
                                gin[ ii ] += wiel * giel[ ir ]
                        else:
                            fin = fin + fiel;
                            for ir in range(len(irange)):
                                ii = irange[ ir ]
                                gin[ ii ] += giel[ ir ]
                    #end isobj
                #end for iiel
            #end if hasattr(pbm,"grelt")
            
            if debug: #D
                print( "ig = ", ig, "  Hinv(nonlinear) = ", Hinv, "  HJv = ", HJv  )#D

            #  Include contribution from the group function.
            
            #  1) the non-TRIVIAL case
            
            if hasattr(pbm,"grftype") and ig < len( pbm.grftype ) and not pbm.grftype[ ig ] is None:
                egname = pbm.grftype[ ig ]; #  the group's ftype
            else:
                egname = "TRIVIAL"
                
            if mode == "Hv":
                if egname == "TRIVIAL":
                    HJv += Hinv / gsc
                else:
                    fa, grada, Hessa = eval('self.'+egname+'( self.pbm, 3, fin, ig )' )
                    sgin = lil_matrix(gin);
                    HJv  += ( ( Hessa * sgin) * (sgin.transpose().dot(v)) + grada * Hinv) / gsc
            elif mode == "HIv":
                if egname == "TRIVIAL":
                    ic  += 1
                    HJv += y[ic] * Hinv / gsc
                else:
                    fa, grada, Hessa = eval('self.'+egname+'( self.pbm, 3, fin, ig )' )
                    sgin = lil_matrix(gin);
                    ic  += 1
                    HJv += y[ic] * ( ( Hessa * sgin) * (sgin.transpose().dot(v)) + grada * Hinv) / gsc
            else:
                ic += 1
                if egname == "TRIVIAL":
                    HJv[ic] = gin.transpose().dot( v ) / gsc
                else:
                    [ fa, grada ] = feval( pbm.name, egname, fin, ig )
                    HJv[ic] = grada * gin.transpose().dot( v ) / gsc
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
            if len( gobjlist ) or hasattr( self.pbm, "H" ):
               Lxy = self.evalgrsum( True, gobjlist, x, 1 )
            else:
               Lxy = 0.
            if len( gconlist ):
                c   = self.evalgrsum( False, gconlist, x, 1 )
                Lxy = Lxy + y.T.dot(c)
            return float(Lxy)
        elif nargout == 2:
            if len( gobjlist ) or hasattr( self.pbm, "H" ):
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
            if len( gobjlist ) or hasattr( self.pbm, "H" ):
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

        if len( gobjlist ):
            LHxyv = self.evalHJv( "Hv", gobjlist, x, v, [] )
        else:
            LHxyv = np.zeros((len(x),1))
        if len( gconlist ):
            for ig in range( len( gconlist ) ):
                LHxyv += self.evalHJv( "HIv", gconlist[ig:ig+1], x, v, y )
           
        return LHxyv

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   Computes the effective index name in List, or add name to List if not in there already. Return
#   the index of name in List and new = 1 if the List has been enlarged or 0 if name was already
#   present in List at the call.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def s2x_ii( name, List ):
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

def s2x_nlx( name, List, pb=None, getxnames=None, xlowdef=None, xuppdef=None, x0def=None ):
    
    iv, List, newvar = s2x_ii( name ,List );
    if( newvar ):
        pb.n = pb.n + 1;
        if getxnames:
            pb.xnames = arrset( pb.xnames, iv, name )
        if hasattr( pb, "xlower" ):
            thelen = len( pb.xlower )
            if ( iv <= thelen ):
                pb.xlower = np.append( pb.xlower, np.full( (iv-thelen+1,1), float(0.0) ), 0 )
            if not xlowdef is None:
                pb.xlower[iv] 
        if hasattr( pb, "xupper" ):
            thelen = len( pb.xupper )
            if ( iv <= thelen ):
                pb.xupper = np.append( pb.xupper, np.full( (iv-thelen+1,1), float('Inf') ), 0 )
            if not xuppdef is None:
                pb.xupper[iv] 
        try:
            pb.xtype  = arrset( pb.xtype, iv, 'r' )
        except:
            pass
        thelen = len( pb.x0 )
        if ( iv <= thelen ):
            pb.x0 = np.append( pb.x0, np.full( (iv-thelen+1,1), 0.0 ), 0 )
        if not x0def is None:
            pb.x0[iv] =  x0def
    return iv, List, pb

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#   An emulation of the Matlab find() function for everything that can ne enumerated
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def find(lst,condition):
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

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
