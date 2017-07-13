VERSION_NAME = 'SAPCOR_v04f'

# MODULE NAME
MODULE_NAME_THIS__FORCE_BLOCK_V_F = VERSION_NAME+'_Force_Block_V_F'


import glob

# -------------------------------------------------------------------
# MODULE NAMES TO BE IMPORTED
MODULE_NAMES_FORCE_BLOCK_V_F = []
MODULE_NAMES_FORCE_BLOCK_V_F.append('Misc')
MODULE_NAMES_FORCE_BLOCK_V_F.append('Input')
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Import Modules
for Module_Name in MODULE_NAMES_FORCE_BLOCK_V_F:
  FileName = VERSION_NAME+'_'+Module_Name+'_Build*.py'
  FList = glob.glob(FileName)
  assert len(FList)!=0,'ERROR: Module not found : '+FileName
  FList.sort()
  FileName = FList[-1][0:-3]
  exec('from '+FileName+' import *')
# -------------------------------------------------------------------


from numpy import loadtxt
from numpy import array,zeros,hstack



#-------------------------------------------------------------------------

def Force_Block_V_F (w,t,Core,K,L,Accel):
  
  """
    Defines the differential equations for the coupled spring-mass system.

    Arguments:
        w :  vector of the state variables at time, t:
                  w = [U,DU,W,DW,R,DR,...]
        t :  time
        Core : Var space
  """

  # VERBOSE
  if '/Force/' in VERBOSE:
    Verbose('<FORCE_BLOCK_V_F>')
    Verbose('time=%e, (K,L)=(%d,%d)'%(t,K,L))
    

  # Column Size
  N = Core['N'][K]

  # Get Block Prop of (K,L)
#  BN = Core['Array'][K][L]
  BN = Core['BTNsKL'][K][L]
  BT = BlockTypes[BN]
  A = BT['a']
  H = BT['h']
  M = BT['M']
  I = BT['I']
  Kv = BT['Kv']
  Cv = BT['Cv']
  mu_s = BT['mu_s']
  mu_k = BT['mu_k']
  d_mu = BT['d_mu']
  xi_F_cr = BT['xi_F_cr']
#  Fixed = BT['Fixed']

  # Get Block Prop of (K,L-1)
  # Subscript, m, means 'Minus 1'
  if L!=1:
#    BN = Core['Array'][K][L-1]
    BN = Core['BTNsKL'][K][L-1]
    BT = BlockTypes[BN]
#    Am = BT['a']
    Hm = BT['h']
    Mm = BT['M']
    Im = BT['I']
#    Kvm = BT['Kv']
#    Cvm = BT['Cv']
#    Fixedm = BT['Fixed']
  else:
    if Core['Flag_CSB']==True: # CSB
      BT = BlockTypes['CSB']
      Hm = 0
      Mm = BT['M']
      Im = BT['I']
    else: # Base
      Hm = 0
      Mm = Im = 1
      
  

  # Get current state of (K,L)
  U,DU,W,DW,R,DR = GetSolFromVect(Core,w,K,L) # U,DU,W,DW,R,DR

  # Get current state of (K,L-1)
  # If L==1:
  #   It means lower component is CSB or Base.
  #   If Flag_CSB==True: Get state of CSB  (which means (K,L)==(1,0))
  #   Else             : Get state of Base (which means (K,L)==(0,0))
  if L!=1:
    Um,DUm,Wm,DWm,Rm,DRm = GetSolFromVect(Core,w,K,L-1) # U,DU,W,DW,R,DR
  else:
    if Core['Flag_CSB']==True:
      Um,DUm,Wm,DWm,Rm,DRm = GetSolFromVect(Core,w,1,0) # CSB
    else:
      Um,DUm,Wm,DWm,Rm,DRm = GetSolFromVect(Core,w,0,0) # Base
    


##############################################################################
#  VERTICAL FORCES
##############################################################################

#(v04f)  if Core['ApplyForces']['Block_V']==True:
# ['Block_V'] is already checked in VectorField()

  # Spring Contraction between (K,L-1)~(K,L)
  if L!=1:
    Gamma_L  = Wm - Hm * (1-cos( Rm )) + A  * sin( Rm )
    Gamma_L -= W  + H  * (1-cos( R  )) + A  * sin( R  )
    Gamma_L *= 0.5
    Gamma_R  = Wm - Hm * (1-cos( Rm )) - A  * sin( Rm )
    Gamma_R -= W  + H  * (1-cos( R  )) - A  * sin( R  )
    Gamma_R *= 0.5

    DGamma_L  = DWm - ( Hm * sin( Rm ) - A  * cos( Rm ) ) * DRm
    DGamma_L -= DW  + ( H  * sin( R  ) + A  * cos( R  ) ) * DR
    DGamma_L *= 0.5
    DGamma_R  = DWm - ( Hm * sin( Rm ) + A  * cos( Rm ) ) * DRm
    DGamma_R -= DW  + ( H  * sin( R  ) - A  * cos( R  ) ) * DR
    DGamma_R *= 0.5
  else: # Interaction between CSB(or Base)~(K,1)
    Gamma_L  = Wm
    Gamma_L -= W  + H  * (1-cos( R  )) + A  * sin( R  )
    Gamma_L *= 0.5
    Gamma_R  = Wm
    Gamma_R -= W  + H  * (1-cos( R  )) - A  * sin( R  )
    Gamma_R *= 0.5

    DGamma_L  = DWm
    DGamma_L -= DW  + ( H  * sin( R  ) + A  * cos( R  ) ) * DR
    DGamma_L *= 0.5
    DGamma_R  = DWm
    DGamma_R -= DW  + ( H  * sin( R  ) - A  * cos( R  ) ) * DR
    DGamma_R *= 0.5

  
  # Forces and Moments
  F_VL,   F_VR = 0,0
  M_VL, M_VR, M_VLm, M_VRm = 0,0,0,0
 
  if Gamma_L>0:
    F_VL  =  Kv * Gamma_L  + Cv * DGamma_L
  if Gamma_R>0:
    F_VR  =  Kv * Gamma_R  + Cv * DGamma_R

  M_VL  =  F_VL * ( A*cos(R ) + H *sin(R ) )
  M_VR  = -F_VR * ( A*cos(R ) - H *sin(R ) )
  M_VLm = -F_VL * ( A*cos(Rm) - Hm*sin(Rm) )
  M_VRm =  F_VR * ( A*cos(Rm) + Hm*sin(Rm) )

  
  #---------------------------
  # ELIMINATE STICKING FORCE
  #---------------------------
  if FLAG_NoStickForce == True:
    if M_VL  < 0. : M_VL  = 0.
    if M_VR  < 0. : M_VR  = 0.
    if M_VLm > 0. : M_VLm = 0.
    if M_VRm > 0. : M_VRm = 0.
  #---------------------------


  # Assemble Forces
  AccelF  =  (F_VL  + F_VR ) / M
  AccelM  =  (M_VL  + M_VR ) / I
  AccelFm = -(F_VL  + F_VR ) / Mm
  AccelMm =  (M_VLm + M_VRm) / Im
  
  # Add Accel on (K,L)
  IndexW = Core['Index'][K][L]*6
  Accel[IndexW+3] += AccelF # DDW
  Accel[IndexW+5] += AccelM # DDR

  # Add Accel on (K,L-1)
  
  if L!=1: # Normal Blocks
    IndexWm = Core['Index'][K][L-1]*6
    Accel[IndexWm+3] += AccelFm # DDW
    Accel[IndexWm+5] += AccelMm # DDR
 
  elif Core['Flag_CSB']==True: # Add Accel on CSB
    pass # if CSB: Do nothing

  else: pass # if Base: Do nothing


  # VERBOSE
  if '/ForceV/' in VERBOSE:
    Verbose('-'*80)
    Verbose('Force_Block_V(), F_V Caculation Verification')
    Verbose('Block (%d,%d)'%(K,L))
    Verbose('U, DU, W, DW, R, DR  = %e, %e, %e, %e, %e, %e'%(U,DU,W,DW,R,DR))
    Verbose('Um,DUm,Wm,DWm,Rm,DRm = %e, %e, %e, %e, %e, %e'%(Um,DUm,Wm,DWm,Rm,DRm))
    Verbose('Gamma  = %e, %e'%(Gamma_L,Gamma_R))
    Verbose('DGamma = %e, %e'%(DGamma_L,DGamma_R))
    Verbose('FV  = %e, %e'%(F_VL,F_VR))
    Verbose('MFV = %e, %e, %e, %e'%(M_VL,M_VR,M_VLm,M_VRm))
    Verbose('Accel(%d,%d) = %e, %e'%(K,L,AccelF,AccelM))
    Verbose('Accel[%d] = %e'%(IndexW+3,Accel[IndexW+3]))
    Verbose('Accel[%d] = %e'%(IndexW+5,Accel[IndexW+5]))
    if L!=1:
      Verbose('Accel(%d,%d) = %e, %e'%(K,L-1,AccelFm,AccelMm))
      Verbose('Accel[%d] = %e'%(IndexWm+3,Accel[IndexWm+3]))
      Verbose('Accel[%d] = %e'%(IndexWm+5,Accel[IndexWm+5]))






##############################################################################
#  FRICTION FORCES
##############################################################################

#(v04f)  if Core['ApplyForces']['Block_F']==True:
  if Core['ApplyForces']['Block_VF']==True:

    # Init
    xi_L=xi_R=0
    F_FL=F_FLm=M_FL=M_FLm=F_FR=F_FRm=M_FR=M_FRm=0
    
    #===========================================================================
    # LEFT start
    
    # Vertical Reaction should be positive
    if F_VL>0: 
      
      # Relative Velocity
      xi_L    = DU  - ( H  * cos( R  ) - A * sin( R  ) ) * DR
      if L!=1:
        xi_L -= DUm + ( Hm * cos( Rm ) + A * sin( Rm ) ) * DRm
      else:
        xi_L -= DUm
      
      # Relative Velocity should be nonzero
      if xi_L!=0:
        
        # Friction Force: F_FL      
        if abs(xi_L) <= xi_F_cr:
          # Viscous Slip Condition
          F_FL = - mu_s * F_VL * xi_L / xi_F_cr
        else:
          # Slip Condition
          mu = mu_k + (mu_s - mu_k) * exp( -d_mu*(abs(xi_L)-xi_F_cr) )
          F_FL = - (+(xi_L>0) or -(xi_L<0)) * mu * F_VL
          
        # Other Forces & Moments: F_FLm, M_FL, M_FLm
        F_FLm = -F_FL
        M_FL    = -F_FL  * ( H  * cos( R  ) - A * sin( R  ) )
        if L!=1:
          M_FLm =  F_FLm * ( Hm * cos( Rm ) + A * sin( Rm ) )
        
      # Check Point: xi==0 -> Do nothing
        
    # Check Point: F_VL<=0 -> Do nothing
    
    
    # LEFT end
    #===========================================================================
    
  
    #===========================================================================
    # RIGHT start
    
    # Vertical Reaction should be positive
    if F_VR>0: 
      
      # Relative Velocity
      xi_R    = DU  - ( H  * cos( R  ) + A * sin( R  ) ) * DR
      if L!=1:
        xi_R -= DUm + ( Hm * cos( Rm ) - A * sin( Rm ) ) * DRm
      else:
        xi_R -= DUm
      
      # Relative Velocity should be nonzero
      if xi_R!=0:
        
        # Friction Force: F_FR      
        if abs(xi_R) <= xi_F_cr:
          # Viscous Slip Condition
          F_FR = - mu_s * F_VR * xi_R / xi_F_cr
        else:
          # Slip Condition
          mu = mu_k + (mu_s - mu_k) * exp( -d_mu*(abs(xi_R)-xi_F_cr) )
          F_FR = - (+(xi_R>0) or -(xi_R<0)) * mu * F_VR
        
        # Other Forces & Moments: F_FRm, M_FR, M_FRm
        F_FRm = -F_FR
        M_FR    = -F_FR  * ( H  * cos( R  ) + A * sin( R  ) )
        if L!=1:
          M_FRm =  F_FRm * ( Hm * cos( Rm ) - A * sin( Rm ) )
        
      # Check Point: xi==0 -> Do nothing
        
    # Check Point: F_VR<=0 -> Do nothing
    
    
    # RIGHT end
    #===========================================================================
  
    # Assemble Accel
    AccelF  =  (F_FL  + F_FR ) / M
    AccelM  =  (M_FL  + M_FR ) / I
    AccelFm = -(F_FL  + F_FR ) / Mm
    AccelMm =  (M_FLm + M_FRm) / Im
    
    # Add Accel on (K,L)
    IndexW = Core['Index'][K][L]*6
    Accel[IndexW+1] += AccelF # DDU
    Accel[IndexW+5] += AccelM # DDR
  
    # Add Accel on (K,L-1)
    
    if L!=1: # Normal Blocks
      IndexWm = Core['Index'][K][L-1]*6
      Accel[IndexWm+1] += AccelFm # DDU
      Accel[IndexWm+5] += AccelMm # DDR
   
    elif Core['Flag_CSB']==True: # Add Accel on CSB
#(v04f) IndexWm = Core['Index'][K][L-1]*6
      IndexWm = Core['Index'][1][0]*6
      Accel[IndexWm+1] += AccelFm # DDU
      pass # No moment on CSB
  
    else: pass # if Base: Do nothing


    # VERBOSE
    if '/ForceVF/' in VERBOSE:
      Verbose('-'*80)
      Verbose('Force_Block_V_F(), Friction Caculation Verification')
      Verbose('Block (%d,%d)'%(K,L))
      Verbose('U, DU, W, DW, R, DR  = %e, %e, %e, %e, %e, %e'%(U,DU,W,DW,R,DR))
      Verbose('Um,DUm,Wm,DWm,Rm,DRm = %e, %e, %e, %e, %e, %e'%(Um,DUm,Wm,DWm,Rm,DRm))
      Verbose('xi  = %e, %e'%(xi_L,xi_R))
      Verbose('FF  = %e, %e'%(F_FL,F_FR))
      Verbose('MF  = %e, %e, %e, %e'%(M_FL,M_FR,M_FLm,M_FRm))
      Verbose('Accel(%d,%d) = %e, %e'%(K,L,AccelF,AccelM))
      Verbose('Accel[%d] = %e'%(IndexW+1,Accel[IndexW+1]))
      Verbose('Accel[%d] = %e'%(IndexW+5,Accel[IndexW+5]))
      if L!=1:
        Verbose('Accel(%d,%d) = %e, %e'%(K,L-1,AccelFm,AccelMm))
        Verbose('Accel[%d] = %e'%(IndexWm+1,Accel[IndexWm+1]))
        Verbose('Accel[%d] = %e'%(IndexWm+5,Accel[IndexWm+5]))



  # END
  return Accel
