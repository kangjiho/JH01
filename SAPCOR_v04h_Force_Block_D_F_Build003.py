VERSION_NAME = 'SAPCOR_v04h'

# MODULE NAME
MODULE_NAME_THIS__FORCE_BLOCK_D = VERSION_NAME+'_Force_Block_D'

###############################################################################
#{    SYSTEM INITIALIZATION

import glob

# -------------------------------------------------------------------
# MODULE NAMES TO BE IMPORTED
MODULE_NAMES_FORCE_BLOCK_D = []
MODULE_NAMES_FORCE_BLOCK_D.append('Misc')
MODULE_NAMES_FORCE_BLOCK_D.append('Input')
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Import Modules
for Module_Name in MODULE_NAMES_FORCE_BLOCK_D:
  FileName = VERSION_NAME+'_'+Module_Name+'_Build*.py'
  FList = glob.glob(FileName)
  assert len(FList)!=0,'ERROR: Module not found : '+FileName
  FList.sort()
  FileName = FList[-1][0:-3]
  exec('from '+FileName+' import *')
# -------------------------------------------------------------------




from numpy import loadtxt
from numpy import array,zeros,hstack

#}#############################################################################



#-------------------------------------------------------------------------
def Force_Block_D_F (w,t,Core,K,L,Accel,POST=False,FO_DIR=FO_DIR):
#{
  

  ##############################################################################
  #{ DOWEL FORCES

  
  #{ VERBOSE
  if '/Force/' in VERBOSE:
    Verbose('<FORCE_BLOCK_D_F>')
    Verbose('time=%e (K,L)=(%d,%d)'%(t,K,L))
  #} END OF VERBOSE
  
  # Column Size
  N = len(Core['Array'][K])

  #{ Get Block Prop of (K,L)
  BN = Core['BTNsKL'][K][L]
  BT = BlockTypes[BN]
  D = BT['d']
  H = BT['h']
  M = BT['M']
  I = BT['I']
  Kd = BT['Kd']
  Cd = BT['Cd']
  mu_s = BT['mu_s']
  mu_k = BT['mu_k']
  d_mu = BT['d_mu']
  xi_F_cr = BT['xi_F_cr']
  omega_F_cr = xi_F_cr  # Use same critical relative velocity (can be changed)
#  Fixed = BT['Fixed']
  #} END OF GET BLOCK PROP (K,L)

  #{ Get Block Prop of (K,L-1)
  # Subscript, m, means 'Minus 1'
  if L!=1:
    BN = Core['BTNsKL'][K][L-1]
    BT = BlockTypes[BN]
    Hm = BT['h']
    Mm = BT['M']
    Im = BT['I']
#    Fixedm = BT['Fixed']
  else:
    if Core['Flag_CSB']==True: # CSB
      BT = BlockTypes['CSB']
      Hm = BT['h']
      Mm = BT['M']
      Im = BT['I']
#      Fixedm = BT['Fixed']
    else: # Base
      Hm = 0
      Mm = Im = 1
#      Fixedm = 'FixedToBase'
  #} END OF GET BLOCK PROP OF (K,L-1)



  #{ Get current state of (K,L)
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
  #} END OF GET CURRENT STATE OF (K,L)  




  #{ Signs
  # (K,L-1)~(K,L)
  #       h1, d1, h2, d2, delta
  SL  = ( +1, +1, -1, +1, +1 )
  SR  = ( +1, -1, -1, -1, -1 )
  SDL = ( +1, +1, -1, +1 )
  SDR = ( +1, -1, -1, -1 )
  #         h,  d
  SML  = ( +1, -1 )
  SMR  = ( +1, +1 )
  SMLm = ( +1, +1 )
  SMRm = ( +1, -1 )
  #} END OF SIGNS



  #{ Spring Contraction between (K,L-1)~(K,L)
  Beta_L  = Um + SL[0] * Hm * sin( Rm ) + SL[1] * D * (1-cos( Rm ))
  Beta_L -= U  + SL[2] * H  * sin( R  ) + SL[3] * D * (1-cos( R  ))
  Beta_L += SL[4] * DeltaD
  Beta_R  = Um + SR[0] * Hm * sin( Rm ) + SR[1] * D * (1-cos( Rm ))
  Beta_R -= U  + SR[2] * H  * sin( R  ) + SR[3] * D * (1-cos( R  ))
  Beta_R += SR[4] * DeltaD

  DBeta_L  = DUm + ( SDL[0]*Hm*cos(Rm) + SDL[1]*D*sin(Rm) ) * DRm
  DBeta_L -= DU  + ( SDL[2]*H *cos(R ) + SDL[3]*D*sin(R ) ) * DR
  DBeta_R  = DUm + ( SDR[0]*Hm*cos(Rm) + SDR[1]*D*sin(Rm) ) * DRm
  DBeta_R -= DU  + ( SDR[2]*H *cos(R ) + SDR[3]*D*sin(R ) ) * DR
  #} END OF SPRING CONTRACTION



  #{ Forces and Moments
  F_DL, F_DR = 0,0
  M_DL, M_DR, M_DLm, M_DRm = 0,0,0,0
  
  if Beta_L<0:  # <0 means contact point is on the left side of dowel
    F_DL  =  Kd * Beta_L  + Cd * DBeta_L

    #-------------------------------------------
    # ELIMINATE STICKING FORCE
    # Left Side
    # Upper Block, Force Direction : <-- (Negative)
    #    
    #    *+++++++++++++++++++++++
    #    *++**********+++++++++++
    #    *++*        *+++++++++++
    #    *++*        *+++++++++++
    #    ****<--F_DL ************
    #
    # Eliminate Positve Force
    #-------------------------------------------
    if FLAG_NoStickForce == True:
      if F_DL > 0. : F_DL = 0.
    #-------------------------------------------

    M_DL  = -F_DL * ( SML [0] *H *cos(R ) + SML [1] *D*sin(R ) )
    M_DLm = -F_DL * ( SMLm[0] *Hm*cos(Rm) + SMLm[1] *D*sin(Rm) )
  if Beta_R>0:
    F_DR  =  Kd * Beta_R  + Cd * DBeta_R

    #-------------------------------------------
    # ELIMINATE STICKING FORCE
    # Right Side
    # Upper Block, Force Direction : --> (Positive)
    #    
    #    +++++++++++++++++++++++*
    #    +++++++++++**********++*
    #    +++++++++++*        *++*
    #    +++++++++++*        *++*
    #    ************ F_DR-->****
    #
    # Eliminate Negative Force
    #-------------------------------------------
    if FLAG_NoStickForce == True:
      if F_DR < 0. : F_DR = 0.
    #-------------------------------------------

    M_DR  = -F_DR * ( SMR [0] *H *cos(R ) + SMR [1] *D*sin(R ) )
    M_DRm = -F_DR * ( SMRm[0] *Hm*cos(Rm) + SMRm[1] *D*sin(Rm) )
  #} END OF FORCES AND MOMENTS



  #{ Assemble Forces
  AccelF  =  (F_DL  + F_DR ) / M
  AccelM  =  (M_DL  + M_DR ) / I
  AccelFm = -(F_DL  + F_DR ) / Mm
  AccelMm =  (M_DLm + M_DRm) / Im
  
  # Add Accel on (K,L)
  IndexW = Core['Index'][K][L]*6
  Accel[IndexW+1] += AccelF # DDU
  Accel[IndexW+5] += AccelM # DDR

  # Add Accel on (K,L-1)
  
  if L!=1: # Normal Blocks
    IndexW = Core['Index'][K][L-1]*6
    Accel[IndexW+1] += AccelFm # DDU
    Accel[IndexW+5] += AccelMm # DDR
 
  elif Core['Flag_CSB']==True: # Add Accel on CSB
    IndexW = 6
    Accel[IndexW+1] += AccelFm # DDU
    # Igonore AccelMm for CSB

  pass # if Base: Do nothing
        
  #} END OF ASSEMBLE FORCES

  #{ POST
  if POST==True:

    # FILENAME
    Filename = FO_DIR+'/(%2d,%2d)_D.csv'%(K,L)

    #{ HEADER
    try: fo = open(Filename,'r')
    except: # File doesn't exist -> This is the first time. Make header.
      temp  = 't,'
      temp += 'BetaL,F_DLk,dBetaL,F_DLc,F_DL,M_DL,M_DLm,'
      temp += 'BetaR,F_DRk,dBetaR,F_DRc,F_DR,M_DR,M_DRm\n'
      fo = open(Filename,'w')
      fo.write(temp)
      fo.close()
    #}

    #{ F_DLk, F_DLc, F_DRk, F_DRc
    F_DLk = F_DLc = 0.
    F_DRk = F_DRc = 0.
    if Beta_L<0:
      F_DLk = Kd * Beta_L
      F_DLc = Cd * DBeta_L
    if Beta_R>0:
      F_DRk = Kd * Beta_R
      F_DRc = Cd * DBeta_R
    #}

    #{ WRITE
    temp  = '%e, '%t
    #        1   2   3   4   5   6   7    : 1       2      3        4      5     6     7
    temp += '%e, %e, %e, %e, %e, %e, %e, '%(Beta_L, F_DLk, DBeta_L, F_DLc, F_DL, M_DL, M_DLm)
    temp += '%e, %e, %e, %e, %e, %e, %e\n'%(Beta_R, F_DRk, DBeta_R, F_DRc, F_DR, M_DR, M_DRm)
    fo = open(Filename,'a')
    fo.write(temp)
    fo.close()
    #}
  #}

  #} END OF DOWEL FORCES
  ##############################################################################


  ##############################################################################
  #{  FRICTION FORCES
  ##############################################################################



  if Core['ApplyForces']['Block_DF']==True:

    # Init
    omega_L=omega_R=0
    F_DFL=F_DFLm=M_DFL=M_DFLm=F_DFR=F_DFRm=M_DFR=M_DFRm=0
    
    #===========================================================================
    #{ LEFT start
    
    # Vertical Reaction should be positive
    if F_DL<0: 
      
      # Relative Velocity
      omega_L    = DW  + (   H  * sin( R  ) + D * cos( R  ) ) * DR
      if L!=1:
        omega_L -= DWm + ( - Hm * sin( Rm ) + D * cos( Rm ) ) * DRm
      else:
        omega_L -= DWm
      
      # Relative Velocity should be nonzero
      if omega_L!=0:
        
        # Friction Force: F_FL      
        if abs(omega_L) <= omega_F_cr:
          # Viscous Slip Condition
          F_DFL = + mu_s * F_DL * omega_L / omega_F_cr
        else:
          # Slip Condition
          mu = mu_k + (mu_s - mu_k) * exp( -d_mu*(abs(omega_L)-omega_F_cr) )
          F_DFL = + (+(omega_L>0) or -(omega_L<0)) * mu * F_DL
          
        # Other Forces & Moments: F_DFLm, M_DFL, M_DFLm
        F_DFLm   = -F_DFL
        M_DFL    = +F_DFL  * ( D * cos( R  ) + H  * sin( R  ) )
        if L!=1:
          M_DFLm = +F_DFLm * ( D * cos( Rm ) - Hm * sin( Rm ) )
        
      # Check Point: omega==0 -> Do nothing
        
    # Check Point: F_VL<=0 -> Do nothing
    
    
    #} LEFT end
    #===========================================================================
    
  
    #===========================================================================
    #{ RIGHT start
    
    # Vertical Reaction should be positive
    if F_DR>0: 
      
      # Relative Velocity
      omega_R    = DW  - ( + H  * sin( R  ) - D * cos( R  ) ) * DR
      if L!=1:
        omega_R -= DWm + ( - Hm * sin( Rm ) - D * cos( Rm ) ) * DRm
      else:
        omega_R -= DWm
      
      # Relative Velocity should be nonzero
      if omega_R!=0:
        
        # Friction Force: F_FR      
        if abs(omega_R) <= omega_F_cr:
          # Viscous Slip Condition
          F_DFR = - mu_s * F_DR * omega_R / omega_F_cr
        else:
          # Slip Condition
          mu = mu_k + (mu_s - mu_k) * exp( -d_mu*(abs(omega_R)-omega_F_cr) )
          F_DFR = - (+(omega_R>0) or -(omega_R<0)) * mu * F_DR
        
        # Other Forces & Moments: F_FRm, M_FR, M_FRm
        F_DFRm   = -F_DFR
        M_DFR    = -F_DFR  * ( D * cos( R  ) - H  * sin( R  ) )
        if L!=1:
          M_DFRm = -F_DFRm * ( D * cos( Rm ) + Hm * sin( Rm ) )
        
      # Check Point: omega==0 -> Do nothing
        
    # Check Point: F_VR<=0 -> Do nothing
    
    
    #} RIGHT end
    #===========================================================================
  
    #{ Assemble Accel
    AccelF  =  (F_DFL  + F_DFR ) / M
    AccelM  =  (M_DFL  + M_DFR ) / I
    AccelFm = -(F_DFL  + F_DFR ) / Mm
    AccelMm =  (M_DFLm + M_DFRm) / Im
    
    # Add Accel on (K,L)
    IndexW = Core['Index'][K][L]*6
    Accel[IndexW+3] += AccelF # DDW
    Accel[IndexW+5] += AccelM # DDR
  
    # Add Accel on (K,L-1)
    
    if L!=1: # Normal Blocks
      IndexWm = Core['Index'][K][L-1]*6
      Accel[IndexWm+3] += AccelFm # DDW
      Accel[IndexWm+5] += AccelMm # DDR
   
    elif Core['Flag_CSB']==True: 
      pass # No vertical accel and moment on CSB
  
    pass # if Base: Do nothing
    #} END OF ASSEMBLE ACCEL
    

    #{ VERBOSE
    if '/ForceDF/' in VERBOSE:
      Verbose('-'*80)
      Verbose('Force_Dowel_F(), Caculation Verification')
      Verbose('Block (%d,%d)'%(K,L))
      Verbose('U, DU, W, DW, R, DR  = %e, %e, %e, %e, %e, %e'%(U,DU,W,DW,R,DR))
      Verbose('Um,DUm,Wm,DWm,Rm,DRm = %e, %e, %e, %e, %e, %e'%(Um,DUm,Wm,DWm,Rm,DRm))
      Verbose('omega  = %e, %e'%(omega_L,omega_R))
      Verbose('FDF  = %e, %e'%(F_DFL,F_DFR))
      Verbose('MDF  = %e, %e, %e, %e'%(M_DFL,M_DFR,M_DFLm,M_DFRm))
      Verbose('Accel(%d,%d) = %e, %e'%(K,L,AccelF,AccelM))
      Verbose('Accel[%d] = %e'%(IndexW+3,Accel[IndexW+3]))
      Verbose('Accel[%d] = %e'%(IndexW+5,Accel[IndexW+5]))
      if L!=1:
        Verbose('Accel(%d,%d) = %e, %e'%(K,L-1,AccelFm,AccelMm))
        Verbose('Accel[%d] = %e'%(IndexWm+3,Accel[IndexWm+3]))
        Verbose('Accel[%d] = %e'%(IndexWm+5,Accel[IndexWm+5]))
    #}
    
    #{ POST
    if POST==True:

      # FILENAME
      Filename = FO_DIR+'/(%2d,%2d)_DF.csv'%(K,L)

      #{ HEADER
      try: fo = open(Filename,'r')
      except: # File doesn't exist -> This is the first time. Make header.
        temp  = 't,'
        temp += 'OmegaL,F_DFL,M_DFL,M_DFLm,'
        temp += 'OmegaR,F_DFR,M_DFR,M_DFRm\n'
        fo = open(Filename,'w')
        fo.write(temp)
        fo.close()
      #}

      #{ WRITE
      temp  = '%e, '%t
      temp += '%e, %e, %e, %e, '%( omega_L, F_DFL, M_DFL, M_DFLm )
      temp += '%e, %e, %e, %e  '%( omega_R, F_DFR, M_DFR, M_DFRm )
      temp += '\n'
      fo = open(Filename,'a')
      fo.write(temp)
      fo.close()
      #}
    #}


  
  #}  END OF FRICTION FORCES
  ##############################################################################
  
  return Accel
#}
#-------------------------------------------------------------------------
