VERSION_NAME = 'SAPCOR_v05a'

# MODULE NAME
MODULE_NAME_THIS__FORCE_BLOCK_H = VERSION_NAME+'_Force_Block_H'

###############################################################################
#region                   SYSTEM INITIALIZATION

import glob

# -------------------------------------------------------------------
# MODULE NAMES TO BE IMPORTED
MODULE_NAMES_FORCE_BLOCK_H = []
MODULE_NAMES_FORCE_BLOCK_H.append('Misc')
MODULE_NAMES_FORCE_BLOCK_H.append('Input')
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Import Modules
for Module_Name in MODULE_NAMES_FORCE_BLOCK_H:
  FileName = VERSION_NAME+'_'+Module_Name+'_Build*.py'
  FList = glob.glob(FileName)
  assert len(FList)!=0,'ERROR: Module not found : '+FileName
  FList.sort()
  FileName = FList[-1][0:-3]
  exec('from '+FileName+' import *')
# -------------------------------------------------------------------




from numpy import loadtxt
from numpy import array,zeros,hstack

#endregion           END OF SYSTEM INITIALIZATION
###############################################################################


pass


###############################################################################
def Force_Block_H (w,t,Core,K,L,Accel,POST=False,FO_DIR=FO_DIR):

  #region Comment of Variables
  #----------------------------------------------------------------------------
  #
  # Z_lt_l : Z coord Left  Top    corner of (K,L)
  # Z_lb_l : Z coord Left  Bottom corner of (K,L)
  # Z_rt_l : Z coord Right Top    corner of (K,L)
  # Z_rb_l : Z coord Right Bottom corner of (K,L)
  #
  # Z_lt_m : Z coord Left  Top    corner of (K+1,M) or (K-1,M)
  # Z_lb_m : Z coord Left  Bottom corner of (K+1,M) or (K-1,M)
  # Z_rt_m : Z coord Right Top    corner of (K+1,M) or (K-1,M)
  # Z_rb_m : Z coord Right Bottom corner of (K+1,M) or (K-1,M)
  #
  # m : 2nd index of (K+1,M) or (K-1,M)
  #
  # Accel[0] : DDU (Fu)
  # Accel[1] : DDW (Fw)
  # Accel[2] : DDR (M)
  #----------------------------------------------------------------------------
  #endregion

  # VERBOSE
  if '/ForceH/' in VERBOSE:
    Verbose('='*80)
    Verbose('<FORCE_BLOCK_H>')
    Verbose('time=%f  K,L=%d,%d'%(t,K,L))

  #region                             INIT

  # Accel[IndexW:IndexW+6]==[DU,DDU,DW,DDW,DR,DDR]
  IndexW = Core['Index'][K][L]*6


  #region Get Block Prop of (K,L)
  BN = Core['Array'][K][L]
  BT = BlockTypes[BN]
  B = BT['b']
  H = BT['h']
  M = BT['M']
  I = BT['I']
  Kh = BT['Kh']
  Ch = BT['Ch']
  Fixed = BT['Fixed']
  #endregion
  
  # Get current state of (K,L)
  U,DU,W,DW,R,DR = GetSolFromVect(Core,w,K,L) # U,DU,W,DW,R,DR
  
  #region Corner Coords of (K,L)
  temp = GetBlockCoords (Core,K,L,U,W,R) # ret: [Center,LT,LB,RT,RB]
  Z_lt_l = temp[1][1] # Z coord of LT (cf. temp[1][0]: X coord)
  Z_lb_l = temp[2][1]
  Z_rt_l = temp[3][1]
  Z_rb_l = temp[4][1]
  #endregion
  
  #endregion                      END OF INIT

  # VERBOSE - (K,L)
  if '/ForceH/' in VERBOSE:
    Verbose('-'*80)
    Verbose('Block = (%d,%d)'%(K,L))
    Verbose('Block Type Name = %s'%BN,NewLine=False)
    Verbose(', H  = %f, B  = %f'%(H,B))
    Verbose('U, DU, R, DR  = %e, %e, %e, %e'%(U,DU,R,DR))



  # ==========================================================================
  #region                         RIGHT SIDE



  # ----------------------------------------------------------------------------
  #region                        COMMENTS


  # K=0:LReflector, K=1~M:Blocks, K=M+1:RReflector
  # Right Side Force is only for K=1~M-1
  # Therefore, need to check K<=M-1


  #endregion                END OF COMMENTS
  # ----------------------------------------------------------------------------
  #region                            K < M



  if K<Core['M']:
    # --------------------------------------------------------------------------
    #region                     (1) UPPER RIGHT



    # --------------------------------------------------------------------------
    #region { Init
    Epsilon = 0.
    DEpsilon = 0.
    Force = 0.
    Moment = 0.
    #endregion }
  
    #region { Start of Find (K+1,m) Block
    FlagFound = False
    len_m = len(Core['Array'][K+1])
    for m in xrange(1,len_m+1): # L starts from 1 to m
  
      # Get current state of (K+1,m)
      Um,DUm,Wm,DWm,Rm,DRm = GetSolFromVect(Core,w,K+1,m) # U,DU,W,DW,R,DR
      
      # Corner Coords of (K+1,m)  
      temp = GetBlockCoords (Core,K+1,m,Um,Wm,Rm) # [Center,LU,LD,RU,RD]
      Z_lt_m = temp[1][1]
      Z_lb_m = temp[2][1]
      Z_rt_m = temp[3][1]
      Z_rb_m = temp[4][1]
      
      # Lambda
      Lambda = ( Z_lt_m - Z_rt_l ) / ( Z_lt_m - Z_lb_m )
      
      # Check 0<=Lamda<1
      if 0 <= Lambda < 1:
        FlagFound = True
        if m == len_m: FlagTopBlock = True # Top Block
        else: FlagTopBlock = False
        break
    #endregion } End of Find (K+1,m) Block
        
    #region { Start of Get Force and Moment
    if FlagFound == True:
  
      # Get Block Prop of (K+1,m)
      BN = Core['Array'][K+1][m]
      BT = BlockTypes[BN]
      Bm = BT['b']
      Hm = BT['h']
    
      # Spring Contraction
      Epsilon  =  U  + (              H *sin(R ) - B *(1-cos(R )) )
      Epsilon += -Um - ( (1-2*Lambda)*Hm*sin(Rm) + Bm*(1-cos(Rm)) )
      if FlagTopBlock == False: Epsilon -= Delta
      else: Epsilon -= DeltaTop
  
      # Calculate Remains Only If Penetration Occurs
      if Epsilon > 0:
  
        # Spring Velocity    
        DEpsilon  =  DU  + (              H *cos(R ) - B *sin(R ) ) * DR
        DEpsilon += -DUm - ( (1-2*Lambda)*Hm*cos(Rm) + Bm*sin(Rm) ) * DRm
      
        # Force and Moment
        Force  = -Kh * Epsilon - Ch * DEpsilon
        Moment =  Force * ( H*cos(R) - B*sin(R) )


        #region ------------- ELIMINATE STICKING FORCE ------------------------------
        # Right Side, Force Direction : <-- (Left)
        # Eliminate Positive Force
        #-------------------------------------------
        if FLAG_NoStickForce == True:
          if Force > 0. : Force = Moment = 0.
        #endregion -------------------------------------------

      
        # Acceleration
        Accel[IndexW+1] += Force/M  # DDU
        Accel[IndexW+5] += Moment/I # DDR
  
    #endregion } End of Get Force and Moment
    
    #region { VERBOSE - UPPER RIGHT
    if '/ForceH/' in VERBOSE:
      ContactConfig = 'Upper & Right'
      
      if  'Right' in ContactConfig: SignK = 1
      elif 'Left' in ContactConfig: SignK = -1
      else: ERROR('Error: ContactConfig("%s") was invalid.'%ContactConfig)
      
      Verbose('-'*80)
      Verbose('Force_Block_V() - %s, Caculation Verification'%ContactConfig)
      if FlagFound == True:
        Verbose('Matching Block = (%d,%d)'%(K+SignK,m))
        Verbose('Block Type Name = %s'%BN,NewLine=False)
        Verbose(', Lambda = %f, Hm = %f, Bm = %f'%(Lambda,Hm,Bm))
        Verbose('Um,DUm,Rm,DRm = %e, %e, %e, %e'%(Um,DUm,Rm,DRm))
        Verbose('Epsilon = %e'%Epsilon,NewLine=False)
        if Epsilon > 0:
          Verbose(', dEpsilon = %e'%DEpsilon)
          Verbose('Force, Moment = %e, %e'%(Force,Moment))
          Verbose('Accel = %e, %e'%(Force/M,Moment/I))
        else:
          Verbose(' <= 0   --->  No Contact')
      else:
        Verbose('  Mathing block was not found.')
    #endregion }

    #region { POST TR
    if POST==True:

      # FILENAME
      Filename = FO_DIR+'/(%2d,%2d)_TR.csv'%(K,L)

      #{ HEADER
      try: fo = open(Filename,'r')
      except: # File doesn't exist -> This is the first time. Make header.
        temp  = 't,'
        temp += 'Epsil,F_TRk,dEpsil,F_TRc,F_TR,M_TR\n'
        fo = open(Filename,'w')
        fo.write(temp)
        fo.close()
      #}

      #{ F_TRk, F_TRc
      F_TRk = F_TRc = 0.
      if Epsilon > 0:
        F_TRk = -Kh *  Epsilon
        F_TRc = -Ch * DEpsilon
      #}

      #{ WRITE
      temp  = '%e, '%t
#              1   2   3   4   5   6    : 1        2      3         4      5      6
      temp += '%e, %e, %e, %e, %e, %e\n'%(Epsilon, F_TRk, DEpsilon, F_TRc, Force, Moment )
      fo = open(Filename,'a')
      fo.write(temp)
      fo.close()
      #}
    #endregion }
    # --------------------------------------------------------------------------



    #endregion                END OF (1) UPPER RIGHT
    # --------------------------------------------------------------------------
    #region                   (2) UPPER-MIDDLE RIGHT



    # --------------------------------------------------------------------------

    #{ Init
    Epsilon = 0.
    DEpsilon = 0.
    Force = 0.
    Moment = 0.
    #}
  
    #{ Start of Find (K+1,m) Block
    FlagFound = False
    len_m = len(Core['Array'][K+1])
    for m in xrange(1,len_m+1): # L starts from 1 to m
  
      # Get current state of (K+1,m)
      Um,DUm,Wm,DWm,Rm,DRm = GetSolFromVect(Core,w,K+1,m) # U,DU,W,DW,R,DR
      
      # Corner Coords of (K+1,m)  
      temp = GetBlockCoords (Core,K+1,m,Um,Wm,Rm) # [Center,LU,LD,RU,RD]
      Z_lt_m = temp[1][1]
      Z_lb_m = temp[2][1]
      Z_rt_m = temp[3][1]
      Z_rb_m = temp[4][1]
   
      # Lambda
      Lambda = ( Z_rt_l - Z_lt_m ) / ( Z_rt_l - Z_rb_l )
      
      # Check 0<Lamda<1
      if 0 < Lambda < 1:
        FlagFound = True
        if m == len_m: FlagTopBlock = True # Top Block
        else: FlagTopBlock = False
        break
    #} End of Find (K+1,m) Block
        
    #{ Start of Get Force and Moment
    if FlagFound == True:
  
      # Get Block Prop of (K+1,m)
      BN = Core['Array'][K+1][m]
      BT = BlockTypes[BN]
      Bm = BT['b']
      Hm = BT['h']
    
      # Spring Contraction
      Epsilon  =  U  + ( (1-2*Lambda)*H *sin(R ) - B *(1-cos(R )) )
      Epsilon += -Um - (              Hm*sin(Rm) + Bm*(1-cos(Rm)) )
      if FlagTopBlock == False: Epsilon -= Delta
      else: Epsilon -= DeltaTop
  
      # Calculate Remains Only If Penetration Occurs
      if Epsilon > 0:
  
        # Spring Velocity    
        DEpsilon  =  DU  + ( (1-2*Lambda)*H *cos(R ) - B *sin(R ) ) * DR
        DEpsilon += -DUm - (              Hm*cos(Rm) + Bm*sin(Rm) ) * DRm
      
        # Force and Moment
        Force  = -Kh * Epsilon - Ch * DEpsilon
        Moment =  Force * ( (1-2*Lambda)*H*cos(R) - B*sin(R) )
      

        #-------------------------------------------
        # ELIMINATE STICKING FORCE
        # Right Side, Force Direction : <-- (Left)
        # Eliminate Positive Force
        #-------------------------------------------
        if FLAG_NoStickForce == True:
          if Force > 0. : Force = Moment = 0.
        #-------------------------------------------


        # Acceleration
        Accel[IndexW+1] += Force/M  # DDU
        Accel[IndexW+5] += Moment/I # DDR
  
    #} End of Get Force and Moment
    
    #{ VERBOSE - UPPER-MIDDLE RIGHT
    if '/ForceH/' in VERBOSE:
      ContactConfig = 'Upper & Middle Right'
      
      if  'Right' in ContactConfig: SignK = 1
      elif 'Left' in ContactConfig: SignK = -1
      else: ERROR('Error: ContactConfig("%s") was invalid.'%ContactConfig)
      
      Verbose('-'*80)
      Verbose('Force_Block_V() - %s, Caculation Verification'%ContactConfig)
      if FlagFound == True:
        Verbose('Matching Block = (%d,%d)'%(K+SignK,m))
        Verbose('Block Type Name = %s'%BN,NewLine=False)
        Verbose(', Hm = %f, Bm = %f'%(Hm,Bm))
        Verbose('Um,DUm,Rm,DRm = %e, %e, %e, %e'%(Um,DUm,Rm,DRm))
        Verbose('Epsilon = %e'%Epsilon,NewLine=False)
        if Epsilon > 0:
          Verbose(', dEpsilon = %e'%DEpsilon)
          Verbose('Force, Moment = %e, %e'%(Force,Moment))
          Verbose('Accel = %e, %e'%(Force/M,Moment/I))
        else:
          Verbose(' <= 0   --->  No Contact')
      else:
        Verbose('  Mathing block was not found.')
    #}

    #{ POST MR
    if POST==True:

      # FILENAME
      Filename = FO_DIR+'/(%2d,%2d)_MR.csv'%(K,L)

      #{ HEADER
      try: fo = open(Filename,'r')
      except: # File doesn't exist -> This is the first time. Make header.
        temp  = 't,'
        temp += 'Epsil,F_MRk,dEpsil,F_MRc,F_MR,M_MR\n'
        fo = open(Filename,'w')
        fo.write(temp)
        fo.close()
      #}

      #{ F_MRk, F_MRc
      F_MRk = F_MRc = 0.
      if Epsilon > 0:
        F_MRk = -Kh *  Epsilon
        F_MRc = -Ch * DEpsilon
      #}

      #{ WRITE
      temp  = '%e, '%t
#              1   2   3   4   5   6    : 1        2      3         4      5      6
      temp += '%e, %e, %e, %e\n'%( Epsilon, F_MRk, DEpsilon, F_MRc, Force, Moment )
      fo = open(Filename,'a')
      fo.write(temp)
      fo.close()
      #}
    #}
    
    # --------------------------------------------------------------------------



    #endregion            END OF (2) UPPER-MIDDLE RIGHT
    # --------------------------------------------------------------------------
    #region                  (3) LOWER-MIDDLE RIGHT



    # --------------------------------------------------------------------------

    #{ Init
    Epsilon = 0.
    DEpsilon = 0.
    Force = 0.
    Moment = 0.
    #}
  
    #{ Start of Find (K+1,m) Block
    FlagFound = False
    len_m = len(Core['Array'][K+1])
    for m in xrange(1,len_m+1): # L starts from 1 to m
  
      # Get current state of (K+1,m)
      Um,DUm,Wm,DWm,Rm,DRm = GetSolFromVect(Core,w,K+1,m) # U,DU,W,DW,R,DR
      
      # Corner Coords of (K+1,m)  
      temp = GetBlockCoords (Core,K+1,m,Um,Wm,Rm) # [Center,LU,LD,RU,RD]
      Z_lt_m = temp[1][1]
      Z_lb_m = temp[2][1]
      Z_rt_m = temp[3][1]
      Z_rb_m = temp[4][1]
   
      # Lambda
      Lambda = ( Z_rt_l - Z_lb_m ) / ( Z_rt_l - Z_rb_l )
      
      # Check 0<Lamda<1
      if 0 < Lambda < 1:
        FlagFound = True
        if m == len_m: FlagTopBlock = True # Top Block
        else: FlagTopBlock = False
        break
    #} End of Find (K+1,m) Block
        
    #{ Start of Get Force and Moment
    if FlagFound == True:

      # Get Block Prop of (K+1,m)
      BN = Core['Array'][K+1][m]
      BT = BlockTypes[BN]
      Bm = BT['b']
      Hm = BT['h']
    
      # Spring Contraction
      Epsilon  =  U  + ( (1-2*Lambda)*H *sin(R ) - B *(1-cos(R )) )
      Epsilon += -Um - (            - Hm*sin(Rm) + Bm*(1-cos(Rm)) )
      if FlagTopBlock == False: Epsilon -= Delta
      else: Epsilon -= DeltaTop
  
      # Calculate Remains Only If Penetration Occurs
      if Epsilon > 0:
  
        # Spring Velocity    
        DEpsilon  =  DU  + ( (1-2*Lambda)*H *cos(R ) - B *sin(R ) ) * DR
        DEpsilon += -DUm - (            - Hm*cos(Rm) + Bm*sin(Rm) ) * DRm
      
        # Force and Moment
        Force  = -Kh * Epsilon - Ch * DEpsilon
        Moment =  Force * ( (1-2*Lambda)*H*cos(R) - B*sin(R) )
      

        #-------------------------------------------
        # ELIMINATE STICKING FORCE
        # Right Side, Force Direction : <-- (Left)
        # Eliminate Positive Force
        #-------------------------------------------
        if FLAG_NoStickForce == True:
          if Force > 0. : Force = Moment = 0.
        #-------------------------------------------


        # Acceleration
        Accel[IndexW+1] += Force/M  # DDU
        Accel[IndexW+5] += Moment/I # DDR
  
    #} End of Get Force and Moment

    #{ VERBOSE - LOWER-MIDDLE RIGHT
    if '/ForceH/' in VERBOSE:
      ContactConfig = 'Lower & Middle Right'
      
      if  'Right' in ContactConfig: SignK = 1
      elif 'Left' in ContactConfig: SignK = -1
      else: ERROR('Error: ContactConfig("%s") was invalid.'%ContactConfig)
      
      Verbose('-'*80)
      Verbose('Force_Block_V() - %s, Caculation Verification'%ContactConfig)
      if FlagFound == True:
        Verbose('Matching Block = (%d,%d)'%(K+SignK,m))
        Verbose('Block Type Name = %s'%BN,NewLine=False)
        Verbose(', Hm = %f, Bm = %f'%(Hm,Bm))
        Verbose('Um,DUm,Rm,DRm = %e, %e, %e, %e'%(Um,DUm,Rm,DRm))
        Verbose('Epsilon = %e'%Epsilon,NewLine=False)
        if Epsilon > 0:
          Verbose(', dEpsilon = %e'%DEpsilon)
          Verbose('Force, Moment = %e, %e'%(Force,Moment))
          Verbose('Accel = %e, %e'%(Force/M,Moment/I))
        else:
          Verbose(' <= 0   --->  No Contact')
      else:
        Verbose('  Mathing block was not found.')
    #}

    #{ POST NR
    if POST==True:

      # FILENAME
      Filename = FO_DIR+'/(%2d,%2d)_NR.csv'%(K,L)

      #{ HEADER
      try: fo = open(Filename,'r')
      except: # File doesn't exist -> This is the first time. Make header.
        temp  = 't,'
        temp += 'Epsil,F_NRk,dEpsil,F_NRc,F_NR,M_NR\n'
        fo = open(Filename,'w')
        fo.write(temp)
        fo.close()
      #}

      #{ F_NRk, F_NRc
      F_NRk = F_NRc = 0.
      if Epsilon > 0:
        F_NRk = -Kh *  Epsilon
        F_NRc = -Ch * DEpsilon
      #}

      #{ WRITE
      temp  = '%e, '%t
#              1   2   3   4   5   6    : 1        2      3         4      5      6
      temp += '%e, %e, %e, %e\n'%( Epsilon, F_NRk, DEpsilon, F_NRc, Force, Moment )
      fo = open(Filename,'a')
      fo.write(temp)
      fo.close()
      #}
    #}
    
    # --------------------------------------------------------------------------



    #endregion           END OF (3) LOWER-MIDDLE RIGHT
    # --------------------------------------------------------------------------
    #region                    (4) LOWER RIGHT



    # --------------------------------------------------------------------------
    #region                         INIT


    Epsilon = 0.
    DEpsilon = 0.
    Force = 0.
    Moment = 0.


    #endregion                    END OF INIT
    # --------------------------------------------------------------------------
    #region                   Find (K+1,m) Block


    FlagFound = False
    len_m = len(Core['Array'][K+1])
    for m in xrange(1,len_m+1): # L starts from 1 to m
  
      # Get current state of (K+1,m)
      Um,DUm,Wm,DWm,Rm,DRm = GetSolFromVect(Core,w,K+1,m) # U,DU,W,DW,R,DR
      
      # Corner Coords of (K+1,m)  
      temp = GetBlockCoords (Core,K+1,m,Um,Wm,Rm) # [Center,LU,LD,RU,RD]
      Z_lt_m = temp[1][1]
      Z_lb_m = temp[2][1]
      Z_rt_m = temp[3][1]
      Z_rb_m = temp[4][1]
   
      # Lambda
      Lambda = ( Z_lt_m - Z_rb_l ) / ( Z_lt_m - Z_lb_m )
      
      # Check 0<Lamda<=1
      if 0 < Lambda <= 1:
        FlagFound = True
        if m == len_m: FlagTopBlock = True # Top Block
        else: FlagTopBlock = False
        break


    #endregion            End of Find (K+1,m) Block
    # --------------------------------------------------------------------------
    #region                  Get Force and Moment



    if FlagFound == True:

      # Get Block Prop of (K+1,m)
      BN = Core['Array'][K+1][m]
      BT = BlockTypes[BN]
      Bm = BT['b']
      Hm = BT['h']
    
      # Spring Contraction
      Epsilon  =  U  + (            - H *sin(R ) - B *(1-cos(R )) )
      Epsilon += -Um - ( (1-2*Lambda)*Hm*sin(Rm) + Bm*(1-cos(Rm)) )
      if FlagTopBlock == False: Epsilon -= Delta
      else: Epsilon -= DeltaTop
  
      # Calculate Remains Only If Penetration Occurs
      if Epsilon > 0:
        #----------------------------------------------------------------------
        #region Spring Velocity
        DEpsilon  =  DU  + (            - H *cos(R ) - B *sin(R ) ) * DR
        DEpsilon += -DUm - ( (1-2*Lambda)*Hm*cos(Rm) + Bm*sin(Rm) ) * DRm
        #endregion
        #----------------------------------------------------------------------
        #region Force and Moment
        Force  = -Kh * Epsilon - Ch * DEpsilon
        Moment =  Force * (-H*cos(R) - B*sin(R) )
        #endregion
        #----------------------------------------------------------------------
        #region ELIMINATE STICKING FORCE
        #-------------------------------------------
        # Right Side, Force Direction : <-- (Left)
        # Eliminate Positive Force
        #-------------------------------------------
        if FLAG_NoStickForce == True:
          if Force > 0. : Force = Moment = 0.
        #endregion
        #----------------------------------------------------------------------
        #region                        Accel


        Accel[IndexW+1] += Force/M  # DDU
        Accel[IndexW+5] += Moment/I # DDR


        #endregion                  END OF ACCEL
        #----------------------------------------------------------------------



    #endregion           End of Get Force and Moment
    # --------------------------------------------------------------------------
    #region                  VERBOSE - LOWER RIGHT


    if '/ForceH/' in VERBOSE:
      ContactConfig = 'Lower & Right'
      
      if  'Right' in ContactConfig: SignK = 1
      elif 'Left' in ContactConfig: SignK = -1
      else: ERROR('Error: ContactConfig("%s") was invalid.'%ContactConfig)
      
      Verbose('-'*80)
      Verbose('Force_Block_V() - %s, Caculation Verification'%ContactConfig)
      if FlagFound == True:
        Verbose('Matching Block = (%d,%d)'%(K+SignK,m))
        Verbose('Block Type Name = %s'%BN,NewLine=False)
        Verbose(', Hm = %f, Bm = %f'%(Hm,Bm))
        Verbose('Um,DUm,Rm,DRm = %e, %e, %e, %e'%(Um,DUm,Rm,DRm))
        Verbose('Epsilon = %e'%Epsilon,NewLine=False)
        if Epsilon > 0:
          Verbose(', dEpsilon = %e'%DEpsilon)
          Verbose('Force, Moment = %e, %e'%(Force,Moment))
          Verbose('Accel = %e, %e'%(Force/M,Moment/I))
        else:
          Verbose(' <= 0   --->  No Contact')
      else:
        Verbose('  Mathing block was not found.')


    #endregion            END OF VERBOSE - LOWER RIGHT
    # --------------------------------------------------------------------------
    #region                           POST BR


    if POST==True:

      # FILENAME
      Filename = FO_DIR+'/(%2d,%2d)_BR.csv'%(K,L)

      #region { HEADER
      try: fo = open(Filename,'r')
      except: # File doesn't exist -> This is the first time. Make header.
        temp  = 't,'
        temp += 'Epsil,F_BRk,dEpsil,F_BRc,F_BR,M_BR\n'
        fo = open(Filename,'w')
        fo.write(temp)
        fo.close()
      #endregion }

      #region { F_BRk, F_BRc
      F_BRk = F_BRc = 0.
      if Epsilon > 0:
        F_BRk = -Kh *  Epsilon
        F_BRc = -Ch * DEpsilon
      #endregion }

      #region { WRITE
      temp  = '%e, '%t
#              1   2   3   4   5   6    : 1        2      3         4      5      6
      temp += '%e, %e, %e, %e\n'%( Epsilon, F_BRk, DEpsilon, F_BRc, Force, Moment )
      fo = open(Filename,'a')
      fo.write(temp)
      fo.close()
      #endregion }
    pass #EOIF


    #endregion                    END OF POST BR
    # --------------------------------------------------------------------------



    #endregion              END OF (4) LOWER RIGHT
    # --------------------------------------------------------------------------



  #endregion                      END OF K < M
  #----------------------------------------------------------------------------
  #region                         K >= M -> SKIP


  else:
    # VERBOSE - LOWER RIGHT
    if '/ForceH/' in VERBOSE:
      Verbose('-'*80)
      Verbose('K(%d)>=M(%d)  --->  No ForceH on Right Side'%(K,Core['M']))


  #endregion                  END OF K >= M -> SKIP
  # ----------------------------------------------------------------------------



  #endregion                  END OF RIGHT SIDE
  # ==========================================================================
  #region                           LEFT SIDE



  # K=0:LReflector, K=1~M:Blocks, K=M+1:RReflector
  # Left Side Force is only for K=2~M
  # Therfore, need to check K>=2
  if K>1:

    # --------------------------------------------------------------------------
    #{ (1) UPPER LEFT
    # --------------------------------------------------------------------------

    #{ Init
    Epsilon = 0.
    DEpsilon = 0.
    Force = 0.
    Moment = 0.
    #}
  
    #{ Start of Find (K-1,m) Block
    FlagFound = False
    len_m = len(Core['Array'][K-1])
    for m in xrange(1,len_m+1): # L starts from 1 to m
  
      # Get current state of (K-1,m)
      Um,DUm,Wm,DWm,Rm,DRm = GetSolFromVect(Core,w,K-1,m) # U,DU,W,DW,R,DR
      
      # Corner Coords of (K-1,m)  
      temp = GetBlockCoords (Core,K-1,m,Um,Wm,Rm) # [Center,LU,LD,RU,RD]
      Z_lt_m = temp[1][1]
      Z_lb_m = temp[2][1]
      Z_rt_m = temp[3][1]
      Z_rb_m = temp[4][1]
      
      # Lambda
      # For Left Case, Lambda needs only to change 
      # subscripts l(Left) and r(Right) with each other from the Right Case.
      # Subscript t,b,l(meand column index,L), and m should not be modified.
      Lambda = ( Z_rt_m - Z_lt_l ) / ( Z_rt_m - Z_rb_m )
      
      # Check 0<=Lamda<1
      if 0 <= Lambda < 1:
        FlagFound = True
        if m == len_m: FlagTopBlock = True # Top Block
        else: FlagTopBlock = False
        break
    #} End of Find (K-1,m) Block
        
    #{ Start of Get Force and Moment
    if FlagFound == True:
  
      # Get Block Prop of (K-1,m)
      BN = Core['Array'][K-1][m]
      BT = BlockTypes[BN]
      Bm = BT['b']
      Hm = BT['h']
    
      # Spring Contraction
      Epsilon  = -U  - (              H *sin(R ) + B *(1-cos(R )) )
      Epsilon +=  Um + ( (1-2*Lambda)*Hm*sin(Rm) - Bm*(1-cos(Rm)) )
      if FlagTopBlock == False: Epsilon -= Delta
      else: Epsilon -= DeltaTop
  
      # Calculate Remains Only If Penetration Occurs
      if Epsilon > 0:
  
        # Spring Velocity    
        DEpsilon  = -DU  - (              H *cos(R ) + B *sin(R ) ) * DR
        DEpsilon +=  DUm + ( (1-2*Lambda)*Hm*cos(Rm) - Bm*sin(Rm) ) * DRm
      
        # Force and Moment
        Force  =  Kh * Epsilon + Ch * DEpsilon
        Moment =  Force * ( H*cos(R) + B*sin(R) )


        #-------------------------------------------
        # ELIMINATE STICKING FORCE
        # Left Side, Force Direction : --> (Right)
        # Eliminate Negative Force
        #-------------------------------------------
        if FLAG_NoStickForce == True:
          if Force < 0. : Force = Moment = 0.
        #-------------------------------------------


        # Acceleration
        Accel[IndexW+1] += Force/M  # DDU
        Accel[IndexW+5] += Moment/I # DDR
  
    #} End of Get Force and Moment

    #{ VERBOSE - UPPER LEFT
    if '/ForceH/' in VERBOSE:
      ContactConfig = 'Upper & Left'
      
      if  'Right' in ContactConfig: SignK = 1
      elif 'Left' in ContactConfig: SignK = -1
      else: ERROR('Error: ContactConfig("%s") was invalid.'%ContactConfig)
      
      Verbose('-'*80)
      Verbose('Force_Block_V() - %s, Caculation Verification'%ContactConfig)
      if FlagFound == True:
        Verbose('Matching Block = (%d,%d)'%(K+SignK,m))
        Verbose('Block Type Name = %s'%BN,NewLine=False)
        Verbose(', Hm = %f, Bm = %f'%(Hm,Bm))
        Verbose('Um,DUm,Rm,DRm = %e, %e, %e, %e'%(Um,DUm,Rm,DRm))
        Verbose('Epsilon = %e'%Epsilon,NewLine=False)
        if Epsilon > 0:
          Verbose(', dEpsilon = %e'%DEpsilon)
          Verbose('Force, Moment = %e, %e'%(Force,Moment))
          Verbose('Accel = %e, %e'%(Force/M,Moment/I))
        else:
          Verbose(' <= 0   --->  No Contact')
      else:
        Verbose('  Mathing block was not found.')
    #}

    #{ POST TL
    if POST==True:

      # FILENAME
      Filename = FO_DIR+'/(%2d,%2d)_TL.csv'%(K,L)

      #{ HEADER
      try: fo = open(Filename,'r')
      except: # File doesn't exist -> This is the first time. Make header.
        temp  = 't,'
        temp += 'Epsil,F_TLk,dEpsil,F_TLc,F_TL,M_TL\n'
        fo = open(Filename,'w')
        fo.write(temp)
        fo.close()
      #}

      #{ F_TLk, F_TLc
      F_TLk = F_TLc = 0.
      if Epsilon > 0:
        F_TLk = +Kh *  Epsilon
        F_TLc = +Ch * DEpsilon
      #}

      #{ WRITE
      temp  = '%e, '%t
#              1   2   3   4   5   6    : 1        2      3         4      5      6
      temp += '%e, %e, %e, %e, %e, %e\n'%(Epsilon, F_TLk, DEpsilon, F_TLc, Force, Moment )
      fo = open(Filename,'a')
      fo.write(temp)
      fo.close()
      #}
    #}
       
    # --------------------------------------------------------------------------
    #} END OF (1) UPPER LEFT
    # --------------------------------------------------------------------------
    
    



    
    # --------------------------------------------------------------------------
    #{ (2) UPPER-MIDDLE LEFT
    # --------------------------------------------------------------------------
  
    #{ Init
    Epsilon = 0.
    DEpsilon = 0.
    Force = 0.
    Moment = 0.
    #}

    #{ Start of Find (K-1,m) Block
    FlagFound = False
    len_m = len(Core['Array'][K-1])
    for m in xrange(1,len_m+1): # L starts from 1 to m
  
      # Get current state of (K-1,m)
      Um,DUm,Wm,DWm,Rm,DRm = GetSolFromVect(Core,w,K-1,m) # U,DU,W,DW,R,DR
      
      # Corner Coords of (K-1,m)  
      temp = GetBlockCoords (Core,K-1,m,Um,Wm,Rm) # [Center,LU,LD,RU,RD]
      Z_lt_m = temp[1][1]
      Z_lb_m = temp[2][1]
      Z_rt_m = temp[3][1]
      Z_rb_m = temp[4][1]
   
      # Lambda
      # For Left Case, Lambda needs only to change 
      # subscripts l(Left) and r(Right) with each other from the Right Case.
      # Subscript t,b,l(meand column index,L), and m should not be modified.
      Lambda = ( Z_lt_l - Z_rt_m ) / ( Z_lt_l - Z_lb_l )
      
      # Check 0<Lamda<1
      if 0 < Lambda < 1:
        FlagFound = True
        if m == len_m: FlagTopBlock = True # Top Block
        else: FlagTopBlock = False
        break
    #} End of Find (K-1,m) Block
        
    #{ Start of Get Force and Moment
    if FlagFound == True:
  
      # Get Block Prop of (K-1,m)
      BN = Core['Array'][K-1][m]
      BT = BlockTypes[BN]
      Bm = BT['b']
      Hm = BT['h']
    
      # Spring Contraction
      Epsilon  = -U  - ( (1-2*Lambda)*H *sin(R ) + B *(1-cos(R )) )
      Epsilon +=  Um + (              Hm*sin(Rm) - Bm*(1-cos(Rm)) )
      if FlagTopBlock == False: Epsilon -= Delta
      else: Epsilon -= DeltaTop
  
      # Calculate Remains Only If Penetration Occurs
      if Epsilon > 0:
  
        # Spring Velocity    
        DEpsilon  = -DU  - ( (1-2*Lambda)*H *cos(R ) + B *sin(R ) ) * DR
        DEpsilon +=  DUm + (              Hm*cos(Rm) - Bm*sin(Rm) ) * DRm
      
        # Force and Moment
        Force  =  Kh * Epsilon + Ch * DEpsilon
        Moment =  Force * ( (1-2*Lambda)*H*cos(R) + B*sin(R) )
      

        #-------------------------------------------
        # ELIMINATE STICKING FORCE
        # Left Side, Force Direction : --> (Right)
        # Eliminate Negative Force
        #-------------------------------------------
        if FLAG_NoStickForce == True:
          if Force < 0. : Force = Moment = 0.
        #-------------------------------------------


        # Acceleration
        Accel[IndexW+1] += Force/M  # DDU
        Accel[IndexW+5] += Moment/I # DDR
  
    #} End of Get Force and Moment

    #{ VERBOSE - UPPER-MIDDLE LEFT
    if '/ForceH/' in VERBOSE:
      ContactConfig = 'Upper & Middle Left'
      
      if  'Right' in ContactConfig: SignK = 1
      elif 'Left' in ContactConfig: SignK = -1
      else: ERROR('Error: ContactConfig("%s") was invalid.'%ContactConfig)
      
      Verbose('-'*80)
      Verbose('Force_Block_V() - %s, Caculation Verification'%ContactConfig)
      if FlagFound == True:
        Verbose('Matching Block = (%d,%d)'%(K+SignK,m))
        Verbose('Block Type Name = %s'%BN,NewLine=False)
        Verbose(', Hm = %f, Bm = %f'%(Hm,Bm))
        Verbose('Um,DUm,Rm,DRm = %e, %e, %e, %e'%(Um,DUm,Rm,DRm))
        Verbose('Epsilon = %e'%Epsilon,NewLine=False)
        if Epsilon > 0:
          Verbose(', dEpsilon = %e'%DEpsilon)
          Verbose('Force, Moment = %e, %e'%(Force,Moment))
          Verbose('Accel = %e, %e'%(Force/M,Moment/I))
        else:
          Verbose(' <= 0   --->  No Contact')
      else:
        Verbose('  Mathing block was not found.')
    #}

    #{ POST ML
    if POST==True:
      
      # FILENAME
      Filename = FO_DIR+'/(%2d,%2d)_ML.csv'%(K,L)

      #{ HEADER
      try: fo = open(Filename,'r')
      except: # File doesn't exist -> This is the first time. Make header.
        temp  = 't,'
        temp += 'Epsil,F_MLk,dEpsil,F_MLc,F_ML,M_ML\n'
        fo = open(Filename,'w')
        fo.write(temp)
        fo.close()
      #}

      #{ F_MLk, F_MLc
      F_MLk = F_MLc = 0.
      if Epsilon > 0:
        F_MLk = +Kh *  Epsilon
        F_MLc = +Ch * DEpsilon
      #}

      #{ WRITE
      temp  = '%e, '%t
#              1   2   3   4   5   6    : 1        2      3         4      5      6
      temp += '%e, %e, %e, %e, %e, %e\n'%(Epsilon, F_MLk, DEpsilon, F_MLc, Force, Moment )
      fo = open(Filename,'a')
      fo.write(temp)
      fo.close()
      #}
    #}
    
    # --------------------------------------------------------------------------
    #} END OF (2) UPPER-MIDDLE LEFT
    # --------------------------------------------------------------------------
  
  
  
  
  
  
    
    # --------------------------------------------------------------------------
    #{ (3) LOWER-MIDDLE LEFT
    # --------------------------------------------------------------------------
  
    #{ Init
    Epsilon = 0.
    DEpsilon = 0.
    Force = 0.
    Moment = 0.
    #}
  
    #{ Start of Find (K-1,m) Block
    FlagFound = False
    len_m = len(Core['Array'][K-1])
    for m in xrange(1,len_m+1): # L starts from 1 to m
  
      # Get current state of (K-1,m)
      Um,DUm,Wm,DWm,Rm,DRm = GetSolFromVect(Core,w,K-1,m) # U,DU,W,DW,R,DR
      
      # Corner Coords of (K-1,m)  
      temp = GetBlockCoords (Core,K-1,m,Um,Wm,Rm) # [Center,LU,LD,RU,RD]
      Z_lt_m = temp[1][1]
      Z_lb_m = temp[2][1]
      Z_rt_m = temp[3][1]
      Z_rb_m = temp[4][1]
   
      # Lambda
      # For Left Case, Lambda needs only to change 
      # subscripts l(Left) and r(Right) with each other from the Right Case.
      # Subscript t,b,l(meand column index,L), and m should not be modified.
      Lambda = ( Z_lt_l - Z_rb_m ) / ( Z_lt_l - Z_lb_l )
      
      # Check 0<Lamda<1
      if 0 < Lambda < 1:
        FlagFound = True
        if m == len_m: FlagTopBlock = True # Top Block
        else: FlagTopBlock = False
        break
    #} End of Find (K-1,m) Block
        
    #{ Start of Get Force and Moment
    if FlagFound == True:
  
      # Get Block Prop of (K-1,m)
      BN = Core['Array'][K-1][m]
      BT = BlockTypes[BN]
      Bm = BT['b']
      Hm = BT['h']
    
      # Spring Contraction
      Epsilon  = -U  - ( (1-2*Lambda)*H *sin(R ) + B *(1-cos(R )) )
      Epsilon +=  Um + (            - Hm*sin(Rm) - Bm*(1-cos(Rm)) )
      if FlagTopBlock == False: Epsilon -= Delta
      else: Epsilon -= DeltaTop
  
      # Calculate Remains Only If Penetration Occurs
      if Epsilon > 0:
  
        # Spring Velocity    
        DEpsilon  = -DU  - ( (1-2*Lambda)*H *cos(R ) + B *sin(R ) ) * DR
        DEpsilon +=  DUm + (            - Hm*cos(Rm) - Bm*sin(Rm) ) * DRm
      
        # Force and Moment
        Force  =  Kh * Epsilon + Ch * DEpsilon
        Moment =  Force * ( (1-2*Lambda)*H*cos(R) + B*sin(R) )
      

        #-------------------------------------------
        # ELIMINATE STICKING FORCE
        # Left Side, Force Direction : --> (Right)
        # Eliminate Negative Force
        #-------------------------------------------
        if FLAG_NoStickForce == True:
          if Force < 0. : Force = Moment = 0.
        #-------------------------------------------


        # Acceleration
        Accel[IndexW+1] += Force/M  # DDU
        Accel[IndexW+5] += Moment/I # DDR
  
    #} End of Get Force and Moment

    #{ VERBOSE - LOWER-MIDDLE LEFT
    if '/ForceH/' in VERBOSE:
      ContactConfig = 'Lower & Middle Left'
      
      if  'Right' in ContactConfig: SignK = 1
      elif 'Left' in ContactConfig: SignK = -1
      else: ERROR('Error: ContactConfig("%s") was invalid.'%ContactConfig)
      
      Verbose('-'*80)
      Verbose('Force_Block_V() - %s, Caculation Verification'%ContactConfig)
      if FlagFound == True:
        Verbose('Matching Block = (%d,%d)'%(K+SignK,m))
        Verbose('Block Type Name = %s'%BN,NewLine=False)
        Verbose(', Hm = %f, Bm = %f'%(Hm,Bm))
        Verbose('Um,DUm,Rm,DRm = %e, %e, %e, %e'%(Um,DUm,Rm,DRm))
        Verbose('Epsilon = %e'%Epsilon,NewLine=False)
        if Epsilon > 0:
          Verbose(', dEpsilon = %e'%DEpsilon)
          Verbose('Force, Moment = %e, %e'%(Force,Moment))
          Verbose('Accel = %e, %e'%(Force/M,Moment/I))
        else:
          Verbose(' <= 0   --->  No Contact')
      else:
        Verbose('  Mathing block was not found.')
    #}

    #{ POST NL
    if POST==True:

      # FILENAME
      Filename = FO_DIR+'/(%2d,%2d)_NL.csv'%(K,L)

      #{ HEADER
      try: fo = open(Filename,'r')
      except: # File doesn't exist -> This is the first time. Make header.
        temp  = 't,'
        temp += 'Epsil,F_NLk,dEpsil,F_NLc,F_NL,M_NL\n'
        fo = open(Filename,'w')
        fo.write(temp)
        fo.close()
      #}

      #{ F_NLk, F_NLc
      F_NLk = F_NLc = 0.
      if Epsilon > 0:
        F_NLk = +Kh *  Epsilon
        F_NLc = +Ch * DEpsilon
      #}

      #{ WRITE
      temp  = '%e, '%t
#              1   2   3   4   5   6    : 1        2      3         4      5      6
      temp += '%e, %e, %e, %e, %e, %e\n'%(Epsilon, F_NLk, DEpsilon, F_NLc, Force, Moment )
      fo = open(Filename,'a')
      fo.write(temp)
      fo.close()
      #}
    #}
    
    # --------------------------------------------------------------------------
    #} END OF (3) LOWER-MIDDLE LEFT
    # --------------------------------------------------------------------------
  
  
  
  
  
  
    # --------------------------------------------------------------------------
    #{ (4) LOWER LEFT
    # --------------------------------------------------------------------------
  
    #{ Init
    Epsilon = 0.
    DEpsilon = 0.
    Force = 0.
    Moment = 0.
    #}
  
    #{ Start of Find (K-1,m) Block
    FlagFound = False
    len_m = len(Core['Array'][K-1])
    for m in xrange(1,len_m+1): # L starts from 1 to m
  
      # Get current state of (K-1,m)
      Um,DUm,Wm,DWm,Rm,DRm = GetSolFromVect(Core,w,K-1,m) # U,DU,W,DW,R,DR
      
      # Corner Coords of (K-1,m)  
      temp = GetBlockCoords (Core,K-1,m,Um,Wm,Rm) # [Center,LU,LD,RU,RD]
      Z_lt_m = temp[1][1]
      Z_lb_m = temp[2][1]
      Z_rt_m = temp[3][1]
      Z_rb_m = temp[4][1]
   
      # Lambda
      # For Left Case, Lambda needs only to change 
      # subscripts l(Left) and r(Right) with each other from the Right Case.
      # Subscript t,b,l(meand column index,L), and m should not be modified.
      Lambda = ( Z_rt_m - Z_lb_l ) / ( Z_rt_m - Z_rb_m )
      
      # Check 0<Lamda<=1
      if 0 < Lambda <= 1:
        FlagFound = True
        if m == len_m: FlagTopBlock = True # Top Block
        else: FlagTopBlock = False
        break
    #} End of Find (K-1,m) Block
        
    #{ Start of Get Force and Moment
    if FlagFound == True:
  
      # Get Block Prop of (K-1,m)
      BN = Core['Array'][K-1][m]
      BT = BlockTypes[BN]
      Bm = BT['b']
      Hm = BT['h']
    
      # Spring Contraction
      Epsilon  = -U  - (            - H *sin(R ) + B *(1-cos(R )) )
      Epsilon +=  Um + ( (1-2*Lambda)*Hm*sin(Rm) - Bm*(1-cos(Rm)) )
      if FlagTopBlock == False: Epsilon -= Delta
      else: Epsilon -= DeltaTop
  
      # Calculate Remains Only If Penetration Occurs
      if Epsilon > 0:
  
        # Spring Velocity    
        DEpsilon  = -DU  - (            - H *cos(R ) + B *sin(R ) ) * DR
        DEpsilon +=  DUm + ( (1-2*Lambda)*Hm*cos(Rm) - Bm*sin(Rm) ) * DRm
      
        # Force and Moment
        Force  =  Kh * Epsilon + Ch * DEpsilon
        Moment =  Force * (-H*cos(R) + B*sin(R) )
      

        #-------------------------------------------
        # ELIMINATE STICKING FORCE
        # Left Side, Force Direction : --> (Right)
        # Eliminate Negative Force
        #-------------------------------------------
        if FLAG_NoStickForce == True:
          if Force < 0. : Force = Moment = 0.
        #-------------------------------------------


        # Acceleration
        Accel[IndexW+1] += Force/M  # DDU
        Accel[IndexW+5] += Moment/I # DDR
  
    #} End of Get Force and Moment

    #{ VERBOSE - LOWER LEFT
    if '/ForceH/' in VERBOSE:
      ContactConfig = 'Lower & Left'
      
      if  'Right' in ContactConfig: SignK = 1
      elif 'Left' in ContactConfig: SignK = -1
      else: ERROR('Error: ContactConfig("%s") was invalid.'%ContactConfig)
      
      Verbose('-'*80)
      Verbose('Force_Block_V() - %s, Caculation Verification'%ContactConfig)
      if FlagFound == True:
        Verbose('Matching Block = (%d,%d)'%(K+SignK,m))
        Verbose('Block Type Name = %s'%BN,NewLine=False)
        Verbose(', Hm = %f, Bm = %f'%(Hm,Bm))
        Verbose('Um,DUm,Rm,DRm = %e, %e, %e, %e'%(Um,DUm,Rm,DRm))
        Verbose('Epsilon = %e'%Epsilon,NewLine=False)
        if Epsilon > 0:
          Verbose(', dEpsilon = %e'%DEpsilon)
          Verbose('Force, Moment = %e, %e'%(Force,Moment))
          Verbose('Accel = %e, %e'%(Force/M,Moment/I))
        else:
          Verbose(' <= 0   --->  No Contact')
      else:
        Verbose('  Mathing block was not found.')
    #}

    #{ POST BL
    if POST==True:

      # FILENAME
      Filename = FO_DIR+'/(%2d,%2d)_BL.csv'%(K,L)

      #{ HEADER
      try: fo = open(Filename,'r')
      except: # File doesn't exist -> This is the first time. Make header.
        temp  = 't,'
        temp += 'Epsil,F_BLk,dEpsil,F_BLc,F_BL,M_BL\n'
        fo = open(Filename,'w')
        fo.write(temp)
        fo.close()
      #}

      #{ F_BLk, F_BLc
      F_BLk = F_BLc = 0.
      if Epsilon > 0:
        F_BLk = +Kh *  Epsilon
        F_BLc = +Ch * DEpsilon
      #}

      #{ WRITE
      temp  = '%e, '%t
#              1   2   3   4   5   6    : 1        2      3         4      5      6
      temp += '%e, %e, %e, %e, %e, %e\n'%(Epsilon, F_BLk, DEpsilon, F_BLc, Force, Moment )
      fo = open(Filename,'a')
      fo.write(temp)
      fo.close()
      #}
    #}
    
    # --------------------------------------------------------------------------
    #} END OF (4) LOWER LEFT
    # --------------------------------------------------------------------------
  
  else:
    # VERBOSE - LOWER RIGHT
    if '/ForceH/' in VERBOSE:
      Verbose('-'*80)
      Verbose('K(%d)<=1  --->  No ForceH on Left Side'%K)



  #endregion                  END OF LEFT SIDE
  # ==========================================================================


  pass #  raw_input('pause FORCE_BLOCK_H')

  # VERBOSE
  if '/ForceH/' in VERBOSE:
    Verbose('-'*80)
    Verbose('Force_Block_V(), Caculation Verification')
    Verbose('DDU, DDR = %e, %e'%(Accel[IndexW+1],Accel[IndexW+5]))
    Verbose('Accel[%d:%d] = '%(IndexW,IndexW+6),NewLine=False)
    Verbose('%e, '%Accel[IndexW  ],NewLine=False)
    Verbose('%e, '%Accel[IndexW+1],NewLine=False)
    Verbose('%e, '%Accel[IndexW+2],NewLine=False)
    Verbose('%e, '%Accel[IndexW+3],NewLine=False)
    Verbose('%e, '%Accel[IndexW+4],NewLine=False)
    Verbose('%e'%Accel[IndexW+5])



  return Accel
####################### END OF FORCE_BLOCK_H() ################################


