# coding=utf-8
VERSION_NAME = 'SAPCOR_v05a'

# MODULE NAME
MODULE_NAME_THIS__FORCE = VERSION_NAME+'_Force'

# -----------------------------------------------------------------------------
#region                       NORMAL IMPORTS


import glob
from numpy import loadtxt,pi,sin,cos
from numpy import array,zeros,hstack
import time


#endregion                 END OF NORMAL IMPORTS
# -----------------------------------------------------------------------------
#region                         USER IMPORTS


# region MODULE NAMES TO BE IMPORTED
MODULE_NAMES_FORCE = []
MODULE_NAMES_FORCE.append('Misc')
MODULE_NAMES_FORCE.append('Input')
MODULE_NAMES_FORCE.append('Force_Block_D_F')
MODULE_NAMES_FORCE.append('Force_Block_H')
#MODULE_NAMES_FORCE.append('Force_Block_V')
MODULE_NAMES_FORCE.append('Force_Block_V_F')
# endregion MODULE NAMES TO BE IMPORTED

# region IMPORT MODULES
for Module_Name in MODULE_NAMES_FORCE:
  FileName = VERSION_NAME+'_'+Module_Name+'_Build*.py'
  FList = glob.glob(FileName)
  assert len(FList)!=0,'ERROR: Module not found : '+FileName
  FList.sort()
  FileName = FList[-1][0:-3]
  exec('from '+FileName+' import *')
# endregion IMPORT MODULES


#endregion                   END OF USER IMPORTS
# -----------------------------------------------------------------------------
# region                       VARS - GLOBAL


Accel = []


# endregion                END OF VARS - GLOBAL
# -----------------------------------------------------------------------------
# region               VARS - PERFORMANCE MEASUREMENT


# For VectorField()
NRun_VectorField = 0
RunTime_VectorField = 0.
NRun_VectorField_Prev = 0
RunTime_VectorField_Prev = 0.


# endregion        END OF VARS - PERFORMANCE MEASUREMENT
# =============================================================================
def MallocAccel( Core ): # MEMORY ALLOCATION OF ACCEL
  global Accel
  LCoreReverseIndex = len(Core['ReverseIndex'])
  for i in range(LCoreReverseIndex * 6):
    Accel.append(0.)
  return
# =============================================================================
def VectorField(w, t, Core, VERBOSE='', POST=False, FO_DIR=''):
  """
    Defines the differential equations for the coupled spring-mass system.

    Arguments:
        w :  vector of the state variables at time, t:
                  w = [U,DU,W,DW,R,DR,...]
        t :  time
        Core : Var space
    
    Return:
        [DU,DDU,DW,DDW,DR,DDR,...]
  """




  # ---------------------------------------------------------------------------
  # region INIT VARS - PERFORMANCE MEASUREMENT
  # Global Vars:
  #   NRun_VectorField : number of runs of VectorField()
  #   AvgRunTime_VectorField : Average Run Time of VectorField()
  global NRun_VectorField
  global RunTime_VectorField
  StartTime=time.time()
  # endregion INIT VARS - PERFORMANCE MEASUREMENT
  # ---------------------------------------------------------------------------
  # region VERBOSE
  if '/Force/' in VERBOSE:
    Verbose('='*80)
    Verbose('<VectorField()>')
  # endregion VERBOSE
  # ---------------------------------------------------------------------------
  # region READ GLOBAL VAR SPACE
  pass #[Deleted] CoreIndex = Core['Index']
  LCoreReverseIndex = len(Core['ReverseIndex'])
  M2 = Core['M'] + 2  # 0: Left restraint, 1~(M-1): Blocks, M+1: Right restraint
  BodyForce = Core['BodyForce']
  # endregion READ GLOBAL VAR SPACE
  # ---------------------------------------------------------------------------
  # region INIT VARS
  global Accel
  for i in range(LCoreReverseIndex*6): Accel[i]=1.
  IndexBlock = 0                      # Index for CoreIndex
  IndexW = 0                          # Index for Solution Vector
  #[Deleted] IndexW = IndexBlock*6
  # endregion INIT LOCAL VARS
  # ---------------------------------------------------------------------------
  # region              FIX "FIXED TO BASE" BLOCKS


  # COPY STATE OF BASE TO "FIXED TO BASE" BLOCKS
  #
  #[Deleted v0f] StateOfBase=GetSolFromVect(Core,w,0,0)
  StateOfBase = w[0:6]
  for (K,L) in Core['KLFixedToBase']:
    PutValToVect(Core,w,K,L,StateOfBase)


  # endregion       END OF FIX "FIXED TO BASE" BLOCKS
  # ---------------------------------------------------------------------------
  #region                   MAIN LOOP FOR CORE


  for IndexBlock in xrange(LCoreReverseIndex):


    # -------------------------------------------------------------------------
    # region                  Index on StateVector


    IndexW = IndexBlock*6


    # endregion               Index on StateVector
    # -------------------------------------------------------------------------
    # region                       Read Block Info


    # Get (K,L) index from sequential index list
    K,L = Core['ReverseIndex'][IndexBlock]
    BTN = Core['BTNs'][IndexBlock]

    # Read 'Fixed' Setting
    if (K,L)==(0,0): # Base
      Fixed = 'Base'
    else: # CSB or normal blocks
      Fixed = BlockTypes[BTN]['Fixed']


    # endregion                   Read Block Info
    # -------------------------------------------------------------------------
    # region                            VERBOSE


    if '/Force/' in VERBOSE:
      Verbose('IndexBlock,IndexW=%d,%d  K,L=%d,%d'%(IndexBlock,IndexW,K,L))


    # endregion                         VERBOSE
    # -------------------------------------------------------------------------
    # region               FORCES ON DIFFERENT BLOCK TYPES



    # -------------------------------------------------------------------------
    #region                           BASE


    if [K,L]==[0,0]:
      # -------------------------------------------------------------------------
      # region Loop for Directions


      for i in (1,3,5):

        # Read Type
        Type = Load[i]['Type'] # 'None', 'Accel', 'Disp'

        # Error Check in Type
        if Type not in ('None', 'Accel', 'Disp'):
          Verbose('='*79)
          Verbose('ERROR!!!')
          Verbose('Invalid Type in External Load Input: '+Type)
          Verbose('Skip it.')
          Verbose('='*79)
          continue

        # None -> continue
        if Type == 'None': continue

        # Read Params
        FnType = Load[i]['FnType'] # 'Sin', 'Cos', or 'Data'

        # Error Check in FnType
        if FnType not in ('Sin', 'Cos', 'Data'):
          Verbose('='*79)
          Verbose('ERROR!!!')
          Verbose('Invalid FnType in External Load Input: '+FnType)
          Verbose('Skip it.')
          Verbose('='*79)
          continue

        # -------------------------------------------------------------------------
        # region          Calculate External Forces on Base


        # ---------------------------------------------------------------------
        #  region IF External Load Data Reading -> Not Yet Supported: Continue


        if FnType == 'Data':
          Verbose('='*79)
          Verbose('WARNING!!!')
          Verbose('External Load Data Reading is not supported yet.')
          Verbose('Data Reading option is ignored and no load is imposed.')
          Verbose('='*79)
          continue


        # endregion IF External Load Data Reading -> Not Yet Supported: Continue
        # ---------------------------------------------------------------------
        # region                ELSE: Analytic Functions


        else:
          # ---------------------------------------------------------------------
          # region                        Read Remains


          Amp   = Load[i]['Amp'] # L/T2 or L
          Freq  = Load[i]['Freq'] # Hz
          Omega = 2*pi*Freq
          Phase = Load[i]['Phase'] # rad


          # endregion                      Read Remains
          # -------------------------------------------------------------------------
          # region                        Cal Value


          if   FnType == 'Sin': Val = Amp*sin(Omega*t-Phase)
          elif FnType == 'Cos': Val = Amp*cos(Omega*t-Phase)
          else: Val = 0. # Logically cannot be reached here


          # endregion                     Cal Value
          # -------------------------------------------------------------------------
          # region Imposing Loads


          # -------------------------------------------------------------------------
          # region                    IF Type == 'Accel':


          if Type == 'Accel':
            # -------------------------------------------------------------------------
            # region         Accel = [DU,DDU,DW,DDW,DR,DDR]


            if i==1: Accel[IndexW+1] = Val # DDU
            if i==3: Accel[IndexW+3] = Val # DDW
            if i==5: Accel[IndexW+5] = Val # DDR


            # endregion                    Accel
            # -------------------------------------------------------------------------
            # region      Need to impose initial velocity if Sin Accel


            if t == 0 and FnType == 'Sin':
              # w = [U,DU,W,DW,R,DR]
              if i==1: w[IndexW+1] = -Amp/Omega # DU
              if i==3: w[IndexW+3] = -Amp/Omega # DW
              if i==5: w[IndexW+5] = -Amp/Omega # DR
            else: pass


            # endregion    Need to impose initial velocity if Sin Accel
            # -------------------------------------------------------------------------


          # endregion                IF Type == 'Accel':
          # -------------------------------------------------------------------------
          # region                  ELIF Type == 'Disp':


          elif Type == 'Disp':
            pass # w = [U,DU,W,DW,R,DR]
            if i==1: w[IndexW+0] = Val # U
            if i==3: w[IndexW+2] = Val # W
            if i==5: w[IndexW+4] = Val # R


          # endregion               ELIF Type == 'Disp':
          # -------------------------------------------------------------------------
          # region         ELSE: Logically cannot be reached here


          else:
            continue # Logically cannot be reached here
          pass


          # endregion     ELSE: Logically cannot be reached here
          # -------------------------------------------------------------------------


          # endregion Imposing Loads
          # -------------------------------------------------------------------------
        pass #EOIF


        # endregion                Analytic Functions
        # -------------------------------------------------------------------------


        #endregion      END OF Calculate External Forces on Base
        # -------------------------------------------------------------------------

      pass #EOFOR

      #endregion           END OF Loop for Directions
      # -------------------------------------------------------------------------


    #endregion                     END OF BASE
    # -------------------------------------------------------------------------
    # region                             CSB


    elif [K,L]==[1,0]:
      Accel = Force_CSB(w, t, Core, K, L, Accel, POST, FO_DIR)


    # endregion                       END OF CSB
    # -------------------------------------------------------------------------
    # region                     CHECK INVALID K,L


    elif K>1 and L==0: ERROR('Invalid Core Index')


    # endregion               END OF CHECK INVALID K,L
    # -------------------------------------------------------------------------
    # region                      NORMAL BLOCKS


    else: # Blocks

      # VERBOSE
      if '/Force/' in VERBOSE: Verbose('Call FORCE_BLOCK')

      # FORCE ON LEFT REFLECTORS
      if K==0: Force_ReflectorL(w,t,Core,K,L)

      # FORCE ON RIGHT REFLECTORS
      elif K==Core['M']+1: Force_ReflectorR(w,t,Core,K,L)

      # FORCE ON INNER BLOCKS
      else:
        pass #[Deleted] Accel=Force_Block(w,t,Core,K,L,Accel)

        # HORIZONTAL
        if Core['ApplyForces']['Block_H']:
          Accel = Force_Block_H(w,t,Core,K,L,Accel,POST,FO_DIR)
          pass

        # VERTICAL & FRICTION
        if Core['ApplyForces']['Block_V']:
          Accel = Force_Block_V_F(w,t,Core,K,L,Accel,POST,FO_DIR)
          pass #[Deleted] Accel = Force_Block_V_F(w,t,Core,K,L,Accel)

        # DOWEL & FRICTION
        if Core['ApplyForces']['Block_D']:
          Accel = Force_Block_D_F(w,t,Core,K,L,Accel,POST,FO_DIR)
          pass

      pass #EOIF

    pass #EOIF


    # endregion                END OF NORMAL BLOCKS
    # -------------------------------------------------------------------------



    #endregion      END OF FORCES ON DIFFERENT BLOCK TYPES
    # -------------------------------------------------------------------------
    # region                     BODY FORCES


    # -------------------------------------------------------------------------


    # -------------------------------------------------------------------------
    #{ BODY FORCES
    #
    # Input: Body Force = (0:DU 1:DDU, 2:DW, 3:DDW, 4:DR, 5:DDR)
    #
    # If BASE           ([K,L]==[0  ,0]) -> NO BODY FORCE
    # If CSB            ([K,L]==[1  ,0]) -> NO BODY FORCE
    # If LEFT RESTRAINT ([K,L]==[0  ,:]) -> NO BODY FORCE
    # If RIGHT RESTRAINT([K,L]==[M+1,:]) -> NO BODY FORCE

    #{ OLD VERSION
    #(v04f)    if [K,L]!=[0,0] and [K,L]!=[1,0] and K!=0 and K!=L+1 and Fixed=='None':
    #(v04f)      # DO SOMETHING HERE
    #(v04f)      Accel[IndexW+1] += BodyForce[1] # DDU
    #(v04f)      Accel[IndexW+3] += BodyForce[3] # DDW
    #(v04f)      Accel[IndexW+5] += BodyForce[5] # DDR
    #(v04f)      pass
    #} END OF OLD VERSION
    # ------------------------------------------------------------------------

    pass #[Deleted v04f BUG] if [K,L]!=[0,0] and [K,L]!=[1,0] and K!=0 and K!=L+1:
    if [K,L]!=[0,0] and [K,L]!=[1,0] and K!=0 and K!=Core['M']+1: #{
      if Fixed=='None': #{
        Accel[IndexW+1] += BodyForce[1] # DDU
        Accel[IndexW+3] += BodyForce[3] # DDW
        Accel[IndexW+5] += BodyForce[5] # DDR
      elif Fixed=='FixedU':
        Accel[IndexW+3] += BodyForce[3] # DDW
        Accel[IndexW+5] += BodyForce[5] # DDR
      elif Fixed=='FixedUR':
        Accel[IndexW+3] += BodyForce[3] # DDW
      pass #} EOIF if Fixed=='None':
    pass #} EOIF if [K,L]!=[0,0] and [K,L]!=[1,0] and K!=0 and K!=Core['M']+1:


    # endregion               END OF BODY FORCES
    # -------------------------------------------------------------------------
    # region                     SENSOR DAMPING


   # ------------------------------------------------------------------------
    # DISP SENSOR DAMPING
    #
    # If BASE           ([0  ,0]) -> NO DISP SENSOR DAMPING
    # If CSB            ([1  ,0]) -> YES DISP SENSOR DAMPING
    # If LEFT RESTRAINT ([0  ,:]) -> NO DISP SENSOR DAMPING
    # If RIGHT RESTRAINT([M+1,:]) -> NO DISP SENSOR DAMPING
    # ------------------------------------------------------------------------

    if [K,L]!=[0,0] and K!=0 and K!=Core['M']+1: #{

      # 1st Order Disp Damping Coefficient
      # C = k1*ABS(U)+k2
      C = DispSensorDamping[0]*abs(w[IndexW]) + DispSensorDamping[1]

      if Fixed=='None':
        Accel[IndexW+1] -= C*w[IndexW+1] # DDU
      elif Fixed=='FixedU':
        pass
      elif Fixed=='FixedUR':
        pass
      pass #EOIF

    pass #EOIF

    # endregion             END OF SENSOR DAMPING
    # -------------------------------------------------------------------------
    # region                 ASSEMBLE FORCES



    # -------------------------------------------------------------------------
    #region                          ACCEL


    Accel[IndexW  ] = w[IndexW+1] # DU
    Accel[IndexW+2] = w[IndexW+3] # DW
    Accel[IndexW+4] = w[IndexW+5] # DR


    #endregion                    END OF ACCEL
    # -------------------------------------------------------------------------
    #region                       Manage Fixed


    if Fixed=='Fixed':
      Accel[IndexW  ] = w[IndexW+1] # DU
      Accel[IndexW+1] = 0 # DDU
      Accel[IndexW+2] = w[IndexW+3] # DW
      Accel[IndexW+3] = 0 # DDW
      Accel[IndexW+4] = w[IndexW+5] # DR
      Accel[IndexW+5] = 0 # DDR
    elif Fixed=='FixedToBase':
      Accel[IndexW  ] = 0 # DU
      Accel[IndexW+1] = 0 # DDU
      Accel[IndexW+2] = 0 # DW
      Accel[IndexW+3] = 0 # DDW
      Accel[IndexW+4] = 0 # DR
      Accel[IndexW+5] = 0 # DDR
    elif Fixed=='FixedU':
      Accel[IndexW  ] = 0 #  DU=0 -> next step:  U=Constant
      Accel[IndexW+1] = 0 # DDU=0 -> next step: DU=Constant
      pass #[deleted]      Accel[IndexW+2] = 0 #  DW=0 -> next step:  W=Constant
      pass #[deleted]      Accel[IndexW+3] = 0 # DDW=0 -> next step: DW=Constant
      pass #[deleted]      Accel[IndexW+4] = 0 #  DR=0 -> next step:  R=Constant
      pass #[deleted]      Accel[IndexW+5] = 0 # DDR=0 -> next step: DR=Constant
    elif Fixed=='FixedUR':
      Accel[IndexW  ] = 0 #  DU=0 -> next step:  U=Contant
      Accel[IndexW+1] = 0 # DDU=0 -> next step: DU=Constant
      pass #[deleted]      Accel[IndexW+2] = 0 #  DW=0 -> next step:  W=Constant
      pass #[deleted]      Accel[IndexW+3] = 0 # DDW=0 -> next step: DW=Constant
      Accel[IndexW+4] = 0 #  DR=0 -> next step:  R=Constant
      Accel[IndexW+5] = 0 # DDR=0 -> next step: DR=Constant
    pass


      #endregion                END OF Manage Fixed
    # -------------------------------------------------------------------------


    # endregion            END OF ASSEMBLE FORCES
    # ------------------------------------------------------------------------
    pass #[Deleted] IndexW += 6
    # -------------------------------------------------------------------------


  #endregion             END OF MAIN LOOP FOR CORE
  # ---------------------------------------------------------------------------
  # region PERFORMANCE MEASUREMENT
  RunTime_VectorField += time.time()-StartTime
  NRun_VectorField += 1
  # endregion PERFORMANCE MEASUREMENT
  # ---------------------------------------------------------------------------




  return Accel         # END OF def VectorField
# =============================================================================
def PerfResult_VectorField():



  """
Return: (a,b,c,d,e)
    a: NRun_VectorField    of this section only
    b: RunTime_VectorField of this section only
    c: Average RunTime     of this section only
    d: NRun_VectorField    of all calculation from beginning
    e: RunTime_VectorField of all calculation from beginning
    f: Average RunTime     of all calculation from beginning
  """
  global NRun_VectorField
  global RunTime_VectorField
  global NRun_VectorField_Prev
  global RunTime_VectorField_Prev
  a = NRun_VectorField    - NRun_VectorField_Prev
  b = RunTime_VectorField - RunTime_VectorField_Prev
  c = b / a
  d = NRun_VectorField
  e = RunTime_VectorField
  f = RunTime_VectorField / NRun_VectorField
  NRun_VectorField_Prev    = NRun_VectorField
  RunTime_VectorField_Prev = RunTime_VectorField



  return a,b,c,d,e,f # END OF DEF - PERFORMANCE RESULT OF VECTOR FIELD
#==============================================================================
def Force_CSB( w, t, Core, K, L, Accel, POST, FO_DIR ):
  return Accel
def Force_ReflectorL (w,t,Core,K,L):
  return
def Force_ReflectorR (w,t,Core,K,L):
  return
#==============================================================================
def ExtractAccelFromXV (Solution,t,Core,VERBOSE=''): #{
# Extract Acceleration from State Vector



  start = time.time()

  AccelTemp = zeros([len(t),len(Core['ReverseIndex'])*3])

  end = time.time()
  print '1:',end-start
  start = end

  print 'len(t)=',len(t)

#  for sol,t_i in Solution,t:
  for i in range(len(t)):
    sol = Solution[i]
    t_i = t[i]

    # f = [ DX,DDX,DW,DDW,DR,DDR of Base,
    #       DX,DDX,DW,DDW,DR,DDR of CSB,
    #       DX,DDX,DW,DDW,DR,DDR of Left Restraint,
    #       DX,DDX,DW,DDW,DR,DDR of First Block,
    #       ... ] at t_i
    f=VectorField(sol,t_i,Core,VERBOSE)

    # acceleration = [ DDX, DDW, DDR,  ... ]
    #              = [ f[1],f[3],f[5], ... ] at t_i
    # Accel.append (acceleration)
#    Accel.append(f[1::2])

    Accel[i]=f[1::2]

  end = time.time()
  print '4:',end-start
  start = end



  return AccelTemp # END OF DEF - Extract Acceleration from State Vector
#==============================================================================
def testfunc():
  print 'testfunc'
  return
#==============================================================================
pass #============================== END OF FILE ==============================