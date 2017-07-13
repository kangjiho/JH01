VERSION_NAME = 'SAPCOR_v04j'

# MODULE NAME
MODULE_NAME_THIS__FORCE = VERSION_NAME+'_Force'


import glob

# -------------------------------------------------------------------
#{ MODULE NAMES TO BE IMPORTED
MODULE_NAMES_FORCE = []
MODULE_NAMES_FORCE.append('Misc')
MODULE_NAMES_FORCE.append('Input')
MODULE_NAMES_FORCE.append('Force_Block_D_F')
MODULE_NAMES_FORCE.append('Force_Block_H')
#MODULE_NAMES_FORCE.append('Force_Block_V')
MODULE_NAMES_FORCE.append('Force_Block_V_F')
#} END OF MODULE NAMES TO BE IMPORTED
# -------------------------------------------------------------------


# -------------------------------------------------------------------
#{ IMPORT MODULES
for Module_Name in MODULE_NAMES_FORCE:
  FileName = VERSION_NAME+'_'+Module_Name+'_Build*.py'
  FList = glob.glob(FileName)
  assert len(FList)!=0,'ERROR: Module not found : '+FileName
  FList.sort()
  FileName = FList[-1][0:-3]
  exec('from '+FileName+' import *')
#} END OF IMPORT MODULES
# -------------------------------------------------------------------



from numpy import loadtxt,pi,sin,cos
from numpy import array,zeros,hstack
import time



# PERFORMANCE MEASUREMENT
# VectorField()
NRun_VectorField = 0
RunTime_VectorField = 0.
NRun_VectorField_Prev = 0
RunTime_VectorField_Prev = 0.




0# === VECTOR FIELD ===========================================================
0

def VectorField(w, t, Core, VERBOSE='', POST=False, FO_DIR=''): #{
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

  1# --------------------------------------------------------------------------
  1# INIT VARS FOR PERFORMANCE MEASUREMENT
  while True:

  # Global Vars:
  #   NRun_VectorField : number of runs of VectorField()
  #   AvgRunTime_VectorField : Average Run Time of VectorField()


  global NRun_VectorField
  global RunTime_VectorField
  StartTime=time.time()


  1# --------------------------------------------------------------------------
  1# VERBOSE
  if '/Force/' in VERBOSE:
    Verbose('='*80)
    Verbose('<VectorField()>')

  1# ------------------------------------------------------------------------
  1# READ GLOBAL VAR SPACE
#  CoreIndex = Core['Index']
  LCoreReverseIndex = len(Core['ReverseIndex'])
  M2 = Core['M']+2 # 0: Left restraint, 1~(M-1): Blocks, M+1: Right restraint
  BodyForce=Core['BodyForce']

  1# ------------------------------------------------------------------------
  1# INIT LOCAL VARS
  Accel = [0]*(LCoreReverseIndex*6)   # InitTemp Memory
  IndexBlock = 0                      # Index for CoreIndex
  IndexW = 0                          # Index for Solution Vector
  # IndexW = IndexBlock*6

  # ------------------------------------------------------------------------
  #{ FIX "FIXED TO BASE" BLOCKS
  # COPY STATE OF BASE TO "FIXED TO BASE" BLOCKS
  #
  #(v0f) StateOfBase=GetSolFromVect(Core,w,0,0)
  StateOfBase = w[0:6]
  for (K,L) in Core['KLFixedToBase']:
    PutValToVect(Core,w,K,L,StateOfBase)
  #} END OF FIX "FIXED TO BASE" BLOCKS
  # ------------------------------------------------------------------------

  # ------------------------------------------------------------------------
  #{ MAIN LOOP FOR CORE
  #
  for IndexBlock in xrange(LCoreReverseIndex):

    # Index on StateVector
    IndexW = IndexBlock*6

    # ------------------------------------------------------------------------
    #{ Read Block Info
    # Get (K,L) index from sequential index list
    K,L = Core['ReverseIndex'][IndexBlock]
    BTN = Core['BTNs'][IndexBlock]

    # Read 'Fixed' Setting
    if (K,L)==(0,0): # Base
      Fixed = 'Base'
    else: # CSB or normal blocks
      Fixed = BlockTypes[BTN]['Fixed']
    #} END OF Read Block Info
    # ------------------------------------------------------------------------

    #{ VERBOSE
    if '/Force/' in VERBOSE:
      Verbose('IndexBlock,IndexW=%d,%d  K,L=%d,%d'%(IndexBlock,IndexW,K,L))
    #} END OF VERBOSE

    # ------------------------------------------------------------------------
    #{ FORCE ON BASE
    if [K,L]==[0,0]:

      #{ Loop for Directions
      for i in (1,3,5):

        # Read Type
        Type = Load[i]['Type'] # 'None', 'Accel', 'Disp'

        #{ Error Check in Type
        if Type not in ('None', 'Accel', 'Disp'):
          Verbose('='*79)
          Verbose('ERROR!!!')
          Verbose('Invalid Type in External Load Input: '+Type)
          Verbose('Skip it.')
          Verbose('='*79)
          continue
        #}

        # None -> continue
        if Type == 'None': continue

        # Read Params
        FnType = Load[i]['FnType'] # 'Sin', 'Cos', or 'Data'

        #{ Error Check in FnType
        if FnType not in ('Sin', 'Cos', 'Data'):
          Verbose('='*79)
          Verbose('ERROR!!!')
          Verbose('Invalid FnType in External Load Input: '+FnType)
          Verbose('Skip it.')
          Verbose('='*79)
          continue
        #}

        #{ External Load Data Reading -> Not Yet Supported -> Continue
        if FnType == 'Data':
          Verbose('='*79)
          Verbose('WARNING!!!')
          Verbose('External Load Data Reading is not supported yet.')
          Verbose('Data Reading option is ignored and no load is imposed.')
          Verbose('='*79)
          continue
        #}

        #{ Analytic Functions
        else:

          #{ Read Remains
          Amp   = Load[i]['Amp'] # L/T2 or L
          Freq  = Load[i]['Freq'] # Hz
          Omega = 2*pi*Freq
          Phase = Load[i]['Phase'] # rad
          #}

          #{ Cal Value
          if   FnType == 'Sin': Val = Amp*sin(Omega*t-Phase)
          elif FnType == 'Cos': Val = Amp*cos(Omega*t-Phase)
          else: Val = 0. # Logically cannot be reached here
          #}

          #{ Imposing Loads
          if Type == 'Accel':

            # Accel = [DU,DDU,DW,DDW,DR,DDR]
            if i==1: Accel[IndexW+1] = Val # DDU
            if i==3: Accel[IndexW+3] = Val # DDW
            if i==5: Accel[IndexW+5] = Val # DDR

            # Need to impose initial velocity if Sin Accel
            if t == 0 and FnType == 'Sin':
              # w = [U,DU,W,DW,R,DR]
              if i==1: w[IndexW+1] = -Amp/Omega # DU
              if i==3: w[IndexW+3] = -Amp/Omega # DW
              if i==5: w[IndexW+5] = -Amp/Omega # DR

          elif Type == 'Disp':
            # w = [U,DU,W,DW,R,DR]
            if i==1: w[IndexW+0] = Val # U
            if i==3: w[IndexW+2] = Val # W
            if i==5: w[IndexW+4] = Val # R
          else:
            continue # Logically cannot be reached here
          #}

        #} END OF Analytic Functions

      #} END OF Loop for Directions


    #} END OF FORCE ON BASE
    # ------------------------------------------------------------------------

    # ------------------------------------------------------------------------
    #{ FORCE ON CSB
    #
    elif [K,L]==[1,0]: Force_CSB(w,t,Core)

      # DO FOLLOWIN IN Force_CSB()
#      if Fixed=='None': Force_CSB(w,t,Core,f)

      # STATE AND FORCE UPDATES ARE DETERMINED BY 'Fix' SETTING
      # WHICH WAS SET IN BlockTypes IN INPUT FILE

      # If not fixed -> Calculate Accel
#      if Fixed=='None': Accel=Force_CSB(w,t,Core)

      # If 'Fixed' -> Make Accel={0} (means maintaining initial state)
#      elif Fixed=='Fixed': Accel=[0.,0.,0.]

      # If 'FixedToBase' -> Copy State of Base to current state
      #                     And make Accel={0}
#      elif Fixed=='FixedToBase':
#        StateOfBase=GetSolFromVect(Core,w,0,0)
#        PutValToVect(Core,w,K,L,StateOfBase)
#        Accel=[0.,0.,0.]

      # INPUT ERROR
#      else: ERROR('Invalid CSB Fix Setting')
    #
    #} END OF FORCE ON CSB
    # ------------------------------------------------------------------------

    # ------------------------------------------------------------------------
    #{ INVAILD K,L
    elif K>1 and L==0: ERROR('Invalid Core Index')
    #} END OF INVAILD K,L
    # ------------------------------------------------------------------------

    # ------------------------------------------------------------------------
    #{ FORCES ON NORMAL BLOCKS
    #
    else: # Blocks


        if '/Force/' in VERBOSE:
          Verbose('Call FORCE_BLOCK')

        # FORCE ON LEFT REFLECTORS
        if K==0: Force_ReflectorL(w,t,Core,K,L)

        # FORCE ON RIGHT REFLECTORS
        elif K==Core['M']+1: Force_ReflectorR(w,t,Core,K,L)

        # FORCE ON INNER BLOCKS
        else:
#          Accel=Force_Block(w,t,Core,K,L,Accel)

          # HORIZONTAL
          if Core['ApplyForces']['Block_H']==True:
            Accel = Force_Block_H(w,t,Core,K,L,Accel,POST,FO_DIR)


          # VERTICAL & FRICTION
          if Core['ApplyForces']['Block_V']==True:
            Accel = Force_Block_V_F(w,t,Core,K,L,Accel,POST,FO_DIR)
#          Accel = Force_Block_V_F(w,t,Core,K,L,Accel)

          # DOWEL & FRICTION
          if Core['ApplyForces']['Block_D']==True:
            Accel = Force_Block_D_F(w,t,Core,K,L,Accel,POST,FO_DIR)

      #{ Seems to be OLD
      # DO FOLLOWING IN SEPARATE FORCE FUNCTIONS
      # STATE AND FORCE UPDATES ARE DETERMINED BY 'Fix' SETTING
      # WHICH WAS SET IN BlockTypes IN INPUT FILE

      # If not fixed -> Calculate Accel
#      if Fixed=='None':


        # FORCE ON LEFT REFLECTORS
#        if K==0: # Left Restraint
#          Accel=Force_RestraintL(w,t,Core,K,L)

        # FORCE ON RIGHT REFLECTORS
#        elif K==Core['M']+1: # Right Restraint
#          Accel=Force_RestraintR(w,t,Core,K,L)

        # FORCE ON INNER BLOCKS
#        else:
#          Accel=Force_Block(w,t,Core,K,L)

      # If 'Fixed' -> Make Accel={0] (means maintaining initial state)
#      elif Fixed=='Fixed': Accel=[0.,0.,0.]

      # If 'FixedToBase' -> Copy State of Base to current state
      #                     And make Accel={0}
#      elif Fixed=='FixedToBase':
#        StateOfBase=GetSolFromVect(Core,w,0,0)
#        PutValToVect(Core,w,K,L,StateOfBase)
#        Accel=[0.,0.,0.]

      # INPUT ERROR
#      else: ERROR('Invalid CSB Fix Setting')
      #} END OF Seems to be OLD

    #
    #} END OF FORCES ON BLOCKS
    # ------------------------------------------------------------------------

    # ------------------------------------------------------------------------
    #{ BODY FORCES
    #
    # Input: Body Force = (0:DU 1:DDU, 2:DW, 3:DDW, 4:DR, 5:DDR)
    #
    # If BASE([0,0]) -> NO BODY FORCE
    # If CSB ([1,0]) -> NO BODY FORCE
    # If LEFT RESTRAINT ([0,:]) -> NO BODY FORCE
    # If RIGHT RESTRAINT([M+1,:]) -> NO BODY FORCE

#{ OLD VERSION
#(v04f)    if [K,L]!=[0,0] and [K,L]!=[1,0] and K!=0 and K!=L+1 and Fixed=='None':
#(v04f)      # DO SOMETHING HERE
#(v04f)      Accel[IndexW+1] += BodyForce[1] # DDU
#(v04f)      Accel[IndexW+3] += BodyForce[3] # DDW
#(v04f)      Accel[IndexW+5] += BodyForce[5] # DDR
#(v04f)      pass
#} END OF OLD VERSION

#(v04f BUG) if [K,L]!=[0,0] and [K,L]!=[1,0] and K!=0 and K!=L+1:
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
      #} EOIF if Fixed=='None':
    #} EOIF if [K,L]!=[0,0] and [K,L]!=[1,0] and K!=0 and K!=Core['M']+1:

    #} END OF BODY FORCES
    # ------------------------------------------------------------------------

    # ------------------------------------------------------------------------
    #{ DISP SENSOR DAMPING
    #
    # If BASE([0,0]) -> NO  DISP SENSOR DAMPING
    # If CSB ([1,0]) -> YES DISP SENSOR DAMPING
    # If LEFT RESTRAINT ([ 0 ,:]) -> NO DISP SENSOR DAMPING
    # If RIGHT RESTRAINT([M+1,:]) -> NO DISP SENSOR DAMPING

    if [K,L]!=[0,0] and K!=0 and K!=Core['M']+1: #{

      # 1th Order Disp Damping Coefficient
      # C = k1*ABS(U)+k2
      C = DispSensorDamping[0]*abs(w[IndexW]) + DispSensorDamping[1]

      if Fixed=='None': #{
        Accel[IndexW+1] -= C*w[IndexW+1] # DDU
      elif Fixed=='FixedU':
        pass
      elif Fixed=='FixedUR':
        pass
      pass #} EOIF

    #} EOIF if [K,L]!=[0,0] and K!=0 and K!=Core['M']+1:

    #     
    #} END OF SYSTEM DAMPING
    # ------------------------------------------------------------------------

    # ------------------------------------------------------------------------
    #{ ASSEMBLE FORCES
    #

    Accel[IndexW  ] = w[IndexW+1] # DU
    Accel[IndexW+2] = w[IndexW+3] # DW
    Accel[IndexW+4] = w[IndexW+5] # DR


    # Manage Fixed
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
      Accel[IndexW  ] = 0 #  DU=0 -> next step:  U=Contant
      Accel[IndexW+1] = 0 # DDU=0 -> next step: DU=Constant
#      Accel[IndexW+2] = 0 #  DW=0 -> next step:  W=Constant
#      Accel[IndexW+3] = 0 # DDW=0 -> next step: DW=Constant
#      Accel[IndexW+4] = 0 #  DR=0 -> next step:  R=Constant
#      Accel[IndexW+5] = 0 # DDR=0 -> next step: DR=Constant
    elif Fixed=='FixedUR':
      Accel[IndexW  ] = 0 #  DU=0 -> next step:  U=Contant
      Accel[IndexW+1] = 0 # DDU=0 -> next step: DU=Constant
#      Accel[IndexW+2] = 0 #  DW=0 -> next step:  W=Constant
#      Accel[IndexW+3] = 0 # DDW=0 -> next step: DW=Constant
      Accel[IndexW+4] = 0 #  DR=0 -> next step:  R=Constant
      Accel[IndexW+5] = 0 # DDR=0 -> next step: DR=Constant


    #
    #} END OF ASSEMBLE FORCES
    # ------------------------------------------------------------------------

#    IndexW += 6
  #
  #} END OF MAIN LOOP FOR CORE
  # ------------------------------------------------------------------------

  # ------------------------------------------------------------------------
  #{ TEST
#  if Core['TEST']!=Core['InitialBlockCenters']:
#    print "Force-Loop2 Core['InitialBlockCenters']"
#    for K in Core['InitialBlockCenters'].keys():
#      Col = Core['InitialBlockCenters'][K]
#      for L in Col.keys():
#        print '(%d,%d)'%(K,L),'=',Col[L]
#    raw_input('pause')
  #} END OF TEST
  # ------------------------------------------------------------------------

  # ------------------------------------------------------------------------
  #{ PERFORMANCE MEASUREMENT
  RunTime_VectorField += time.time()-StartTime
  NRun_VectorField += 1
  #} END OF PERFORMANCE MEASUREMENT
  # ------------------------------------------------------------------------

  return Accel





#==== PERFORMANCE RESULT OF VECTOR FIELD =================================

def PerfResult_VectorField(): #{
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
  return a,b,c,d,e,f
#}=========================================================================




#=========================================================================
#{ BLOCK FORCES




#-------------------------------------------------------------------------

def Force_CSB (w,t,Core):
  return

#-------------------------------------------------------------------------


def Force_ReflectorL (w,t,Core,K,L):
  return

#-------------------------------------------------------------------------

def Force_ReflectorR (w,t,Core,K,L):
  return

#-------------------------------------------------------------------------

def OLD_Force_Block (w,t,Core,K,L,Accel):

  # HORIZONTAL
  if Core['ApplyForces']['Block_H']==True:
    if '/Force/' in VERBOSE:
      Verbose('-'*80)
      Verbose('Call FORCE_BLOCK_H')
    Accel = Force_Block_H(w,t,Core,K,L)


  # VERTICAL & FRICTION
#  if Core['ApplyForces']['Block_V']==True:
  if Core['ApplyForces']['Block_V_F']==True:
    if '/Force/' in VERBOSE:
      Verbose('-'*80)
#      Verbose('Call FORCE_BLOCK_V')
      Verbose('Call FORCE_BLOCK_V_F')
#    Accel = Force_Block_V(w,t,Core,K,L,Accel)
    Accel = Force_Block_V_F(w,t,Core,K,L,Accel)

  # DOWEL
  if Core['ApplyForces']['Block_D']==True:
    if '/Force/' in VERBOSE:
      Verbose('-'*80)
      Verbose('Call FORCE_BLOCK_D')
    Accel = Force_Block_D(w,t,Core,K,L,Accel)

  return Accel


#=========================================================================
# END OF BLOCK FORCES
#}=========================================================================






#=========================================================================
# Extract Acceleration from State Vector
#=========================================================================

def OLD_ExtractAccelFromXV (Solution,t,Core,VERBOSE=''): #{
  Accel = []
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
    Accel.append(f[1::2])
#    Accel=numpy.array(Accel)
  return Accel
#}

def ExtractAccelFromXV (Solution,t,Core,VERBOSE=''): #{

  start = time.time()

  Accel = zeros([len(t),len(Core['ReverseIndex'])*3])

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

  return Accel
#}

#=========================================================================
