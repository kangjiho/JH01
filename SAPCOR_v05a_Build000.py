#
# SAPCOR odeint Version.py
#
"""Seismic Analysis for a Prismatic CORe of a HTGR
"""

VERSION_NAME = 'SAPCOR_v05a'


#region #################### SYSTEM INITIALIZATION ####################

#{ IMPORT =====================================================================
from scipy.integrate import odeint

from numpy import loadtxt,savetxt
from numpy import array,zeros
from pylab import figure, plot, xlabel, grid, hold, legend, title, savefig,ion,ioff,draw
import sys,numpy,time,shutil,glob,os


# -------------------------------------------------------------------
# Module Names
MODULE_NAMES = []
MODULE_NAMES.append('Misc')
MODULE_NAMES.append('Force')
MODULE_NAMES.append('Force_Block_D_F')
MODULE_NAMES.append('Force_Block_H')
MODULE_NAMES.append('Force_Block_V_F')
MODULE_NAMES.append('Input')
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Import Modules
for Module_Name in MODULE_NAMES:
  FileName = VERSION_NAME+'_'+Module_Name+'_Build*.py'
  FList = glob.glob(FileName)
  assert len(FList)!=0,'ERROR: .py not found : '+FileName
  FList.sort()
  FileName = FList[-1][0:-3]
  exec('from '+FileName+' import *')
# -------------------------------------------------------------------

#} IMPORT =====================================================================

#{ Make Output Folder =========================================================
start_time = time.time()
temp = time.localtime(start_time)
year = temp.tm_year
month = temp.tm_mon
day = temp.tm_mday
hour = temp.tm_hour
minute = temp.tm_min
second = temp.tm_sec
FO_DIR = 'Output_%4d-%02d-%02d_%02dh%02dm%02ds'%(year,month,day,hour,minute,second)
os.mkdir(FO_DIR)
VerboseInit(FO_DIR)
#} Make Output Folder =========================================================

Verbose('//INIT//')

# endregion #################### SYSTEM INITIALIZATION ####################


# region ##################### LOCAL FUNCTIONS ########################
def CheckIndex (Core):
  Verbose(' Index Check')
  Verbose('-------------')
  RIndex = Core['ReverseIndex']
  Index = Core['Index']
  NError = 0
  for KL in RIndex:
    i = RIndex.index(KL)
    j = Index[KL[0]][KL[1]]
#    print '(K,L)=',KL,'ReverseIndex=',i,'Index=',j,
    temp='  (K,L)=%s'%KL+' ReverseIndex=%d'%i+' Index=%d'%j
    Verbose(temp,False)
    if i!=j:
#      print '[ERROR] Index Mistach !'
      Verbose(' --> [ERROR] Index Mismatach !')
      NError+=1
    else:
#      print
      Verbose()
  if NError>0:
    ERROR('ERROR: Index Mismatch')
  return
# endregion ##################### LOCAL FUNCTIONS ########################


# region ######################## INITIALIZATION ##############################

Verbose('Started at '+time.ctime(start_time))
Verbose('='*40)

#{ Import Modules - LOG ---------------------------------------------
Verbose(' Imported Modules')
Verbose('------------------')
for Module_Name in MODULE_NAMES:
  FileName = VERSION_NAME+'_'+Module_Name+'_Build*.py'
  FList = glob.glob(FileName)
  FList.sort()
  FileName = FList[-1][0:-3]
  Verbose('  '+FileName)
Verbose()
#} Import Modules - LOG ----------------------------------------------


#----------------------------------------------------------------------------
#{ MAKE BACKUP OF INPUT FILE
FList = glob.glob(VERSION_NAME+'_Input_Build*.py')
FList.sort()
FileName = FList[-1]
shutil.copyfile(FileName,FO_DIR+'/'+FN_OUT+'.in')
#}----------------------------------------------------------------------------


Core={}


#-------------------------------------------------------------------------
#{ Index ('ReverseIndex')
CoreArray2={}
i_Index=1
Index2={}
Index=[]
Index.append([0,0]) # Base
K=0
Core['Flag_CSB']=False
Core['N']=[]


for Column in CoreArray: # Read a column (could be 'CSB')

  # If it is CSB
  if Column=='CSB':
    i_Index+=1 # If CSB exists, reserve Index[1][0]
    Index.append([1,0]) # CSB
    Core['Flag_CSB']=True

  # If it is Reflectors or Blocks
  else:
    CoreArray2[K]={}
    Index2[K]=[-1]
    Core['N'].append(len(Column))
    for i in range(len(Column)):
      L=i+1
      Index2[K].append(i_Index)
      i_Index+=1
      Index.append([K,L])
      CoreArray2[K][L]=Column[i]
    K+=1



# K : 0, 1, 2, ..., M, M+1, M+2
Core['M']=K-2 # It is not K-1 beacause last K+=1 in above loop

Index2[0][0]=0 # (0,0) ; Base ; Index=0

if Core['Flag_CSB']==True:
#{
  Index2[1][0]=1 # (1,0) ; CSB ; Index=1 if exists
#}

#CoreArray=CoreArray2

#} -------------------------------------------------------------------------


#-------------------------------------------------------------------------
#{ BTNsKL
BTNsKL={}
K=0

for Column in CoreArray: # Read a column (could be 'CSB')

  # If it is CSB: Process it later
  if Column=='CSB': continue

  BTNsKL[K]={}

  for i in xrange(len(Column)):
    L=i+1
    BTNsKL[K][L]=Column[i]
    
  K+=1

BTNsKL[0][0] = 'Base'

if CoreArray[0]=='CSB': BTNsKL[1][0] = 'CSB'
#} -------------------------------------------------------------------------



#-------------------------------------------------------------------------
#{ BTNs
BTNs=[]
for index in Index:
  if index==[0,0]: # BTN='Base' which is not in CoreArray
    BTN = 'Base'
  elif index==[1,0]: # BTN='CSB' which is not in CoreArray
    BTN = 'CSB'
  else:
    BTN = CoreArray2[index[0]][index[1]]
  BTNs.append(BTN)
#} -------------------------------------------------------------------------


#-------------------------------------------------------------------------
#{ Index for FixedToBase
IndexFixedToBase=[]
for i in xrange(len(Index)):
  if i==0: # Base
    continue
  BTN = BTNs[i]
  Fixed = BlockTypes[BTN]['Fixed']
  if Fixed=='FixedToBase':
    IndexFixedToBase.append(i)
  
#}-------------------------------------------------------------------------


#-------------------------------------------------------------------------
#{ (K,L) for FixedToBase
KLFixedToBase=[]
for i in IndexFixedToBase:
  (K,L) = Index[i]
  KLFixedToBase.append([K,L])
  
#}-------------------------------------------------------------------------


#-------------------------------------------------------------------------
#{ Index for Fixed
IndexFixed=[]
for i in xrange(len(Index)):
  if i==0: # Base
    continue
  BTN = BTNs[i]
  Fixed = BlockTypes[BTN]['Fixed']
  if Fixed=='Fixed':
    IndexFixed.append(i)
  
#}-------------------------------------------------------------------------


#-------------------------------------------------------------------------
#{ (K,L) for Fixed
KLFixed=[]
for i in IndexFixed:
  (K,L) = Index[i]
  KLFixed.append([K,L])
  
#}-------------------------------------------------------------------------


#-------------------------------------------------------------------------
#{ Core
Core['BlockTypes']=BlockTypes
Core['MatProp']=MatProp
Core['Array']=CoreArray2
Core['ReverseIndex']=Index
Core['Index']=Index2
Core['BTNs']=BTNs
Core['BTNsKL']=BTNsKL
if len(Core['ReverseIndex'])!=len(Core['BTNs']):
  ERROR("ERROR: len(Core['ReverseIndex']) != len(Core['BTNs'])")
Core['IndexFixedToBase']=IndexFixedToBase
Core['KLFixedToBase']=KLFixedToBase
Core['IndexFixed']=IndexFixed
Core['KLFixed']=KLFixed
# Core['Flag_CSB'] is assigned in Index
Core['BodyForce']=BodyForce
Core['Load']=Load
Core['ApplyForces']=ApplyForces
#MakeConnectivity(Core,VERBOSE=VERBOSE)
GetInitBlockCenters(Core,VERBOSE=VERBOSE)

MakeTxtCoreArray (Core) # Make text form of core array -> Core['TxtArray']

CheckIndex(Core)
#}-------------------------------------------------------------------------



#-------------------------------------------------------------------------
#{ LOG - INITIALIZATION
if '/Init/' in VERBOSE:
  Verbose('='*80)
  Verbose('<Initialization>')

  Verbose('-'*80)
  Verbose('Core Size')
  Verbose('  M=%2d'%Core['M'])
  for i in range(Core['M']+2):
    Verbose('  Column %2d, N=%2d'%(i,Core['N'][i]))

  Verbose('-'*80)
  Verbose('Rearranged CoreArray')
  for i in xrange(len(Core['Array'])):
    Verbose('%s'%CoreArray2[i])

  Verbose('-'*80)
  Verbose('Index, [K,L], BlockTypeName')
  tempi = len(Core['ReverseIndex'])
  for i in xrange(tempi):
    temps = '%d, %s, %s'%(i, Core['ReverseIndex'][i], Core['BTNs'][i])
    Verbose(temps)

  Verbose('-'*80)
  Verbose('BTNsKL')
  for K in Core['BTNsKL'].keys():
    for L in Core['BTNsKL'][K].keys():
      Verbose('(%d,%d)=%s'%(K,L,Core['BTNsKL'][K][L]))

  Verbose('-'*80)
  Verbose('FixedToBase')
  if len(Core['IndexFixedToBase'])==0:
    Verbose('  None')
  for i in xrange(len(Core['IndexFixedToBase'])):
    j = Core['IndexFixedToBase'][i]
    temps = '%d, %s, %s'%(i, Core['KLFixedToBase'][i], Core['BTNs'][j])
    Verbose(temps)
    
  Verbose('-'*80)
  Verbose('Fixed')
  if len(Core['IndexFixed'])==0:
    Verbose('  None')
  for i in xrange(len(Core['IndexFixed'])):
    j = Core['IndexFixed'][i]
    temps = '%d, %s, %s'%(i, Core['KLFixed'][i], Core['BTNs'][j])
    Verbose(temps)

  Verbose('-'*80)
  Verbose("Core['InitialBlockCenters']")
  for K in Core['InitialBlockCenters'].keys():
    Col = Core['InitialBlockCenters'][K]
    for L in Col.keys():
      Verbose('(%d,%d)=%s'%(K,L,Col[L]))
    
  Verbose('\n%s'%Core['TxtArray'])
  PrintIndex(Core)

  Verbose('='*80)
#}-------------------------------------------------------------------------


CurrentTime=0.



#----------------------------------------------------
#{ Initial Conditions
# InitialConditions.append([K,L, U,DU, W,DW, R,DR])

StartValues = [0. for i in xrange(6*len(Core['ReverseIndex']))]
for i in xrange(len(InitialConditions)):
  # InitialConditions[i]==[K,L, U,DU, W,DW, R,DR]
#  KL = InitialConditions[i][:2]
#  index = Core['ReverseIndex'].index(KL)
  K,L = InitialConditions[i][:2]
  index = Core['Index'][K][L]
  for j in xrange(6):
    StartValues[6*index+j] = InitialConditions[i][j+2]

if '/Init/' in VERBOSE:
  Verbose('\n<Initial Conditions>')
  for i in xrange(len(InitialConditions)):
    Verbose('(K,L)=%s'%InitialConditions[i][:2],False)
#    print ', Index=',6*Core['Index'].index(InitialConditions[i][:2]),
    K,L = InitialConditions[i][:2]
    Verbose(', Index=%d'%(6*Core['Index'][K][L]),False)
    Verbose(', (U,DU,W,DW,R,DR)=%s'%InitialConditions[i][2:])
  for i in xrange(len(StartValues)):
    if StartValues[i]!=0.:
      Verbose('StartValues[%d]=%f'%(i,StartValues[i]))
#}----------------------------------------------------
      
#{ SAVE INITIAL T,STATEVECTOR      
fo = open(FO_DIR+'/Time.csv','w')
savetxt(fo,[0.])
fo.close()

fo = open(FO_DIR+'/StateVector.csv','w')
savetxt(fo,[StartValues],delimiter=',')
fo.close()
#}





#-------------------------------------------------------------------------
# region INITIALIZATION FOR PLOTTING
# Plotting Setting
FlagColorToggle=True
FlagColorToggleVal=0

# Output Setting
FlagHeader=True
FlagHeader_Shape=True

# Initialize Postprocess-Core Shape
Post_CoreShapeInit (Core)

# Save Initial Core Shape, Solution
#Post_CoreShape([0],[StartValues],Core,FO_DIR)
#Accel = ExtractAccelFromXV([StartValues],[0],Core)
#Post_Solution([0],[StartValues],Accel,Core,FO_DIR)

# endregion INITIALIZATION FOR PLOTTING ---------------------------------------





#------------------------------------------------------------------------------
# region MEMORY ALLOCATION OF ACCEL
MallocAccel(Core)
# endregion MEMORY ALLOCATION OF ACCEL ----------------------------------------





#------------------------------------------------------------------------------
Verbose('\n\n=== INITIALIZATION ENDED ===\n\n')
raw_input('PRESS ENTER TO CONTINUE!')

# endregion ######################## INITIALIZATION ##############################





# region ======================= TEST FORCES ONLY =============================
#
# Objective: For module verification, calculate initial acceleration only.
#
# Usage:
#   Set FLAG_TestForce=True in Input file to invoke this test.
#   Set FLAG_TestForce=False to skip this test.
#   
# ============================================================================

if( FLAG_TestForce==True ):
#{  
  Accel = VectorField (StartValues,0,Core)
  print '='*80
  for IndexBlock in range(len(Core['ReverseIndex'])):
    K,L=Core['ReverseIndex'][IndexBlock]
    IndexW = IndexBlock*6
    AccelBl = Accel[IndexW:IndexW+6]
    Verbose('(K,L)=(%d,%d) : Accel=%s'%(K,L,AccelBl))
  #print 'in main : Accel=',Accel

  print '='*80
  
  print '( K- L), U       , DU      , DDU     , W       , DW      , DDW     , R        , DR      , DDR'
  for i in xrange(len(Core['ReverseIndex'])):
    (K,L)=Core['ReverseIndex'][i]
    BlockState=GetSolFromVect (Core,StartValues,K,L)
    BlockAccel=GetSolFromVect (Core,Accel,K,L)
    Msg  = '(%2d-%2d),'%(K,L)
    Msg += '%9.2e,%9.2e,          %9.2e,%9.2e,           %9.2e,%9.2e\n'%(BlockState)
    Msg += '       ,'
    Msg += '          %9.2e,%9.2e,          %9.2e,%9.2e,           %9.2e,%9.2e'%BlockAccel
    print Msg
  
  Verbose('='*80)
  
  Verbose('( K- L), DU        , DDU       , DW        , DDW       , DR        , DDR')
  for i in xrange(len(Core['ReverseIndex'])):
    (K,L)=Core['ReverseIndex'][i]
    BlockAccel=GetSolFromVect (Core,Accel,K,L)
    Msg  = '(%2d-%2d),'%(K,L)
    for j in range(6):
      if BlockAccel[j]==0.0:
        Msg += ' '*10+'0,'
      else:
        Msg += '%11.4e,'%BlockAccel[j]
#    Msg += '%11.4e,%11.4e,%11.4e,%11.4e,%11.4e,%11.4e'%BlockAccel
    Verbose(Msg)
  
  Msg  = '=== Performance Measurement ===\n'
  Msg += 'VectorField() : NRun=%d, RunTime=%f, AvgTime=%f'%PerfResult_VectorField()
  print Msg


  raw_input('\n\nTEST FORCES ENDED. PRESS ENTER TO QUIT!')
  sys.exit(-1)
#}

# endregion ======================= TEST FORCES ONLY =============================





# region ######################### M A I N   L O O P ##########################

#==============================================================================
#{ WHILE CURRENTTIME<TOTALTIME

while CurrentTime<TotalTime:

  Msg  = 'Current Time = %f'%CurrentTime
  Verbose(Msg)
  
  # Prepare Integration #{
  # Reset SectionTime
  SectionTime=SectionTime0
  Tol=Tol0
  dt=TimeIncr # from INPUT
  #}

  #-----------------------------------------------------------------------------
  #  CONVERGE CONTROL LOOP #{
  while True:

#    print '  SectionTime=%f'%SectionTime

    # Integration
#    t=[SectionTime*float(i)/(NPoints-1)+CurrentTime for i in range(NPoints)]
#    print '  t[0]=%f, t[-1]=%f, len(t)=%d'%(t[0],t[-1],len(t))

    # PREPARE T,SOLUTION #{
    temp=SectionTime/dt
    NPoints=int(temp)
    if NPoints<temp: NPoints+=1
    t=[dt*i+CurrentTime for i in range(NPoints+1)]
    #}

    # SOLVE #{
    Solution,Status = odeint(VectorField, StartValues, t, args=(Core,VERBOSE), atol=Tol, rtol=Tol, full_output=True)
    #}

    # If odeint succeeded, break and go to next CurrentTime #{
    if 'successful' in Status['message']: break
    #}

    #---------------------------------------------------------------------------
    # If odeint failed, #{
    
    Verbose('CONVERGENCE FAILED')      
    
    if SectionTime/2. >= SectionTimeMin:
      
    # SectionTime Control (Halfing)
      SectionTime/=2.
      Verbose('Half the SectionTime = %f'%SectionTime)
    else:
      
    # Tolerance Control (Doubling)
      while True:
#        Tol*=2.
#        Verbose('Double the Tol = %e'%Tol)
        Tol*=10.
        Verbose('10 times the Tol = %e'%Tol)

        # Bound Check
        if Tol>TolMax: # Error
          print '='*70
          print 'ERR: Tol(%e) > TolMax(%e)'%(Tol,TolMax)
          print 'STOP'
          print '='*70
          sys.exit(-1)
        
        Solution,Status = odeint(VectorField, StartValues, t, args=(Core,VERBOSE),
            atol=Tol, rtol=Tol, full_output=True)
      
        if 'successful' in Status['message']:
        
          # Try average of previous Tol and current Tol
          Tol*=0.75
          Solution2,Status2 = odeint(VectorField, StartValues, t, args=(Core,VERBOSE),
              atol=Tol, rtol=Tol, full_output=True)
          if 'successful' in Status2['message']:
            Solution = Solution2
          break
    # END OF If odeint failed
    #}--------------------------------------------------------------------------

    
  #  END OF INIFINITE LOOP  
  #}----------------------------------------------------------------------------



  #-----------------------------------------------------------------------------
  #  PERFORMANCE MEASUREMENT #{
  a,b,c,d,e,f = PerfResult_VectorField()
  Msg  = '  === Performance Measurement ===\n'
  Msg += '  VectorField() of this section: NRun=%d, RunTime=%.1f, Avg=%.1fms\n'%(a,b,c*1000)
  Msg += '  VectorField() of total run   : NRun=%d, RunTime=%.1f, Avg=%.1fms'%(d,e,f*1000)
  Verbose(Msg)
  #}----------------------------------------------------------------------------



  #-----------------------------------------------------------------------------
  #   Save Results #{

  if OP_CoreShape_TimeFreq != 0:
    print '  Post_CoreShape',
    Post_CoreShape(t,Solution,Core,FO_DIR)
    print 'Finished'

#  Accel = ExtractAccelFromXV(Solution,t,Core)
#  print '3'

  if OP_Solution_TimeFreq != 0:
    print '  Post_Solution',
#    Post_Solution(t,Solution,Accel,Core,FO_DIR)
    Post_Solution(t,Solution,None,Core,FO_DIR)
    print 'Finished'

  if OP_Block_TimeFreq != 0:
    print '  Post_Block',
#    Post_Block(t,Solution,Accel,Core,FO_DIR)
    Post_Block(t,Solution,None,Core,FO_DIR)
    print 'Finished'
  
  #{ SAVE STATEVECTOR (DEFAULT)
  fo = open(FO_DIR+'/Time.csv','a')
  savetxt(fo,t[1:])
  fo.close()

  fo = open(FO_DIR+'/StateVector.csv','a')
  savetxt(fo,Solution[1:],delimiter=',')
  fo.close()
  #}
  
  #}----------------------------------------------------------------------------



  # DELETED #{
  # PLOT CURRENT INTEGRATION RESULT  
  # PLOT CURRENT INTEGRATION RESULT  
 
  # Toggle color for each section
#  if FlagColorToggle==True:
#    if FlagColorToggleVal==0:
#      c1='b-'
#      c2='g-'
#      FlagColorToggleVal=1
#    else:
#      c1='r-'
#      c2='c-'
#      FlagColorToggleVal=0

  # Plot
#  x1=Solution[:,18] # 1st block: X(18),DX(19),W(20),DW(21),R(22),DR(23)
#  x2=Solution[:,24] # 2nd block: X(24),DX(25),W(26),DW(27),R(28),DR(29)
#  x3=Solution[:,32] # 3rd block: X(30),DX(31),W(32),DW(33),R(34),DR(35)
#  w1=Solution[:,20] # 1st block: X(18),DX(19),W(20),DW(21),R(22),DR(23)
#  w2=Solution[:,26] # 2nd block: X(24),DX(25),W(26),DW(27),R(28),DR(29)
#  w3=Solution[:,32] # 3rd block: X(30),DX(31),W(32),DW(33),R(34),DR(35)
#  r1=Solution[:,22] # 1st block: X(18),DX(19),W(20),DW(21),R(22),DR(23)
#  r2=Solution[:,28] # 2nd block: X(24),DX(25),W(26),DW(27),R(28),DR(29)
#  r3=Solution[:,34] # 3rd block: X(30),DX(31),W(32),DW(33),R(34),DR(35)
#  v1=Solution[:,19] # 1st block: X(18),DX(19),W(20),DW(21),R(22),DR(23)
#  v2=Solution[:,25] # 2nd block: X(24),DX(25),W(26),DW(27),R(28),DR(29)
#  v3=Solution[:,33] # 3rd block: X(30),DX(31),W(32),DW(33),R(34),DR(35)
#  a1=Accel[:,9] # 1st block: DDX(9),DDW(10),DDR(11)
#  a2=Accel[:,12] # 2nd block: DDX(12),DDW(13),DDR(14)
#  a3=Accel[:,15] # 3rd block: DDX(15),DDW(16),DDR(17)



  
  # Block 1 Shape
#  if FlagHeader_Shape==True:
#    text = 't, LD x, RD x, RU x, LU x, LD x, LD z, RD z, RU z, LU z, LD z\n'
#    fo_shape.write(text)
#    FlagHeader_Shape=False
#  NDiv = len(x1)/10
#  NDiv = len(x1)
#  for i in range(0,len(x1),NDiv):
#    T = t[i]
#    U,W,R=x1[i],w1[i],r1[i]
#    Cent,LU,LD,RU,RD = GetBlockCoords(Core,1,1,U,W,R,OP_CoreShape_Scale)
#    text = '%e,'%T
#    text += '%e,%e,%e,%e,%e,'%(LD[0],RD[0],RU[0],LU[0],LD[0])
#    text += '%e,%e,%e,%e,%e\n'%(LD[1],RD[1],RU[1],LU[1],LD[1])
#    fo_shape.write(text)
   
    
    
  
  
#  ion()
#  figure(1, figsize=(20, 3))
#  plot(t, x1, c1, linewidth=1)
#  savefig(FN_FIG+'1_%e.png'%t[-1])
#  figure(2, figsize=(20, 3))
#  plot(t, r1, c2, linewidth=1)
#  savefig(FN_FIG+'2_%e.png'%t[-1])
#  figure(3, figsize=(20, 3))
#  plot(t, a1, c1, linewidth=1)
#  savefig(FN_FIG+'3_%e.png'%t[-1])
#  figure(4, figsize=(20, 3))
#  plot(t, a2, c1, linewidth=1)
#  savefig(FN_FIG+'4_%e.png'%t[-1])






  # WRITE OUTPUT  
  # WRITE OUTPUT  
  # WRITE OUTPUT  
  
#  if FlagHeader==True:
##    text = 't, x1, w1, r1, x2, w2, r2, v1, v2, v3, a1, a2, a3\n'
#    text='t,'
#    for (K,L) in Core['ReverseIndex']:
#      text+='(%d-%d) U,'%(K,L)
#      text+='(%d-%d) DU,'%(K,L)
#      text+='(%d-%d) W,'%(K,L)
#      text+='(%d-%d) DW,'%(K,L)
#      text+='(%d-%d) R,'%(K,L)
#      text+='(%d-%d) DR,'%(K,L)
#    text+='\n'
#    fo.write(text)
#    FlagHeader=False

  
#  for i in xrange(len(x1)):
#    text  = '%e,%e,%e,%e,' %(t[i],x1[i],w1[i],r1[i])
#    text +=   '%e,%e,%e,' %(     x2[i],w2[i],r2[i])
#    text +=   '%e,%e,%e,' %(     v1[i],v2[i],v3[i])
#    text +=   '%e,%e,%e\n'%(     a1[i],a2[i],a3[i])

#  for i in xrange(len(t)):
#    text='%e,'%t[i]
#    for Sol in Solution[i]:
#      text+='%e,'%Sol
#    text+='\n'
#    fo.write(text)

  #}

  # Next Step #{
  CurrentTime=t[-1]
  StartValues=Solution[-1]
  #}
  
#} END OF WHILE CURRENTTIME<TOTALTIME
#==============================================================================


# endregion ######################### M A I N   L O O P ##########################






# region ======================= TIME MEASUREMENT =============================

end_time = time.time()
run_time_min = (end_time-start_time) / 60
run_time_sec = int((end_time-start_time) % 60)
Verbose('='*40)
Verbose('Finished at '+time.ctime(end_time))
Verbose('Run Time : %d min %d sec'%(run_time_min,run_time_sec))

# endregion ======================= TIME MEASUREMENT =============================
