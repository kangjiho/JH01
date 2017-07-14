#
# SAPCOR-POST.py
#
"""Seismic Analysis POST for a Prismatic CORe of a HTGR
"""

VERSION_NAME = 'SAPCOR_v04h'

###############################################################################
#{    SYSTEM INITIALIZATION

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
FO_DIR = 'POST_%4d-%02d-%02d_%02dh%02dm%02ds'%(year,month,day,hour,minute,second)
os.mkdir(FO_DIR)
VerboseInit(FO_DIR)
#} Make Output Folder =========================================================

Verbose('//INIT//')

#}    SYSTEM INITIALIZATION
###############################################################################


###############################################################################
#{                               LOCAL FUNCTIONS                              #
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

def ERROR (ERR_MSG_LIST,Align='Left'):
#{
  Msg  = '\n'*3+'='*79+'\n'
  Msg += '='+' '*34+'E R R O R'+' '*34+'=\n'
  Msg += '='*79+'\n'
  Msg += '='+' '*77+'=\n'

  if Align == 'Left':
    # LEFT-ALIGNED FORMAT
    for ERR_MSG in ERR_MSG_LIST:
      Length = len(ERR_MSG)
      NBlankR = 75 - Length
      Msg += '='+'  '+ERR_MSG+' '*NBlankR+'=\n'
  else:
    # CENTERED FORMAT
    for ERR_MSG in ERR_MSG_LIST:
      Length = len(ERR_MSG)
      NBlankL = NBlankR = 38 - Length/2
      if Length%2 == 0: NBlankR += 1
      Msg += '='+' '*NBlankL+ERR_MSG+' '*NBlankR+'=\n'

  Msg += '='+' '*77+'=\n'
  Msg += '='*79+'\n'
  Msg += '\n'*3
  Verbose(Msg)  
  raise
  return
#}

#}                               LOCAL FUNCTIONS                              #
###############################################################################


###############################################################################
#{                               INITIALIZATION                               #
#                                                                             #
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
#{
  #{ If it is CSB
  if Column=='CSB':
    i_Index+=1 # If CSB exists, reserve Index[1][0]
    Index.append([1,0]) # CSB
    Core['Flag_CSB']=True
  #}

  #{ If it is Reflectors or Blocks
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
  #}
#}

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



#-------------------------------------------------------------------------
#{ INITIALIZATION FOR PLOTTING
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

#}-------------------------------------------------------------------------



Verbose('\n\n=== INITIALIZATION ENDED ===\n\n')
raw_input('PRESS ENTER TO CONTINUE!')

#                                                                             #
#}                               INITIALIZATION                               #
###############################################################################



###############################################################################
#{                              M A I N   L O O P                             #
#                                                                             #
#                                                                             #

#-----------------------------------------------------------------------------
#{ READ STATE VECTOR
try:
  Ts=loadtxt('Time.csv')
  StateVectors=loadtxt('StateVector.csv',delimiter=',')
except:
  ERROR(["Cannot find 'Time.csv' and 'StateVector.csv'.",
         "Please copy both files to this directory."])
temp = zip(Ts,StateVectors)
#}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#{ LOOP
for t,StateVector in temp:

  Msg = 'Time = %g'%t
  Verbose(Msg)  
  Accel = VectorField(StateVector, t, Core, VERBOSE='', POST=True, FO_DIR=FO_DIR)
  

#} END OF LOOP
#-----------------------------------------------------------------------------
 



#-----------------------------------------------------------------------------
#  PERFORMANCE MEASUREMENT #{
a,b,c,d,e,f = PerfResult_VectorField()
Msg  = '=== Performance Measurement ===\n'
Msg += 'VectorField() of this section: NRun=%d, RunTime=%.1f, Avg=%.1fms\n'%(a,b,c*1000)
Msg += 'VectorField() of total run   : NRun=%d, RunTime=%.1f, Avg=%.1fms'%(d,e,f*1000)
Verbose(Msg)
#}----------------------------------------------------------------------------


#                                                                             #
#                                                                             #
#}                              M A I N   L O O P                             #
###############################################################################





# =============================================================================
#{    TIME MEASUREMENT

end_time = time.time()
run_time_min = (end_time-start_time) / 60
run_time_sec = int((end_time-start_time) % 60)
Verbose('='*40)
Verbose('Finished at '+time.ctime(end_time))
Verbose('Run Time : %d min %d sec'%(run_time_min,run_time_sec))

#}    TIME MEASUREMENT
# =============================================================================
