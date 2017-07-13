VERSION_NAME = 'SAPCOR_v04j'

# MODULE NAME - MISC
MODULE_NAME_THIS_MISC = VERSION_NAME+'_Misc'


from numpy import array,zeros,savetxt,column_stack
import glob
from pylab import figure,plot,axis,savefig,clf

# -------------------------------------------------------------------
#{ MODULE NAMES TO BE IMPORTED
MODULE_NAMES_MISC = []
MODULE_NAMES_MISC.append('Input')
#} -------------------------------------------------------------------

# -------------------------------------------------------------------
#{ Import Modules
for Module_Name in MODULE_NAMES_MISC:
  FileName = VERSION_NAME+'_'+Module_Name+'_Build*.py'
  FList = glob.glob(FileName)
  assert len(FList)!=0,'ERROR: Module not found : '+FileName
  FList.sort()
  FileName = FList[-1][0:-3]
  exec('from '+FileName+' import *')
#} -------------------------------------------------------------------





import sys
from math import sin,cos,pi,exp


# ------------------------------------------------------------------------
#{  GLOBAL PARAMETERS

# Output Dir
FO_DIR = ''
 
# Incresing number of CoreShape figures
NPost_CoreShape = 0

# CoreShape output time pointer
OP_CoreShape_TimePointer = 0.

# Incresing number of Solution Files
NPost_Solution = 0

# CoreShape output time pointer
OP_Solution_TimePointer = 0.

# Block output time pointer
OP_Block_TimePointer = 0.

#} ------------------------------------------------------------------------



#-------------------------------------------------------------------------
def ERROR (msg,code=-1): #{
  Verbose(msg)
  sys.exit(code)

#}-------------------------------------------------------------------------
def VerboseInit (OUTPUT_DIR): #{
  global FO_DIR
  FO_DIR=OUTPUT_DIR
  return
 
#}-------------------------------------------------------------------------
def Verbose (Msg='',NewLine=True): #{
  
#(v04f)  FN_VERBOSE = FN_OUT+'.log'
  global FO_DIR
  FN_VERBOSE = FO_DIR+'/'+FN_OUT+'.log'

  if Msg=='//INIT//':
    fp=file(FN_VERBOSE,'w',0)
    fp.close()
    return

  if NewLine==True: Msg+='\n'
    
  print Msg,
  
  fp=file(FN_VERBOSE,'a',0)
  fp.write(Msg)
  fp.close()
  
  return

#}-------------------------------------------------------------------------




#=========================================================================
#      POST - CORE SHAPE
#=========================================================================
def Post_CoreShape (t,Solution,Core,FO_DIR,VERBOSE=''): #{

  global NPost_CoreShape
  global OP_CoreShape_TimePointer
  
  # Determine time points for output
  IndicesTime = []
  temp = range(len(t))
  for i in temp:
    if OP_CoreShape_TimePointer <= t[i]:
      IndicesTime.append(i)
#      print 'time=%f, index=%d, t[index]=%f'%(OP_CoreShape_TimePointer,i,t[i])
      OP_CoreShape_TimePointer += OP_CoreShape_TimeFreq

  # Prepare Figure Window  
  figure(100,(10,10))

  # Main Loop
  for IndexTime in IndicesTime:

    # Get Final Results    
    Ti = t[IndexTime]
    Si = Solution[IndexTime]
      
    # Clear Figure
    clf()
      
    # Clear CSV
    ShapeCSV = 'x,y\n'
    
    # CSB
    if Core['Flag_CSB']==True:
      IndexBlock = 1
      IndexW = IndexBlock*6
  
      U=Si[IndexW+0]*OP_CoreShape_Scale
      W=Si[IndexW+2]*OP_CoreShape_Scale
      R=Si[IndexW+4]
  
      Coords = GetBlockCoords(Core,0,1,U,W,R)
      LU = [Coords[2][0],Coords[2][1]]
      LD = [Coords[2][0],Coords[2][1]-BlockTypes['CSB']['h']]
      Coords = GetBlockCoords(Core,2,1,U,W,R)
      RU = [Coords[4][0],Coords[4][1]]
      RD = [Coords[4][0],Coords[4][1]-BlockTypes['CSB']['h']]
      ShapeCSV+='%f,%f,"CSB"\n'%(LU[0],LU[1])
      ShapeCSV+='%f,%f\n'%(RU[0],RU[1])
      ShapeCSV+='%f,%f\n'%(RD[0],RD[1])
      ShapeCSV+='%f,%f\n'%(LD[0],LD[1])
      ShapeCSV+='%f,%f\n'%(LU[0],LU[1])
      ShapeCSV+=',\n'
        
      x=[LU[0],RU[0],RD[0],LD[0],LU[0]]
      y=[LU[1],RU[1],RD[1],LD[1],LU[1]]
      plot(x,y)
      
    # Blocks
    LCoreReverseIndex = len(Core['ReverseIndex'])
    for IndexBlock in xrange(LCoreReverseIndex):
        
      K,L = Core['ReverseIndex'][IndexBlock]
      IndexW = IndexBlock*6 
  
      # Skip Condition
      if (K,L)==(0,0): continue # Base
      if (K,L)==(1,0): continue # CSB
  
      U=Si[IndexW+0]*OP_CoreShape_Scale
      W=Si[IndexW+2]*OP_CoreShape_Scale
      R=Si[IndexW+4]
        
      Coords = GetBlockCoords(Core,K,L,U,W,R)
      LU = Coords[1]
      LD = Coords[2]
      RU = Coords[3]
      RD = Coords[4]
      ShapeCSV+='%f,%f,"(%d,%d)"\n'%(LU[0],LU[1],K,L)
      ShapeCSV+='%f,%f\n'%(RU[0],RU[1])
      ShapeCSV+='%f,%f\n'%(RD[0],RD[1])
      ShapeCSV+='%f,%f\n'%(LD[0],LD[1])
      ShapeCSV+='%f,%f\n'%(LU[0],LU[1])
      ShapeCSV+=',\n'
  
      x=[LU[0],RU[0],RD[0],LD[0],LU[0]]
      y=[LU[1],RU[1],RD[1],LD[1],LU[1]]
      plot(x,y)
   
    
    # Save Core Shape in CSV format
    FnTimeFormat = '%05d_%10.4e'%(NPost_CoreShape,Ti)
  #  fo = open(FN_CoreShape+'_t%e.csv'%Ti,'w')
    fo = open(FO_DIR+'/'+FN_CoreShape+'_'+FnTimeFormat+'.csv','w')
    fo.write(ShapeCSV)
    fo.close()  # Accelerations
    
    # Save Core Shape in PNG format
    axis('image')
    axis(Core['Post_CoreShapeAxis'])
    savefig(FO_DIR+'/'+FN_CoreShape+'_'+FnTimeFormat+'.png')
  
    NPost_CoreShape += 1
    
    
  return
#}========================================================================


# ===========================================================================
# Get Initial Core Shape and Axis Info
# ===========================================================================
def Post_CoreShapeInit (Core,VERBOSE=''): #{

  figure(100,(10,10))
  clf()
  
  # CSB
  if Core['Flag_CSB']==True:
    U,W,R=0,0,0
    Coords = GetBlockCoords(Core,0,1,U,W,R)
    LU = [Coords[2][0],Coords[2][1]]
    LD = [Coords[2][0],Coords[2][1]-BlockTypes['CSB']['h']]
    Coords = GetBlockCoords(Core,Core['M'],1,U,W,R)
    RU = [Coords[4][0],Coords[4][1]]
    RD = [Coords[4][0],Coords[4][1]-BlockTypes['CSB']['h']]
    
    x=[LU[0],RU[0],RD[0],LD[0],LU[0]]
    y=[LU[1],RU[1],RD[1],LD[1],LU[1]]
    plot(x,y)
  
  # Blocks
  LCoreReverseIndex = len(Core['ReverseIndex'])
  for IndexBlock in xrange(LCoreReverseIndex):
    
    K,L = Core['ReverseIndex'][IndexBlock]

    # Skip Condition
    if (K,L)==(0,0): continue # Base
    if (K,L)==(1,0): continue # CSB

    U,W,R=0,0,0
    Coords = GetBlockCoords(Core,K,L,U,W,R)
    LU = Coords[1]
    LD = Coords[2]
    RU = Coords[3]
    RD = Coords[4]

    x=[LU[0],RU[0],RD[0],LD[0],LU[0]]
    y=[LU[1],RU[1],RD[1],LD[1],LU[1]]
    plot(x,y)
  
  axis('image')
  temp=axis()

# X & Y Margin : Change from scales of figure size to absolute length
#  scale=OP_CoreShapeScaleMargin
#  temp=(temp[0]-(temp[1]-temp[0])*scale,temp[1]+(temp[1]-temp[0])*scale,temp[2]-(temp[3]-temp[2])*scale,temp[3]+(temp[3]-temp[2])*scale)
  temp=( temp[0] - OP_CoreShape_MarginX,
         temp[1] + OP_CoreShape_MarginX,
         temp[2] - OP_CoreShape_MarginY,
         temp[3] + OP_CoreShape_MarginY)
  
  axis(temp)
  
  Core['Post_CoreShapeAxis']=temp
  
  return

#}===========================================================================



# ===========================================================================
def OLD_Post_CoreShape_FnTimeFormat (Increment,t): #{

  # ----------------------------------------------------------------------
  # Time Format in File Name  

  # Left of point
  temp = TotalTime
  LeftOfPoint = 0
  if TotalTime >= 1:
    while temp >= 1:
      temp /= 10.
      LeftOfPoint += 1 
  else: LeftOfPoint = 1
  
  # Right of point
  if Increment >= len(t): # if OP_CoreShapeNFreq==1 -> DeltaTi=SectionLength
    DeltaTi = t[-1]+t[1]-2*t[0] # SectionLength=t[-1]-t[0]+(t[1]-t[0])
  else:
    DeltaTi = t[Increment]-t[0]
  RightOfPoint = 0
  if DeltaTi>=1: RightOfPoint = 0
  if DeltaTi<1:
    while DeltaTi<1:
      DeltaTi *= 10.
      RightOfPoint += 1
  RightOfPoint += 1

  # Time format
  FnTimeFormat = '%%%d.%df'%(LeftOfPoint+RightOfPoint+1,RightOfPoint)

  # End of Time Format in File Name  
  # ----------------------------------------------------------------------

  return FnTimeFormat

#      END OF POST - CORE SHAPE
#}========================================================================





#=========================================================================
#      POST - SOLUTION AT LAST OF SECTION
#=========================================================================
def Post_Solution (t,Solution,Accel,Core,FO_DIR,VERBOSE=''): #{

  global NPost_Solution
  global OP_Solution_TimePointer
  
  # Determine time points for output
  IndicesTime = []
  temp = range(len(t))
  for i in temp:
    if OP_Solution_TimePointer <= t[i]:
      IndicesTime.append(i)
      OP_Solution_TimePointer += OP_Solution_TimeFreq

  # Main Loop
  for IndexTime in IndicesTime:
  
    Ti = t[IndexTime]
    Si = Solution[IndexTime]
    if Accel!=None: Ai = Accel[IndexTime]
      
    # Clear CSV
    if Accel!=None: SolutionCSV = 'K,L,U,W,R,DU,DW,DR,DDU,DDW,DDR\n'
    else:  SolutionCSV = 'K,L,U,W,R,DU,DW,DR\n'
    
    # Loop for Blocks
    LCoreReverseIndex = len(Core['ReverseIndex'])
    for IndexBlock in xrange(LCoreReverseIndex):
        
      K,L = Core['ReverseIndex'][IndexBlock]
      IndexW = IndexBlock*6
      IndexA = IndexBlock*3
  
      # K,L
      SolutionCSV+='%d, %d, '%(K,L)
        
      # Displacement
      SolutionCSV+='%e, %e, %e, '%(Si[IndexW],Si[IndexW+2],Si[IndexW+4])
  
      # Velocity
      SolutionCSV+='%e, %e, %e, '%(Si[IndexW+1],Si[IndexW+3],Si[IndexW+5])
        
      # Acceleration
      if Accel!=None: 
        SolutionCSV+='%e, %e, %e\n'%(Ai[IndexA],Ai[IndexA+1],Ai[IndexA+2])
      else: SolutionCSV+='\n'
  
    # Save Solution
    FnTimeFormat = '%05d_%10.4e'%(NPost_Solution,Ti)
    fo = open(FO_DIR+'/'+FN_Solution+'_'+FnTimeFormat+'.csv','w')
    fo.write(SolutionCSV)
    fo.close()
    
    NPost_Solution += 1
    
  return

#      END OF POST - SOLUTION
#}========================================================================



#=========================================================================
#      POST - DISP,ACCEL BY BLOCK
#=========================================================================
def Post_Block (t,Solution,Accel,Core,FO_DIR,VERBOSE=''): #{

  global OP_Block_TimePointer
  
  # Determine time points for output
  IndicesTime = []
  temp = range(len(t))
  for i in temp:
    if OP_Block_TimePointer <= t[i]:
      IndicesTime.append(i)
      OP_Block_TimePointer += OP_Block_TimeFreq

  # Main Loop
  # Loop for Blocks
  LCoreReverseIndex = len(Core['ReverseIndex'])
  for IndexBlock in xrange(LCoreReverseIndex):
  
    K,L = Core['ReverseIndex'][IndexBlock]
    IndexW = IndexBlock*6
    IndexA = IndexBlock*3
    
    # Clear CSV
    BlockCSV = ''

    # Loop for Time    
#    for IndexTime in xrange(len(t)):
    for IndexTime in IndicesTime:
      
      Ti = t[IndexTime]
      Si = Solution[IndexTime]
      if Accel!=None: Ai = Accel[IndexTime]
      LocationBase = (Si[0],Si[2])

      # Get Location
      if K==0 and L==0: # Base
        Location = LocationBase
      else:
        U=Si[IndexW+0]
        W=Si[IndexW+2]
        R=Si[IndexW+4]
        Coords = GetBlockCoords(Core,K,L,U,W,R)
        Location = Coords[0]
        Location[0] -= LocationBase[0]
        Location[1] -= LocationBase[1]
      
      # Time
      BlockCSV+='%12.5e, '%Ti
      
      # Displacement
      BlockCSV+='%12.5e, %12.5e, %12.5e, '%(Si[IndexW],Si[IndexW+2],Si[IndexW+4])

      # Location
      BlockCSV+='%12.5e, %12.5e, '%(Location[0],Location[1])
      
      # Acceleration
      if Accel!=None:
        BlockCSV+='%12.5e, %12.5e, %12.5e\n'%(Ai[IndexA],Ai[IndexA+1],Ai[IndexA+2])
      else: BlockCSV+='\n'

    # Save
    fo = open(FO_DIR+'/(%2d,%2d).csv'%(K,L),'a')
    fo.write(BlockCSV)
    fo.close()
  
  return

#      END OF POST - DISP,ACCEL BY BLOCK
#}========================================================================




#=========================================================================
#{        INITIALIZE
#=========================================================================

#-------------------------------------------------------------------------
#
# DO NOT USE !!!  Concept of INDEX was changed in v0.4.
# DO NOT USE !!!  Concept of INDEX was changed in v0.4.
# DO NOT USE !!!  Concept of INDEX was changed in v0.4.
# DO NOT USE !!!  Concept of INDEX was changed in v0.4.
#
# Index and Rearrange Array
# Result: Make Core['Index']
#         Make Core['Flag_CSB']
#         Make Core['Array']
def MakeIndex_RearrangeArray (Core,CoreArray,VERBOSE=''): #{
  
  CoreArray2={}
  Index=[]
  Index.append([0,0]) # Base
  K=0
  Core['Flag_CSB']=False
  for Column in CoreArray:
    if Column=='CSB':
      Index.append([1,0]) # CSB
      Core['Flag_CSB']=True
    else:
  #    CoreArray2[K]=[None,]
      CoreArray2[K]={}
      for i in range(len(Column)):
        L=i+1
        Index.append([K,L])
  #      CoreArray2[K].append(Column[i])
        CoreArray2[K][L]=Column[i]
      K+=1
  Core['M']=K-2
  
  if '/Init/' in VERBOSE:
    print '<Rearranged CoreArray>'
    for i in range(len(CoreArray2)):
      print CoreArray2[i]
  
  Core['Index']=Index
  Core['Array']=CoreArray2

#} END OF Index and Rearrange Array
#-------------------------------------------------------------------------




#-------------------------------------------------------------------------
# Connectivity
# Tested: 140303_Connectivity Test #1.ppt
# Result: Make Core['Connectivity']
def MakeConnectivity (Core,VERBOSE=''): #{
  
  # Var Spaces
  BlockTypes = Core['BlockTypes']
  CoreArray = Core['Array']
  
  # Init
  Connectivity=[]
 
  # Loop to explore columns
  for K in range(len(CoreArray)):
    
    # Init
    Connectivity.append([])
    Z_KL_LU = 0.
  
    # Loop to explore blocks
    for L in range(len(CoreArray[K])):
      
      # Allocate storage
      Connectivity[K].append(None)
      
      # Zero-th row is reserved for special elements (Base,CSB)
      # They do not need Connectivity
      # Fill None in Connectivity for them
      if L==0: continue
    
      # BlockType of block (K,L)
      BlockType = CoreArray[K][L]
      
      # Z coords on LD,LU corner of block (K,L)
      Z_KL_LD = Z_KL_LU
      Z_KL_LU += BlockTypes[BlockType]['h']*2
      Z_KL_RD = Z_KL_LD
      Z_KL_RU = Z_KL_LU
      
      # Init Temp Storages
      ConnLU,ConnLD,ConnRU,ConnRD=None,None,None,None
  
      # Check Left
      # Check Left
      if K>=1:  # If K==0 -> Left Restraint -> Skip
  
        # Init Z coord of lowest block in (K-1) column
        Z_Km1l_RU = 0
       
        # Flag for Conn for LD found
        Flag_LD_Found=False
        
        # Loop to explore blocks in (K-1) column
        for l in range(1,len(CoreArray[K-1])):
  
          # BlockType of block (K-1,l)
          BlockType = CoreArray[K-1][l]
        
          # Update Z coord on RU corner of block (K-1,l) 
          Z_Km1l_RU += BlockTypes[BlockType]['h']*2
          
          # If Conn for LD not found yet -> Check for LD only
          if Flag_LD_Found==False:
            if Z_Km1l_RU > Z_KL_LD:
              Flag_LD_Found=True # Toggle Flag to begin check for LU 
              ConnLD=l
              # DO NOT CONTINUE OR BREAK HERE!!!
          
          # If Conn for LD was found already -> Check for LU only
          if Flag_LD_Found==True:
            if Z_Km1l_RU >= Z_KL_LU:
              ConnLU=l
              break # To be here, both ConnLD, ConnLU were found -> break
  
          # Next upper block in K-1 column
          l+=1
        # End of Check Left
        # End of Check Left
        
      # Check Right
      # Check Right
      if K<=Core['M']:  # If K==M+1 -> Right Restraint -> Skip
  
        # Init Z coord of lowest block in (K+1) column
        Z_Kp1l_LU = 0
       
        # Flag for Conn for RD found
        Flag_RD_Found=False
        
        # Loop to explore blocks in (K+1) column
        for l in range(1,len(CoreArray[K+1])):
  
          # BlockType of block (K+1,l)
          BlockType = CoreArray[K+1][l]
        
          # Update Z coord on LU corner of block (K-1,l) 
          Z_Kp1l_LU += BlockTypes[BlockType]['h']*2
          
          # If Conn for RD not found yet -> Check for RD only
          if Flag_RD_Found==False:
            if Z_Kp1l_LU > Z_KL_RD:
              Flag_RD_Found=True # Toggle Flag to begin check for RU 
              ConnRD=l
              # DO NOT CONTINUE OR BREAK HERE!!!
          
          # If Conn for RD was found already -> Check for RU only
          if Flag_RD_Found==True:
            if Z_Kp1l_LU >= Z_KL_RU:
              ConnRU=l
              break # To be here, both ConnLD, ConnLU were found -> break
  
          # Next upper block in K-1 column
          l+=1
        # End of Check Right
        # End of Check Right
  
      Connectivity[K][L]=[ConnLU,ConnLD,ConnRU,ConnRD]
  
  if '/Init/' in VERBOSE:
    print '<Connectivity>'
    for K in range(Core['M']+2):
      for L in range(1,len(CoreArray[K])):
        print '(%d,%d)='%(K,L),Connectivity[K][L]

  Core['Connectivity']=Connectivity

#} END OF Connectivity
#-------------------------------------------------------------------------





#-------------------------------------------------------------------------
# Get Initial Block Center Coords
# Origin at LD corner of first lowest left restratin
# Results: Core['InitialBlockCenters'][K][L]=[x,z]
def GetInitBlockCenters (Core,VERBOSE=''): #{
  
  # Var Space
  BlockTypes=Core['BlockTypes']
  Array=Core['Array']
  Flag_CSB=Core['Flag_CSB']
  
  # Init Storage
  BlockCenters={}
  
  # No of columns
  M2 = len(Array)
  
  # Init X
  X = 0.
  
  # Loop for columns
  for K in range(M2):

    # Init Storage    
    BlockCenters[K]={}
    
    # BlockType of (K,1)
    BlockType = Array[K][1]
    
    # Block half width
    Width2 = BlockTypes[BlockType]['b']
    
    # Add Block half width
    X += Width2
    
    # Add Gap
    X += Delta
    
    # No of blocks in column
    N = len(Array[K])
    
    # Init Z
    Z = 0.
    
    # Loop for blocks
    for L in range(1,N+1):
    #{
      # BlockType
      BlockType = Array[K][L]
      
      # Add lower half
      Z+=BlockTypes[BlockType]['h']
      
      # Store coords
      BlockCenters[K][L]=[X,Z]
     
      # Add upper half
      Z+=BlockTypes[BlockType]['h']
    #}
      
    # Add remained half width of block
    X += Width2
    
  # Base and CSB
  BlockCenters[0][0]=(0.,0.)
  BlockCenters[1][0]=(0.,0.)
  
  Core['InitialBlockCenters']=BlockCenters

  if '/Init/' in VERBOSE:
    print '<Initial Block Centers>'
    for K in range(M2):
      for L in range(1,len(Array[K])+1):
        print '(%d,%d)='%(K,L),BlockCenters[K][L]

#} END OF GetInitBlockCenters
#-------------------------------------------------------------------------

#} END OF INITIALIZE
#=========================================================================




#=========================================================================
#{        STATUS
#=========================================================================



#-------------------------------------------------------------------------
# Get Current Coords of (K,L) Block Center and Corners
# Origin at LD corner of first lowest left restratin
# Return: [Center,LU,LD,RU,RD]

def GetBlockCoords (Core,K,L,U,W,R,ShapeScale=1.,VERBOSE=''):
  
  # Var Space
  BlockTypes=Core['BlockTypes']
  Array=Core['Array']
  InitialCenters = Core['InitialBlockCenters']

# TEST #
#  if Core['TEST']!=Core['InitialBlockCenters']:
#    print "GBC1 Core['InitialBlockCenters']"
#    for K in Core['InitialBlockCenters'].keys():
#      Col = Core['InitialBlockCenters'][K]
#      for L in Col.keys():
#        print '(%d,%d)'%(K,L),'=',Col[L]
#    raw_input('pause')


  # Block type & dimension
  if K==0 and L==0:
    BlockType = 'Base'
    h = 0
    b = 0
  elif K==1 and L==0:
    BlockType = 'CSB'
    BTS = BlockTypes[BlockType]
    h = BTS['h']
    b = BTS['b']
  else:
    BlockType = Array[K][L]
    BTS = BlockTypes[BlockType]
    h = BTS['h']
    b = BTS['b']

# TEST #
#  if Core['TEST']!=Core['InitialBlockCenters']:
#    print "GBC2 Core['InitialBlockCenters']"
#    for K in Core['InitialBlockCenters'].keys():
#      Col = Core['InitialBlockCenters'][K]
#      for L in Col.keys():
#        print '(%d,%d)'%(K,L),'=',Col[L]
#    raw_input('pause')


  
  # Center
  Center = [InitialCenters[K][L][0],InitialCenters[K][L][1]]
  Center[0] += U*ShapeScale
  Center[1] += W*ShapeScale

# TEST #
#  if Core['TEST']!=Core['InitialBlockCenters']:
#    print "GBC3 Core['InitialBlockCenters']"
#    for K in Core['InitialBlockCenters'].keys():
#      Col = Core['InitialBlockCenters'][K]
#      for L in Col.keys():
#        print '(%d,%d)'%(K,L),'=',Col[L]
#    raw_input('pause')


  
  # LU
  LU = [Center[0]-b,Center[1]+h]
  LU[0] +=  h*   sin(R) +b*(1-cos(R))
  LU[1] += -h*(1-cos(R))+b*   sin(R)  

# TEST #
#  if Core['TEST']!=Core['InitialBlockCenters']:
#    print "GBC4 Core['InitialBlockCenters']"
#    for K in Core['InitialBlockCenters'].keys():
#      Col = Core['InitialBlockCenters'][K]
#      for L in Col.keys():
#        print '(%d,%d)'%(K,L),'=',Col[L]
#    raw_input('pause')


  
  # LD
  LD = [Center[0]-b,Center[1]-h]
  LD[0] += -h*   sin(R) +b*(1-cos(R))
  LD[1] +=  h*(1-cos(R))+b*   sin(R)
  

# TEST #
#  if Core['TEST']!=Core['InitialBlockCenters']:
#    print "GBC5 Core['InitialBlockCenters']"
#    for K in Core['InitialBlockCenters'].keys():
#      Col = Core['InitialBlockCenters'][K]
#      for L in Col.keys():
#        print '(%d,%d)'%(K,L),'=',Col[L]
#    raw_input('pause')


  # RU
  RU = [Center[0]+b,Center[1]+h]
  RU[0] +=  h*   sin(R) -b*(1-cos(R))
  RU[1] += -h*(1-cos(R))-b*   sin(R)
  

# TEST #
#  if Core['TEST']!=Core['InitialBlockCenters']:
#    print "GBC6 Core['InitialBlockCenters']"
#    for K in Core['InitialBlockCenters'].keys():
#      Col = Core['InitialBlockCenters'][K]
#      for L in Col.keys():
#        print '(%d,%d)'%(K,L),'=',Col[L]
#    raw_input('pause')


  # RD
  RD = [Center[0]+b,Center[1]-h]
  RD[0] += -h*   sin(R) -b*(1-cos(R))
  RD[1] +=  h*(1-cos(R))-b*   sin(R)
    

# TEST #
#  if Core['TEST']!=Core['InitialBlockCenters']:
#    print "GBC7 Core['InitialBlockCenters']"
#    for K in Core['InitialBlockCenters'].keys():
#      Col = Core['InitialBlockCenters'][K]
#      for L in Col.keys():
#        print '(%d,%d)'%(K,L),'=',Col[L]
#    raw_input('pause')


  if '/Coords/' in VERBOSE:
    print '<(%d,%d) Block Coords>'%(K,L)
    print 'U,W,R=',U,W,R
    print 'Center,LU,LD,RU,RD'
    print Center[0]
    print Center[1]
    print LU[0]
    print LU[1]
    print LD[0]
    print LD[1]
    print RU[0]
    print RU[1]
    print RD[0]
    print RD[1]


# TEST #
#  if Core['TEST']!=Core['InitialBlockCenters']:
#    print "GBC8 Core['InitialBlockCenters']"
#    for K in Core['InitialBlockCenters'].keys():
#      Col = Core['InitialBlockCenters'][K]
#      for L in Col.keys():
#        print '(%d,%d)'%(K,L),'=',Col[L]
#    raw_input('pause')


  return [Center,LU,LD,RU,RD]
# END OF Get Current Coords of (K,L) Block Center and Corners
#-------------------------------------------------------------------------






#-------------------------------------------------------------------------
# Print Core Array in a Vertical Arrangement

def PrintCoreArray (Core):
  print '<Print Array>'  
  
  M2=len(Core['Array'])
  
  # Find max column
  Nmax=0
  for K in range(M2):
    N=len(Core['Array'][K])
    if N>Nmax: Nmax=N
  
  # Fill out blanks
  Array=['' for i in range(M2)]
  for K in range(M2):
    for L in range(1,Nmax+1):
      N=len(Core['Array'][K])
      if L>N:
        Array[K] += ' '
      else:
        Array[K] += Core['Array'][K][L]
        
  # CSB
  StringCSB = 'CSB'
  if 2*M2-1>3: StringCSB += '='*(2*M2-1-3)
  
  # Print
  for L in range(Nmax-1,-1,-1):
    for K in range(M2):
      print Array[K][L],
    print
  if Core['Flag_CSB']==True:
    print StringCSB

  print
# END OF Print Core Array in a Vertical Arrangement
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# Make Text Core Array in a Vertical Arrangement

def MakeTxtCoreArray (Core):
  
  M2=len(Core['Array'])
  
  # Find max column
  Nmax=0
  for K in range(M2):
    N=len(Core['Array'][K])
    if N>Nmax: Nmax=N
  
  # Fill out blanks
  Array=['' for i in range(M2)]
  for K in range(M2):
    for L in range(1,Nmax+1):
      N=len(Core['Array'][K])
      if L>N:
        Array[K] += ' '
      else:
        Array[K] += Core['Array'][K][L]
        
  # CSB
  StringCSB = 'CSB'
  if 2*M2-1>3: StringCSB += '='*(M2-3)

  # Make string
  lines=''
  for L in range(Nmax-1,-1,-1):
    for K in range(M2):
      lines += Array[K][L]
    lines += '\n'
  if Core['Flag_CSB']==True:
    lines += StringCSB+'\n'

  Core['TxtArray']=lines
# END OF Make Text Core Array in a Vertical Arrangement
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Print index

def PrintIndex (Core):
  
  len_index=len(Core['ReverseIndex'])
  for i in range(len_index):
    [K,L]=Core['ReverseIndex'][i]
    if [K,L]==[0,0]: BTN='Base'
    elif [K,L]==[1,0]: BTN='CSB'
    else: BTN=Core['Array'][K][L]
    
    Verbose('index=%d, [%d,%d]=%s'%(i,K,L,BTN))

# END OF Make Text Core Array in a Vertical Arrangement
#-------------------------------------------------------------------------





#-------------------------------------------------------------------------
# GET SOLUTION OF (K,L) BLOCK FROM STATE VECTOR
# Input: Core Var Space, Solution Vector, K, L
# Output: [U,DU,W,DW,R,DR]

def GetSolFromVect (Core,Vect,K,L):
#  index = Core['Index'].index([K,L])
  index = Core['Index'][K][L]
  retval=tuple([ Vect[6*index+i] for i in range(6) ])
  return retval

def GetSolFromVect_OLD (Core,Vect,K,L):
#  index = Core['Index'].index([K,L])
  index = Core['Index'][K][L]
  retval=[]
  for i in range(6):
    retval.append(Vect[6*index+i])
  return retval

#-------------------------------------------------------------------------






#-------------------------------------------------------------------------
# PUT VALUES FOR (K,L) BLOCK TO STATE VECTOR
# Input: Core Var Space, Solution Vector, K, L, [U,DU,W,DW,R,DR]
# Output: None

def PutValToVect (Core,Vect,K,L,Values):
#  index = Core['Index'].index([K,L])
  index = Core['Index'][K][L]
  for i in range(6):
    Vect[6*index+i] = Values[i]

#-------------------------------------------------------------------------


#} END OF STATUS
#=========================================================================



#=========================================================================
#{      TEST
#=========================================================================




if __name__=='__main__':
  
  # -------------------------------------------------------------------
  # IMPORT MODULE - INPUT
  MODULE_NAME_THIS = VERSION_NAME+'_Input'
  exec('from '+MODULE_NAME_THIS+' import *')
  # -------------------------------------------------------------------

  Core={}
  Core['BlockTypes']=BlockTypes
  Core['MatProp']=MatProp
  Core['BodyForce']=BodyForce
  Core['Load']=Load
  
  MakeIndex_RearrangeArray(Core,CoreArray,VERBOSE='/Init/')
  MakeConnectivity(Core,VERBOSE='/Init/')
  GetInitBlockCenters(Core,VERBOSE='/Init/')
  PrintCoreArray(Core)
  U,W,R = 0,0,0
  K,L = 1,1
  GetBlockCoords(Core,K,L,U,W,R,VERBOSE='/Coords/')
  U,W,R = 1,2,15*pi/180
  K,L = 1,1
  GetBlockCoords(Core,K,L,U,W,R,VERBOSE='/Coords/')
  print 'TEST COMPLETED'
  
#} END OF TEST
#=========================================================================
