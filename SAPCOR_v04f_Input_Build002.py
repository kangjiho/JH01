VERSION_NAME = 'SAPCOR_v04f'





###############################################################################
# USER OUTPUT FILENAME (WITHOUT EXTENSION)
FN_OUT = 'Verif_05_Single'
FN_FIG = 'Fig'
FN_CoreShape = 'CoreShape'
FN_Solution = 'Solution'
###############################################################################



###############################################################################
# USER CONTROL

# Eliminate sticking force
FLAG_NoStickForce = False

# Test Forces Only
# Will not run MAIN LOOP
# Will calculate acceleration using intial conditions
FLAG_TestForce = False

# END OF USER CONTROL
###############################################################################




###############################################################################
# NUMERICAL INTEGRATION SETTING
# Error Tolerance
Tol0 = 1e-9
TolMax = 1e-1

# ODEINT Initial Tolerance Criteria

# Lower Limit of Section Time
# Initial Section Time
SectionTime0 = 0.001
# If SectionTime<SectionTimeMin: Tolerance Control Starts
SectionTimeMin = 1e-4


# Number of Sampling Points in a Integration Section
#NPoints = 32000
#NPoints = 32000
# END OF NUMERICAL INTEGRATION SETTING
###############################################################################



###############################################################################
# VERBOSE CONTROL
# VERBOSER FLAG
VERBOSE = ''
VERBOSE += '/Init/'
#VERBOSE += '/Force/'
#VERBOSE += '/ForceDF/'
#VERBOSE += '/ForceH/'
#VERBOSE += '/ForceV/'
#VERBOSE += '/ForceF/'
#VERBOSE += '/Coords/'
#VERBOSE += '/RK/'
#VERBOSE += '/Error/'
#VERBOSE += '/Adaptive/'
#VERBOSE += '/Output/'
# END OF VERBOSE CONTROL
###############################################################################













###############################################################################
# USER INPUT
M = 0.00628
C = 1.77
K = 2.5e4
a = 2.5
b = 7.23
b2= 1.
h = 14.23
h2 = 2.
d = 5.2
I = 0.524
Kh = K
Ch = C
Kv = K
Cv = C
Kd = K
Cd = C
mu_s = 0.2
mu_k = 0.2
d_mu = 100
xi_F_cr = 0.01

# Gap
Delta = 0.1     # General Block Gap 
DeltaTop = 0.05 # Top Block Gap
DeltaD = 0.025    # Dowel Gap
#DeltaL = 0.05   # Dowel Gap Left
#DeltaR = 0.05   # Dowel Gap Right

# Total time for analysis
# Real analysis time >= TotalTime
TotalTime = 0.0005
TotalTime = 1.001
#TotalTime = 2.001


# Toggle Forces On/Off for Test
ApplyForces={}
ApplyForces['Block_H']=True
ApplyForces['Block_D']=True
ApplyForces['Block_DF']=True
ApplyForces['Block_V']=True
ApplyForces['Block_VF']=True


# Block Types
#
# 'Fixed':'None' if the block is not fixed
#         'Fixed' if want to maintain initial U,DU,W,DW,R,DR
#         'FixedToBase' if U,W,R=U,W,R of Base
#         'FixedU' if U=Fixed to Initial,DU=0
#         'FixedUR' if  U=Fixed to Initial, DU=0
#                   and R=Fixed to Initial, DR=0
#
BlockTypes={}
BlockTypes['W']={'a':a,'b':b,'h':1e-10,  'd':d,'M':M,'I':I,'Kh':Kh,'Ch':Ch,'Kv':Kv,'Cv':Cv,'Kd':Kd,'Cd':Cd,'mu_s':mu_s,'mu_k':mu_k,'d_mu':d_mu,'xi_F_cr':xi_F_cr,'Fixed':'None'}
BlockTypes['X']={'a':a,'b':b,'h':h*2,'d':d,'M':M,'I':I,'Kh':Kh,'Ch':Ch,'Kv':Kv,'Cv':Cv,'Kd':Kd,'Cd':Cd,'mu_s':mu_s,'mu_k':mu_k,'d_mu':d_mu,'xi_F_cr':xi_F_cr,'Fixed':'None'}
BlockTypes['Y']={'a':a,'b':b,'h':h*3,'d':d,'M':M,'I':I,'Kh':Kh,'Ch':Ch,'Kv':Kv,'Cv':Cv,'Kd':Kd,'Cd':Cd,'mu_s':mu_s,'mu_k':mu_k,'d_mu':d_mu,'xi_F_cr':xi_F_cr,'Fixed':'None'}
BlockTypes['Z']={'a':a,'b':b,'h':h*4,'d':d,'M':M,'I':I,'Kh':Kh,'Ch':Ch,'Kv':Kv,'Cv':Cv,'Kd':Kd,'Cd':Cd,'mu_s':mu_s,'mu_k':mu_k,'d_mu':d_mu,'xi_F_cr':xi_F_cr,'Fixed':'None'}
BlockTypes['B']={'a':a,'b':b,'h':h,  'd':d,'M':M,'I':I,'Kh':Kh,'Ch':Ch,'Kv':Kv,'Cv':Cv,'Kd':Kd,'Cd':Cd,'mu_s':mu_s,'mu_k':mu_k,'d_mu':d_mu,'xi_F_cr':xi_F_cr,'Fixed':'None'}
BlockTypes['C']={'a':a,'b':b,'h':h,  'd':d,'M':M,'I':I,'Kh':Kh,'Ch':Ch,'Kv':Kv,'Cv':Cv,'Kd':Kd,'Cd':Cd,'mu_s':mu_s,'mu_k':mu_k,'d_mu':d_mu,'xi_F_cr':xi_F_cr,'Fixed':'FixedUR'}
BlockTypes['D']={'a':a,'b':b,'h':h*3,'d':d,'M':M,'I':I,'Kh':Kh,'Ch':Ch,'Kv':Kv,'Cv':Cv,'Kd':Kd,'Cd':Cd,'mu_s':mu_s,'mu_k':mu_k,'d_mu':d_mu,'xi_F_cr':xi_F_cr,'Fixed':'None'}
BlockTypes['E']={'a':a,'b':b,'h':h*4,'d':d,'M':M,'I':I,'Kh':Kh,'Ch':Ch,'Kv':Kv,'Cv':Cv,'Kd':Kd,'Cd':Cd,'mu_s':mu_s,'mu_k':mu_k,'d_mu':d_mu,'xi_F_cr':xi_F_cr,'Fixed':'None'}
BlockTypes['CSB']={'a':0,'b':0,'h':h/2,'d':d,'M':M,'I':I,'Kh':Kh,'Ch':Ch,'Kv':Kv,'Cv':Cv,'Kd':Kd,'Cd':Cd,'mu_s':mu_s,'mu_k':mu_k,'d_mu':d_mu,'xi_F_cr':xi_F_cr,'Fixed':'FixedToBase'}

# Material Properties
MatProp={}
MatProp['Block']={}
MatProp['Block']['K']=K
MatProp['Block']['C']=C
MatProp['Restraint']={}
MatProp['Restraint']['K']=K
MatProp['Restraint']['C']=C

# Core Array
CoreArray = []
CoreArray.append('CSB')
CoreArray.append('W')
#CoreArray.append('BBBBBBBBBBBBB')
CoreArray.append('BB')
CoreArray.append('W')

# Initial Conditions
# InitialConditions.append([K,L, U,DU, W,DW, R,DR])
InitialConditions=[]
#InitialConditions.append([1,1, 0.,0., 0.,0., 0.,0.])
#
#                          K  L  U       DU  W       DW  R       DR 
InitialConditions.append([ 1, 1, 1.2677, 0., 0.5760, 0., 0.0873, 0. ])
InitialConditions.append([ 1, 2, 1.2677, 0., 1.7280, 0., -0.0873, 0. ])

# Body Force (acceleration)
BodyForce = {}
# Unit : cm/s2
BodyForce[1] = 0.
BodyForce[3] = -981.0
BodyForce[5] = 0.

# Loads on support frame (F_SupportFrame()) (acceleration)
Load = {}; Load[1]={}; Load[3]={}; Load[5]={}
# Load[i] : i is axis direction
#   i=1 : Coord=X, Disp=U
#   i=3 : Coord=Z, Disp=W
#   i=5 : Coord=Y, Rot=R
#   i=2,4,6 are not used in 2D
# Load[i]['Type'] = 'None', 'Accel', 'Disp', or 'Data'
# if Load[i]['Type'] != 'None
#   Load[i]['FnType'] = 'Sin', 'Cos', or 'Data
#   if Load[i]['FnType'] == 'Sin' or 'Cos'
#     Load[i]['Amp'] : Max amplitude (cm/s2 or cm)
#     Load[i]['Freq'] : Frequency (fixed during calculation) (Hz)
#   if Load[i]['FnType'] == 'Data'
#     Load[i]['FileName'] = '<<Data File Name>>'
#     ! Data File Format
#     ! 1  : Dummy Header Line (eg. t, Accel)
#     ! 2~ : t, Acceleration (cm/s2) or Displacement (cm)
Load[1]['Type']='None'
Load[1]['FnType']='Cos'
Load[1]['Amp']=500.
Load[1]['Freq']=10.
Load[3]['Type']='None'
Load[5]['Type']='None'

# Time Frequency for CoreShape Output
OP_CoreShape_TimeFreq = 0.02

# Scale Factor for Deformation
OP_CoreShape_Scale = 1.

# CoreShape Figure Boarder Margins
OP_CoreShape_MarginX = 10
OP_CoreShape_MarginY = 10

# Time Frequency for Solution Output
# What is "Solution" : Core State (All Block States) at a Specified Time Point
OP_Solution_TimeFreq = 0.02

# Block Data Sampling Frequency
# What is "Block Data" : All time response from beginning of a block
OP_Block_TimeFreq = 0.001

# Filename for Verbose
FN_Verbose = 'Verbose.txt'


# END OF USER INPUT
###############################################################################

