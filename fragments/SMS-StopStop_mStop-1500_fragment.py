import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import *
from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *

import math

baseSLHATable="""
BLOCK MASS  # Mass Spectrum
# PDG code           mass       particle
   1000001     1.00000000E+05   # ~d_L
   2000001     1.00000000E+05   # ~d_R
   1000002     1.00000000E+05   # ~u_L
   2000002     1.00000000E+05   # ~u_R
   1000003     1.00000000E+05   # ~s_L
   2000003     1.00000000E+05   # ~s_R
   1000004     1.00000000E+05   # ~c_L
   2000004     1.00000000E+05   # ~c_R
   1000005     1.00000000E+05   # ~b_1
   2000005     1.00000000E+05   # ~b_2
   1000006     %MSTOP%          # ~t_1
   2000006     1.00000000E+05   # ~t_2
   1000011     1.00000000E+05   # ~e_L
   2000011     1.00000000E+05   # ~e_R
   1000012     1.00000000E+05   # ~nu_eL
   1000013     1.00000000E+05   # ~mu_L
   2000013     1.00000000E+05   # ~mu_R
   1000014     1.00000000E+05   # ~nu_muL
   1000015     1.00000000E+05   # ~tau_1
   2000015     1.00000000E+05   # ~tau_2
   1000016     1.00000000E+05   # ~nu_tauL
   1000021     1.00000000E+05    # ~g
   1000022     %MLSP%           # ~chi_10
   1000023     1.00000000E+05   # ~chi_20
   1000025     1.00000000E+05   # ~chi_30
   1000035     1.00000000E+05   # ~chi_40
   1000024     1.00000000E+05   # ~chi_1+
   1000037     1.00000000E+05   # ~chi_2+

# DECAY TABLE
#         PDG            Width
DECAY   1000001     0.00000000E+00   # sdown_L decays
DECAY   2000001     0.00000000E+00   # sdown_R decays
DECAY   1000002     0.00000000E+00   # sup_L decays
DECAY   2000002     0.00000000E+00   # sup_R decays
DECAY   1000003     0.00000000E+00   # sstrange_L decays
DECAY   2000003     0.00000000E+00   # sstrange_R decays
DECAY   1000004     0.00000000E+00   # scharm_L decays
DECAY   2000004     0.00000000E+00   # scharm_R decays
DECAY   1000005     0.00000000E+00   # sbottom1 decays
DECAY   2000005     0.00000000E+00   # sbottom2 decays
DECAY   1000006     1.00000000E+00   # stop1 decays
    0.00000000E+00    3    1000022      5     24  # dummy allowed decay, in order to turn on off-shell decays
    1.00000000E+00    2    1000022      6
DECAY   2000006     0.00000000E+00   # stop2 decays

DECAY   1000011     0.00000000E+00   # selectron_L decays
DECAY   2000011     0.00000000E+00   # selectron_R decays
DECAY   1000012     0.00000000E+00   # snu_elL decays
DECAY   1000013     0.00000000E+00   # smuon_L decays
DECAY   2000013     0.00000000E+00   # smuon_R decays
DECAY   1000014     0.00000000E+00   # snu_muL decays
DECAY   1000015     0.00000000E+00   # stau_1 decays
DECAY   2000015     0.00000000E+00   # stau_2 decays
DECAY   1000016     0.00000000E+00   # snu_tauL decays
DECAY   1000021     0.00000000E+00   # gluino decays
DECAY   1000022     0.00000000E+00   # neutralino1 decays
DECAY   1000023     0.00000000E+00   # neutralino2 decays
DECAY   1000024     0.00000000E+00   # chargino1+ decays
DECAY   1000025     0.00000000E+00   # neutralino3 decays
DECAY   1000035     0.00000000E+00   # neutralino4 decays
DECAY   1000037     0.00000000E+00   # chargino2+ decays
"""

generator = cms.EDFilter("Pythia8GeneratorFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13600.),
    RandomizedParameters = cms.VPSet(),
)

model = "T2tt"
batch = 1

# weighted average of matching efficiencies for the full scan
# must equal the number entered in McM generator params
mcm_eff = 0.248

# Number of events: min(goalLumi*xsec, maxEvents) (always in thousands)
goalLumi = 400
minLumi = 1e-40 # Skip minimum lumi
minEvents, maxEvents = 20, 1000
diagStep, bandStep = 25, 50
midDM, maxDM = 87, 100
addDiag = [183, 167] # DeltaM for additional diagonal lines to be added

class gridBlock:
  def __init__(self, xmin, xmax, xstep, ystep):
    self.xmin = xmin
    self.xmax = xmax
    self.xstep = xstep
    self.ystep = ystep

scanBlocks = []
if (batch==1):
  scanBlocks.append(gridBlock(1500,  1501, 100, 100))
minDM = 87
ymin, ymed, ymax = 0, 0, 2500

def matchParams(mass):
  if mass < 649: return 62., 0.274
  elif mass < 699: return 64., 0.269
  elif mass < 749: return 64., 0.269
  elif mass < 799: return 66., 0.259
  elif mass < 849: return 66., 0.261
  elif mass < 899: return 68., 0.257
  elif mass < 949: return 68., 0.252
  elif mass < 999: return 70., 0.250
  elif mass < 1049: return 70., 0.248
  elif mass < 1099: return 70., 0.248
  elif mass < 1149: return 70., 0.249
  elif mass < 1199: return 70., 0.242
  elif mass < 1249: return 70., 0.239
  elif mass < 1299: return 70., 0.242
  elif mass < 1349: return 70., 0.241
  elif mass < 1399: return 70., 0.237
  elif mass < 1449: return 70., 0.239
  elif mass < 1499: return 70., 0.241
  elif mass < 1549: return 70., 0.235
  elif mass < 1599: return 70., 0.239
  elif mass < 1649: return 70., 0.239
  elif mass < 1699: return 70., 0.237
  elif mass < 1749: return 70., 0.241
  elif mass < 1799: return 70., 0.237
  elif mass < 1849: return 70., 0.237
  elif mass < 1899: return 70., 0.240
  elif mass < 1949: return 70., 0.241
  elif mass < 1999: return 70., 0.244
  elif mass < 2049: return 70., 0.246
  elif mass < 2099: return 70., 0.249
  elif mass < 2149: return 70., 0.246
  elif mass < 2199: return 70., 0.246
  elif mass < 2249: return 70., 0.251
  elif mass < 2299: return 70., 0.249
  elif mass < 2349: return 70., 0.257
  elif mass < 2399: return 70., 0.257
  elif mass < 2449: return 70., 0.261
  elif mass < 2499: return 70., 0.264
  elif mass < 2549: return 70., 0.266
  ### Just for testing
  else: return 70., 0.243

def xsec(mass):
  if mass < 300: return 319925471928717.38*math.pow(mass, -4.10396285974583*math.exp(mass*0.0001317804474363))
  else: return 4855957031250000*math.pow(mass, -4.71716128867804*math.exp(mass*6.175277146619076e-05))

# Number of events for mass point, in thousands
def events(mass):
  xs = xsec(mass)
  nev = min(goalLumi*xs, maxEvents*1000)
  if nev < xs*minLumi: nev = xs*minLumi
  nev = max(nev/1000, minEvents)
  return math.ceil(nev) # Rounds up

cols = []
Nevents = []
xmin, xmax = 9999, 0
for block in scanBlocks:
  Nbulk, Ndiag = 0, 0
  for mx in range(block.xmin, block.xmax, min(bandStep, diagStep)):
    if (batch==2 and mx==block.xmin): continue
    xmin = min(xmin, block.xmin)
    xmax = max(xmax, block.xmax)
    col = []
    my = 0
    begBand = min(max(ymed, mx-maxDM), mx-minDM)
    begDiag = min(max(ymed, mx-midDM), mx-minDM)
    # Adding bulk points
    if (mx-block.xmin)%block.xstep == 0 :
      for my in range(ymin, begBand, block.ystep):
        if my > ymax: continue
        # Adding extra diagonals to the bulk
        nev = events(mx)
        col.append([mx,my, nev])
        Nbulk += nev
    
    cols.append(col)

  Nevents.append([Nbulk, Ndiag])

mpoints = []
for col in cols: mpoints.extend(col)

for point in mpoints:
    mstop, mlsp = point[0], point[1]
    qcut, tru_eff = matchParams(mstop)
    wgt = point[2]*(mcm_eff/tru_eff)
    print(mstop, mlsp, wgt)
    
    if mlsp==0: mlsp = 1
    slhatable = baseSLHATable.replace('%MSTOP%','%e' % mstop)
    slhatable = slhatable.replace('%MLSP%','%e' % mlsp)
################################################

    basePythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        pythia8PSweightsSettingsBlock,
        processParameters = cms.vstring(
            'JetMatching:setMad = off',
            'JetMatching:scheme = 1',
            'JetMatching:merge = on',
            'JetMatching:jetAlgorithm = 2',
            'JetMatching:etaJetMax = 5.',
            'JetMatching:coneRadius = 1.',
            'JetMatching:slowJetPower = 1',
            'JetMatching:qCut = %.0f' % qcut, #this is the actual merging scale
            'JetMatching:nQmatch = 5', #4 corresponds to 4-flavour scheme (no matching of b-quarks), 5 for 5-flavour scheme
            'JetMatching:nJetMax = 2', #number of partons in born matrix element for highest multiplicity
            'JetMatching:doShowerKt = off', #off for MLM matching, turn on for shower-kT matching
            '6:m0 = 172.5',
            'Check:abortIfVeto = on',
        ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'pythia8PSweightsSettings',
                                    'processParameters'
                                    )
    )
    generator.RandomizedParameters.append(
        cms.PSet(
            ConfigWeight = cms.double(wgt),
            GridpackPath =  cms.string('/afs/cern.ch/user/t/taiwoo/public/Run3_NPS_StopStop/gridpacks/SMS-StopStop_mStop-%i_slc7_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz' % mstop),
            ConfigDescription = cms.string('%s_%i_%i' % (model, mstop, mlsp)),
            SLHATableForPythia8 = cms.string('%s' % slhatable),
            PythiaParameters = basePythiaParameters,
        ),
    )
    ProductionFilterSequence = cms.Sequence(generator)