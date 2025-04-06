import math

model = "T2tt"
batch = 1
# batch 1 contains 33,504,000 events in 37 points
# batch 2 contains 31,675,000 events in 368 points

# weighted average of matching efficiencies for the full scan
# must equal the number entered in McM generator params
mcm_eff = 0.295
if batch==2: mcm_eff = 0.261

# Number of events: min(goalLumi*xsec, maxEvents) (always in thousands)
goalLumi = 400
minLumi = 1e-40 # Skip minimum lumi
minEvents, maxEvents = 20, 1000
diagStep, bandStep = 25, 50
midDM, maxDM = 300, 700
addDiag = [183, 167] # DeltaM for additional diagonal lines to be added

# Parameters that define the grid in the bulk and diagonal
class gridBlock:
  def __init__(self, xmin, xmax, xstep, ystep):
    self.xmin = xmin
    self.xmax = xmax
    self.xstep = xstep
    self.ystep = ystep

scanBlocks = []
if (batch==1): scanBlocks.append(gridBlock(1200,  3001, 300, 100))
elif (batch==2): scanBlocks.append(gridBlock(400,  1201, 50, 50))
minDM = 87
ymin, ymed, ymax = 0, 0, 500

def matchParams(mass):
  if mass>99 and mass<199: return 62., 0.498
  elif mass<299: return 62., 0.361
  elif mass<399: return 62., 0.302
  elif mass<499: return 64., 0.275
  elif mass<599: return 64., 0.254
  elif mass<1299: return 68., 0.237
  elif mass<1801: return 70., 0.243
  ### Just for testing
  else: return 70., 0.243

def xsec(mass):
  if mass < 300: return 319925471928717.38*math.pow(mass, -4.10396285974583*math.exp(mass*0.0001317804474363))
  else: return 6953884830281245*math.pow(mass, -4.7171617288678069*math.exp(mass*6.1752771466190749e-05))

# Number of events for mass point, in thousands
def events(mass):
  xs = xsec(mass)
  nev = min(goalLumi*xs, maxEvents*1000)
  if nev < xs*minLumi: nev = xs*minLumi
  nev = max(nev/1000, minEvents)
  return math.ceil(nev) # Rounds up

# -------------------------------
#    Constructing grid

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
        for dm in addDiag:
          if(len(cols)==0 and batch==1): continue # Don't add point before the beginning
          dm_before = mx-block.xstep -my
          dm_after = mx - my
          if(dm>dm_before and dm<dm_after):
            nev = events(my+dm)
            col.append([my+dm,my, nev])
            Nbulk += nev
        nev = events(mx)
        col.append([mx,my, nev])
        Nbulk += nev
    # Adding diagonal points in inside band
    if (mx-block.xmin)%bandStep == 0 :
      for my in range(begBand, mx-midDM, bandStep):
        if my > ymax: continue
        # Adding extra diagonals to the band
        for dm in addDiag:
          if(len(cols)==0 and batch==1): continue # Don't add point before the beginning
          dm_before = mx-bandStep -my
          dm_after = mx - my
          if(dm>dm_before and dm<dm_after):
            nev = events(my+dm)
            col.append([my+dm,my, nev])
            Ndiag += nev
        # Adding standard diagonal points
        nev = events(mx)
        col.append([mx,my, nev])
        Ndiag += nev
    # Adding diagonal points in band closest to outer diagonal
    for my in range(begDiag, mx-minDM+1, diagStep):
      if my > ymax: continue
      # Adding extra diagonals to the band
      for dm in addDiag:
        if(len(cols)==0 and batch==1): continue # Don't add point before the beginning
        dm_before = mx-diagStep -my
        dm_after = mx - my
        if(dm>dm_before and dm<dm_after):
          nev = events(my+dm)
          col.append([my+dm,my, nev])
          Ndiag += nev
      nev = events(mx)
      col.append([mx,my, nev])
      Ndiag += nev
    if(my !=  mx-minDM and mx-minDM <= ymax):
      my = mx-minDM
      nev = events(mx)
      col.append([mx,my, nev])
      Ndiag += nev
    cols.append(col)
  Nevents.append([Nbulk, Ndiag])

mpoints = []
for col in cols: mpoints.extend(col)

for point in mpoints:
    mstop, mlsp = point[0], point[1]
    qcut, tru_eff = matchParams(mstop)
    wgt = point[2]*(mcm_eff/tru_eff)
    print(mstop, mlsp, wgt)