def lineConfig():

#	Take the definition of the line and clean it up
	length = len(lineDef)
	for pnt, val in enumerate(lineDef):
		if pnt == 2:
			if val == "=": 
				continue
			else:
				print ("Unexpected structure of the line definition")
				exit()
		if (pnt == 3):
			val = val[1:-1]
			lineDesc.append(val)			
		if pnt > 3 and pnt < length-1:
			val = val[0:-1]
			lineDesc.append(val)
		if pnt == length-1:
			val = val[0:-2]
			lineDesc.append(val)

def beamLineCreate():

	outLine = "Stage,Section,Element,Type,Parameter,Value,Unit,Comment\n"
	outFile.write(outLine)

	for pnt, element in enumerate(lineDesc):
		if element in driftDict:
			driftSpaceWrite(element)
		elif element in sbendDict:
			sbendWrite(element)
		elif element in quadDict:
			quadWrite(element)
		else:
			print ("no matching entry")
			exit()


def driftSpaceWrite(element):

	props = driftDict[element]
	elType = props[0]
	elLength = props[1]
	eqPos = elLength.find("=")
	starPos = elLength.find("*")
	length = float(elLength[eqPos+1:starPos])

	units = elLength[starPos+1:-1]
	if units == "mm":
		length = length/1000.0
	elif units =="cm":
		length = length/100.0
	elif units =="m":
		length = length
	else:
		print ("unrecongised length units ", units)
		exit()

	outString = "1,ProductionStraight,Drift,,Length," + str(length) + ',m,"' + element + '"\n'
	outFile.write(outString)

def sbendWrite(element):

#s_rm_dp_2: sbend, l=250*cm, angle=0.09496865082224923;
#1,Arc,Dipole,"Sector (Length, angle)",Length,0.8,m,
#1,Arc,Dipole,"Sector (Length, angle)",Angle,45,degrees,

	props = sbendDict[element]
	eqPos = props[2].find("=")
	elAngle = props[2][eqPos+1:-1]

	eqPos = props[1].find("=")
	starPos = props[1].find("*")
	elLength = props[1][eqPos+1:starPos]
	units = props[1][starPos+1:-1]
	if (units == "cm"):
		length = float(elLength)/100.0
	else:
		print ("unrecognised units in sbend")
		exit()

	outString0 = '1,ProductionStraight,Dipole,"Sector (Length, angle)",Length,' + str(length) + ',m,"' + element + '"\n'
	outString1 = '1,ProductionStraight,Dipole,"Sector (Length, angle)",Angle,' + elAngle + ',degrees,\n'
	outFile.write(outString0)	
	outFile.write(outString1)

def quadWrite(element):

#1,ProductionStraight,Dquad,,Length,1.0,m,"bdsim q3"
#1,ProductionStraight,Dquad,,kq,0.359935,m-2,
#q1: quadrupole, l=1*m, k1=-0.253169, magnetGeometryType="polesfacet", aper1=300*mm;

	props = quadDict[element]
	field = float(props[2][3:-1])
	elLength = props[1]
	eqPos = elLength.find("=")
	starPos = elLength.find("*")

	length = elLength[eqPos+1:starPos]

	if (field < 0.0):
		field = -field
		outString0 = "1,ProductionStraight,Dquad,,Length," + length + ',m,"' + element + '"\n' 
		outString1 = "1,ProductionStraight,Dquad,,kq," + str(field) + ',m-2,\n'
		outFile.write(outString0)
		outFile.write(outString1)
	else:
		outString0 = "1,ProductionStraight,Fquad,,Length," + length + ',m,"' + element + '"\n' 
		outString1 = "1,ProductionStraight,Fquad,,kq," + str(field) + ',m-2,\n'
		outFile.write(outString0)
		outFile.write(outString1)

#============================================
 
# read in the whole bdsim line definition file
inFileName = 'nuSTORM_PionBeamline_ProductionStraight_5GeV.gmad'
print ("reading bdSim file ", inFileName)

file1 = open(inFileName, 'r')
Lines = file1.readlines()
file1.close()

#	Split the lines into quadrupoles, sbends, drifts and the definition of 
#	the beam line.
quadDict={}
driftDict={}
lineDef=[]
lineDesc=[]
sbendDict={}
# Strips the newline character
for line in Lines:
	quad=[]
	drift=[]
	sbend=[]
# Strips the newline character
	bits = line.strip().split()
	if (len(bits) > 1):
		if (bits[0][0:1] == "!"):
			continue
		if (bits[1] == "drift,"):
			driftDict[bits[0][:-1]]=bits[1:]
		elif (bits[1] == "quadrupole,"):
			quadDict[bits[0][:-1]]=bits[1:]
		elif (bits[1] == "sbend,"):
			sbendDict[bits[0][:-1]]=bits[1:]
		elif (bits[1] == "line"):
			lineDef = bits

lineConfig()
outFileName = "optics.csv"
outFile = open(outFileName, 'w')
beamLineCreate()
outFile.close()

print ("Created Linear Optics file ", outFileName)
print ("Warning ***** apertures not included")

exit()
