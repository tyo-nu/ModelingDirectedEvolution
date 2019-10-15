
# This program will take stimulated and unstimulated .fcs files and produce the following:
#  1. an overlay plot of the stimulated and unstimulated files
#  2. a breakdown of the library into True Positives, Low True Positives, High True Positives,
#			Constitutive Low and Constitutive High, with proportions and bin frequencies in a
#			.txt file named Populations.txt
#  3. an overlay plot of the True Positive, Constitutive Low, Constitutive High and 
#			Unstimulated True Positive populations
#  4. predictions of the enrichment of True Positives by sorting threshold for Stimulated cells
#			below the threshold and Unstimulated cells above the threshold in a .txt file named
#			Thresholds.txt
#  5. a plot of the post-sort enrichment predictions 
#
# Script should be run in console directly as 'python FCStoThresholdPrediction.py'
#
# NOTE: The file Template_FCStoThresholdPrediction.txt must be in the same folder.
#
# Template file 'Template_FCStoThresholdPrediction.txt' should be a plain text file containing:
# 1. The absolute path to directory containing the unstimulated .fcs file
# 2. The absolute path to directory containing the stimulated .fcs file
# 3.  the channel name (choose from Cy3-A, DAPI-A, PE-Texas Red-A, CFP-A, Alexa Fluor
#       488-A, etc.) The channel names are listed in the FCS file header when it is
#       opened as a text file. No quotes. 
# 4.  the log(FITC) cutoff between Low and High True Positives
# 5.  the log(FITC) cutoff between Low and High Constitutive


from __future__ import division
import numpy as np
import matplotlib.pylab as plt
from FlowCytometryTools import FCMeasurement, FCPlate, PolyGate

__author__ = 'Sarah Stainbrook'

def threshold_finder(array,thresholdvalue):
	for item in range(1,len(array)+1):
		total = reduce(lambda x, y: x+y, array[0:item])
		if total >= thresholdvalue:
			return item
	return len(array)


# Parse template file.
with open('Template_FCStoThresholdPrediction.txt','r') as f:
	NumFiles = 2
	FileLocs = [0 for x in range(NumFiles)]
	FileNames = [0 for x in range(NumFiles)]
	for file in range(NumFiles):
		FileLoc = f.readline()
		FileLocs[file] = FileLoc.strip()
	channel = f.readline()
	channel = str(channel.strip())
	TPLowVsHighCutoff=f.readline()
	TPLowVsHighCutoff = float(TPLowVsHighCutoff.strip())
	ConstitLowVsHighCutoff = f.readline()
	ConstitLowVsHighCutoff = float(ConstitLowVsHighCutoff.strip())

#Retrieve data from file and sort into bins
Frequencies = []
X = np.linspace(0, 13.5, 1350)
counter=0
for file in range(NumFiles):
	data = FCMeasurement(ID = 'Data', datafile = FileLocs[file])
	FP_data = data[channel]
	log_data = [np.log(x) for x in FP_data if x>1]
	scaledlog_data = [100*x for x in log_data]
	freqs = (np.bincount(scaledlog_data))/100               # one bin per every 0.1 log fluorescence unit
	scaledfreqs = [a/(len(scaledlog_data)/100) for a in freqs]  # scales the data by the total number of points such that files of different lengths can be compared.
	Frequencies.append(scaledfreqs)
	if len(scaledfreqs) < len(X):
		for counter in range(0,(len(X)-len(scaledfreqs))):
			scaledfreqs.append(0.0)
			counter += 1

# Plot the original .fcs files
Unstimulated, = plt.plot(X,Frequencies[0],lw=2,color='dimgray')#linestyle=":",marker='.',markersize=4,
Stimulated, = plt.plot(X,Frequencies[1],lw=2,color = 'darkorange')#linestyle=':',marker='.',markersize=4
plt.ylim(-0.002,0.0121)
plt.xlim(2,10)
plt.xlabel('log GFP expression', size = 15)
plt.ylabel('Population Proportion', size = 15)
plt.legend([Unstimulated, Stimulated],['Unstimulated','Stimulated'],loc='upper left')
plt.title('Unstimulated vs Stimulated Cells', size = 20)
plt.show()



# Subtract the Stimulated from Unstimulated histograms to derive populations
SubtractedHistograms = []
UnstimTruePositives = []
TruePositives = []
LowTruePositives = []
HighTruePositives = []
ConstitutiveLow = []
ConstitutiveHigh = []

for item in range(0,len(X)):
	Subtracted = Frequencies[1][item] - Frequencies[0][item]
	SubtractedHistograms.append(Subtracted)
	if X[item]<= TPLowVsHighCutoff:
		HighTruePositives.append(0)
		if SubtractedHistograms[item] >= 0:
			LowTruePositives.append(SubtractedHistograms[item])
			UnstimTruePositives.append(0)
			TruePositives.append(SubtractedHistograms[item])
		else:
			TruePositives.append(0)
			LowTruePositives.append(0)
			UnstimTruePositives.append(SubtractedHistograms[item])
	if X[item]>TPLowVsHighCutoff:
		LowTruePositives.append(0)	
		if SubtractedHistograms[item] >= 0:
			HighTruePositives.append(SubtractedHistograms[item])
			TruePositives.append(SubtractedHistograms[item])
			UnstimTruePositives.append(0)
		else:
			HighTruePositives.append(0)
			TruePositives.append(0)
			UnstimTruePositives.append(SubtractedHistograms[item])
	if X[item]<=ConstitLowVsHighCutoff:
		ConstitutiveHigh.append(0)
		ConstitutiveLow.append(Frequencies[0][item]+UnstimTruePositives[item])	
	if X[item]>ConstitLowVsHighCutoff:
		ConstitutiveLow.append(0)
		ConstitutiveHigh.append(Frequencies[0][item]+UnstimTruePositives[item])

ProportionLowTruePositives=sum(LowTruePositives)
ProportionHighTruePositives=sum(HighTruePositives)
ProportionTruePositives = sum(LowTruePositives)+sum(HighTruePositives) #+ ProportionInducers
ProportionConstitutiveLow = sum(ConstitutiveLow)
ProportionConstitutiveHigh = sum(ConstitutiveHigh)


# Write file containing population information
f = open('Populations.txt', 'w')
f.write('\n%.3f'%(ProportionTruePositives))
f.write('\n%3f'%(ProportionLowTruePositives))
f.write('\n%3f'%(ProportionHighTruePositives))
f.write('\n%.3f'%(ProportionConstitutiveLow))
f.write('\n%.3f'%(ProportionConstitutiveHigh))
f.write('\n')
f.write(str(TruePositives))
f.write('\n')
f.write(str(LowTruePositives))
f.write('\n')
f.write(str(HighTruePositives))
f.write('\n')
f.write(str(ConstitutiveLow))
f.write('\n')
f.write(str(ConstitutiveHigh))
f.write('\n')
f.write(str(UnstimTruePositives))
f.close()

# Plot populations
TruePos, = plt.plot(X,TruePositives,color='blue')
ConstLow, = plt.plot(X,ConstitutiveLow,color='darkorchid')
ConstHigh, = plt.plot(X,ConstitutiveHigh,color='olivedrab')
UnstimTruePos,=plt.plot(X,UnstimTruePositives,color='deepskyblue')
#plt.plot((LowVsHighCutoff,LowVsHighCutoff),(-0.1,0.1),'k',linestyle=':',lw=1)
plt.ylim(-0.002,0.0121)
plt.xlim(2,10)
plt.xlabel('log GFP expression', size = 15)
plt.ylabel('Population Proportion', size = 15)
plt.legend([TruePos,ConstLow,ConstHigh,UnstimTruePos],['True Positives','Constitutive Low','Constitutive High','Unstimulated True Positives'],loc='upper left')
plt.title('Populations', size = 20)
plt.show()



# Predict enrichment for sorting stimulated population, keeping cells below the threshold
Thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]
LowOutputProportions = []
LowTruePosPostSort = []
for item in range(0,len(Thresholds)):
	cutoff = threshold_finder(Frequencies[1],Thresholds[item]/100)
	LowTruePosKept = 0
	HighTruePosKept = 0
	ConstitutiveHighKept = 0
	ConstitutiveLowKept = 0
	for item in range(0, len(X)):
		if item <= cutoff:
			LowTruePosKept += LowTruePositives[item]*(0.003)# + 0.02*HighTruePosValues[item]*Proportions[1]
			HighTruePosKept += HighTruePositives[item]*0.98
			ConstitutiveHighKept += ConstitutiveHigh[item]*0.98
			ConstitutiveLowKept += 0.04*ConstitutiveLow[item]+((HighTruePositives[item]+ConstitutiveHigh[item])*0.02)

	TruePosKept = LowTruePosKept+HighTruePosKept
	TotalKept = LowTruePosKept + HighTruePosKept + ConstitutiveHighKept + ConstitutiveLowKept
	Proportions = [TruePosKept/TotalKept, LowTruePosKept/TotalKept,HighTruePosKept/TotalKept,ConstitutiveLowKept/TotalKept,ConstitutiveHighKept/TotalKept]
	LowTruePosPostSort.append(TruePosKept/TotalKept)
	LowOutputProportions.append(Proportions)

# Plot Below predictions 
Below, = plt.plot(Thresholds,LowTruePosPostSort,color='k')

BelowDataX = [1.1333,2.666666667,10.54666667,17.03]
BelowDataY = [0.072,0.073666667,0.218333333,0.174666667]
BelowStDevs = [0.013747727,0.006506407,0.029143324,0.008504901]
BelowXStDevs = [0.115470054,0.057735027,0.427239199,0.530282943]
plt.ylim(0.0,0.4)
plt.xlim(0.0,50)

plt.errorbar(BelowDataX,BelowDataY,yerr=BelowStDevs,xerr=BelowXStDevs,color='k',lw=0, elinewidth=1,marker='o')
plt.xlabel('Percent Stimulated Cells Below Threshold', size = 15)
plt.ylabel('Post-Sort TP Proportion', size = 15)
plt.text(9.9,0.25,'**')
plt.text(16.4,0.185,'**')
plt.text(1.13,0.04,'^',size=12, stretch='ultra-condensed',weight=1000, ha='center')
plt.title('Sorting of Stimulated Library', size = 20)
plt.show()





# Predict enrichment for sorting unstimulated population, keeping cells above the threshold
HighOutputProportions = []
HighTruePosPostSort = []
for item in range(0,len(Thresholds)):
	cutoff = threshold_finder(Frequencies[0],1-(Thresholds[item]/100))
	TruePosKept = 0
	ConstitutiveHighKept = 0
	ConstitutiveLowKept = 0
	for item in range(0, len(X)):
		if item >= cutoff:
			TruePosKept += -UnstimTruePositives[item]*0.98
	 		ConstitutiveHighKept += ConstitutiveHigh[item]*0.98
	 		ConstitutiveLowKept += ConstitutiveLow[item]*0.02+(LowTruePositives[item]+HighTruePositives[item]+ConstitutiveHigh[item])*0.02

	TotalKept = TruePosKept + ConstitutiveHighKept + ConstitutiveLowKept
	Proportions = [TruePosKept/TotalKept, ConstitutiveLowKept/TotalKept,ConstitutiveHighKept/TotalKept]
	HighTruePosPostSort.append(TruePosKept/TotalKept)
	HighOutputProportions.append(Proportions)



# Plot Below predictions 
Above, = plt.plot(Thresholds,HighTruePosPostSort,color='k')

AboveDataX = [0.064733333,3.510133333,14.50666667]
AboveDataY = [0.175666667,0.145666667,0.138333333]
AboveStDevs = [0.004163332,0.021079216,0.017156146]
AboveXStDevs = [0.034999619,0.139344944,1.029190621]

plt.errorbar(AboveDataX,AboveDataY,yerr=AboveStDevs,xerr=AboveXStDevs,color='k',lw=0, elinewidth=1,marker='o')
plt.xlabel('Percent Unstimulated Cells Above Threshold', size = 15)
plt.ylabel('Post-Sort TP Proportion', size = 15)
plt.title('Sorting of Unstimulated Library', size = 20)
plt.ylim(0.0,0.4)
plt.xlim(0,50)
plt.text(3.51,0.105,'^',size=12, stretch='ultra-condensed',weight=1000, ha='center')
plt.show()


# Write file containing predictions
f = open('Thresholds.txt', 'w')
f.write('In-Gate %')
f.write('     TP from Stim Below')
f.write('     TP from Unstim Above')
for item in range(0,len(Thresholds)):
	f.write('\n%.3f'%Thresholds[item])
	f.write('     ')
	f.write('%.3f'%LowTruePosPostSort[item])
	f.write('     ')
	f.write('%.3f'%HighTruePosPostSort[item])
f.close()




