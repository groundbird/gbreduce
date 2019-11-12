# Quick script to read in the calibrated temperature files, and plot them.
import numpy as np
import matplotlib.pyplot as plt
import os

indir = "/Volumes/iOmega/GroundBIRD/data/logdata/thermo/2019/09/"

# Get the available dates to plot
infiles = os.listdir(indir)
prefixes = []
for line in infiles:
	print(line[0:8])
	if line[0:8] not in prefixes:
		prefixes.append(line[0:8])

# Settings for the different postfixes
postfixes = ['_detector.cal', '_shield.cal']#, '_he10.cal'
dtypes = [{'names':('date','unix','t1','t2','t3'),'formats':('S1','d','f','f','f')},\
{'names':('date','unix','t1','t2','t3','t4','t5','t6','t7','t8'),'formats':('S1','d','f','f','f','f','f','f','f','f')}]#\
#{'names':('date','unix','string','t1','t2','t3','string2','t4','t5','t6','t7','t8','t9','t10'),'formats':('S1','d','S1','f','f','f','S1','f','f','f','f','f','f','f')},]
numtemps=[3,8]#,10

for prefix in prefixes:
	plt.figure(figsize=(8.0, 5.0), dpi=100)
	for postfix_count in range(len(postfixes)):
		print(indir+prefix+postfixes[postfix_count])
		try:
			data = np.loadtxt(indir+prefix+postfixes[postfix_count],unpack=False,dtype=dtypes[postfix_count])
		except:
			continue
		print(data)
		print(data['unix'])
		for i in range(1,numtemps[postfix_count]+1):
			plt.plot(data['unix'],data['t'+str(i)],label=postfixes[postfix_count].replace('.cal','').replace('_','') + '_t'+str(i),linewidth=0.3)
	l = plt.legend(prop={'size':6},loc=1)
	l.set_zorder(20)
	plt.savefig(indir+prefix+'_plot.png', dpi=100)
	plt.ylim(0.2,0.4)
	plt.savefig(indir+prefix+'_plot_low.png', dpi=100)
	plt.ylim(30,35)
	plt.savefig(indir+prefix+'_plot_mid.png', dpi=100)
