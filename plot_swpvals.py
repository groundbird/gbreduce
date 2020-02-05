import numpy as np
import matplotlib.pyplot as plt
import os

outputdir='/Volumes/Toshiba5TB2/GroundBIRD/analysis/'

folderlist = os.listdir(outputdir)
todolist = []
for folder in folderlist:
	if '.' not in folder:
		subfolderlist = os.listdir(outputdir+folder)
		search = 'log.txt'
		test = [f for f in subfolderlist if search in f]
		# print(test)
		if test != []:
			todolist.append(outputdir+folder+'/'+test[0])
print(todolist)

i = 0
arrvals = ['freq', 'arga', 'absa', 'tau', 'fr', 'Qr', 'Qc', 'phi0', 'c']
values = []
for file in todolist:
	infile = open(file,"r+")
	for line in infile:
		if "[{'arga" in line:
			vals = line.split('}, {')
			for val in vals:
				arr = np.zeros(9)
				split = val.split(',')
				for s in split:
					arr[arrvals.index(s.split("'")[1])] = float(s.split(':')[1].replace('}]','').strip())
				arr[0] = np.round(arr[4]*1e-7)
				# int(str(arr[4])[0:3])
				print(val)
				print(arr)
				values.append(arr)
				i += 1
	infile.close()
	# print(values)
	# if len(values) > 50:
	# 	break
print('hi')
values = np.asarray(values)
print(np.shape(values))
# print(values[:,0])

plt.style.use('classic')
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)

plt.plot(values[:,4],np.ones(len(values[:,4])),'x')
plt.savefig(outputdir+'_freqs.pdf')
# exit()

print('Freqs:')
freqs = np.unique(np.asarray(values[:,0],dtype=int))
print(freqs)
# exit()
print('Doing plots')
for i in range(0,len(freqs)):
	fig, axs = plt.subplots(4,2)
	axs[0,0].plot(values[values[:,0]==freqs[i],1])
	axs[0,0].set_title(arrvals[1])
	axs[0,1].plot(values[values[:,0]==freqs[i],2])
	axs[0,1].set_title(arrvals[2])
	axs[1,0].plot(values[values[:,0]==freqs[i],3])
	axs[1,0].set_title(arrvals[3])
	axs[1,1].plot(values[values[:,0]==freqs[i],4]*1e-6)
	axs[1,1].set_title(arrvals[4] + ' x1e6')
	axs[2,0].plot(values[values[:,0]==freqs[i],5])
	axs[2,0].set_title(arrvals[5])
	axs[2,1].plot(values[values[:,0]==freqs[i],6])
	axs[2,1].set_title(arrvals[6])
	axs[3,0].plot(values[values[:,0]==freqs[i],7])
	axs[3,0].set_title(arrvals[7])
	axs[3,1].plot(values[values[:,0]==freqs[i],8]*1e9)
	axs[3,1].set_title(arrvals[8]+' x1e-9')

	# for ax in axs.flat:
	    # ax.set(xlabel='Count', ylabel='Val')

	plt.setp(axs[0,0].get_xticklabels(), visible=False)
	plt.setp(axs[0,1].get_xticklabels(), visible=False)
	plt.setp(axs[1,0].get_xticklabels(), visible=False)
	plt.setp(axs[1,1].get_xticklabels(), visible=False)
	plt.setp(axs[2,0].get_xticklabels(), visible=False)
	plt.setp(axs[2,1].get_xticklabels(), visible=False)
	# plt.setp(axs[3,0].get_xticklabels(), visible=False)
	# plt.setp(axs[3,1].get_xticklabels(), visible=False)

	# # Hide x labels and tick labels for top plots and y ticks for right plots.
	# for ax in axs.flat:
	#     ax.label_outer()
	plt.savefig(outputdir+'freq_'+str(freqs[i])+'.png')
	plt.clf()
print(np.unique(freqs))


# plt.title('Positions of max values in data minus positions of moon')
# plt.ylabel('Distance in elevation')
# plt.xlabel('Distance in azimuth')
# plt.savefig(outputdir+'pixel_positions_from_moon.pdf')
# plt.xlim([-15,15])
# plt.ylim([-15,15])
# plt.savefig(outputdir+'pixel_positions_from_moon_zoom.pdf')

# plt.xlim([-5,5])
# plt.ylim([-10,10])
# plt.savefig(outputdir+'pixel_positions_from_moon_zoom_comp.pdf')
