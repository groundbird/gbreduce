import numpy as np
import matplotlib.pyplot as plt
import os

outputdir='/Volumes/Toshiba5TB2/GroundBIRD/analysis_az42/'
# infile = 'gb_pixinfo_20190919.txt'
infile = 'gb_pixinfo_20191217.txt'
npix,mod,mod_pix,x,y,theta,phi,omt1,omt2, pixelfreq = np.loadtxt(infile,unpack=True)
pixelfreq2 = pixelfreq.tolist()
focallength = 500.0
x = (x/focallength)*180.0/np.pi
y = (y/focallength)*180.0/np.pi
# x[0:14] *= (220.0/150.0)
# y[0:14] *= (220.0/150.0)
print(x)
print(y)
plt.xlim([-6,2])
plt.xlabel('Azimuth offset')
plt.ylabel('Elevation offset')
plt.plot(x,y,'o')
plt.savefig(outputdir+'pixel_expected_pos.pdf')
plt.clf()

angle = 180.0-84.0
new_x = x*np.cos(angle*np.pi/180.0) - y*np.sin(angle*np.pi/180.0)
new_y = x*np.sin(angle*np.pi/180.0) + y*np.cos(angle*np.pi/180.0)
plt.xlim([-6,2])
plt.xlabel('Azimuth offset')
plt.ylabel('Elevation offset')
plt.plot(new_x,new_y,'o')
for i, txt in enumerate(mod_pix):
    plt.annotate(txt, (new_x[i], new_y[i]))
plt.savefig(outputdir+'pixel_expected_pos_rot.pdf')
plt.clf()

docorr=True

# exit()

moonlist = ['20200107_data_211314_swp_newaz42_GB01','20200107_data_211314_swp_newaz42_GB02','20200107_data_232558_swp_newaz42_GB01','20200108_data_214017_swp_newaz42_GB01','20200108_data_222408_swp_newaz42_GB02','20200108_data_235207_swp_newaz42_GB01','20200108_data_235207_swp_newaz42_GB02','20200109_data_005248_swp_newaz42_GB01','20200109_data_005248_swp_newaz42_GB02','20200109_data_223525_swp_newaz42_GB01','20200109_data_223525_swp_newaz42_GB02','20200109_data_231913_swp_newaz42_GB02','20200110_data_005219_swp_newaz42_GB02','20200110_data_014538_swp_newaz42_GB01','20200110_data_014538_swp_newaz42_GB02','20200110_data_230732_swp_newaz42_GB01','20200110_data_230732_swp_newaz42_GB02','20200110_data_235634_swp_newaz42_GB02','20200111_data_005133_swp_newaz42_GB02','20200111_data_011920_swp_newaz42_GB02','20200111_data_021950_swp_newaz42_GB01','20200111_data_030412_swp_newaz42_GB01','20200111_data_030412_swp_newaz42_GB02','20200113_data_015007_swp_newaz42_GB01','20200113_data_015007_swp_newaz42_GB02','20200113_data_025032_swp_newaz42_GB02','20200113_data_035056_swp_newaz42_GB01','20200113_data_035056_swp_newaz42_GB02']

# moonlist = ['20200107_data_211314_swp_newaz42_GB01','20200113_data_035056_swp_newaz42_GB01']


folderlist = os.listdir(outputdir)
todolist = []
for folder in folderlist:
	# if '.' not in folder:
	if folder in moonlist:
		subfolderlist = os.listdir(outputdir+folder)
		search = 'log.txt'
		test = [f for f in subfolderlist if search in f]
		# print(test)
		if test != []:
			todolist.append(outputdir+folder+'/'+test[0])
	# 		print('No data')
	# 		continue
	# 	test = [f for f in subfolderlist if 'tod' in f]
	# 	if test == []:
	# 		print('No data')
	# 		continue
	# 	subfolderlist = [f for f in subfolderlist if search not in f]
	# 	if len(subfolderlist) > 2:
	# 		todolist.append(subdir+folder)
print(todolist)

if docorr == False:
	plt.plot(new_x,new_y,'o')
	for i, txt in enumerate(mod_pix):
	    plt.annotate(txt, (new_x[i], new_y[i]))

i = 0
arrvals = ['freq', 'arga', 'absa', 'tau', 'fr', 'Qr', 'Qc', 'phi0', 'c']
outfile = open(outputdir+'poslist.txt', "w+")

for i in range(len(new_x)):
	outfile.write(str(new_x[i]) + " " + str(new_y[i]) + " " + str(i) + "\n")
outfile.write('\n')

for file in todolist:
	print(file)
	infile = open(file,"r+")
	az = []
	el = []
	freq = []
	values = []
	j = 0
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
		if 'Difference is' in line:
			line = line.replace('Difference is','').replace('in azimuth and',',').replace('in elevation.','').strip()
			line = line.split(',')
			print(i)
			if j % 4 == 3:
				azval = -(float(line[0])-0.75)*np.cos(70.0*np.pi/180.0)
				elval = -float(line[1])#+0.5
				if elval < -9.0:
					elval += 10.0

				print(j)
				print((j-3)/2)
				freqval = values[int((j-3)/2)][4]

				if docorr==True:
					matchfreq = int(freqval/1e6)
					print(pixelfreq2)
					print(matchfreq)
					try:
						match = pixelfreq2.index(float(matchfreq))
						if match >= 0:
							azval -= new_x[match]
							elval -= new_y[match]
					except:
						null = 0
				az.append(azval)
				el.append(elval)
				freq.append(freqval)
				outfile.write(str(az[-1]) + " " + str(el[-1]) + " " + str(freq[-1])+'\n')
				print(str(az[-1]) + " " + str(el[-1]) + " " + str(freq[-1])+'\n')
				# input('continue?')
			i+=1
			j+=1
	# if freq[0] > 6.0e9:
	plt.plot(az,el,'x')
	# for k, txt in enumerate(freq):
	    # plt.annotate(int(txt/1e6), (az[k], el[k]),color='g')
	infile.close()
outfile.close()
if docorr == False:
	plt.title('Positions of max values in data minus positions of moon')
	plt.ylabel('Distance in elevation')
	plt.xlabel('Distance in azimuth')
	plt.savefig(outputdir+'pixel_positions_from_moon.pdf')
	plt.xlim([-15,15])
	plt.ylim([-15,15])
	plt.savefig(outputdir+'pixel_positions_from_moon_zoom.pdf')

	plt.xlim([-5,5])
	plt.ylim([-10,10])
	plt.savefig(outputdir+'pixel_positions_from_moon_zoom_comp.pdf')

	plt.xlim([0,1])
	plt.ylim([-0.5,1])
	plt.savefig(outputdir+'pixel_positions_from_moon_zoom_comp_sron.pdf')

	plt.xlim([-2,1])
	plt.ylim([5.0,7])
	plt.savefig(outputdir+'pixel_positions_from_moon_zoom_comp_riken.pdf')
else:
	plt.title('Positions of max values in data minus positions of moon, with position corrections')
	plt.ylabel('Distance in elevation')
	plt.xlabel('Distance in azimuth')
	plt.savefig(outputdir+'pixel_positions_corrected.pdf')
	plt.xlim([-1.0, 1.0])
	plt.ylim([-1.0, 1.0])
	plt.savefig(outputdir+'pixel_positions_corrected_zoom.pdf')
