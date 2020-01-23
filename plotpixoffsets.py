import numpy as np
import matplotlib.pyplot as plt
import os

outputdir='/Volumes/iOmega/GroundBIRD/analysis/'
outputdir=''
# infile = 'gb_pixinfo_20190919.txt'
infile = 'gb_pixinfo_20191217.txt'
npix,mod,mod_pix,x,y,theta,phi,omt1,omt2 = np.loadtxt(infile,unpack=True)
focallength = 500.0
x = (x/focallength)*180.0/np.pi
y = (y/focallength)*180.0/np.pi
x[0:14] *= (220.0/150.0)
y[0:14] *= (220.0/150.0)
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
plt.savefig(outputdir+'pixel_expected_pos_rot.pdf')
plt.clf()

exit()

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

plt.plot(new_x,new_y,'o')

i = 0
for file in todolist:
	infile = open(file,"r+")
	az = []
	el = []
	for line in infile:
		if 'Difference is' in line:
			line = line.replace('Difference is','').replace('in azimuth and',',').replace('in elevation.','').strip()
			line = line.split(',')
			if i % 4 == 3:
				az.append(-float(line[0]))
				el.append(-float(line[1]))
			i+=1
	plt.plot(az,el,'x')
	infile.close()
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
