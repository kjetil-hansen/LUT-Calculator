##LUT Calculation Program
##KJ Hansen
##Created: 09/05/2019
#Last edited: 26/09/2019

##Import required modules
import os
import sys
import glob
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from scipy import interpolate

###Required Functions
#Query continue
def query_yes_no(question, default="yes"):
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    prompt = " [Y/n] "
    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

#Query which charge?
def query_charge(question, default='2'):
    valid = {"1": 1, "2": 2, "3": 3, "4": 4, "all": 0}                           ##valid options should be unique peptide.z values
    prompt = ("Please choose: 1, 2, 3, 4 or all ")
    while True:
        sys.stdout.write(question + prompt)
        charge = input().lower()
        if charge == '':
            return valid[default]
        elif charge in valid:
            return valid[charge]
        else:
            sys.stdout.write("Please respond with '1', '2', '3', '4' or 'all' \n")

#Query colour by intensity?
def color_int(question, default="yes"):
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    prompt = ("[Y/n]")
    while True:
        sys.stdout.write(question + prompt)
        cint = input().lower()
        if cint == '':
            return valid[default]
        elif cint in valid:
            return valid[cint]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

##Intro Text
print('\n\n\n--------------------LUT Calculator--------------------',
    '\nThis python script takes PLGS ion accounting files and creates a ',
    '\nlook up table (LUT) containing optimised transfer collision energy values',
    '\nbased on the mobility of your peptide ions for your Synapt G2Si system.\n')

##Input options
parser = argparse.ArgumentParser()
parser.add_argument('-f', help='Filepath to folders containing PLGS outputs (No slashes) (Default = Current Directory)')
parser.add_argument('-s', help='Name of protein for filename output (optional)')
args = parser.parse_args()

##Set paths to PLGS files
if args.f:
    path = os.getcwd() + '/' + args.f
    path2 = path + "/**/*peptide.csv"
    print('\nLooking in filepath: ' + path + '\n')
else:
    path = os.getcwd()
    path2 = path + "/**/*peptide.csv"
    print('\nSearching in current directory and subdirectories: ' + path + '\n')


##Find PLGS files and print
filenames = glob.glob(path2, recursive=True)
print('\n' + str(len(filenames)) + ' PLGS files found: \n', *filenames, sep='\n')

#Query continue?
answer = query_yes_no('\nCreate CE LUT using these files?')
if answer == True:
    print('\nCombining reference files...\n')
else:
    sys.exit('\nExiting...\n')

##Load PLGS files into dataframes (encoding due to degree symbol not recognised by utf-8)
dfs = []
for filename in filenames:
    dfs.append(pd.read_csv(filename, encoding = "latin-1"))

##Concatenate all data into one DataFrame
mt1 = pd.concat(dfs, ignore_index=True, sort=True)

##Create ID column
mt1['peptide.id'] = mt1[['peptide.seq', 'peptide.modification']].apply(lambda x: '-'.join(x.map(str)), axis=1)

#Only green filtered peptides
mt1 = mt1[mt1['peptide.AutoCurate'] == 'Green']

##Create table with only mob, m/z, z and ID
mt1 = mt1[['peptide.id', 'precursor.z', 'precursor.mz', 'precursor.Mobility', 'precursor.inten']].copy()

##Find duplicate entries and average
mt1 = mt1.groupby('peptide.id').mean().reset_index()

##User defines charge to build LUT
charge = query_charge('Which charge state should be used to create the LUT? Default = 2+ (recommended)\n')
if charge == 0:
    print('Creating LUT based on all charge states... \n')
else:
    print('\nCreating LUT based on ' + str(charge) + '+ charge peptides...\n')
    mt1 = mt1.loc[mt1['precursor.z'] == charge]

##If no charge state found, exit
if mt1.empty:
     print("No peptides of " + str(charge) + "+ charge state found.")
     sys.exit('\nExiting...\n')

##Order by mobility
mt1 = mt1.sort_values(by=['precursor.Mobility'])

##Calculate Transfer CE
mt1['transfer.ce'] = (mt1['precursor.mz']*0.034)+8                               #EDIT FORMULA HERE FOR OTHER SYSTEMS!

##Plot scatter of bins vs CE
mt1.dropna(inplace=True)                                                         #removes nan values due to csv format
x = mt1['precursor.Mobility']
y = mt1['transfer.ce']
z = mt1['precursor.inten']
z2 = mt1['precursor.z']
p = np.polyfit(x, y, 1)

##Colour by intensity?
if charge != 0:
    int = color_int("\nColour regression by precursor signal intensity? ")
    if int == True and charge !=0:                                                  #only one that works
        print('\nColouring by intensity... \n')
        zx = z
        nx = colors.LogNorm(vmin=z.min(), vmax=z.max())                             #normalize z logarithmically
        alx = 0.8
        cm = 'Reds'
    elif int == False and charge !=0:
        print('\nDefault colours used... \n')
        zx = None
        nx = None
        alx = 0.8
        cm = None
else:
    zx = z2
    nx = None
    alx = 0.8
    cm = 'Dark2'

f1 = plt.figure(1)
plt.scatter(x, y, c=zx, edgecolors = 'none', s=80, alpha = alx, norm=nx, cmap = cm)

##Plot line of best fit
plt.plot(x, np.poly1d(np.polyfit(x, y, 1))(x), c = '#7bb8ed')
sl = p[0]
int = p[1]
f1.suptitle('LUT Regression', fontsize=16)
plt.xlabel('Mobility (Bins)', fontsize=12)
plt.ylabel('Theoretical Collision Energy (V)', fontsize=12)

##Calculate LUT values
lutdata = [[0, 0], [19, 0], [20, (sl*20)+int], [120, (sl*120)+int], [200, ((sl*120)+int)+5]]

##Store lut values
lut = pd.DataFrame(lutdata, columns=['Bin', 'CE'])

##Create LUT gradient
x = lut['Bin']
y = lut['CE']
f = interpolate.interp1d(x, y)
xnew = np.linspace(x[0],x[4],num=201)

##Create gradient figure
f2 = plt.figure(2)
plt.plot(xnew, f(xnew))
plt.plot(x, y,'o')
f2.suptitle('LUT Gradient', fontsize=16)
plt.xlabel('Mobility (Bins)', fontsize=12)
plt.ylabel('Optimised Collision Energy (V)', fontsize=12)

##Variable for full gradient
lut2 = pd.DataFrame({'Bin':xnew, 'CE':f(xnew)})

##Define output folder and filename
if args.s:
    path3 = os.getcwd() + '/Output/' + args.s + r'/LUTCalc_Output/'
else:
    path3 = path + r'/LUTCalc_Output/'

os.makedirs(path3, exist_ok=True)

if args.s:
    pathout = path3 + args.s + '_LUT_' + str(charge) + 'Charge'
else:
    pathout = path3 + 'LUT_' + str(charge) + 'Charge'

##Export peptides CSV and save figures
mt1.to_csv(pathout + '_peptides.csv', index=None, header=True)
grad_out = f2.savefig(pathout + '_gradient.png')
reg_out = f1.savefig(pathout + '_regression.png')

##Export LUT to CSV
lut.to_csv(pathout + '_LUT.csv', index=None, header=True)

##Endscreen
print("LUT succesfully calculated:\n" + str(lut))
print("\n---------------Additional Information----------------\n\n"
    "\nThe following 3 additional files were created: \n\n" + pathout + "_peptides.csv\n"
    + pathout + "_regression.png\n" + pathout + "_gradient.png\n")
