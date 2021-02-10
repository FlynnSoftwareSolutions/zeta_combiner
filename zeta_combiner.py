#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''ELSCompiler combines multiple zeta potential scans into one Excel workbook

ELSCompiler combines multiple electrophoretic mobility / zeta potential scans
performed by a Malvern NanoBrook Omni and exported as separate csv files,
because the Malvern database software cannot export all repeat scans of the
same sample into one file except as a PDF report.
'''
__author__ = "Michael Flynn"
__date__ = "20210210"

from copy import deepcopy
from os import path, listdir
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

# Get directory containing this script
try:
    currentFile = path.realpath(__file__).replace('\\','/')
except:
    currentFile = input(
        ('Python needs your help determining the location of this script.\n'
        'Please copy the full filepath of this Python script from above,\n'
        'then paste it here using right-click. Then press Enter to begin.\n'
        )
        )
currentFile = currentFile.replace('\\','/')
scriptdir = currentFile[:currentFile.rfind('/')+1]

binSize = 1 # bin size in mV if rebinning, 0 to not rebin
# Get first csv files of each experiment in script directory
fileGroupPaths = [f'{scriptdir}{f.replace("1.csv","")}'
                for f in listdir(scriptdir) if '-1.csv' in f or '_1.csv' in f]
# If any valid csv files found
if len(fileGroupPaths):
    # Process one group at a time
    for fileGroupPath in fileGroupPaths:
        # Get experiment name without counter or file extension
        expName = fileGroupPath.replace(scriptdir,'')[:-1]
        # Count how many files have same base experiment name
        numFiles = 0
        while path.isfile(f'{fileGroupPath}{numFiles + 1}.csv'):
            numFiles += 1
        if numFiles:
            # Count through files in group, all of which are repeat
            # measurements of the same experiment
            for i in range(numFiles):
                # Load csv to pandas dataframe
                df = pd.read_csv(
                        f'{fileGroupPath}{i + 1}.csv',
                        sep=',',
                        )
                # Process columns of pandas dataframe
                rawZeta = df['Zeta Potential']
                rawPower = df['Power']
                # Determine number of zeta potentials on x axis
                # This number should be constant between runs
                nPoints = len(rawZeta)
                # If first file
                if i == 0:
                    # Create arrays to store multiple files of data in columns
                    rawZetas = np.zeros((nPoints, numFiles), dtype=np.float64)
                    rawPowers = np.zeros((nPoints, numFiles), dtype=np.float64)
                # Store current file data in appropriate column
                rawZetas[:,i] = rawZeta
                rawPowers[:,i] = rawPower
            # Sort all files' zeta potential axes and preserve pairing of
            # relative power values with their corresponding zeta potential
            # values
            zetasSortedTuple, powerSortedTuple = zip(*sorted(zip(
                rawZetas.flatten(),
                rawPowers.flatten())))
            if binSize:
                # Binned array mode
                zetas = np.arange(
                        zetasSortedTuple[0],
                        zetasSortedTuple[-1],
                        binSize,
                        dtype=np.float64,
                        )
            else:
                # All unique zeta values mode
                # find all unique zeta potential values collected
                # (they vary between runs)
                zetas = np.unique(zetasSortedTuple)
            # Calculate number of bins or number of
            # unique zeta potential values in data set
            nBins = len(zetas)
            powers = np.zeros([nBins], dtype=np.float64)
            powerStd = deepcopy(powers)
            powerN = deepcopy(powers)
            indsTable = [[]]*nBins
            if binSize:
                for j,b in enumerate(zetas):
                    indsTable[j] = np.where(np.all((
                        zetasSortedTuple < b + 0.5*binSize,
                        zetasSortedTuple >= b - 0.5*binSize,
                        ),axis=0))[0]
            else:
                for j,zeta in enumerate(zetas):
                    indsTable[j] = np.where(zetasSortedTuple==zeta)[0]
            for j,inds in enumerate(indsTable):
                powers[j] = np.mean(np.take(
                    powerSortedTuple,inds)) if len(inds) else 0
                powerStd[j] = np.std(np.take(
                    powerSortedTuple,inds)) if len(inds) > 1 else 0
                powerN[j] = len(inds)
            # Smoothing nSmooths times using moving
            # average 1% of data set in size
            nSmooths = 2
            powerSmoothed = deepcopy(powers)
            w = int(nBins * 0.01)
            if w:
                if binSize:
                    box = np.ones(w, dtype=np.float64)
                else:
                    box = np.ones(w, dtype=np.float64)
                for i in range(nSmooths):
                    powerSmoothed = np.convolve(
                            powerSmoothed,
                            box,
                            mode='same',
                            ) / int(0.01 * nBins)
            powerSmoothed = powerSmoothed/np.max(powerSmoothed)
            # Plotting smoothed data and saving as svg
            if binSize:
                plt.plot(zetas,powerSmoothed)
            else:
                plt.plot(zetas,powerSmoothed)
            plt.xlabel('Zeta Potential [mV]')
            plt.ylabel('Relative Power')
            plt.title(expName)
            plt.savefig(expName+'.svg')
            plt.cla()
           
            # Saving smoothed and raw data as xlsx
            headers = [expName+str(i) for i in range(1,numFiles+1)]
            writer = pd.ExcelWriter(
                scriptdir+expName+'.xlsx', engine='xlsxwriter',)
            df = pd.DataFrame(np.vstack((
                zetas,
                powerSmoothed,
                powerStd,
                powerN,
                )).T)
            df.to_excel(writer,
                sheet_name='Summary',
                header=['Zeta Potential [mV]','Relative Power','Stdev','N'],
                index=False,
                )
            df = pd.DataFrame([''])
            df.to_excel(writer,
                sheet_name='RawData',
                header=['Zeta Potential [mV]'],
                index=False,
                )
            df2 = pd.DataFrame(rawZetas)
            df2.to_excel(writer,
                sheet_name='RawData',
                startcol=1,
                header=headers,
                index=False,
                )
            df3 = pd.DataFrame([''])
            df3.to_excel(writer,
                sheet_name='RawData',
                startcol=numFiles+1,
                header=['Relative Power'],
                index=False,
                )
            df4 = pd.DataFrame(rawPowers)
            df4.to_excel(writer,
                sheet_name='RawData',
                startcol=numFiles+2,
                header=headers,
                index=False,
                )
           
            # Producing summary line plot in excel
            fontsize = 20
            chart = writer.book.add_chart({'type': 'scatter',
                                           'subtype':'straight',})
            chart.add_series({'name': expName,
                              'categories': ['Summary', 1, 0, nBins, 0],
                              'values': ['Summary', 1, 1, nBins, 1],
                              'line': {'color': 'red'},
                              })
            chart.set_title({'name':expName,'name_font':{'size':20},})
            chart.set_x_axis({'name': 'Zeta Potential [mV]',
                              'name_font': {'size': fontsize},
                              'num_font': {'size': fontsize},
                              'major_gridlines': {'visible': False},
                              'min': -150,'max': 150, 'major_unit':25,
                              'label_position': 'low',
                              })
            chart.set_y_axis({'name': 'Relative Power',
                              'name_font': {'size': fontsize},
                              'num_font': {'size': fontsize},
                              'major_gridlines': {'visible': False},
                              'min': 0,'label_position': 'low',
                              })
            chart.set_legend({'position': 'none'})
            chart.set_size({'width': 1280, 'height': 720})
            writer.sheets['Summary'].insert_chart('F2', chart)

            # Producing all lines on one plot in excel
            chart = writer.book.add_chart({
                'type': 'scatter',
                'subtype':'straight',
                })
            for i in range(numFiles):
                chart.add_series({
                    'name': expName+'-'+str(i),
                    'categories': ['RawData', 1, i, nBins, i],
                    'values': ['RawData', 1, numFiles+i+1,
                               nBins, numFiles+i+1],
                    'line': {'color': 'blue'},
                    })
            chart.add_series({
                'name': expName,
                'categories': ['Summary', 1, 0, nBins, 0],
                'values': ['Summary', 1, 1, nBins, 1],
                'line': {'color': 'red'},
                })
            chart.set_title({
                'name':expName,
                'name_font':{'size':20},
                })
            chart.set_x_axis({
                'name': 'Zeta Potential [mV]',
                'name_font': {'size': fontsize},
                'num_font': {'size': fontsize},
                'major_gridlines': {'visible': False},
                'min': -150, 'max': 150,
                'major_unit': 25,
                'label_position': 'low',
                })
            chart.set_y_axis({
                'name': 'Relative Power',
                'name_font': {'size': fontsize},
                'num_font': {'size': fontsize},
                'major_gridlines': {'visible': False},
                'min': 0,
                'label_position': 'low',
                })
            chart.set_legend({
                'position': 'right',
                'name_font': {'size': fontsize},
                })
            chart.set_size({
                'width': 1280,
                'height': 720,
                })
            writer.sheets['RawData'].insert_chart('D2', chart)
            writer.close()
