#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''ELSCompiler combines multiple zeta potential scans into one Excel workbook

ELSCompiler combines multiple electrophoretic mobility / zeta potential scans
performed by a Malvern NanoBrook Omni and exported as separate csv files,
because the Malvern database software cannot export all repeat scans of the
same sample into one file except as a PDF report.
'''
__author__ = "Michael Flynn"
__date__ = "20210308"

import os
import sys
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

def processFolder(d):
    # Remove trailing slash on directory d
    d = d[:-1] if d[-1] == '/' else d
    # Isolate name of deepest folder in full path directory d
    endFolderName = d[d.rfind('/') + 1:]
    # Get first csv files of each experiment in script directory
    fileGroupPaths = [f'{d}/{f.replace("1.csv","")}'
                for f in os.listdir(scriptdir) if '-1.csv' in f or '_1.csv' in f]
    nFileGroups = len(fileGroupPaths)
    # If any valid csv files found
    if nFileGroups:
        # Get experiment name without counter or file extension
        expNames = [fileGroupPath.replace(f'{d}/','')[:-1]
                    for fileGroupPath in fileGroupPaths]
        writer = pd.ExcelWriter(
                f'{d}/{endFolderName}_Summary.xlsx',
                engine='xlsxwriter',
                )
        # Label every 4th column of Summary sheet with the name of experiment
        # whose data will be placed there
        spacedNames = [''] * (4 * nFileGroups)
        for z,expName in zip(range(0, 4 * nFileGroups, 4), expNames):
            spacedNames[z] = expName
        df = pd.DataFrame([spacedNames])
        df.to_excel(writer,
                    sheet_name='Summary',
                    startcol=22,
                    header=False,
                    index=False,
                    )
        df = pd.DataFrame(["Sample", "Mode Zeta Potential [mV]"]).T
        df.to_excel(writer,
                    sheet_name='Summary',
                    header=False,
                    index=False,
                    )

        # Produce summary line plot in Summary sheet
        summaryChart = writer.book.add_chart({
            'type': 'scatter',
            'subtype':'straight',
            })
        summaryChart.set_title({
            'name': endFolderName,
            'name_font': {'size': 20},
            })
        summaryChart.set_x_axis({
            'name': 'Zeta Potential [mV]',
            'name_font': {'size': fontsize},
            'num_font': {'size': fontsize},
            'major_gridlines': {'visible': False},
            'min': -150,
            'max': 150,
            'major_unit': 25,
            'label_position': 'low',
            })
        summaryChart.set_y_axis({
            'name': 'Relative Power',
            'name_font': {'size': fontsize},
            'num_font': {'size': fontsize},
            'major_gridlines': {'visible': False},
            'min': 0,
            'label_position': 'low',
            })
        summaryChart.set_legend({
            'position': 'right',
            'name_font': {'size': fontsize},
            })
        summaryChart.set_size({
            'width': 1280,
            'height': 720,
            })
        writer.sheets['Summary'].insert_chart('C2', summaryChart)

        # Process each group of CSV files
        for fileGroupNumber, expName, fileGroupPath in \
                zip(range(nFileGroups), expNames, fileGroupPaths):
            n = processFileGroup(writer, fileGroupNumber,
                                 expName, fileGroupPath)
            # Produce line on chart in Summary sheet
            summaryChart.add_series({
                'name': ['Summary', 0, 4 * fileGroupNumber + 22,
                         0, 4 * fileGroupNumber + 22],
                'categories': ['Summary', 2, 4 * fileGroupNumber + 22,
                               n + 1, 4 * fileGroupNumber + 22],
                'values': ['Summary', 2, 4 * fileGroupNumber + 23,
                           n + 1, 4 * fileGroupNumber + 23]})
        writer.close()
    else:
        print('No CSV file groups found')

def processFileGroup(writer, fileGroupNumber, expName, fileGroupPath):
    # Count how many files have same base experiment name
    nFiles = 0
    while os.path.isfile(f'{fileGroupPath}{nFiles+1}.csv'):
        nFiles += 1
    if nFiles:
        # Load csv to pandas dataframe
        df = pd.read_csv(f'{fileGroupPath}1.csv', sep=',')
        # Process columns of pandas dataframe
        rawZeta = df['Zeta Potential'].values
        rawPower = df['Power'].values
        # Determine number of zeta potentials on x axis
        # This number should be constant between runs
        nPoints = len(rawZeta)
        # Create arrays to store multiple files of data in columns
        rawZetas = np.zeros((nFiles, nPoints), dtype=np.float64)
        rawPowers = np.zeros((nFiles, nPoints), dtype=np.float64)
        rawZetas[0] = rawZeta
        rawPowers[0] = rawPower
        # Count through files in group, all of which are repeat
        # measurements of the same experiment
        for i, rawZeta, rawPower in zip(
                range(2, nFiles + 1), rawZetas[1:], rawPowers[1:]):
            # Load csv to pandas dataframe
            df = pd.read_csv(f'{fileGroupPath}{i}.csv', sep=',')
            # Process columns of pandas dataframe
            rawZeta[:] = df['Zeta Potential'].values
            rawPower[:] = df['Power'].values
        # Sort all files' zeta potential axes and preserve pairing of 
        # relative power values with their corresponding zeta potential values
        z = rawZetas.flatten()
        p = rawPowers.flatten()
        indOrder = z.argsort()
        sortedZetas = z[indOrder]
        sortedPowers = p[indOrder]
        del z, p
        if binSize:
            # Binned array mode
            zetas = np.arange(
                    sortedZetas[0],
                    sortedZetas[nFiles * nPoints - 1],
                    binSize,
                    dtype=np.float64,
                    )
        else:
            # All unique zeta values mode
            # Find all unique zeta potential values collected
            # (they vary between runs)
            zetas = np.unique(sortedZetas)
        # Calculate number of bins or number of
        # unique zeta potential values in data set
        nBins = len(zetas)
        powers = np.zeros(nBins, dtype=np.float64)
        powerStdevs = np.zeros(nBins, dtype=np.float64)
        powerNs = np.zeros(nBins, dtype=np.float64)
        halfBinSize = 0.5 * binSize
        if binSize:
            inds = [np.argmax(sortedZetas >= zeta - halfBinSize)
                    for zeta in zetas]
            inds += [nFiles * nPoints]
        else:
            # Untested
            inds = [np.argmax(sortedZetas == zeta) for zeta in zetas]
        
        for i, startInd, stopInd in \
                zip(range(nBins), inds[:nFiles * nPoints - 1], inds[1:]):
            powerN = stopInd - startInd
            powerNs[i] = powerN
            yVals = sortedPowers[startInd:stopInd]
            powers[i] = np.mean(yVals) if powerN else 0
            powerStdevs[i] = np.std(yVals) if powerN > 1 else 0
        # Smooth nSmooths times using moving average 1% of data set in size
        nSmooths = 2
        powerSmoothed = np.copy(powers)
        w = int(nBins * 0.01)
        if w:
            box = np.ones(w, dtype=np.float64)
            for _ in range(nSmooths):
                powerSmoothed = np.convolve(
                        powerSmoothed,
                        box,
                        mode='same',
                        ) / w
        powerSmoothed /= np.max(powerSmoothed)

        # Plot smoothed data and save as svg
        if binSize:
            plt.plot(zetas, powerSmoothed)
        else:
            plt.plot(zetas, powerSmoothed)
        plt.xlabel('Zeta Potential [mV]')
        plt.ylabel('Relative Power')
        plt.title(expName)
        plt.savefig(f'{fileGroupPath}.svg')
        plt.cla()
        
        # Save raw data to Excel sheet
        headers = [f'{expName}-{i}' for i in range(1, nFiles+1)]
        df = pd.DataFrame([expName, 'Zeta Potential [mV]'])
        df.to_excel(writer,
                    sheet_name=f"{fileGroupNumber + 1}",
                    header=False,
                    index=False,
                    )
        df2 = pd.DataFrame(rawZetas.T)
        df2.to_excel(writer,
                     sheet_name=f"{fileGroupNumber + 1}",
                     startrow=1,
                     startcol=1,
                     header=headers,
                     index=False,)
        df3 = pd.DataFrame(['Relative Power'])
        df3.to_excel(writer,
                     sheet_name=f"{fileGroupNumber + 1}",
                     startcol=nFiles + 1,
                     header=False,
                     index=False,
                     )
        df4 = pd.DataFrame(rawPowers.T)
        df4.to_excel(writer,
                     sheet_name=f"{fileGroupNumber + 1}",
                     startrow=1,
                     startcol=nFiles + 2,
                     header=headers,
                     index=False,
                     )
        
        # Produce all lines on chart in Sheet numbered for
        # the current file group in Excel
        chart = writer.book.add_chart({
            'type': 'scatter',
            'subtype': 'straight',
            })
        for i in range(nFiles):
            chart.add_series({
                'name': f'{expName}-{i}',
                'categories': [f"{fileGroupNumber + 1}", 2,
                               i, nBins + 1, i],
                'values': [f"{fileGroupNumber + 1}", 2, nFiles + i + 1,
                           nBins + 1, nFiles + i + 1],
                'line': {'color': 'blue'},
                })
        chart.add_series({
            'name': expName,
            'categories': ['Summary', 2, 4 * fileGroupNumber + 22,
                           nBins + 1, 4 * fileGroupNumber + 22],
            'values': ['Summary', 2, 4 * fileGroupNumber + 23,
                       nBins + 1, 4 * fileGroupNumber + 23],
            'line': {'color': 'red'},
            })
        chart.set_title({'name': expName, 'name_font': {'size':fontsize},})
        chart.set_x_axis({
            'name': 'Zeta Potential [mV]',
            'name_font': {'size': fontsize},
            'num_font': {'size': fontsize},
            'major_gridlines': {'visible': False},
            'min': -150,
            'max': 150,
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
        writer.sheets[f"{fileGroupNumber + 1}"].insert_chart('B2', chart)

        # Fill out four column area in Summary sheet with zeta values,
        # smoothed average power distribution, stdev, and N
        df = pd.DataFrame(np.vstack((
            zetas,
            powerSmoothed,
            powerStdevs,
            powerNs
            )).T)
        df.to_excel(writer,
                    sheet_name='Summary',
                    startrow=1,
                    startcol=4 * fileGroupNumber + 22,
                    header=['Zeta Potential [mV]',
                            'Relative Power',
                            'Stdev',
                            'N'
                            ],
                    index = False,
                    )

        # Fill out table of mode zeta values
        df = pd.DataFrame([
            expName,
            np.mean(zetas[np.where(powerSmoothed == 1.0)[0]])
            ]).T
        df.to_excel(writer,
                    sheet_name='Summary',
                    startrow=fileGroupNumber + 1,
                    header=False,
                    index=False,
                    )
        return nBins

if __name__ == "__main__":
    if len(sys.argv) == 2:
        scriptdir = sys.argv[1]
        scriptdir = scriptdir.replace('\\','/')
        scriptdir = scriptdir if scriptdir[-1] == '/' else f'{scriptdir}/'
    elif len(sys.argv) == 1:
        # Get directory containing this script
        try:
            currentFile = os.path.realpath(__file__).replace('\\','/')
        except:
            currentFile = input(
                ('Python needs help determining the location of this script\n'
                'Please copy the full filepath of this Python script\n'
                'then paste it here using right-click\n'
                'Then press Enter to begin.\n'))
        currentFile = currentFile.replace('\\','/')
        scriptdir = currentFile[:currentFile.rfind('/')]
    else:
        print('Only one command line argument may be provided: '
              'the full path to the directory containing the'
              'CSV files to be processed')

    # Required globals
    # Font size on graphs
    fontsize = 20
    # Set binSize in mV if rebinning is desired, 0 to disable rebinning
    binSize = 1

    processFolder(scriptdir)
