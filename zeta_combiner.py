#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''ELSCompiler combines multiple zeta potential scans into one Excel workbook

ELSCompiler combines multiple electrophoretic mobility / zeta potential scans
performed by a Malvern NanoBrook Omni and exported as separate csv files,
because the Malvern database software cannot export all repeat scans of the
same sample into one file except as a PDF report.
'''
__author__ = "Michael Flynn"
__date__ = "20220712"

import os
import sys
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

def processFolder(d: str):
    # Remove trailing slash on directory d
    d = d[:-1] if d[-1] == '/' else d
    # Isolate name of deepest folder in full path directory d
    endFolderName = d[d.rfind('/') + 1:]
    # Get first data files of each experiment in script directory
    fileGroupPaths = [f'{d}/{f.replace("1" + extension, "")}'
                for f in os.listdir(scriptdir)
                if f.endswith(f'1{extension}')
                ]
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
                    startcol=27,
                    header=False,
                    index=False,
                    )
        df = pd.DataFrame(["Sample",
                           "Mode Zeta Potential [mV]",
                           "Mean of non-zero modes [mV]",
                           "Stdev of non-zero modes [mV]",
                           "N",
                           ]).T
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
            'min': minZeta,
            'max': maxZeta,
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
        writer.sheets['Summary'].insert_chart('G2', summaryChart)

        # Process each group of CSV files
        totalFilesProcessed = 0
        for fileGroupNumber, expName, fileGroupPath in \
                zip(range(nFileGroups), expNames, fileGroupPaths):
            nFiles, nBins = processFileGroup(
                writer, fileGroupNumber, totalFilesProcessed,
                nFileGroups, expName, fileGroupPath)
            # Produce line on chart in Summary sheet
            summaryChart.add_series({
                'name': ['Summary', 0, 4 * fileGroupNumber + 27,
                         0, 4 * fileGroupNumber + 27],
                'categories': ['Summary', 2, 4 * fileGroupNumber + 27,
                               nBins + 1, 4 * fileGroupNumber + 27],
                'values': ['Summary', 2, 4 * fileGroupNumber + 28,
                           nBins + 1, 4 * fileGroupNumber + 28]})
            # Update totalFilesProcessed for calculating row for
            # individual file statistics
            totalFilesProcessed += nFiles
        writer.close()
    else:
        print(f'No {extension} file groups found')

def processFileGroup(
        writer: pd.ExcelWriter, fileGroupNumber: int,
        totalFilesProcessed: int, nFileGroups: int,
        expName: str, fileGroupPath: str):
    # Count how many files have same base experiment name
    nFiles = 0
    while os.path.isfile(f'{fileGroupPath}{nFiles+1}{extension}'):
        nFiles += 1
    if nFiles:
        # Load data file to pandas dataframe
        if extension.startswith('.xls'):
            df = pd.read_excel(f'{fileGroupPath}1{extension}')
        elif extension == '.csv':
            df = pd.read_csv(f'{fileGroupPath}1{extension}', sep=',')
        # Process column of pandas dataframe
        rawZeta = df['Zeta Potential'].values
        # Determine number of zeta potentials on x axis
        # This number should be constant between runs
        nPoints = len(rawZeta)
        # Create array of x values to interpolate each run onto
        zetas = np.arange(
                minZeta,
                maxZeta,
                binSize,
                dtype=np.float64,
                )
        nBins = len(zetas)
        # Create arrays to store multiple files of data in columns
        rawZetas = np.empty((nFiles, nPoints), dtype=np.float64)
        rawPowers = np.empty((nFiles, nPoints), dtype=np.float64)
        powers = np.empty((nFiles, nBins), dtype=np.float64)
        powerNs = np.empty((nFiles, nBins), dtype=np.float64)
        modes = np.empty((nFiles), dtype=np.float64)
        rawZetas[0] = rawZeta
        # Count through files in group, all of which are repeat
        # measurements of the same sample
        for i, rawZeta, rawPower, power, powerN in zip(
                range(nFiles), rawZetas, rawPowers, powers, powerNs):
            # Load data file to pandas dataframe
            if extension.startswith('.xls'):
                df = pd.read_excel(f'{fileGroupPath}{i + 1}{extension}')
            elif extension == '.csv':
                df = pd.read_csv(f'{fileGroupPath}{i + 1}{extension}', sep=',')
            # Process columns of pandas dataframe
            rawZeta[:] = df['Zeta Potential'].values
            rawPower[:] = df['Power'].values
            # For each file, find zeta potential value where power is highest
            modes[i] = rawZeta[np.argmax(rawPower)]
            # Interpolate power values onto same x axis (zetas)
            power[:], powerN[:] = interpolate(rawZeta, rawPower, zetas)
        powerMeans = np.mean(powers, axis=0)
        powerStdevs = np.std(powers, ddof=1, axis=0)
        powerNs = np.sum(powerNs, axis=0)
        powerMeans /= np.max(np.nan_to_num(powerMeans))

        # Plot smoothed data and save as svg
        plt.plot(zetas, powerMeans)
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
                'name': [f"{fileGroupNumber + 1}", 1,
                         nFiles + i + 2, 1, nFiles + i + 2],
                'categories': [f"{fileGroupNumber + 1}", 2,
                               i, nBins + 1, i],
                'values': [f"{fileGroupNumber + 1}", 2, nFiles + i + 1,
                           nBins + 1, nFiles + i + 1],
                'line': {'color': 'blue'},
                })
        chart.add_series({
            'name': ['Summary', 1, 4 * fileGroupNumber + 27, 1, 4 * fileGroupNumber + 27],
            'categories': ['Summary', 2, 4 * fileGroupNumber + 27,
                           nBins + 1, 4 * fileGroupNumber + 27],
            'values': ['Summary', 2, 4 * fileGroupNumber + 28,
                       nBins + 1, 4 * fileGroupNumber + 28],
            'line': {'color': 'red'},
            })
        chart.set_title({'name': expName, 'name_font': {'size':fontsize},})
        chart.set_x_axis({
            'name': 'Zeta Potential [mV]',
            'name_font': {'size': fontsize},
            'num_font': {'size': fontsize},
            'major_gridlines': {'visible': False},
            'min': minZeta,
            'max': maxZeta,
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
            powerMeans,
            powerStdevs,
            powerNs
            )).T)
        df.to_excel(writer,
                    sheet_name='Summary',
                    startrow=1,
                    startcol=4 * fileGroupNumber + 27,
                    header=['Zeta Potential [mV]',
                            'Relative Power',
                            'Stdev',
                            'N'
                            ],
                    index = False,
                    )

        # Fill out table of mode zeta values for each file group
        df = pd.DataFrame([
            expName,
            np.mean(zetas[np.where(powerMeans == 1.0)[0]])
            ]).T
        df.to_excel(writer,
                    sheet_name='Summary',
                    startrow=fileGroupNumber + 1,
                    header=False,
                    index=False,
                    )
        
        # Calculate average and standard deviation of all
        # files with non-zero modes in each file group
        includedModes = modes[np.where(modes!=0.0)]
        df = pd.DataFrame([
            np.mean(includedModes),
            np.std(includedModes, ddof=1),
            len(includedModes)
            ]).T
        df.to_excel(writer,
                    sheet_name='Summary',
                    startrow=fileGroupNumber + 1,
                    startcol=2,
                    header=False,
                    index=False,
                    )
        # Fill out table of mode zeta values from each individual file
        df = pd.DataFrame([
            [f"{expName}-{i + 1}" for i in range(nFiles)],
            modes,
            ]).T
        df.to_excel(writer,
                    sheet_name='Summary',
                    startrow=nFileGroups + totalFilesProcessed + 1,
                    header=False,
                    index=False,
                    )
        return nFiles, nBins

def interpolate(xs: np.ndarray, ys: np.ndarray, zetas: np.ndarray):
    nBins = len(zetas)
    powers = np.empty(nBins, dtype=np.float64)
    powerNs = np.empty(nBins, dtype=np.float64)
    minRawZeta = np.min(xs)
    maxRawZeta = np.max(xs)
    for i, zeta in enumerate(zetas):
        if zeta >= maxRawZeta:
            powers[i] = ys[-1]
            powerNs[i] = 0
        elif zeta <= minRawZeta:
            powers[i] = ys[0]
            powerNs[i] = 0
        else:
            # Linear interpolation
            xi = np.argmax(xs >= zeta)
            powers[i] = (ys[xi] - ys[xi - 1]) \
                / (xs[xi] - xs[xi - 1]) \
                * (zeta - xs[xi - 1]) \
                + ys[xi - 1]
        powerNs[i] = np.count_nonzero(
                (xs >= zeta - 0.5 * binSize) \
                & (xs <= zeta + 0.5 * binSize)
                )
    return powers, powerNs

if __name__ == "__main__":
    if len(sys.argv) == 2:
        scriptdir = sys.argv[1].replace('\\','/')
        scriptdir = scriptdir if scriptdir[-1] == '/' else f'{scriptdir}/'
    elif len(sys.argv) == 1:
        # Get directory containing this script
        currentFile = os.path.realpath(__file__).replace('\\','/')
        scriptdir = currentFile[:currentFile.rfind('/')]
    else:
        print('Only one command line argument may be provided: '
              'the full path to the directory containing the'
              f'{extension} files to be processed')

    # Required globals
    # Font size on graphs
    fontsize = 20
    # Set binSize in mV
    binSize = 1
    # Set min and max zeta values to include in the exported file
    minZeta = -90
    maxZeta = 90
    extension = '.csv'
    processFolder(scriptdir)
