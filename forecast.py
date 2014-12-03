# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 14:08:30 2014
Forecasting using JMA 1.0
@author: phillip
"""
import os
import glob
import datetime
import numpy as np
from interpolate import alphanum
from interpolate import surgetotide_dict
from interpolate import interpolator_generator


def readXTide(tidePath):
    """
    Get tide values from Wxtide.
    One value coressponds to an entry on the list.
    Times are converted to PST
    entry format : [value(m), HH:MM, time_format, MM/DD/YYYY]
    
    tideStationDict={tide_station1:{time1:value1,
                                    time2:value2,
                                    time3:value3
                                    ...},
                     tide_station2:{time1:value1,
                                    time2:value2,
                                    time3:value3
                                    ...},
                     ...}
    """
    print 'Reading tide data...\n'
    xtideFileList = glob.glob(os.path.join(tidePath, '*.txt'))
    tideStationDict = {}
    for filename in sorted(xtideFileList):
        print 'Processing', os.path.basename(filename)
        with open(filename, 'r') as tideFile:
            tideData = tideFile.readlines()
        tideStation = alphanum(tideData[0])
        tideData = tideData[3:]
        cleanedTideData = []
        for i in range(len(tideData)):
            tideData[i] = tideData[i].replace('\r\n', '').split()
            if tideData[i]:
                cleanedTideData.append(tideData[i])
        for i in range(len(cleanedTideData)):
            if len(cleanedTideData[i]) < 4:
                cleanedTideData[i].append(cleanedTideData[i - 1][3])
        tideValueDict = {}
        for i in range(len(cleanedTideData)):
            if len(cleanedTideData[i][1]) < 5:
                cleanedTideData[i][1] = '0' + cleanedTideData[i][1]
            dataTime = datetime.datetime.strptime(cleanedTideData[i][3] +
                       cleanedTideData[i][1], '%Y-%m-%d%H:%M') + (
                       datetime.timedelta(hours=8))
            tideValueDict[dataTime] = float(cleanedTideData[i][0])
        tideValues = [tideValueDict[t] for t in tideValueDict.keys()]
        meanTide = np.mean(tideValues)
        for t in tideValueDict.keys():
            tideValueDict[t] = tideValueDict[t] - meanTide
        tideStationDict[tideStation] = tideValueDict
    print '\nDone reading tide data...\n'
    return tideStationDict


def readJmaOutput(jmapath):
    """
    Get surge values from JMA output
    One value corresponds to an entry on the list.
    entry format : [YYYY/MM/DD, HH:MM, value(cm)]
    surgeStationDict={surge_station1:{time1:value1,
                                      time2:value2,
                                      time3:value3
                                      ...},
                      surge_station2:{time1:value1,
                                      time2:value2,
                                      time3:value3
                                      ...},
                      ...}
    """
    jmaFileList = glob.glob(os.path.join(jmapath, '[0-9]*.txt'))
    try:
        os.mkdir(os.path.join(jmapath, 'total'))
        os.mkdir(os.path.join(jmapath, 'maximum'))
    except OSError:
        pass

    surgeStationDict = {}
    for filename in sorted(jmaFileList):
        print filename
        with open(filename, 'r') as surgeFile:
            surgeData = surgeFile.readlines()
        header = surgeData[0]
        surgeStation = alphanum(header[21:59])
        surgeData = surgeData[1:]
        surgeData = [data.replace('\n', '').split()[1:4] for data in surgeData]
        surgeValueDict = {}
        for line in surgeData:
            try:
                dataTime = datetime.datetime.strptime(line[0] +
                           line[1], '%Y/%m/%d%H:%M') + (
                           datetime.timedelta(hours=8))
            except:
                continue
            try:
                surgeValueDict[dataTime] = float(line[2]) / 100
            except:
                continue
        surgeStationDict[surgeStation] = surgeValueDict
    return surgeStationDict
    
    
def peaksOver(ts, over):
    peaks = [(i, ts[i]) for i in range(1, len(ts) - 1) if ts[i - 1] <= ts[i] >=
             ts[i + 1] and ts[i] >= over]
    return peaks      
    

def forecast(tidePath, inputPath, outputPath):
    TIDEPATH = tidePath
    JMA1INPUTPATH = inputPath
    JMA1OUTPUTPATH = outputPath
    
    tideStationDict = readXTide(TIDEPATH)
    surgeStationDict = readJmaOutput(JMA1INPUTPATH)

    firstMax = True
    firstlist = True
    for station in sorted(surgeStationDict.keys()):
        surgeValueDict = surgeStationDict[station]
#        ## Temp add
#        surgeTimeList, surgeValueList = zip(*surgeValueDict.items())
#        maxSurge = max(surgeValueList)
#        maxTime = surgeTimeList[surgeValueList.index(maxSurge)]
#        ##

        for longStation in surgetotide_dict.keys():
            if station in longStation:
                matchStation = longStation
        vertex1 = alphanum(surgetotide_dict[matchStation]['V1'])
        vertex2 = alphanum(surgetotide_dict[matchStation]['V2'])
        vertex3 = alphanum(surgetotide_dict[matchStation]['V3'])
        tideValue1Dict = tideStationDict[vertex1]
        tideValue2Dict = tideStationDict[vertex2]
        tideValue3Dict = tideStationDict[vertex3]
        interpolate = interpolator_generator(
                          surgetotide_dict[matchStation]['V1_LONG'],
                          surgetotide_dict[matchStation]['V1_LAT'],
                          surgetotide_dict[matchStation]['V2_LONG'],
                          surgetotide_dict[matchStation]['V2_LAT'],
                          surgetotide_dict[matchStation]['V3_LONG'],
                          surgetotide_dict[matchStation]['V3_LAT'],
                          surgetotide_dict[matchStation]['LONG'],
                          surgetotide_dict[matchStation]['LAT']
                      )
    
        interpolatedTideDict = {}
        for dataTime in tideValue1Dict.keys():
            interpolatedTideDict[dataTime] = interpolate(
                                           tideValue1Dict[dataTime],
                                           tideValue2Dict[dataTime],
                                           tideValue3Dict[dataTime]
                                           )
    
        interpolatedTideTimeList, interpolatedTideValueList = zip(*interpolatedTideDict.items())
    
        surgeTimeList, surgeValueList = zip(*surgeValueDict.items())
        """
        # Uncomment this section if you want to align the maximums
        maxInterpolatedTide = max(interpolatedTideValueList)
        maxSurge = max(surgeValueList)
        maxTideTime = interpolatedTideTimeList[interpolatedTideValueList.index(maxInterpolatedTide)]
        maxSurgeTime = surgeTimeList[surgeValueList.index(maxSurge)]
        newsurgevalue_dict={}
        for dataTime in surgeValueDict.keys():
            newsurgevalue_dict[dataTime-maxSurgeTime+maxTideTime]=surgeValueDict[dataTime]
        del surgeValueDict
        surgeValueDict=newsurgevalue_dict
        del newsurgevalue_dict
        """
        commonTime = list(set(interpolatedTideDict.keys())&set(surgeValueDict.keys()))
        commonTime.sort()
        print matchStation
    
        filename = matchStation + '.csv'
        first = True
    
        stormTide = [interpolatedTideDict[dataTime] + surgeValueDict[dataTime]
                     for dataTime in commonTime]
        peaks = peaksOver(stormTide, 2.0)
        indeces = [i[0] for i in peaks]
        stormTidePeaks = [i[1] for i in peaks]
        maxTimes = [commonTime[i] for i in indeces]
        for dataTime, stormTideValue in zip(commonTime, stormTide):
            dataLine = dataTime.strftime('%m-%d-%Y %H:%M ') + '{:.2f}'.format(stormTideValue)
            if first:
                outputFile = open(os.path.join(JMA1OUTPUTPATH, 'total', #dataTime.strftime('%y%m%d_') + 
                                  filename), 'w')
                outputFile.write(dataLine + '\n')
                outputFile.close()
                first = False
            else:
                outputFile = open(os.path.join(JMA1OUTPUTPATH, 'total', #dataTime.strftime('%y%m%d_') + 
                                  filename), 'a')
                outputFile.write(dataLine + '\n')
                outputFile.close()
        
        if firstlist:
            surgelist = open(os.path.join(JMA1OUTPUTPATH, 'list.txt'), 'w')
            surgelist.write(matchStation + '\n')
            surgelist.close()
            firstlist = False
        else:
            surgelist = open(os.path.join(JMA1OUTPUTPATH, 'list.txt'), 'a')
            surgelist.write(matchStation + '\n')
            surgelist.close()
        
        if firstMax and stormTidePeaks != []:
            maxFile = open(os.path.join(JMA1OUTPUTPATH, 'maximum/maxStormTide.csv'), 'w')
            maxFile.write('Station Name\tStorm Surge (m)\tAstronomical Tide (m)\tSurge + Tide (m)\tEstimated Dates and Times of Peaks\n')
            maxFile.write(surgetotide_dict[matchStation]['VERBOSE_NAME'].replace(',',' '))
            for maxTime, maxStormTide in zip(maxTimes, stormTidePeaks):  
                line = ('\t' +
                        '{:.2f}'.format(surgeValueDict[maxTime] - 0.5) + ' - ' +
                        '{:.2f}'.format(surgeValueDict[maxTime] + 0.5) + '\t' +
                        '{:.2f}'.format(interpolatedTideDict[maxTime]) + '\t' +
                        '{:.2f}'.format(maxStormTide - 0.5) + ' - ' + 
                        '{:.2f}'.format(maxStormTide + 0.5) + '\t' + 
                        (maxTime - datetime.timedelta(minutes=15)).strftime('%m/%d/%Y %H:%M - ') +
                        (maxTime + datetime.timedelta(minutes=15)).strftime('%m/%d/%Y %H:%M')
                       )
                maxFile.write(line + '\n')
            maxFile.close()
            firstMax = False
        elif stormTidePeaks != []:
            maxFile = open(os.path.join(JMA1OUTPUTPATH, 'maximum/maxStormTide.csv'), 'a')
            maxFile.write(surgetotide_dict[matchStation]['VERBOSE_NAME'].replace(',',' '))
            for maxTime, maxStormTide in zip(maxTimes, stormTidePeaks):            
                line = ('\t' +
                        '{:.2f}'.format(surgeValueDict[maxTime] - 0.5) + ' - ' + 
                        '{:.2f}'.format(surgeValueDict[maxTime] + 0.5) + '\t' + 
                        '{:.2f}'.format(interpolatedTideDict[maxTime]) + '\t' +
                        '{:.2f}'.format(maxStormTide - 0.5) + ' - ' + 
                        '{:.2f}'.format(maxStormTide + 0.5) + '\t' + 
                        (maxTime - datetime.timedelta(minutes=15)).strftime('%m/%d/%Y %H:%M - ') +
                        (maxTime + datetime.timedelta(minutes=15)).strftime('%m/%d/%Y %H:%M')
                       )
                maxFile.write(line + '\n')
            maxFile.close()
        else:
            continue
                
            
if __name__ == '__main__':
    tidePath = r'C:\Users\user\Dropbox\storm surge\Tides\Tides_September2013'
    
    inputPath = r'C:\Users\user\Documents\timeseries\timeseries'
    outputPath = r'C:\Users\user\Documents\timeseries\timeseries'
    
    forecast(tidePath, inputPath, outputPath)