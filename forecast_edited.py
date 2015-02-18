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
    
    
def ssaLevel(x):
    if x <= 2.0:
        return '1'
    elif 2.0 < x <= 3.0:
        return '2'
    elif 3.0 < x <= 4.0:
        return '3'
    else:
        return '4'

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
        peaks = peaksOver(stormTide, max(2.0, max(stormTide)))
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
            maxFile.write('Station Name\tStorm Surge (m)\tAstronomical Tide (m)\tSurge + Tide (m)\tEstimated Dates and Times of Peaks\t SSA\n')
#            maxFile.write(surgetotide_dict[matchStation]['VERBOSE_NAME'].replace(',',' '))
            for maxTime, maxStormTide in zip(maxTimes, stormTidePeaks):
                name = (' ').join(surgetotide_dict[matchStation]['VERBOSE_NAME'].split(',')[-2:])
                line = (name.replace(name[-2:], provCode[name[-2:]])+ '\t' +
                        '{:.1f}'.format(surgeValueDict[maxTime] - 0.5) + ' - ' +
                        '{:.1f}'.format(surgeValueDict[maxTime] + 0.5) + '\t' +
                        '{:.1f}'.format(interpolatedTideDict[maxTime]) + '\t' +
                        '{:.1f}'.format(maxStormTide - 0.5) + ' - ' + 
                        '{:.1f}'.format(maxStormTide + 0.5) + '\t' + 
                        (maxTime - datetime.timedelta(minutes=15)).strftime('%m/%d/%Y %H:%M - ') +
                        (maxTime + datetime.timedelta(minutes=15)).strftime('%m/%d/%Y %H:%M') + '\t' +
                        ssaLevel(np.round(maxStormTide + 0.5))
                       )
                maxFile.write(line + '\n')
            maxFile.close()
            firstMax = False
        elif stormTidePeaks != []:
            maxFile = open(os.path.join(JMA1OUTPUTPATH, 'maximum/maxStormTide.csv'), 'a')
#            maxFile.write(surgetotide_dict[matchStation]['VERBOSE_NAME'].replace(',',' '))
            for maxTime, maxStormTide in zip(maxTimes, stormTidePeaks):  
                name = (' ').join(surgetotide_dict[matchStation]['VERBOSE_NAME'].split(',')[-2:])
                line = (name.replace(name[-2:], provCode[name[-2:]])+ '\t' +
                        '{:.1f}'.format(surgeValueDict[maxTime] - 0.5) + ' - ' + 
                        '{:.1f}'.format(surgeValueDict[maxTime] + 0.5) + '\t' + 
                        '{:.1f}'.format(interpolatedTideDict[maxTime]) + '\t' +
                        '{:.1f}'.format(maxStormTide - 0.5) + ' - ' + 
                        '{:.1f}'.format(maxStormTide + 0.5) + '\t' + 
                        (maxTime - datetime.timedelta(minutes=15)).strftime('%m/%d/%Y %H:%M - ') +
                        (maxTime + datetime.timedelta(minutes=15)).strftime('%m/%d/%Y %H:%M') + '\t' +
                         ssaLevel(np.round(maxStormTide + 0.5))
                       )
                maxFile.write(line + '\n')
            maxFile.close()
        else:
            continue
provCode = {'DO': 'Davao Oriental', 'BG': 'Benguet', 'BA': 'Bataan', 'DI': 'Dinagat Islands', 'TR': 'Tarlac', 'BN': 'Batanes', 'BO': 'Bohol', 'BI': 'Biliran', 'BK': 'Bukidnon', 'BT': 'Batangas', 'BU': 'Bulacan', 'SS': 'Surigao del Sur',
            'BS': 'Basilan', 'DV': 'Davao del Norte', 'ST': 'Surigao del Norte', 'DR': 'Davao del Sur', 'LG': 'Laguna', 'IB': 'Isabela', 'MT': 'Mountain', 'LN': 'Lanao del Norte', 'TT': 'Tawi-Tawi',
            'NC': 'Cotabato', 'ND': 'Negros Occidental', 'NE': 'Nueva Ecija', 'ZS': 'Zamboanga del Sur', 'QZ': 'Quezon', 'SK': 'Sultan Kudarat', 'LU': 'La Union', 'ZM': 'Zambales', 'LS': 'Lanao del Sur', 'NR': 'Negros Oriental', 'IF': 'Ifugao',
            'PN': 'Pangasinan', 'IS': 'Ilocos Sur', 'NV': 'Nueva Vizcaya', 'PM': 'Pampanga', 'RO': 'Romblon', 'GU': 'Guimaras', 'AB': 'Abra', 'QR': 'Quirino', 'CN': 'Camarines Norte',
            'CM': 'Camiguin', 'CL': 'Compostela Valley', 'ZN': 'Zamboanga del Norte', 'CB': 'Cebu', 'AK': 'Aklan', 'SL': 'Southern Leyte', 'CG': 'Cagayan', 'AL': 'Albay', 'AN': 'Agusan del Norte',
            'AQ': 'Antique', 'AP': 'Apayao', 'AS': 'Agusan del Sur', 'LE': 'Leyte', 'AU': 'Aurora', 'ZY': 'Zamboanga-Sibugay', 'IN': 'Ilocos Norte', 'CS': 'Camarines Sur', 'MD': 'Misamis Occidental', 'CP': 'Capiz', 'CV': 'Cavite', 'ES': 'Eastern Samar',
            'CT': 'Catanduanes', 'II': 'Iloilo', 'KA': 'Kalinga', 'NS': 'Northern Samar', 'MC': 'Occidental Mindoro', 'MB': 'Masbate', 'DC': 'Davao Occidental', 'SR': 'Sorsogon', 'SQ': 'Siquijor',
            'MN': 'Misamis Oriental', 'SU': 'Sulu', 'RI': 'Rizal', 'MG': 'Maguindanao', 'MM': 'Metropolitan Manila', 'MQ': 'Marinduque', 'PW': 'Palawan', 'SM': 'Samar', 'MR': 'Oriental Mindoro', 'SC': 'South Cotabato', 'SG': 'Sarangani'}       
            
            
def rankOutput(inputPath, outputPath):
    with open(os.path.join(inputPath, 'maximum', 'maxStormTide.csv'), 'r') as csv:
        data = csv.readlines()
        header = data[0]
        data = data[1:]
        values = []
        for i in xrange(len(data)):
            line = data[i].strip().split('\t')
            values.append((line[0], line[1], line[2], line[3], line[4], line[5]))
        dtype = [('name', 'S50'), ('surge', 'S50'), ('tide', 'S50'), ('stormtide', 'S50'), ('time', 'S50'), ('ssa', 'S50')]
        array = np.array(values, dtype=dtype)
        array = np.sort(array, order='stormtide')[::-1]
    with open(os.path.join(inputPath, 'maximum', 'maxStormTide_ranked.csv'), 'w') as csv:
        csv.write(header)
        written = []
        for i in array:
            if i[0] in written:
                continue
            else:
                csv.write(('\t').join(i))
                csv.write('\n')
                written.append(i[0])
        
        
if __name__ == '__main__':
    tidePath = r'C:\Users\user\Dropbox\storm surge\Tides\Tides_January2015'
    inputPath = r'C:\Users\user\Google Drive\SS_Forecast\Amang_Mekkhala\20150116\21UTC\TIMESERIES'
    outputPath = r'C:\Users\user\Google Drive\SS_Forecast\Amang_Mekkhala\20150116\21UTC\TIMESERIES'
    
    forecast(tidePath, inputPath, outputPath)
    rankOutput(inputPath, outputPath)