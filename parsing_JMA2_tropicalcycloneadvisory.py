#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Irene Crisologo, John Phillip Lapidez
#
# Created:     04/02/2014
# Copyright:   (c) NOAH Storm Surge
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import pandas as pd
import numpy as np
import datetime
import os

def check_latlon(latlon):
    """
    Checks the direction of latitude or longitude. If latitude/longitude is N/E,
    the function returns a positive value. If latitude/longitude is S/W, the
    function returns a negative value.
    """
    if latlon.endswith('N'):
        return float(latlon[:len(latlon)-1])
    elif latlon.endswith('S'):
        return -1.*float(latlon[:len(latlon)-1])
    elif latlon.endswith('E'):
        return float(latlon[:len(latlon)-1])
    elif latlon.endswith('W'):
        return -1.*float(latlon[:len(latlon)-1])
    else:
        return float(latlon)

def date_to_string(date, fmt='%y%m%d%H'):
    """
    Converts datetime object to string with the given format.
    """
    return datetime.datetime.strftime(date, fmt)

def string_to_date(s, fmt='%y%m%d%H'):
    """
    Converts string to datetime object with the given format.
    """
    return datetime.datetime.strptime(s, fmt)

def typhoon_class(mxwd, units='kt'):
    """
    Returns JMA typhoon class based on maximum sustained windspeed.
    Options for the units of mxwd are 'kt' for knots, 'mps' for m/s,
    'kph' for km/h.
    """
    if mxwd < 0:
        raise Exception('mxwd must be a positive value.')
    valrange = {'kt': [(0, 33), (34, 47), (48, 63)],
                'mps': [(0, 17), (18, 24), (25, 32)],
                'kph': [(0, 62), (63, 88), (89, 118)]
                }

    for i, (low, high) in enumerate(valrange[units]):
        if low <= mxwd <= high:
            return i + 2
    return i + 3


if __name__ == '__main__':

    # set path where files will be saved
    path = r'C:\Users\user\Google Drive\SS_Forecast\Ruby_Hagupit\downloads'

    # filename to contain the latest tropical advisory from JMA
    filename = 'latestTropicalCycloneAdvisory.txt'

    # gets date today                 adjust to UTC                 temporary (for testing)
    dtnow = datetime.datetime.now() - datetime.timedelta(hours=8)# - datetime.timedelta(days=3)

    # read file
    f = open(os.path.join(os.path.join(path, filename)))
    lines = f.readlines()
    f.close()
    lines = [l.strip() for l in lines]

    # set-up dictionary to contain advisory data
    data = {}

    # set up header information
    advisorytime = lines[0].split(' ')[2]
    data['Time'] = datetime.datetime(dtnow.year, dtnow.month, int(advisorytime[0:2]), int(advisorytime[2:4]), int(advisorytime[4:6]))
    data['Type'] = lines[1]
    cyclone_name = lines[2].split('  ')[1]
    try:
        cyclone_name = cyclone_name.split('UPGRADED FROM')[0].strip()
    except:
        pass
    data['Name'] = cyclone_name

    # separate file contents to analysis and forecast
    idx_analysis = lines.index('ANALYSIS')
    idx_forecast = lines.index('FORECAST')
    list_analysis = lines[idx_analysis : idx_forecast]
    list_forecast = lines[idx_forecast + 1:]
    
    # create dictionary for analysis data
    analysis = {}
    anal_pstn = list_analysis[1].split('  ')

    analysis['Latitude'] = check_latlon(anal_pstn[1].split(' ')[1])
    analysis['Longitude'] = check_latlon(anal_pstn[1].split(' ')[2])

    for line in list_analysis[1:]:
        line_ = line.split('  ')
        if line_[0].endswith('KT'):
            analysis['GUST'] = analysis['GUST'] + ' ' + line_[0] + ' ' + line_[1]
        else:
            analysis[line_[0]] = line_[1]

    data['Analysis'] = analysis
    # create dictionary for forecast data

    # find forecast times and positions
    forecast_times = []
    for i in xrange(len(lines)):
        if 'HF' in lines[i] and not lines[i].endswith('='):
            forecast_times.append(lines[i].split()[0])

    forecast = {}
    for line in list_forecast:
        # if line starts with 'XXHF' but ends with '=', then line is ignored (because it does not contain forecast info)
        try:        
            if line.split('  ')[0].endswith('HF') and not line.endswith('='):
                current_hf = line.split('  ')[0]
            else:
                if line.split('  ')[0].endswith('KT'):
                    forecast[current_hf]['GUST'] = forecast[current_hf]['GUST'] + ' ' + line.split('  ')[0] + ' ' + line.split('  ')[1]
                else:
                    forecast[current_hf][line.split('  ')[0]] = line.split('  ')[1]
            if line.startswith(current_hf) and not line.endswith('='):
                forecast[current_hf] = {}
                forecast[current_hf]['Latitude'] = check_latlon(line.split('  ')[1].split(' ')[1])
                forecast[current_hf]['Longitude'] = check_latlon(line.split('  ')[1].split(' ')[2])
                fcast_timedelta = datetime.timedelta(hours=int(current_hf.strip('HF')))
                forecast[current_hf]['Time'] =  data['Time'] + fcast_timedelta
                forecast[current_hf]['70radius'] = float(line.split()[-2].strip('NM'))
        except IndexError:
            if line.split()[0].endswith('HF') and not line.endswith('='):
                current_hf = line.split()[0]
            else:
                if line.split('  ')[0].endswith('KT'):
                    forecast[current_hf]['GUST'] = forecast[current_hf]['GUST'] + ' ' + line.split('  ')[0] + ' ' + line.split('  ')[1]
                else:
                    forecast[current_hf][line.split('  ')[0]] = line.split('  ')[1]
            if line.startswith(current_hf) and not line.endswith('='):
                forecast[current_hf] = {}
                forecast[current_hf]['Latitude'] = check_latlon(line.split()[2])
                forecast[current_hf]['Longitude'] = check_latlon(line.split()[2])
                fcast_timedelta = datetime.timedelta(hours=int(current_hf.strip('HF')))
                forecast[current_hf]['Time'] =  data['Time'] + fcast_timedelta
                forecast[current_hf]['70radius'] = float(line.split()[-2].strip('NM'))

    copy_forecast_times = [t for t in forecast_times]
    
    for hf in copy_forecast_times:
        if 'PRES' not in forecast[hf].keys():
            forecast_times.remove(hf)
            del forecast[hf]
    data['Forecast'] = forecast
    data['Forecast Times'] = forecast_times
    # check if outdata.txt exists
    # outdata.txt contains the dictionary of the parsed data
    outdata_file = os.path.join(path,'outdata.txt')
    if os.path.exists(outdata_file):
        # if file exists, read it
        f = open(outdata_file)
        d = f.read()
        f.close()
        outdata = eval(d)
    else:
        # if file does not exists, make it
        outdata = {}
        outdata['Name'] = ''
        outdata['Start'] = data['Time']
        outdata['End'] = data['Time']
        # setup columns
        outdata['Dates'] = []
        outdata['Pcenter'] = []
        outdata['lon'] = []
        outdata['lat'] = []
        outdata['r0'] = []
        outdata['Coef'] = []
        outdata['Pfar'] = []
        outdata['Mxwd'] = []

    # update the outdata dictionary

    # if Tropical depression, keep name as TD (this avoids including other notes such as "DOWNGRADED TO ..."
    #   from the name
    if data['Name'].startswith('TD'):
        outdata['Name'] = 'TD'
    else:
        outdata['Name'] = data['Name']

    # update End time value to the last forecast value
    outdata_end = data['Forecast'][data['Forecast Times'][-1]]['Time']
    if outdata['End'] < outdata_end:
        outdata['End'] = outdata_end
    else:
        pass

    if outdata['Dates'] != []:
        for idx_to_update, date in enumerate(outdata['Dates']):
            if data['Time'] <= string_to_date(date):
                break
        
        outdata['Dates'] = outdata['Dates'][:idx_to_update]
        outdata['lat'] = outdata['lat'][:idx_to_update]
        outdata['lon'] = outdata['lon'][:idx_to_update]
        outdata['Mxwd'] = outdata['Mxwd'][:idx_to_update]
        outdata['Pfar'] = outdata['Pfar'][:idx_to_update]
        outdata['Pcenter'] = outdata['Pcenter'][:idx_to_update]
        outdata['r0'] = outdata['r0'][:idx_to_update]
        outdata['Coef'] = outdata['Coef'][:idx_to_update]
    
    # Append Analysis Data
    outdata['Dates'].append(date_to_string(data['Time']))
    outdata['Pcenter'].append(float(data['Analysis']['PRES'].strip('HPA')))
    outdata['Mxwd'].append(float(data['Analysis']['MXWD'].strip('KT')))
    outdata['lon'].append(data['Analysis']['Longitude'])
    outdata['lat'].append(data['Analysis']['Latitude'])
    outdata['r0'].append(60.0)
    outdata['Coef'].append(0.70)
    outdata['Pfar'].append(1012.00)
    
    outdata['70radius'] = [0 for d in outdata['Dates']]
    
    # Append Forecast data
    for ftime in data['Forecast Times']:
        ftime_dict = data['Forecast'][ftime]
        outdata['Dates'].append(date_to_string(ftime_dict['Time']))
        outdata['Pcenter'].append(float(ftime_dict['PRES'].strip('HPA')))
        outdata['Mxwd'].append(float(ftime_dict['MXWD'].strip('KT')))
        outdata['lat'].append(ftime_dict['Latitude'])
        outdata['lon'].append(ftime_dict['Longitude'])
        outdata['r0'].append(60)
        outdata['Coef'].append(0.70)
        outdata['Pfar'].append(1012.00)
        outdata['70radius'].append(ftime_dict['70radius'])
        

    # sort data according to Dates
    array = zip(outdata['Dates'],
                outdata['Pcenter'],
                outdata['Mxwd'],
                outdata['lon'],
                outdata['lat'],
                outdata['r0'],
                outdata['Coef'],
                outdata['Pfar'],
                outdata['70radius'])
    sorted_array = sorted(array)
    Dates_sorted = [point[0] for point in sorted_array]
    Pcenter_sorted = [point[1] for point in sorted_array]
    Mxwd_sorted = [point[2] for point in sorted_array]
    lon_sorted = [point[3] for point in sorted_array]
    lat_sorted = [point[4] for point in sorted_array]
    r0_sorted = [point[5] for point in sorted_array]
    Coef_sorted = [point[6] for point in sorted_array]
    Pfar_sorted = [point[7] for point in sorted_array]
    Probradius_sorted = [point[8] for point in sorted_array]

    # interpolate values to 6 hour interval
    dateObjects = [string_to_date(d) for d in Dates_sorted]
    dataArray = np.array([Pcenter_sorted, Mxwd_sorted, lon_sorted, lat_sorted,
                          r0_sorted, Coef_sorted, Pfar_sorted, Probradius_sorted]).transpose()
    df = pd.DataFrame(dataArray, index=dateObjects,
                      columns=['Pcenter', 'Mxwd', 'lon', 'lat', 'r0', 'Coef',
                               'Pfar', '70radius'])
    rng = pd.date_range(dateObjects[0], dateObjects[-1], freq='6H')
    if rng[-1] != dateObjects[-1]:
        rng = pd.date_range(dateObjects[0], dateObjects[-1], freq='6H', closed='left')
    dummyseries = pd.DataFrame([np.NaN] * len(rng), index=rng)

    df = df.join(dummyseries, how='outer')
#    print df
    df = df.interpolate(method='time')
    df = np.round(df, 2)
    df = df.reindex(rng)
    print df
    
    # update value lists
    dtrng = [d for d in pd.to_datetime(rng)]
    strrng = [date_to_string(d) for d in dtrng]
    outdata['Dates'] = strrng
    Dates_sorted = strrng
    Pcenter_sorted = df['Pcenter'].values
    Mxwd_sorted = df['Mxwd'].values
    lon_sorted = df['lon'].values
    lat_sorted = df['lat'].values
    r0_sorted = df['r0'].values
    Coef_sorted = df['Coef'].values
    Pfar_sorted = df['Pfar'].values

    # other info, compute according to new time range
    hours_calc = dtrng[-1] - dtrng[0]
    hours_calc = hours_calc.total_seconds()/3600.
    outdata['Hours Calculation'] = hours_calc

    # write to outdata.txt
    strdata = str(outdata)
    f = open(os.path.join(path,'outdata.txt'),'w')
    f.write(strdata)
    f.close()

    # write to text file
    f = open(os.path.join(path,outdata['Name'].replace(' ','-')+'.txt'),'w')
    f.write('Cyclone '+outdata['Name']+'\t')
    f.write('Start= '+ strrng[0] + ' ')
    f.write('End= ' + strrng[-1] + '\n')
    f.write(str(int(outdata['Hours Calculation'])) + ' hours calculation\n')
    f.write('6 hourly data\n')
    f.write('date' + '\t' + 'Pcenter' + '\t' + 'lon' + '\t' + 'lat' + '\t' + 'r0\tCoef\tPfar\n')

    for i in xrange(len(outdata['Dates'])):
        f.write(' ' + str(Dates_sorted[i]) +
                ' ' * (11 - len('{:.2f}'.format(Pcenter_sorted[i]))) +
                '{:.2f}'.format(Pcenter_sorted[i]) +
                ' ' * (10 - len('{:.2f}'.format(lon_sorted[i]))) +
                '{:.2f}'.format(lon_sorted[i]) +
                ' ' * (9 - len('{:.2f}'.format(lat_sorted[i]))) +
                '{:.2f}'.format(lat_sorted[i]) +
                ' ' * (10 - len('{:.2f}'.format(r0_sorted[i]))) +
                '{:.2f}'.format(r0_sorted[i]) +
                ' ' * (8 - len('{:.2f}'.format(Coef_sorted[i]))) +
                '{:.2f}'.format(Coef_sorted[i]) +
                ' ' * (11 - len('{:.2f}'.format(Pfar_sorted[i]))) +
                '{:.2f}'.format(Pfar_sorted[i]) + '\n')
    f.close()

    # write to text file (JMA1.0 best track)
    year = datetime.datetime.now().year
    f = open(os.path.join(path, 'bst' + str(year) + '.txt'), 'w')
    f.write('66666 ' + str(year)[-2:] + '01   ' + str(len(outdata['Dates'])) +
            '      ' + str(year)[-2:] + '01 0 6 ' +
            outdata['Name'].split()[2] + '\n'  )
    for i in xrange(len(outdata['Dates'])):
        f.write(str(Dates_sorted[i]) + ' 002 ' +
                '{:} '.format(typhoon_class(Mxwd_sorted[i])) +
                '{:03d} '.format(int(lat_sorted[i] * 10)) +
                '{:03d} '.format(int(lon_sorted[i] * 10)) +
                ' ' * (4 - len('{:}'.format(int(Pcenter_sorted[i])))) +
                '{:}'.format(int(Pcenter_sorted[i])) +
                ' ' * (8 - len('{:03d}'.format(int(Mxwd_sorted[i])))) +
                '{:03d}'.format(int(Mxwd_sorted[i])) + '\n')
    f.close()
