#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      irene
#
# Created:     07/02/2014
# Copyright:   (c) irene 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import urllib2
import datetime
import os
from bs4 import BeautifulSoup


def get_filenames(url):
    """
    Gets filenames from the given url using BeautifulSoup
    Parameters
    ----------
    url : string
        URL containing filenames
    """
    open_url = urllib2.urlopen(url)
    soup = BeautifulSoup(open_url)
    fnamelist = []
    for fname in soup.find_all('a'):
        fnamelist.append(str(fname.get('href')))
    fnamelist.sort(cmp)
    fnamelist = fnamelist[1:]
    return fnamelist

def get_timefolder(dtnow):
    if dtnow > datetime.datetime(dtnow.year, dtnow.month, dtnow.day, 21): # 21:00
        timefolder = '210000/'
    elif dtnow > datetime.datetime(dtnow.year, dtnow.month, dtnow.day, 18): # 18:00
        timefolder = '180000/'
    elif dtnow > datetime.datetime(dtnow.year, dtnow.month, dtnow.day, 15): # 15:00
        timefolder = '150000/'
    elif dtnow > datetime.datetime(dtnow.year, dtnow.month, dtnow.day, 12): # 12:00
        timefolder = '120000/'
    elif dtnow > datetime.datetime(dtnow.year, dtnow.month, dtnow.day, 9): # 09:00
        timefolder = '090000/'
    elif dtnow > datetime.datetime(dtnow.year, dtnow.month, dtnow.day, 6): # 06:00
        timefolder = '060000/'
    elif dtnow > datetime.datetime(dtnow.year, dtnow.month, dtnow.day, 3): # 03:00
        timefolder = '030000/'
    elif dtnow > datetime.datetime(dtnow.year, dtnow.month, dtnow.day, 0): # 00:00
        timefolder = '000000/'
    return timefolder

if __name__ == '__main__':


    URLbase = 'http://www.wis-jma.go.jp/d/o/RJTD/Alphanumeric/Warning/Tropical_cyclone/'
    path = r'C:\Users\Windows User\Dropbox\storm surge\Ruby_Hagupit\downloads DEC4 8PM'
    try:
        os.mkdir(path)
    except WindowsError:
        pass


    # gets date today                 adjust to UTC                  temporary (for testing)
    dtnow = datetime.datetime.now() - datetime.timedelta(hours=8)# - datetime.timedelta(hours=3) 

    datefolder = datetime.datetime.strftime(dtnow, '%Y%m%d')
    URLdate= URLbase + datefolder + '/' + get_timefolder(dtnow)
    try:
        filelist = get_filenames(URLdate)
    except:
        print 'Tropical Cyclone Advisory for ' + datefolder+'/' + get_timefolder(dtnow) + ' is not available as of now'
    else:
        
        ta_fname = [f for f in sorted(filelist) if 'WTPQ' in f][-1]
        URLfile = URLdate + ta_fname
        print URLfile
        datfile = urllib2.urlopen(URLfile)
        output = open(os.path.join(path,'latestTropicalCycloneAdvisory.txt'),'w')
        output.write(datfile.read())
        output.close()
