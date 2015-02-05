import numpy as np
import matplotlib.mlab as ml
import os
import string

def alphanum(s):
    letter_set = frozenset(string.ascii_lowercase+string.ascii_uppercase+string.digits)
    tab = string.maketrans(string.ascii_lowercase + string.ascii_uppercase,
                       string.ascii_lowercase * 2)
    deletions = ''.join(ch for ch in map(chr,range(256)) if ch not in letter_set)
    return string.translate(s, tab, deletions)

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def xyz2grd(x,y,z,quality=100):
    #~ X,Y=np.meshgrid(x,y)
    xi=np.linspace(min(x), max(x), quality)
    yi=np.linspace(min(y), max(y), quality)
    zi=ml.griddata(x,y,z,xi,yi)
    return xi,yi,zi

def frange(start,stop,step):
    l=[]
    while start<=stop:
        l.append(round(start,2))
        start+=step
    return l

def getTerminalSize():
    env = os.environ
    def ioctl_GWINSZ(fd):
        try:
            import fcntl, termios, struct, os
            wh = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
        except:
            return
        return wh
    wh = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not wh:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            wh = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not wh:
        wh = (env.get('LINES', 25), env.get('COLUMNS', 80))
        ### Use get(key[, default]) instead of a try/catch
        #try:
        #    wh = (env['LINES'], env['COLUMNS'])
        #except:
        #    wh = (25, 80)
    return int(wh[1]), int(wh[0])

def printInABox(string,width=getTerminalSize()[0],padded_top=True,padded_bottom=True):
    if padded_top:
        print ''
    print '*'*width
    print '*'+' '*(width-2)+'*'
    word_list=string.split()
    lines=[]; line=[]
    strlen=0; count=1
    for item in word_list:
        if strlen+(len(item)+1)<(width-2):
            line.append(item)
            strlen+=(len(item)+1)
            if count==len(word_list):
                lines.append(line)
            count+=1
        else:
            strlen=0
            lines.append(line)
            line=[]
            line.append(item)
            strlen+=(len(item)+1)
            count+=1
    for i in range(len(lines)):
        lines[i]=' '.join(lines[i])
        odd=False
        if (len(lines[i])+width)%2!=0:
            odd=True
        if odd:
            print '*'+' '*((width-2-len(lines[i]))/2)+lines[i]+' '*(1+((width-2-len(lines[i]))/2))+'*'
        else:
            print '*'+' '*((width-2-len(lines[i]))/2)+lines[i]+' '*((width-2-len(lines[i]))/2)+'*'
    print '*'+' '*(width-2)+'*'
    print '*'*width
    if padded_bottom:
        print ''

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def movingAverage(timeseries, radius):
    newTimeseries = []
    for i in range(len(timeseries)):
        if i-radius < 0:
            start = 0
        else:
            start = i-radius
        end = i+radius+1
        segment = np.array(timeseries[start:end])
        newTimeseries.append(segment.mean())
    return np.array(newTimeseries)
