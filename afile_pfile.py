# coding: utf-8
import re
import math
import csv
import sys
import os
import struct
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from obspy import UTCDateTime, Stream, Trace
import argparse

def unpackAfile(infile):

    # == opening Afile ==
    b = os.path.getsize(infile)
    FH = open(infile, 'rb')
    line = FH.read(b)
    fileHeader = struct.unpack("<4s3h6bh6s", line[0:24])

    fileLength = fileHeader[3]
    port = fileHeader[10]
    # FirstStn = fileHeader[11][0:4].decode('ASCII').rstrip()
    # =================================Header=================================

    portHeader = []
    for i in range(24, port * 32, 32):
        port_data = struct.unpack("<4s4s3sbh2b4s12b", line[i:i+32])
        portHeader.append(port_data)

    # =================================Data===================================

    dataStartByte = 24+int(port)*32
    dataPoint = 3*int(port)*int(fileLength)*100
    times = int(port)*3*4
    data = []

    data = struct.unpack("<%di" % dataPoint, line[dataStartByte:dataStartByte + dataPoint*4])


    portHeader = np.array(portHeader)
    data = np.array(data)
    idata = data.reshape((3,port,fileLength*100),order='F')

    #== write to obspy Stream --
    sttime = UTCDateTime(fileHeader[1], fileHeader[4], fileHeader[5], fileHeader[6], fileHeader[7], fileHeader[8], fileHeader[2])
    npts = fileHeader[3]*fileHeader[9]
    samp = fileHeader[9]
    afst = Stream()
    
    for stc in range(fileHeader[10]):
        stn = portHeader[stc][0].decode('ASCII').rstrip()
        instrument = portHeader[stc][1].decode('ASCII').rstrip()
        loc = '0'+str(portHeader[stc][6].decode('ASCII'))
        net = str(portHeader[stc][7].decode('ASCII')).rstrip()
        GPS = int(portHeader[stc][3])
        
        # remove GPS unlock or broken station
        if ( GPS == 1 or GPS == 2 ):
            chc = 0
            if instrument == 'FBA':
                chc = 1
            elif instrument == 'SP':
                chc = 4
            elif instrument == 'BB':
                chc = 7
            
            for ch in range(3):
                chn = 'Ch'+str(chc+ch)
                
                stats = {'network': net, 'station': stn, 'location': loc,
                        'channel': chn, 'npts': npts, 'sampling_rate': samp,
                        'starttime': sttime}
                
                data = np.array(idata[ch][stc], dtype=float)
                sttmp = Stream([Trace(data=data, header=stats)])
                afst += sttmp

    return afst

def unpackPfile(infile):
    
    with open(infile) as f:
        lines = f.readlines()
    
    tmp = lines[0]
    year = int(tmp[1:5])
    month = int(tmp[5:7])
    day = int(tmp[7:9])
    hour = int(tmp[9:11])
    minute = int(tmp[11:13])
    sec = float(tmp[13:19])

    dt = datetime(year,month,day,hour,minute,int(sec//1),int(sec%1 * 1000000))
    mag = float(tmp[40:44])

    pfile_info = {}
    pfile_info["ori_time"] = dt
    pfile_info["mag"] = mag

    intensity = {}
    arrival_time = {}
    weighting = {}
    pga = {}
    for i in lines[1:]:
        sta = i[:5].strip() # strip 去掉左右空格
        weighting[sta] = int(float(i[35:39]))
        intensity[sta] = int(i[76:77])
        pga[sta] = float(i[78:83])
        arrival_time[sta] = pfile_info["ori_time"].replace(minute=int(i[21:23]),second=0,microsecond=0) + timedelta(seconds=float(i[23:29]))
    pfile_info["intensity"] = intensity
    pfile_info["arrival_time"] = arrival_time
    pfile_info["weighting"] = weighting
    pfile_info["pga"] = pga
    
    return pfile_info

def process_command():
    parser = argparse.ArgumentParser()
    # parser.add_argument('--foo', help='foo help')
    parser.add_argument('--a_filename', '-a', type=str, required=True, help='Input a_file name')
    parser.add_argument('--p_filename', '-p', type=str, required=True, help='Input p_file name')
    return parser.parse_args()

def get_StationInfo():
    s = []
    with open('nsta24.dat') as f:
        for line in f.readlines():
            l = line.strip().split()
            # print(l)
            d = {}
            if l[-1] == '99999999':
                d['station'] = l[0]
                d['lon'] = float(l[1])
                d['lat'] = float(l[2])
                d['location'] = '0' + l[5]
                d['institute'] = l[6]
                d['network'] = l[7]
                d['instrument'] = l[8]
                d['E'] = float(l[9])
                d['N'] = float(l[10])
                d['Z'] = float(l[11])
                s.append(d)
    return s

def mulFactor(tr1):
    stationInfo = get_StationInfo()
    for idx, s in enumerate(stationInfo):
        instrument = ''
        if tr1.stats.channel == 'Ch1' or tr1.stats.channel == 'Ch2' or tr1.stats.channel == 'Ch3':
            instrument = 'FBA'
        elif tr1.stats.channel == 'Ch4' or tr1.stats.channel == 'Ch5' or tr1.stats.channel == 'Ch6':
            instrument = 'SP'
        elif tr1.stats.channel == 'Ch7' or tr1.stats.channel == 'Ch8' or tr1.stats.channel == 'Ch9':
            instrument = 'BB'

        if (tr1.stats.station == s['station'] and instrument == s['instrument'] and 
        tr1.stats.network == s['network'] and tr1.stats.location == s['location']):
            if tr1.stats.channel == 'Ch1':
                tr1.data = tr1.data * s['Z']
                
            elif tr1.stats.channel == 'Ch2':
                tr1.data = tr1.data * s['N']
                
            elif tr1.stats.channel == 'Ch3':
                tr1.data = tr1.data * s['E']
                
def calaulateSNR(dd_before, dd_after):
    s = np.sum(np.square(dd_after))
    n = np.sum(np.square(dd_before))

    snr = s / n
    return snr

def zeroCrossing(dd):
    zero_cross = 0
    for i in range(0, len(dd)-1):
        if dd[i]*dd[i+1] < 0:
            zero_cross += 1
    return zero_cross

def checkSignal(dd):
    for i in range(len(dd)-9):
        part = dd[i:i+10]
        count = 0
        for p in part:
            if p == 0.0:
                count += 1
        if count == 10:
            return False   
    return True

def writeOut(filename, pfile_name, name, snr, mag, pga):
    with open(filename, 'a') as f:
        f.write(pfile_name)
        f.write(' ')
        f.write(name)
        f.write(' ')
        f.write(str(snr))
        f.write(' ')
        f.write(str(mag))
        f.write(' ')
        f.write(str(pga))
        f.write('\n')
                         

def run_5stas():
    afile = args.a_filename
    pfile = args.p_filename

    st = unpackAfile(afile)
    pfile_info = unpackPfile(pfile)

    stalist = []
    for i in pfile_info['arrival_time']:
        stalist.append(i[:4].strip())
    
    print(stalist)
    
    for tr in st:
        # 加速度 地表
        # FBA:Ch1~3, SP:Ch4~6, BB:Ch7~9
        if (tr.stats.station in stalist and tr.stats.location=="01" and ((tr.stats.station=="ETLH" and tr.stats.network=="BH" ) or 
        (tr.stats.station=="ETM" and tr.stats.network=="SMT") or (tr.stats.station=="ESL" and tr.stats.network=="SMT") or
        (tr.stats.station=="SHUL" and tr.stats.network=="BH") or (tr.stats.station=="ENA" and tr.stats.network=="SMT"))):

            if tr.stats.channel=="Ch1" or tr.stats.channel=="Ch2" or tr.stats.channel=="Ch3":
                print("Running "+tr.stats.station, tr.stats.channel)
                tr1 = tr.copy()

                # 判斷紀錄有無被picking
                pfile_arrival = int((pfile_info["arrival_time"][tr.stats.station] - tr.stats.starttime.datetime) / timedelta(microseconds=10000))
                if pfile_arrival != 0:
                    # 100 是1秒 
                    time_before = 500 
                    time_after = 3000        

                    st = pfile_arrival - time_before
                    ed = pfile_arrival + time_after

                    dd_count = tr1.data[st:ed+1]

                    # 減平均 去offset
                    tr1.data = tr1.data - np.mean(tr1.data)

                    # 乘係數
                    mulFactor(tr1)
                    dd_fac = tr1.data[st:ed+1]
                    dd_fac_before = tr1.data[st:pfile_arrival]
                    dd_fac_after = tr1.data[pfile_arrival:ed+1]
                    
                    # 判斷有無斷訊
                    if checkSignal(dd_count):
                        # 計算SNR
                        snr = float('%.3f' % calaulateSNR(dd_fac_before, dd_fac_after))
                        # PGA Magnitude
                        pga = float('%.1f' % max(tr1.data[pfile_arrival:]))
                        mag = pfile_info["mag"]
                        # 算通過0的次數
                        zero_cross = zeroCrossing(dd_fac)
                        if zero_cross > 40:
                
                            fig = plt.figure(figsize=(12, 6))
                            ax = fig.add_subplot(1, 1, 1)

                            tr_times = tr.times("matplotlib")
                            ax.plot(tr_times, tr1.data, "b-")  #plot wave
                            ymin, ymax = ax.get_ylim()
                            
                            ax.vlines(tr_times[pfile_arrival], ymin, ymax, color='r', linewidth=3, label="pfile")
                
                            plt.xlim(tr_times[pfile_arrival-time_before], tr_times[pfile_arrival+time_after])  #放大

                            title = "_".join([pfile[-12:],tr.stats.station,tr.stats.channel,tr.stats.network,tr.stats.location
                                ,str(pfile_info["intensity"][tr.stats.station])
                                ,str(pga)
                                ,str(pfile_info["mag"])
                                ,str(pfile_info["weighting"][tr.stats.station])
                                ,str(snr)
                                ])
                        
                            plt.title(title)
                            ax.xaxis_date()
                            fig.autofmt_xdate()
                            plt.legend()

                            plt.savefig("D:\\00000\\" + title + ".png")
                            plt.close()

                            str1 ="D:\\00000\\" + title + ".txt"
                            with open(str1, 'a') as f:
                                for j in dd_fac:
                                    f.write(str(j)+"\n")
                            
                            #name = "_".join([tr.stats.station,tr.stats.channel,tr.stats.network,tr.stats.location])
                            #filename = 'D:\\123.txt'
                            #writeOut(filename, pfile[-12:], name, snr, mag, pga)

args = process_command()
run_5stas()
