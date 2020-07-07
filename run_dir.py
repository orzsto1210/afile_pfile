# coding: utf-8
import os
import subprocess, time

def run():
    for directory, sub_directory, files in os.walk(r"E:\\CWB_data"):
        directory = directory+"\\"
        print("directory:   "+directory)
        # print(sub_directory)
        for i in files:
            if i.find(".A1")>-1:
                #print(i)
                hh = i[:-4]
                ee = i[-2:]
                pfile = hh+".P"+ee
                afile = hh+".A"+ee
                
                p = subprocess.Popen(['python', 'afile_pfile.py', '--a_filename', directory+afile, '--p_filename',  directory+pfile])
                p.communicate()
                p.wait()
                print("do it!")

run()

# count = 0
# for directory, sub_directory, files in os.walk(r"E:\CWB_data\2019\felt"):
#     directory = directory+"\\"
#     #print("directory:   "+directory)
#     # print(sub_directory)
#     for i in files:
#         if i.find(".A1")>-1:
#             print(i)
#             count+=1
# print(count)
