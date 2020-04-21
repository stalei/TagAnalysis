#  © Shahram Talei @ 2020 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#read a sage files

# $python ReadSage.py SageModelFile
# $python ReadSage.py /media/shahram/SD/717/model_z0.000_0

import argparse
import codecs

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("SageFile", type=str)
    args = parser.parse_args()
    #f=open(args.SageFile,"rb")
    #line=f.readline()
    #while line!='':
    #    print(line)
    #    line=f.readline()
    ########
    #t=next(f)
    #t2=next(f)
    #print(t2)
    #num=list(f.read())
    #print (num)
    ########
    with open(args.SageFile, mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        print(fileContent)
        base64_data = codecs.encode(fileContent, 'base64x')
        print(base64_data)
        file.close()
