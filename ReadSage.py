#  © Shahram Talei @ 2020 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#read a sage files

# $python ReadSage.py SageModelFile
# $python ReadSage.py


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("SageFile", type=str)
    
