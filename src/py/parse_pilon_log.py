import sys
import os
import re
import subprocess
import argparse


#This Python script parses a Pilon standard output
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('--pilon_out', type=file, dest='pilon_out', required=True, help='Pilon standard out')
args = parser.parse_args()
print "Test arg parse", args.pilon_out

for line in args.pilon_out.readlines():
    line.strip('\n')
    print line

