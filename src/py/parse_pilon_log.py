import os
import re
import subprocess
import argparse

# March 14, 2018
# This Python script parses a Pilon standard output
# This script is not yet finished

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('--pilon_out', type=file, dest='pilon_out', required=True, help='Pilon standard out')
parser.add_argument('--pilon_summary', type=file, dest='pilon_summary', required=True, help='Pilon Summary of standard out')
args = parser.parse_args()

# Open output file
pilon_summary_file = open(args.pilon_summary, 'w+')

for line in args.pilon_out.readlines():
    if 'Scanning BAMs' in line:
        pilon_summary_file.write(line)
        print line

# grep -w 'Scanning BAMs' -A 6 pilon_log
