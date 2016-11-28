import os
import sys

"""
python ./data/bin/calc_string_length.py file_name column_number
"""

file_name = sys.argv[1]
column_id = int(sys.argv[2])


with open(file_name, 'r') as fh:
    for buf in fh:
        records = buf.split('\t')
        print(len(records[column_id]))





