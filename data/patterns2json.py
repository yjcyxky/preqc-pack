#!/usr/bin/python3

import csv
import json, bson
from os import sep

file = open('./pattern_file.tsv')
csvreader = csv.reader(file, delimiter="\t")

def write2json(data):
  with open("patterns.json", "w") as f:
    f.write(json.dumps(data))
    
def write2bson(data):
  with open("patterns.bson", "wb") as bf:
    bf.write(bson.dumps(data))    

data = {}
for row in csvreader:
  for field in row:
    if len(row) > 0:
      values = field.split('|')
      print(values)
      # values[0] --> Index, values[1] --> Alt/Ref
      data[values[2]] = [int(values[0]), int(values[1])]

write2json(data)
write2bson(data)
    