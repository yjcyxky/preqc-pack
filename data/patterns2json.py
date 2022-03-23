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
idx_values = []
for row in csvreader:
  for field in row:
    if len(row) > 0:
      values = field.split('|')
      print(values)
      # values[0] --> Index, values[1] --> Alt/Ref
      idx_values.append(int(values[0]))
      data[values[2]] = [int(values[0]), int(values[1])]

results = {
  "data": data,
  "indexes": list(set(idx_values)),
  "count": max(idx_values) - min(idx_values) + 1
}

write2json(results)
write2bson(results)
    