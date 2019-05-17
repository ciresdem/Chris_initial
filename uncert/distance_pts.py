import sys
import itertools
import math
from math import cos, asin, sqrt
import numpy as np
#Haversine formula


def distance(latlon1, latlon2):
	lat1=float(latlon1[0])
 	lon1=float(latlon1[1])
 	lat2=float(latlon2[0])
 	lon2=float(latlon2[1])
 	radius = 637100 # m
 	dlat = math.radians(lat2-lat1)
 	dlon = math.radians(lon2-lon1)
 	a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
 	c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
 	d = radius * c
 	return d

def unique(seq, idfun=None): 
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result


#Loop through point
#Find the closest point and get the distance
#sum all minimum distances
#I want to maximum the sum of these minimum distances
data=(sys.argv[1])
with open(data) as f:
    lis=[line.split() for line in f]
    #print lis 

latlon = []
f = open( data, 'rU' ) #open the file in read universal mode
for line in f:
    cells = line.split( "," )
    latlon.append((cells[0], cells[1], cells[2])) #since we want the first, second column
f.close()

print latlon



ids0=[]
ids1=[]
results=[]

for p0, p1 in itertools.combinations(latlon, 2):
	p0_id=p0[0]
	p1_id=p1[0]
	p0_latlon=p0[1:3]
	p1_latlon=p1[1:3]
	dist=distance(p0_latlon,p1_latlon)
	ids0.append(p0_id)
	ids1.append(p1_id)
	results.append(dist)

rows = zip(ids0,ids1,results)
#print "ROWS ARE", rows

all_ids=ids0+ids1
unique_ids=unique(all_ids)
min_distances=[]

for tile in unique_ids:
	tile_distances=[]
	for tuple_ in rows:
		if tuple_[0] == tile or tuple_[1] == tile:
			#print "tile is", tile
			#print "distance is", tuple_[2]
			tile_distances.append(tuple_[2])
	min_dist=min(tile_distances)
	min_distances.append(min_dist)

#print min_distances
print sum(min_distances)
