import sys
import fiona
import shapely
from shapely import geometry
from shapely.geometry import shape
from shapely.geometry import mapping


a = fiona.open(sys.argv[1])
a_pol = a.next()
a_geom = shape(a_pol['geometry'])
a_poly_data = a_pol["geometry"]["coordinates"][0]
inner_poly = geometry.Polygon(a_poly_data)  

b = fiona.open(sys.argv[2])
b_pol = b.next()
b_geom = shape(b_pol['geometry'])
b_poly_data = b_pol["geometry"]["coordinates"][0]
outer_poly = geometry.Polygon(b_poly_data)   

ring_poly=sys.argv[3]
output=inner_poly.symmetric_difference(outer_poly)

schema = {
    'geometry': 'Polygon',
    'properties': {'id': 'int'},
}

# Write a new Shapefile
with fiona.open(ring_poly, 'w', 'ESRI Shapefile', schema) as c:
    ## If there are multiple geometries, put the "for" loop here
    c.write({
        'geometry': mapping(output),
        'properties': {'id': 123},
    })




