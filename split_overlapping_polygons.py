
import fiona
import folium
import geopandas as gpd
import numpy as np
from shapely.affinity import scale
from shapely.geometry import (GeometryCollection, LineString, MultiPoint,
                              Point, Polygon)
from shapely import wkt
from shapely.ops import snap, split
from scipy.spatial import Voronoi, voronoi_plot_2d
import networkx as nx

fiona.drvsupport.supported_drivers['KML'] = 'rw'

def split_polygons_at_intesection(a, b, tolerance):
    a1 = a
    b1 = b

    # Get the points of intersection between the two polygons
    p1 = a1.boundary.intersection(b1.boundary)
    l = LineString(p1.geoms) # Convert to LineString

    # Get the intersection of the polygon ie the overlapping area
    c1 = a1.intersection(b1)
    c = split(c1, l) # Split this at the intersection points

    # Remove the overlap polygon from the original polygons
    a2 = a1.difference(c1)
    b2 = b1.difference(c1)

    # Merge the polygon with the nearest half overlapping polygon
    idx = np.argmin(np.array([a2.centroid.distance(i) for i in c.geoms]))
    a3 = a2.union(c.geoms[idx])

    idx = np.argmin(np.array([b2.centroid.distance(i) for i in c.geoms]))
    b3 = b2.union(c.geoms[idx])

    # Due to tolerance issues, remove points from each polygon that may still
    # be contained in the opposing polygon

    pts = MultiPoint(list(a3.exterior.coords))
    a5 = Polygon([pp for pp in pts.geoms if b3.buffer(-tolerance).contains(pp)==False])

    pts = MultiPoint(list(b3.exterior.coords))
    b5 = Polygon([pp for pp in pts.geoms if a3.buffer(-tolerance).contains(pp)==False])

    return a5, b5

def split_polygons(a, b):
    p1 = a.boundary.intersection(b.boundary)
    l = LineString(p1.geoms)

    c = a.union(b)
    return split(c, l)


def wkt_loads(x):
    try:
        return wkt.loads(x)
    except Exception:
        return None


def wkt_dumps(x):
    try:
        return wkt.dumps(x, output_dimension=2)
    except Exception:
        return None


def split_polygon_voronoi(a, b):
    a1, b1 = a, b

    # Union of both polygons
    d1 = a1.union(b1)

    # Points of intersection
    p1 = a1.boundary.intersection(b1.boundary)

    # Overlap polygon
    c1 = a1.intersection(b1)

    if c1.type=='MultiPolygon':
        c2 = []
        ci = c1.geoms[0]
        for ci in c1.geoms:
            pp = []
            ci_coords = list(ci.exterior.coords)
            pp.append(ci)
            for i in p1.geoms:
                if list(i.coords)[0] in ci_coords:
                    pp.append(Point(i))
            c2.append(pp)
    else:
        c2 = [(c1, p1.geoms[0], p1.geoms[1])]
    # GeometryCollection(c2[0])

    d2 = d1
    for ci, pp, pq in c2:
        cp = list(ci.exterior.coords)

        # Get Voronoi of overlapping polygon
        v = Voronoi(cp)
        v1 = MultiPoint(v.vertices)
        # voronoi_plot_2d(v)

        # Get the vornoi vertices indices that are within the overlapping polygon
        vp = [idx for idx, i in enumerate(v1.geoms) if ci.contains(i)]

        if len(vp)>1:
            # get the index of the vor vertices that is closest to each point of intersection
            v2a = np.array([(i.distance(pp), idx) for idx, i in enumerate(v1.geoms) if (idx in vp)])
            v2a_idx = v2a[np.argmin(v2a, axis=0)[0]][1]

            v2b = np.array([(i.distance(pq), idx) for idx, i in enumerate(v1.geoms) if (idx in vp)])
            v2b_idx = v2b[np.argmin(v2b, axis=0)[0]][1]

            # Get the ridge points for the vertices within the overlap polygon
            ridge_points = [(m, n) for (m, n) in v.ridge_vertices if ((m in vp) & (n in vp))]

            # Create a graph of edges so that we can find the path from the start to end
            g = nx.from_edgelist(ridge_points)

            path = nx.shortest_path(g, source=v2a_idx, target=v2b_idx)
            path = [int(i) for i in path]

        else:
            path = vp

        # Add the two points of intersection at the start and end
        l1 = LineString(list(pp.coords) + list(v.vertices[path]) + list(pq.coords))

        # Take the difference of the union shape and the line with a buffer.
        # This is better than using "split" as for shape with multiple overlapping regions,
        # it will split part a polygon.
        d2 = d2.difference(l1.buffer(1e-6))

    return d2


data = gpd.read_file("./Polygon split.kml", driver='KML')
data = data.to_crs("EPSG:27700")
data = data.set_index('Name')
data['geometry'] = data['geometry'].apply(wkt_dumps)
data['geometry'] = data['geometry'].apply(wkt_loads)
data = gpd.GeoDataFrame(data, geometry='geometry', crs='EPSG:27700')
data

A = data.loc['A', 'geometry']
B = data.loc['B', 'geometry']
C = data.loc['C', 'geometry']
D = data.loc['D', 'geometry']

# these are the overlapping shapes
# DB, AB, AC


a1 = A
b1 = C
GeometryCollection([a1, b1])

split_polygon_voronoi(A, B)
split_polygon_voronoi(A, C)
split_polygon_voronoi(D, B)

##############################

split_polygons(a1, b1)

##############################

