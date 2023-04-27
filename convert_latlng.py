import pandas as pd
import geopandas as gpd
import collections
import re
import time

import itertools

import numpy as np


from shapely import wkb, wkt
from shapely.geometry import Polygon, Point, MultiPolygon

from pyproj import Transformer


def wkt_loads(x):
    try:
        return wkt.loads(x)
    except Exception:
        return None


def wkt_dumps(x):
    try:
        return wkt.dumps(x)
    except Exception:
        return None


# class Geometry(sqlalchemy.types.UserDefinedType):
#     # custom sqlalchemy usertype to convert wkt-column to geometry type

#     def __init__(self, srid: int = 4326):
#         self.srid = srid

#     def get_col_spec(self):
#         return "GEOMETRY"

#     def bind_expression(self, bindvalue):
#         return sqlalchemy.text(f'geometry::STGeomFromText(:{bindvalue.key},{self.srid})').bindparams(bindvalue)


def convert_latlng_to_xy(geom):
    """Convert a Shapley polygon from long lat to x y.

    Args:
        geom (Shapely Polygon): Polygon object to convert

    Returns:
        Polygon: New Polygon with converted coordinates
    """

    transformer = Transformer.from_crs("WGS84", "EPSG:27700", always_xy=False)
    # For some reason, the first time the transformer is used to convert latlngs,
    # the result is inf values - kind of like it needs initialising before use...
    # Running the transformer over the first couple of coordinates, then doing the
    # actual conversion seems to work though.
    coords = np.array(list(geom.exterior.coords))
    lng, lat = list(coords[:, 0]), list(coords[:, 1])
    x, y = transformer.transform(lat[0:1], lng[0:1])
    x, y = transformer.transform(lat, lng)
    xy = list(zip(x, y))

    return Polygon(xy)


def convert_xy_to_latlng(geom):
    """Convert a Shapley polygon from long lat to x y.

    Args:
        geom (Shapely Polygon): Polygon object to convert

    Returns:
        Polygon: New Polygon with converted coordinates
    """

    transformer = Transformer.from_crs('epsg:27700', 'WGS84') # epsg:4326
    if isinstance(geom, Polygon):
        geom = MultiPolygon([geom])

    new_geoms = []
    for g in geom.geoms:
        coords = np.array(list(g.exterior.coords))
        ll = np.array([transformer.transform(x, y) for (x, y) in coords])
        ll = np.array(list(zip(ll[:, 1], ll[:, 0])))
        ll = np.array(ll)
        ll = np.round(ll, 6)
        new_geoms.append(Polygon(ll))


    if len(new_geoms) == 1:
        return Polygon(new_geoms[0])

    return MultiPolygon(new_geoms)


