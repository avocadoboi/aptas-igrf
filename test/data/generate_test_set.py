import ppigrf
from datetime import datetime
from itertools import product
from numpy import linspace

def decimal_year(date: datetime) -> float:
    year_start = datetime(date.year, 1, 1)
    year_end = datetime(date.year + 1, 1, 1)
    year_length = (year_end - year_start).total_seconds()
    return date.year + (date - year_start).total_seconds()/year_length

EARTH_RADIUS = 6371.2
longitudes = linspace(0, 360, 5, endpoint=False)
co_latitudes = linspace(5, 175, 5)
altitudes = linspace(0, 1000, 5)
dates = linspace(datetime(2025, 1, 1), datetime(2030, 1, 1), 5)
with open("test_set", "w") as file:
    for lon, co_lat, altitude, date in product(longitudes, co_latitudes, altitudes, dates):
        Br, Btheta, Bphi = ppigrf.igrf_gc(EARTH_RADIUS + altitude, co_lat, lon, date)
        # print(lon, co_lat, altitude, 1970 + date.timestamp()/60**2/24/365, B_east[0], B_north[0], B_up[0], file = file)
        print(90 - co_lat, lon, altitude, decimal_year(date), Bphi[0], -Btheta[0], Br[0], file = file)
