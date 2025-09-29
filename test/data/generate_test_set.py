import ppigrf
from datetime import datetime
from itertools import product
from numpy import linspace

def decimal_year(date: datetime) -> float:
    year_start = datetime(date.year, 1, 1)
    year_end = datetime(date.year + 1, 1, 1)
    year_length = (year_end - year_start).total_seconds()
    return date.year + (date - year_start).total_seconds()/year_length

longitudes = linspace(0, 360, 5)
co_latitudes = linspace(1, 179, 5)
radii = 6378 + linspace(0, 1000, 5)
dates = linspace(datetime(2020, 1, 1), datetime(2030, 1, 1), 5)
with open("test_set", "w") as file:
    for lon, co_lat, r, date in product(longitudes, co_latitudes, radii, dates):
        B_east, B_north, B_up = ppigrf.igrf_gc(r, co_lat, lon, date)
        # print(lon, co_lat, r, 1970 + date.timestamp()/60**2/24/365, B_east[0], B_north[0], B_up[0], file = file)
        print(lon, co_lat, r, decimal_year(date), B_east[0], B_north[0], B_up[0], file = file)
