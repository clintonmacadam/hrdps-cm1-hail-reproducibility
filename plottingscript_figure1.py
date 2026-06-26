import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Rectangle
import cartopy.io.img_tiles as cimgt
import rasterio
from matplotlib.colors import PowerNorm
import numpy as np
from matplotlib.colors import LightSource
from rasterio.windows import from_bounds
from rasterio.windows import from_bounds, bounds as window_bounds
from pyproj import Transformer
from pyproj import CRS

# -----------------------------
# Map projections
# -----------------------------
proj = ccrs.LambertConformal(
    central_longitude=-100,
    central_latitude=45,
    standard_parallels=(33, 60)
)

pc = ccrs.PlateCarree()

fig = plt.figure(figsize=(13, 10))

pc = ccrs.PlateCarree()

fig = plt.figure(figsize=(10, 12))

# -----------------------------
# Main Alberta map
# -----------------------------
ax = fig.add_axes([0.08, 0.05, 0.84, 0.88], projection=pc)

alberta_extent = [-120.5, -109, 49, 57]
ax.set_extent(alberta_extent, crs=pc)

ax.add_feature(cfeature.BORDERS, linewidth=0.8)
ax.add_feature(cfeature.STATES, linewidth=0.6, edgecolor="0.4")
ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
#ax.add_feature(cfeature.LAKES, edgecolor="0.4", facecolor="none", linewidth=0.5)

import rasterio
from rasterio.windows import from_bounds

west, south, east, north = -120.5, 49, -109, 57

dem_crs_pyproj = "EPSG:3978"   # NAD83 / Canada Atlas Lambert
dem_crs_cartopy = ccrs.epsg(3978)

with rasterio.open("cdem-canada-dem-clipped.tif") as src:
    terrain = src.read(1)
    bounds = src.bounds

terrain_10 = terrain[::10, ::10]

ax.imshow(
    terrain_10,
    origin="upper",
    extent=[bounds.left, bounds.right, bounds.bottom, bounds.top],
    transform=ccrs.PlateCarree(),
    cmap="terrain",
    vmin=300,
    vmax=2500,
)

ax.plot(-114.07, 51.04, marker='o', markersize=10, color='black', transform=ccrs.PlateCarree())
ax.text(-113.95, 51.04, 'Calgary', fontsize=12, transform=ccrs.PlateCarree())
ax.plot(-113.81, 52.27, marker='o', markersize=10, color='black', transform=ccrs.PlateCarree())
ax.text(-113.69, 52.27, 'Red Deer', fontsize=12, transform=ccrs.PlateCarree())
ax.plot(-113.5, 53.54, marker='o', markersize=10, color='black', transform=ccrs.PlateCarree())
ax.text(-113.38, 53.54, 'Edmonton', fontsize=12, transform=ccrs.PlateCarree())
ax.plot(-114.1180, 51.2123, marker='o', markersize=10, color='r', transform=ccrs.PlateCarree())
ax.plot(-114.2874, 52.1809, marker='o', markersize=10, color='r', transform=ccrs.PlateCarree())
ax.plot(-114.135, 49.883, marker='o', markersize=10, color='b', transform=ccrs.PlateCarree())
ax.plot(-110.386, 51.98, marker='o', markersize=10, color='white', transform=ccrs.PlateCarree())
ax.plot(-114.85, 52.46, marker='o', markersize=10, color='white', transform=ccrs.PlateCarree())
ax.plot(-115.54, 52.49, marker='o', markersize=10, color='b', transform=ccrs.PlateCarree())
ax.plot(-114.89, 52.57, marker='o', markersize=10, color='r', transform=ccrs.PlateCarree())
ax.plot(-115.89, 53.10, marker='o', markersize=10, color='r', transform=ccrs.PlateCarree())
ax.plot(-115.22, 53.27, marker='o', markersize=10, color='white', transform=ccrs.PlateCarree())

#ax.set_title("Alberta", fontsize=16)

# -----------------------------
# North America inset map
# -----------------------------
inset_ax = fig.add_axes([0.59, 0.52, 0.32, 0.32], projection=pc)

inset_ax.set_extent([-170, -50, 15, 75], crs=pc)

inset_ax.add_feature(cfeature.LAND, facecolor="0.9")
inset_ax.add_feature(cfeature.OCEAN, facecolor="0.8")
inset_ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
inset_ax.add_feature(cfeature.BORDERS, linewidth=0.4)

# Alberta locator box
rect = Rectangle(
    (-120.5, 48.5),
    11.5,
    8.5,
    linewidth=2,
    edgecolor="red",
    facecolor="none",
    transform=pc,
    zorder=10
)

inset_ax.add_patch(rect)
#inset_ax.set_title("Location in North America", fontsize=10)

plt.savefig("alberta_with_north_america_inset.png", dpi=300, bbox_inches="tight")
plt.show()
