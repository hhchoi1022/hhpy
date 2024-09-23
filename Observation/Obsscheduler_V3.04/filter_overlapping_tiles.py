
#%%
from astropy.coordinates import SkyCoord
import astropy.units as u
import shapely.geometry as geom
from astropy.io import ascii
from astropy.table import Table
#%%
def find_overlapping_tiles(tbl1, tbl2, 
                            FOV_RA_tbl1=None, 
                            FOV_Dec_tbl1=None, 
                            FOV_RA_tbl2=None, 
                            FOV_Dec_tbl2=None, 
                            overlap_threshold=0.5,
                            find_non_overlap : bool = True):
    """
    Filters tiles from the second table that overlap more than a threshold percentage 
    with any tiles from the first table, using spherical geometry for accurate corner calculation.
    
    Parameters:
    - tbl1: astropy.Table containing ra1, ra2, ra3, ra4, dec1, dec2, dec3, dec4 or central ra, dec for telescope 1
    - tbl2: astropy.Table containing ra1, ra2, ra3, ra4, dec1, dec2, dec3, dec4 or central ra, dec for telescope 2
    - FOV_RA_tbl1, FOV_Dec_tbl1: Field of View (RA, Dec) for telescope 1 if corners are not present
    - FOV_RA_tbl2, FOV_Dec_tbl2: Field of View (RA, Dec) for telescope 2 if corners are not present
    - overlap_threshold: float, overlap percentage threshold (default is 0.5 for 50%)
    
    Returns:
    - filtered_tbl2: astropy.Table containing the filtered tiles from tbl2
    """
    
    # Helper function to create corners using spherical geometry
    def create_corners_from_center(ra, dec, FOV_RA, FOV_Dec):
        # Create SkyCoord for the center of the tile
        center = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
        
        # Calculate the four corners of the field of view using offsets
        corners = []
        offsets = [(-FOV_RA/2, FOV_Dec/2),  # Top-left
                   (FOV_RA/2, FOV_Dec/2),   # Top-right
                   (FOV_RA/2, -FOV_Dec/2),  # Bottom-right
                   (-FOV_RA/2, -FOV_Dec/2)] # Bottom-left
        
        for offset in offsets:
            corner = center.spherical_offsets_by(offset[0]*u.deg, offset[1]*u.deg)
            corners.append((corner.ra.deg, corner.dec.deg))
        
        return corners
    
    # Process tbl1, check if it contains corner coordinates or needs to generate from center
    def create_polygons(table, FOV_RA=None, FOV_Dec=None):
        tiles = []
        if 'ra1' in table.colnames and 'dec1' in table.colnames:  # Check if corners are present
            tiles = [
                [(table['ra1'][i], table['dec1'][i]), 
                 (table['ra2'][i], table['dec2'][i]), 
                 (table['ra3'][i], table['dec3'][i]), 
                 (table['ra4'][i], table['dec4'][i])]
                for i in range(len(table))
            ]
        elif 'ra' in table.colnames and 'dec' in table.colnames and FOV_RA is not None and FOV_Dec is not None:
            # If no corner data, generate corners from center ra/dec using FOV_RA and FOV_Dec
            tiles = [
                create_corners_from_center(table['ra'][i], table['dec'][i], FOV_RA, FOV_Dec)
                for i in range(len(table))
            ]
        else:
            raise ValueError("Table must contain either corner coordinates or central (ra, dec) with FOV info.")
        
        return [geom.Polygon(corners) for corners in tiles]

    # Create polygons for each tile in tbl1 and tbl2
    telescope1_polygons = create_polygons(tbl1, FOV_RA_tbl1, FOV_Dec_tbl1)
    telescope2_polygons = create_polygons(tbl2, FOV_RA_tbl2, FOV_Dec_tbl2)

    # List to store indices of tiles from tbl2 that are not overlapping more than the threshold
    non_overlap_indices = []
    overlap_indices = []

    for i, tile2_poly in enumerate(telescope2_polygons):
        overlap = False
        for tile1_poly in telescope1_polygons:
            # Check if polygons intersect
            if tile1_poly.intersects(tile2_poly):
                intersection_area = tile1_poly.intersection(tile2_poly).area
                tile2_area = tile2_poly.area
                # Calculate overlap percentage
                overlap_percentage = intersection_area / tile2_area
                if overlap_percentage > overlap_threshold:
                    overlap = True
                    break
        if not overlap:
            non_overlap_indices.append(i)
        else:
            overlap_indices.append(i)

    # Create filtered table by selecting only valid indices
    overlap_tbl = tbl2[overlap_indices]
    non_overlap_tbl = tbl2[non_overlap_indices]
    
    if find_non_overlap:
        return non_overlap_tbl
    else:
        return overlap_tbl
#%%
# Example usage
def main():
    # Read the astropy tables
    tbl1 = ascii.read('/home/hhchoi1022/S240910ci_PRELIMINARY/SkyGridCatalog_RASA36_90.csv')
    tbl1.remove_columns(['ra1', 'dec1', 'ra2', 'dec2', 'ra3', 'dec3', 'ra4', 'dec4'])
    tbl2 = ascii.read('/home/hhchoi1022/S240910ci_PRELIMINARY/SkyGridCatalog_7DT_90.csv')
    tbl2.remove_columns(['ra1', 'dec1', 'ra2', 'dec2', 'ra3', 'dec3', 'ra4', 'dec4'])

    # Specify FOV values for the two telescopes (in degrees)
    FOV_RA_tbl1 = 2.685  # Example FOV_RA for tbl1
    FOV_Dec_tbl1 = 2.688  # Example FOV_Dec for
    FOV_RA_tbl2 = 1.414  # Example FOV_RA for tbl1
    FOV_Dec_tbl2 = 0.915  # Example FOV_Dec for

    # Filter tiles from tbl2 based on overlap with tbl1
    filtered_tbl2 = find_overlapping_tiles(tbl1, tbl2, 
                                           FOV_RA_tbl1, 
                                           FOV_Dec_tbl1, 
                                           FOV_RA_tbl2, 
                                           FOV_Dec_tbl2,
                                           find_non_overlap=False)

# %%


from astropy.coordinates import SkyCoord
import astropy.units as u
import shapely.geometry as geom
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
from descartes import PolygonPatch

#%%
def find_overlapping_tiles(tbl1, tbl2, 
                            FOV_RA_tbl1=None, 
                            FOV_Dec_tbl1=None, 
                            FOV_RA_tbl2=None, 
                            FOV_Dec_tbl2=None, 
                            overlap_threshold=0.5,
                            find_non_overlap : bool = True):
    """
    Filters tiles from the second table that overlap more than a threshold percentage 
    with any tiles from the first table, using spherical geometry for accurate corner calculation.
    
    Parameters:
    - tbl1: astropy.Table containing ra1, ra2, ra3, ra4, dec1, dec2, dec3, dec4 or central ra, dec for telescope 1
    - tbl2: astropy.Table containing ra1, ra2, ra3, ra4, dec1, dec2, dec3, dec4 or central ra, dec for telescope 2
    - FOV_RA_tbl1, FOV_Dec_tbl1: Field of View (RA, Dec) for telescope 1 if corners are not present
    - FOV_RA_tbl2, FOV_Dec_tbl2: Field of View (RA, Dec) for telescope 2 if corners are not present
    - overlap_threshold: float, overlap percentage threshold (default is 0.5 for 50%)
    
    Returns:
    - filtered_tbl2: astropy.Table containing the filtered tiles from tbl2
    """
    
    # Helper function to create corners using spherical geometry
    def create_corners_from_center(ra, dec, FOV_RA, FOV_Dec):
        # Create SkyCoord for the center of the tile
        center = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
        
        # Calculate the four corners of the field of view using offsets
        corners = []
        offsets = [(-FOV_RA/2, FOV_Dec/2),  # Top-left
                   (FOV_RA/2, FOV_Dec/2),   # Top-right
                   (FOV_RA/2, -FOV_Dec/2),  # Bottom-right
                   (-FOV_RA/2, -FOV_Dec/2)] # Bottom-left
        
        for offset in offsets:
            corner = center.spherical_offsets_by(offset[0]*u.deg, offset[1]*u.deg)
            corners.append((corner.ra.deg, corner.dec.deg))
        
        return corners
    
    # Process tbl1, check if it contains corner coordinates or needs to generate from center
    def create_polygons(table, FOV_RA=None, FOV_Dec=None):
        tiles = []
        if 'ra1' in table.colnames and 'dec1' in table.colnames:  # Check if corners are present
            tiles = [
                [(table['ra1'][i], table['dec1'][i]), 
                 (table['ra2'][i], table['dec2'][i]), 
                 (table['ra3'][i], table['dec3'][i]), 
                 (table['ra4'][i], table['dec4'][i])]
                for i in range(len(table))
            ]
        elif 'ra' in table.colnames and 'dec' in table.colnames and FOV_RA is not None and FOV_Dec is not None:
            # If no corner data, generate corners from center ra/dec using FOV_RA and FOV_Dec
            tiles = [
                create_corners_from_center(table['ra'][i], table['dec'][i], FOV_RA, FOV_Dec)
                for i in range(len(table))
            ]
        else:
            raise ValueError("Table must contain either corner coordinates or central (ra, dec) with FOV info.")
        
        return [geom.Polygon(corners) for corners in tiles]

    # Create polygons for each tile in tbl1 and tbl2
    telescope1_polygons = create_polygons(tbl1, FOV_RA_tbl1, FOV_Dec_tbl1)
    telescope2_polygons = create_polygons(tbl2, FOV_RA_tbl2, FOV_Dec_tbl2)

    # List to store indices of tiles from tbl2 that are not overlapping more than the threshold
    non_overlap_indices = []
    overlap_indices = []

    overlap_polygons = []

    for i, tile2_poly in enumerate(telescope2_polygons):
        overlap = False
        for tile1_poly in telescope1_polygons:
            # Check if polygons intersect
            if tile1_poly.intersects(tile2_poly):
                intersection_area = tile1_poly.intersection(tile2_poly).area
                tile2_area = tile2_poly.area
                # Calculate overlap percentage
                overlap_percentage = intersection_area / tile2_area
                if overlap_percentage > overlap_threshold:
                    overlap = True
                    overlap_polygons.append(tile2_poly)
                    break
        if not overlap:
            non_overlap_indices.append(i)
        else:
            overlap_indices.append(i)

    # Create filtered table by selecting only valid indices
    overlap_tbl = tbl2[overlap_indices]
    non_overlap_tbl = tbl2[non_overlap_indices]
        # Plot the tiles using matplotlib polygons
    fig, ax = plt.subplots(figsize=(10, 8))

    # Helper function to plot polygons using matplotlib
    def plot_polygon(ax, polygon, edgecolor='black', facecolor='none', lw=1, linestyle='-'):
        x, y = polygon.exterior.xy
        ax.plot(x, y, color=edgecolor, lw=lw, linestyle=linestyle)
        ax.fill(x, y, color=facecolor, alpha=0.5 if facecolor != 'none' else 0)

    # Plot tbl1 tiles (black)
    for poly in telescope1_polygons:
        plot_polygon(ax, poly, edgecolor='black')

    # Plot tbl2 tiles (black dashed)
    for poly in telescope2_polygons:
        plot_polygon(ax, poly, edgecolor='blue')

    # Plot overlapping tiles (red)
    for poly in overlap_polygons:
        plot_polygon(ax, poly, edgecolor='red', facecolor='red')

    #ax.set_xlim(0, 360)
    #ax.set_ylim(-90, 90)
    ax.set_xlabel('Right Ascension (RA)')
    ax.set_ylabel('Declination (Dec)')
    ax.set_title('Sky Tiles: Black = All tiles, Red = Overlapping tiles')

    plt.grid(True)
    plt.show()
    if find_non_overlap:
        return non_overlap_tbl
    else:
        return overlap_tbl
#%%
# Example usage
def main():
    # Read the astropy tables
    tbl1 = ascii.read('/home/hhchoi1022/Desktop/Gitrepo/Observation/Obsscheduler_V3.04/7DT_scheduled.ascii_fixed_width', format = 'fixed_width')
    tbl2 = ascii.read('/home/hhchoi1022/S240910ci_UPDATE/SkyGridCatalog_KMTNet_90.csv')
    #tbl2.remove_columns(['ra1', 'dec1', 'ra2', 'dec2', 'ra3', 'dec3', 'ra4', 'dec4'])    
    #tbl2.remove_columns(['ra1', 'dec1', 'ra2', 'dec2', 'ra3', 'dec3', 'ra4', 'dec4'])

    # Specify FOV values for the two telescopes (in degrees)
    FOV_RA_tbl1 = 2.685  # RASA36
    FOV_Dec_tbl1 = 2.688  # RASA36
    FOV_RA_tbl2 = 1.414  # 7DT
    FOV_Dec_tbl2 = 0.915  # 7DT

    # Filter tiles from tbl2 based on overlap with tbl1
    filtered_tbl2 = find_overlapping_tiles(tbl1, tbl2, 
                                           FOV_RA_tbl1, 
                                           FOV_Dec_tbl1, 
                                           FOV_RA_tbl2, 
                                           FOV_Dec_tbl2,
                                           overlap_threshold= 0.01,
                                           find_non_overlap=False)
    filtered_tbl2.write('/home/hhchoi1022/Desktop/Gitrepo/Observation/Obsscheduler_V3.04/RASA36_scheduled.ascii_fixed_width', format = 'ascii.fixed_width', overwrite = True)
    

if __name__ == "__main__":
    main()


# %%
