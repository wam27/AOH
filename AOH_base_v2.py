# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 16:54:04 2022

@author: wam27
"""
#Birds Distribution Map
#This tool allows the creation of species distribution maps based on presence 
#   and absence points for 93 threatened bird species in Ecuador, and generates a
#   table with habitat distribution areas, as well as percentage of distribution 
#   area within protected reserves.

#The tool also have the option to generate a total richness distribution map
#   with a polygon user’s option to select an area in the map to see which and 
#   how many species are contributing to the total distribution.

#Steps
#1.Set data to variables and projections. Data inputs are Presence points 
#(ocurrence of birds which might be either in a series of csv files -one file 
#per species-, or downloading data directly from a website such as e-bird or 
#GBIF. For this case, I will use a data folder in which I stored 93 files after
#debugging information from GBIF databases), Absence points (a shapefile with
#hotspots of birds from e-bird), Elevation raster (dem 30m from ...), Habitat 
#raster (2016 forest cover raster from Environment Governmental Agency in Ecuador).


import os, arcpy
import arcpy.sa
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio import mask
import numpy as np


def AOH(df, n, occ, hull, dem, dem_proj, tc, g, 
        AOH_output, temp_output,se):
    
    epsg_WGS84= 4326
    m2= 1/10**6
    src_dem = rasterio.open(dem)
    src_dem_proj = rasterio.open(dem_proj)
    src_forest = rasterio.open(tc)
    cell_size=src_dem.transform[0]
    cell_size_proj=src_dem_proj.transform[0]
    
    if os.path.isfile(se):    
        se_df=gpd.read_file(se)
        Checklists=True
    else:
        Checklists=False
        
    report_Habitat=[]
    report_Errors=[]
    
    errors= AOH_output + "\Report_Error_Areas.csv"
    
    for b in  range(len(df)): #len(df)
        try:
            c= df.loc[b][0]
            print(str(b+1)+" " + c)
            
            occ_path=  occ + c + "_pres_eg.shp"
            hull_path= hull + c + "_hull2.shp"

            output1= AOH_output + '\{}_AOH{}.tif'.format(c,n)
            AOH_shp= temp_output + '\{}_AOH{}.shp'.format(c,n)
            Non_detection= temp_output +'\{}_Non_Detection{}.shp'.format(c,n)
            Cckl_spp_dir= temp_output +'\Cckl_spp'+n+'.shp'
            Cckl_spp_hab_dir= temp_output +'\Cckl_spp_hab'+n+'.shp'
            Absences= temp_output +'\{}_Absences{}.shp'.format(c,n)
            Chcklists_pixel= temp_output +'\{}_Absences_pix{}.shp'.format(c,n)
            Pres_Abs_path= temp_output +'\{}_Pres_Abs{}.shp'.format(c,n)
            idw_spp= temp_output + '\{}_IDW{}.tif'.format(c,n)
            Potential_AOH_spp= AOH_output + '\{}_PO_AOH{}.tif'.format(c,n)
            Unoccupied_AOH_spp= AOH_output + '\{}_AU_AOH{}.tif'.format(c,n)
            
            Report= AOH_output + "\Report_Habitat{}2.csv".format(n)
    
            if os.path.isfile(occ_path) and os.path.isfile(hull_path):
                
                presences = gpd.read_file(occ_path)
                alphahull = gpd.read_file(hull_path)
                
                GD=presences.COUNTRY.unique()
                GD=', '.join(GD)
                presences.crs = epsg_WGS84
                number_presences= len(presences)
                
                ######################################################################
                ######################## HABITAT #####################################
                ######################################################################
                print ('Refining Habitat')
                ##### TREE COVER #####################################################
                
                print ('1. Tree Cover')
                coords = [(x,y) for x, y in zip(presences.LONGITUDE, presences.LATITUDE)] # Read points from shapefile
                presences['Tree_cover'] = [x[0] for x in src_forest.sample(coords)] # Sample the raster at every point location and store values in DataFrame
                
                presences_above_2000_TC= presences[(presences.YEAR>1999) & (presences.Tree_cover>0)]
                
                if not presences_above_2000_TC.empty:
                    TC_threshold=int(round(pd.Series.quantile(presences_above_2000_TC['Tree_cover'], 0.25)))
                    if TC_threshold > 50:
                        Forest_Species= 'Yes'
                    else:
                        Forest_Species= 'No'
                else:
                    TC_threshold=0
                    Forest_Species= 'No'
                    
                
                #mask forest raster by species alphahull
                
                TC_in_hull=mask.mask(src_forest,alphahull.geometry,crop=True)
                array= TC_in_hull[0] #get the array, a matrix of raster values
                affine= TC_in_hull[1] #get info from the raster
                #Forest trimmed by a threshold within the alphahull
                TC_by_spp=(array>= TC_threshold).astype(int) #subset the array after a condition
                
                
                
                print ('2. Elevation')
                ##### ELEVATION #####################################################
                PresElevation=presences.ELEVATION
                
                E_threshold= pd.Series.quantile(PresElevation, [0,1])
                LE=int(round(E_threshold[0]))
                UE=int(round(E_threshold[1]))
                
                Min=LE
                Max=UE
                  
                  
                #Set elevations to be used
                if df['elevation_upper'][b] == 0 and df['elevation_lower'][b] == 0:
                    Min=LE
                    Max=UE 
               
                if ~df['Min extreme altitude (m)'].isna()[b]:
                    BLmin=int(round(df['Min extreme altitude (m)'][b]))
                    Min= min(BLmin, LE)  
                else:
                    BLmin=int(round(df['elevation_lower'][b]))
                    Min= min(BLmin, LE)
                if  ~df['Max extreme altitude (m)'].isna()[b]:
                    BLmax=int(round(df['Max extreme altitude (m)'][b]))
                    Max= max(BLmax, UE)
                else:
                    BLmax=int(round(df['elevation_upper'][b]))
                    Max= max(BLmax, UE)
                
                
                #mask elevation raster by species alphahull
                E_in_hull=mask.mask(src_dem,alphahull.geometry,crop=True)
                array_E= E_in_hull[0].astype(int) #get the array, a matrix of raster values
                #Elevation trimmed by a threshold within the alphahull
                if np.isnan(Min):
                    E_by_spp=(array_E<=Max).astype(int) #subset the array after a condition
                if np.isnan(Max):
                    E_by_spp=(array_E<=Min).astype(int) #subset the array after a condition
                if np.isnan(Min) and np.isnan(Max):
                    E_by_spp=array_E
                if ~np.isnan(Min) and ~np.isnan(Max):
                    E_by_spp= ((array_E<=Max) & (array_E>=Min)).astype(int)
                    
                    
                #Habitat 2000
                AOH_path=output1
                AOH= TC_by_spp*E_by_spp
                with rasterio.open(
                        AOH_path,
                        'w', driver='GTiff',
                        height=AOH.shape[1],width=AOH.shape[2],
                        count=1,dtype=str(AOH.dtype),
                        crs=epsg_WGS84,transform=affine, nodata=0
                        ) as new_dataset:
                    new_dataset.write(AOH)
                new_dataset.close()
                
                size_raster=len(AOH[AOH==1])
                AOH_area=round((size_raster*cell_size_proj**2)*m2) 
                #########
                
                
                arcpy.env.extent = hull_path
                arcpy.RasterToPolygon_conversion(AOH_path, AOH_shp, "NO_SIMPLIFY")
            
                ######################################################################
                ######################## ABSENCES ####################################
                ######################################################################
                if Checklists:
                    print ('Subsetting Checklists')
                    Cckl_spp=gpd.sjoin(se_df,alphahull,'inner','intersects').reset_index(drop=True) #Get Checklists within alpha hull
                    
                    if not Cckl_spp.empty:
                        Cckl_spp=Cckl_spp.drop('index_right', axis=1) #Get rid of undesired column
                        Cckl_spp.to_file(Cckl_spp_dir, driver='ESRI Shapefile', encoding = 'utf-8') #Save new ChecklistsCckl_spp_hab=  #Checklists in habitat path file
                        arcpy.SpatialJoin_analysis(Cckl_spp_dir, AOH_shp, Cckl_spp_hab_dir, 
                                                   '', 'KEEP_COMMON', '', 'INTERSECT', float(cell_size)) # Intersecting checklists within a radius distance of the cell size of habitat
                        Cckl_spp_hab2=gpd.read_file(Cckl_spp_hab_dir) #Checklists intersected
                        
                        #Create grid feature with Habitat raster as extent
                        Cckl_spp_grid=gpd.sjoin(g,Cckl_spp_hab2,'inner','intersects').reset_index(drop=True)
                        Cckl_spp_grid['Count'] = Cckl_spp_grid.groupby('PageName')['PageName'].transform('count')
                        Cckl_spp_grid25=Cckl_spp_grid[Cckl_spp_grid.Count>=25]
                        Cckl_spp_grid25=Cckl_spp_grid25.drop('index_right', axis=1) #Get rid of undesired column
                        Cckl_25inpres= gpd.sjoin(Cckl_spp_grid25,presences,'inner','intersects').reset_index(drop=True)
                        Cckl_25outpres=(Cckl_spp_grid25[~Cckl_spp_grid25.PageName.isin(Cckl_25inpres.PageName)])
                        
                        if not Cckl_25outpres.empty:
                            Cckl_25outpres.to_file(Non_detection, driver='ESRI Shapefile', encoding = 'utf-8')
                            Cckl_25outpres_copy=Cckl_25outpres.copy()
                            Cckl_25outpres_copy.geometry=Cckl_25outpres.geometry.centroid
                            Cckl_25outpres_copy.to_file(Absences, driver='ESRI Shapefile', encoding = 'utf-8')
                            
                            #Get Checklists from Non_detection pixels
                            Select_abs=arcpy.SelectLayerByLocation_management(Cckl_spp_hab_dir, 'INTERSECT', Non_detection)
                            arcpy.CopyFeatures_management(Select_abs, Chcklists_pixel)
                            
                            #Reading absences
                            Abs=Cckl_25outpres_copy[['PageName','PageNumber','ID','Count','geometry']]
                            Abs=Abs.drop_duplicates()
                            Abs_r=len(Abs)
                            
                if not Checklists or Cckl_spp.empty or Cckl_25outpres.empty:
                    print('No Absences')
                    Abs_r=None
                    PO_AOH_area=None
                    AU_AOH_area=None
                    
            
                else:
                    presences['ID']=0
                    Pres=presences[['ID','geometry']]
                    Pres_Abs=pd.concat([Abs,Pres])
                    Pres_Abs.to_file(Pres_Abs_path, driver='ESRI Shapefile', encoding = 'latin1')
                    
                    ######################################################################
                    ######################## INTERPOLATION ###############################
                    ######################################################################
                    
                    print ('Interpolating')
                    arcpy.env.extent = AOH_path
                    IDW=arcpy.Idw_3d(Pres_Abs_path, 'ID', idw_spp, AOH_path, 1, 'VARIABLE 1')
                    
                    IDW_array=rasterio.open(IDW[0]).read().astype(int)
                   
                    #Apparently Unoccupied AOH early
                    Unoccupied_AOH= IDW_array*AOH
                    with rasterio.open(
                            Unoccupied_AOH_spp,
                            'w', driver='GTiff',
                            height=Unoccupied_AOH.shape[1],width=Unoccupied_AOH.shape[2],
                            count=1,dtype=str(Unoccupied_AOH.dtype),
                            crs=epsg_WGS84,transform=affine, nodata=0
                            ) as new_dataset:
                        new_dataset.write(Unoccupied_AOH)
                    new_dataset.close()
                    
            
                    size_AU_AOH= len(Unoccupied_AOH[Unoccupied_AOH==1])
                    AU_AOH_area=round((size_AU_AOH*cell_size_proj**2)*m2)
                    
                    #Potentially Occupied AOH Early
                    Potential_AOH= AOH-Unoccupied_AOH
                    with rasterio.open(
                            Potential_AOH_spp,
                            'w', driver='GTiff',
                            height=Potential_AOH.shape[1],width=Potential_AOH.shape[2],
                            count=1,dtype=str(Potential_AOH.dtype),
                            crs=epsg_WGS84,transform=affine, nodata=0
                            ) as new_dataset:
                        new_dataset.write(Potential_AOH)
                    new_dataset.close()
                    
                    size_PO_AOH= len(Potential_AOH[Potential_AOH==1])
                    PO_AOH_area=round((size_PO_AOH*cell_size_proj**2)*m2)
                    
                ######################################################################
                ######################## END 1 #######################################
                ######################################################################
            else:
                print('No Presences')
                AOH_area= None
                TC_threshold= None
                Forest_Species= None
                GD= None
                Abs_r= None
                PO_AOH_area=  None
                AU_AOH_area=None
                number_presences= None
                Min= None
                Max= None
            #SAVING DATA TO A REPORT TABLE AND EXPORT TO CSV
            report_Habitat.append([c,
                                   AOH_area, 
                                   PO_AOH_area,
                                   AU_AOH_area,
                                   TC_threshold,Forest_Species,GD,Abs_r,
                                   number_presences, Min,Max])
            Report_Table=pd.DataFrame(report_Habitat,columns=['Scientific Name',
                                                              'Area of Habitat (Km\u00b2)',
                                                              'Potential AOH (Km\u00b2)',
                                                              'Unnocupied AOH (Km\u00b2)',
                                                              'Forest Threshold', 'Forest Species',
                                                              'Geographical Distribution','Abs',
                                                              'Presence Records', 'Lower Elev', 'Upper Elev'])
            Report_Table.to_csv(Report,encoding='latin-1', index=False)
            print('******* ' + c +' ******* done!\n')
            arcpy.ResetEnvironments()
            
        
                
            
        except ValueError as e:
            report_Errors.append([c,e])
            RErrors=pd.DataFrame(report_Errors,columns=['Scientific Name','Error'])
            RErrors.to_csv(errors,encoding='latin-1', index=False)
    print ('Complete!')
    return Report_Table
    
        
    

