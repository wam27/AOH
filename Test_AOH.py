# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:31:17 2022

@author: wam27
"""
import os, arcpy
from arcpy import env
import arcpy.sa
import pandas as pd
import geopandas as gpd
import warnings
from AOH_base_v2 import AOH

#======================================================================================================
#//////////////////////////////////////////////////////////////////////////////////////////////////////

################# ENVIRONMENT #########################################################################
#======================================================================================================

#Set work environment
warnings.filterwarnings("ignore")
#Chek out ArcGis Analyst and allow overwrite files
arcpy.CheckOutExtension('Spatial')
arcpy.env.overwriteOutput = True
#Set relative paths
#scriptWS= os.path.basename(sys.argv[0])
rootWS= os.getcwd()
dataWS=os.path.join(rootWS,'Data' ) 
docsWS=os.path.join(rootWS,'Docs' ) 
scratchWS=os.path.join(rootWS, 'Scratch_Occurrences\Bird movements')
#Set work environment
env.workspace= dataWS
#Call variables
os.chdir(rootWS)
#======================================================================================================
#//////////////////////////////////////////////////////////////////////////////////////////////////////

################# VARIABLES ###########################################################################
#======================================================================================================

df=pd.read_csv(docsWS + "\spp_NA3.csv", encoding='latin-1')


dem = dataWS + r"\DEM\DEM_NorthAndes.tif"
dem_proj = dataWS + r"\DEM\DEM_NorthAndes_proj.tif"
se=dataWS + "\eBird\\CL\\Checklists_dc.shp"
tc = dataWS + r"\TreeCover\TC_2000_90m.tif"

GREAT_GRID= dataWS + '\eBird\GREAT_GRID.shp'
g= gpd.read_file(GREAT_GRID)
NCckl=25
n='_early'

occ= scratchWS + '\Pres\\'
hull= scratchWS + '\Hull\\'

AOH_output= scratchWS + '\Areas\\'
temp_output= scratchWS + '\Temp\\'

        
AOH1=AOH(df, n, occ, hull, dem, dem_proj, tc, g, 
        AOH_output, temp_output, se, NCckl)
            
    