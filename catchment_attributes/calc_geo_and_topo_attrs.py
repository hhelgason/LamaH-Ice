### OBJECT
# Fill empty Polygon-shapefiles with attribute fields and
# calculate some topographic and geologic attributes

### INFO
# TODO: the geological glim classes _______ are not present in the area of interest, and therefore aren´t reflected in the script
# calculation may takes a while, depending on the number of basins
# the code presuppose that all imported shapefiles are already in the same coordinate system

### AUTHOR
# This code accompanies the paper "LamaH-Ice | Large-Sample Data for Hydrology and Environmental Sciences for Iceland" submitted to the journal Earth Syst. Sci. Data (ESSD), 2022
# The code was originally written by Christoph Klingler, Institute for Hydrology and Water Management, University of Natural Resources and Life Sciences, Vienna, 11 June 2021, v1.0
# Originally written by Christoph Klingler, Institute for Hydrology and Water Management, University of Natural Resources and Life Sciences, Vienna, 11 June 2021, v1.0, accompanying the paper "LamaH-CE | Large-Sample Data for Hydrology and Environmental Sciences for Central Europe" published in the journal Earth Syst. Sci. Data (ESSD), 2021
# The original code has been was adapted to QGIS v. 3.14 and updated to create the LamaH-Ice dataset, by Hordur Helgason, University of Washington, Seattle, June 30th 2022
# This version includes aggregation of CORINE land cover classes to basin shapefiles, as well as more detailed geological information.

##########
### IMPORT
import math
import processing
from qgis.core import QgsVectorLayer
import numpy as np

from qgis.core import (
  QgsApplication,
  QgsDataSourceUri,
  QgsCategorizedSymbolRenderer,
  QgsClassificationRange,
  QgsPointXY,
  QgsProject,
  QgsExpression,
  QgsField,
  QgsFields,
  QgsFeature,
  QgsFeatureRequest,
  QgsFeatureRenderer,
  QgsGeometry,
  QgsGraduatedSymbolRenderer,
  QgsMarkerSymbol,
  QgsMessageLog,
  QgsRectangle,
  QgsRendererCategory,
  QgsRendererRange,
  QgsSymbol,
  QgsVectorDataProvider,
  QgsVectorLayer,
  QgsVectorFileWriter,
  QgsWkbTypes,
  QgsSpatialIndex,
  QgsVectorLayerUtils
)

# load shapefiles and add them to the map/registry
basins = QgsVectorLayer(r"C:\Users\hordurbhe\Dropbox\UW\lamah_ice\GIS\watersheds\final_watersheds\final_watersheds\Basins_A.shp","basins")

QgsProject.instance().addMapLayer(basins, True)
gauges = QgsVectorLayer(r"C:\Users\hordurbhe\Dropbox\UW\lamah_ice\GIS\watersheds\final_watersheds\gauges_with_splitted_included.shp","gauges")
QgsProject.instance().addMapLayer(gauges, True)
strmnet = QgsVectorLayer(r"C:\Users\hordurbhe\Dropbox\UW\lamah_ice\GIS\EU-Hydro\River_Net_I_isn_93.shp", "EU-Hydro_network")
QgsProject.instance().addMapLayer(strmnet, True)
glim = QgsVectorLayer(r"C:\Users\hordurbhe\Dropbox\UW\lamah_ice\GIS\glim_geology\glim_isn93_clipped_gpds\glim_isn93_clipped_gpds.shp", "GLiM")
QgsProject.instance().addMapLayer(glim, True)
NI_geo = QgsVectorLayer(r"C:\Users\hordurbhe\Dropbox\UW\lamah_ice\GIS\glim_geology\NI_geology_isn_93.shp", "NI_geo")
QgsProject.instance().addMapLayer(NI_geo, True)
corine = QgsVectorLayer(r"C:\Users\hordurbhe\Dropbox\UW\lamah_ice\GIS\CORINE\CLC18_isn93.shp", "corine")
QgsProject.instance().addMapLayer(corine, True)

#################
### PREPROCESSING

# set crs to crs of layer "basin"
crs = basins.crs().toWkt()

# create list from attribute "ID" of layer basins
ID_list = [i.attribute("id") for i in basins.getFeatures()]

# order list by ascending value
ID_list.sort() 

# create empty lists
filter_list = []

# create a list with geological attribute classes from GLiM (Layer 1)
code_xx = ['ig', 'mt', 'pa', 'pb', 'pi', 'py', 'sc', 'sm', 'ss', 'su', 'va', 'vb', 'wb'] # Here we include some classes that are not present in Iceland
# create a list with more specific geological classes from GLiM (Layers 2 and 3)
code_glim_litho = ['pa__vr','pb____','vapy__','vb____','vb__sr','vbpy__'] # Here we only include the 6 classes present in Iceland
# create lists with geological attribute classes from the Icelandic NI dataset
code_NIgeo_category = ['621', '701', '743', 'binn', 'bnew', 'bold', 'gnew', 'gold', 'hraun', 'mob', 'sgos', 'sinn', 'sn','snew','sold']
# 
corine_categories = ['111','112','121','122','123','124','131','132','133','141','142','211','212','213','221','222','223','231','241','242','243','244','311','312','313','321','322','323','324','331','332','333','334','335','411','412','421','422','423','511','512','521','522','523']

# loop for joining "ID=" to ID_list to filter layers
for i in ID_list:
    string = "id" + "=" + str(i)
    filter_list.append(string)


########################################################
### LOOP THROUGH ALL BASINS TO CALCULATE SOME ATTRIBUTES

# start loop
for idx, j in enumerate(filter_list):

    # print loop status
    print(j)
    
    # filter vector layer and select them
    basins.setSubsetString(j) # filter 
    basins.selectByExpression(j) # selection
    gauges.setSubsetString(j) # filter
    gauges.selectByExpression(j) # selection
        
    ## Calculate some topographic attributes

    # create list with coordinates of x/y from all vertices of layer basin
    vertl = basins.selectedFeatures()[0].geometry().asMultiPolygon()
    vertl = [i for i in vertl[0][0]]

    # get coordinates of gauge
    gc = gauges.selectedFeatures()[0].geometry().asPoint()

    # calculate distance from every single vertex to the gauge in [m]
    dist = [((i.x()-gc.x())**2 + (i.y()-gc.y())**2)**0.5 for i in vertl]

    # get index of vertex, which is most distant to the gauge 
    vertidx = dist.index(max(dist))

    # get maximum vertex distance to the gauge
    L = dist[vertidx]

    # calculate angle from most distant vertex [deg]
    # flow direction from north to south: 180 deg
    # flow direction from east to west: 270 deg
    angle = math.degrees(math.atan(abs(gc.x()-vertl[vertidx].x()) / abs(gc.y()-vertl[vertidx].y()) ))
    if gc.x()-vertl[vertidx].x() > 0 and gc.y()-vertl[vertidx].y() > 0:
        angle = angle
    elif gc.x()-vertl[vertidx].x() > 0 and gc.y()-vertl[vertidx].y() < 0:
        angle = 180 - angle 
    elif gc.x()-vertl[vertidx].x() < 0 and gc.y()-vertl[vertidx].y() < 0:
        angle = 180 + angle 
    elif gc.x()-vertl[vertidx].x() < 0 and gc.y()-vertl[vertidx].y() > 0:
        angle = 360 - angle
        
    # get area of basin in [m2]
    area = basins.selectedFeatures()[0].geometry().area()
            
    # calculate elongation ratio after Schumm
    Re = 1/L * (4*area/math.pi)**0.5

    ## Calculate stream density of basin
    
    # clip river network to filtered basin
    strmclip = processing.run("native:clip", {'INPUT':strmnet,'OVERLAY':QgsProcessingFeatureSourceDefinition(basins.id(), True),'OUTPUT':'memory:'})["OUTPUT"]

    # create list with stream-lengths [m] of layer strmclip
    streamlengths = [i.geometry().length() for i in strmclip.getFeatures()]

    # sum up the stream segment lengths [m]
    streamlen = sum(streamlengths)

    # calculate stream density [m/m2]
    D = streamlen/area
    
    ## Calculate geologic attributes
    
    # create empty lists
    filter_list_glim = []
    glim_area = []
    
    # intersect glim to filtered basin
    glimisec = processing.run("native:intersection", {'INPUT':glim,'OVERLAY':basins,'INPUT_FIELDS':['xx'],'OVERLAY_FIELDS':['ID'],'OUTPUT':'memory:'})["OUTPUT"]
    
    # loop for joining area to list
    for i in code_xx:
        glimisec.selectByExpression('"xx"=\'%s\'' % i)
        if len(glimisec.selectedFeatures()) == 0:
            glim_area_i = 0
        else:
            glim_area_sec = []
            for feat in glimisec.selectedFeatures():
                glim_area_sec.append(feat.geometry().area())
            glim_area_i = sum(glim_area_sec)
        glim_area.append(glim_area_i)

    # get value and index of list, where area is maximum
    glim_area_max = max(glim_area)
    glim_area_max_idx = glim_area.index(glim_area_max)

    # get geologic class, with highest area share
    glim_area_max_xx = code_xx[glim_area_max_idx]

    # get fractions for all geological classes
    fra_gc_ig = glim_area[0]/area
    fra_gc_mt = glim_area[1]/area
    fra_gc_pa = glim_area[2]/area
    fra_gc_pb = glim_area[3]/area
    fra_gc_pi = glim_area[4]/area
    fra_gc_py = glim_area[5]/area
    fra_gc_sc = glim_area[6]/area
    fra_gc_sm = glim_area[7]/area
    fra_gc_ss = glim_area[8]/area
    fra_gc_su = glim_area[9]/area
    fra_gc_va = glim_area[10]/area
    fra_gc_vb = glim_area[11]/area
    fra_gc_wb = glim_area[12]/area
    
    ## Calculate the more specific attributes from GLiM
    # create empty lists
    filter_list_glim_litho = []
    glim_area_litho = []
    
    # intersect glim to filtered basin
    glimisec_litho = processing.run("native:intersection", {'INPUT':glim,'OVERLAY':basins,'INPUT_FIELDS':['Litho'],'OVERLAY_FIELDS':['ID'],'OUTPUT':'memory:'})["OUTPUT"]
    
    # loop for joining area to list
    for i in code_glim_litho:
        glimisec_litho.selectByExpression('"Litho"=\'%s\'' % i)
        if len(glimisec_litho.selectedFeatures()) == 0:
            glim_area_i = 0
        else:
            glim_area_sec_litho = []
            for feat in glimisec_litho.selectedFeatures():
                glim_area_sec_litho.append(feat.geometry().area())
            glim_area_i = sum(glim_area_sec_litho)
        glim_area_litho.append(glim_area_i)

    # get value and index of list, where area is maximum
    glim_litho_area_max = max(glim_area_litho)
    glim_litho_area_max_idx = glim_area_litho.index(glim_litho_area_max)

    # get geologic class, with highest area share
    glim_litho_area_max_xx = code_glim_litho[glim_litho_area_max_idx]

    # get fractions for all geological classes
    fra_litho_pavr = glim_area_litho[0]/area
    fra_litho_pb = glim_area_litho[1]/area
    fra_litho_vapy = glim_area_litho[2]/area
    fra_litho_vb = glim_area_litho[3]/area
    fra_litho_vbsr = glim_area_litho[4]/area
    fra_litho_vbpy = glim_area_litho[5]/area
 
    ## Calculate geologic attributes from the NI Geology dataset
    
    # create empty lists
    filter_list_NIgeo = []
    NIgeo_area = []
    
    # intersect glim to filtered basin
    NIgeoisec = processing.run("native:intersection", {'INPUT':NI_geo,'OVERLAY':basins,'INPUT_FIELDS':['FLOKKUR_1'],'OVERLAY_FIELDS':['ID'],'OUTPUT':'memory:'})["OUTPUT"]
    
    # loop for joining area to list
    for i in code_NIgeo_category:
        NIgeoisec.selectByExpression('"FLOKKUR_1"=\'%s\'' % i)
        if len(NIgeoisec.selectedFeatures()) == 0:
            NIgeo_area_i = 0
        else:
            NIgeo_area_sec = []
            for feat in NIgeoisec.selectedFeatures():
                NIgeo_area_sec.append(feat.geometry().area())
            NIgeo_area_i = sum(NIgeo_area_sec)
        NIgeo_area.append(NIgeo_area_i)

    # get value and index of list, where area is maximum
    NIgeo_area_max = max(NIgeo_area)
    NIgeo_area_max_idx = NIgeo_area.index(NIgeo_area_max)

    # get geologic class, with highest area share
    NIgeo_area_max_category = code_NIgeo_category[NIgeo_area_max_idx]

    # get fractions for all geological classes
    fra_gc_621 = NIgeo_area[0]/area
    fra_gc_701 = NIgeo_area[1]/area
    fra_gc_743 = NIgeo_area[2]/area
    fra_gc_binn = NIgeo_area[3]/area
    fra_gc_bnew = NIgeo_area[4]/area
    fra_gc_bold = NIgeo_area[5]/area
    fra_gc_gnew = NIgeo_area[6]/area
    fra_gc_gold = NIgeo_area[7]/area
    fra_gc_hraun = NIgeo_area[8]/area
    fra_gc_mob = NIgeo_area[9]/area
    fra_gc_sgos = NIgeo_area[10]/area
    fra_gc_sinn = NIgeo_area[11]/area
    fra_gc_sn = NIgeo_area[12]/area
    fra_gc_snew = NIgeo_area[13]/area
    fra_gc_sold = NIgeo_area[14]/area
    
    ## Calculate land cover attributes from the CORINE Land Cover (CLC) dataset
    # create empty lists
    filter_list_clc = []
    clc_area = []
    
    # intersect CLC to filtered basin
    clc_isec = processing.run("native:intersection", {'INPUT':corine,'OVERLAY':basins,'INPUT_FIELDS':['CODE_18'],'OVERLAY_FIELDS':['ID'],'OUTPUT':'memory:'})["OUTPUT"]
    
    # loop for joining area to list
    for i in corine_categories:
        clc_isec.selectByExpression('"CODE_18"=\'%s\'' % i)
        if len(clc_isec.selectedFeatures()) == 0:
            clc_area_i = 0
        else:
            clc_area_sec = []
            for feat in clc_isec.selectedFeatures():
                clc_area_sec.append(feat.geometry().area())
            clc_area_i = sum(clc_area_sec)
        clc_area.append(clc_area_i)

    # get value and index of list, where area is maximum
    clc_area_max = max(clc_area)
    clc_area_max_idx = clc_area.index(clc_area_max)

    # get corine class, with highest area share
    clc_area_max_category = int(corine_categories[clc_area_max_idx])

    # get fractions for specified classes
     
    # fraction of agricultural areas (all categories starting with 2)
    # We first define the indices of ag codes in the corine_categories list
    ag_code_idx = [corine_categories.index(i) for i in corine_categories if i.startswith('2')]
    agr_fra = sum([clc_area[i] for i in ag_code_idx])/area
    # Fraction of bare areas (classes 332, 333)
    bare_fra = sum([clc_area[corine_categories.index(i)] for i in ['332','333']])/area
    #bare_fra = np.round(bare_fra,3)
    # Fraction of forest areas (classes 311,312,313)
    forest_fra = sum([clc_area[corine_categories.index(i)] for i in ['311','312','313']])/area
    # Fraction of glaciers (class 335)
    glac_fra = clc_area[corine_categories.index('335')]/area
    # Fraction of water bodies (class 512)
    lake_fra = clc_area[corine_categories.index('512')]/area
    # Fraction of urban areas (classes 111,112,121,122,123,124)
    urban_fra = sum([clc_area[corine_categories.index(i)] for i in ['111','112','121','122','123','124']])/area
    
##########################
### CHANGE ATTRIBUTE TABLE
    
    # start editing of layer
    basins.startEditing()
    
    
    ## Add fields to attribute table
    
    # append fields in first loop step
    if idx == 0:
        
        # topographic attributes (10)
        basins.addAttribute(QgsField("area_calc", QVariant.Double, len=10, prec=3)) # km2
        basins.addAttribute(QgsField("elev_mean", QVariant.Int, len=4, prec=0))
        basins.addAttribute(QgsField("elev_med", QVariant.Int, len=4, prec=0))
        basins.addAttribute(QgsField("elev_std", QVariant.Int, len=4, prec=0))
        basins.addAttribute(QgsField("elev_ran", QVariant.Int, len=4, prec=0))
        basins.addAttribute(QgsField("slope_mean", QVariant.Int, len=4, prec=0)) # m/km
        basins.addAttribute(QgsField("mvert_dist", QVariant.Double, len=6, prec=1)) # km
        basins.addAttribute(QgsField("mvert_ang", QVariant.Int, len=3, prec=0)) # deg
        basins.addAttribute(QgsField("elon_ratio", QVariant.Double, len=5, prec=3)) # -
        basins.addAttribute(QgsField("strm_dens", QVariant.Double, len=4, prec=2))# km/km2
        
        # climate attributes (11)
        basins.addAttribute(QgsField("p_mean", QVariant.Double, len=4, prec=2))
        basins.addAttribute(QgsField("et0_mean", QVariant.Double, len=4, prec=2))
        basins.addAttribute(QgsField("aridity", QVariant.Double, len=4, prec=2))
        basins.addAttribute(QgsField("p_season", QVariant.Double, len=5, prec=2))
        basins.addAttribute(QgsField("frac_snow", QVariant.Double, len=4, prec=2))
        basins.addAttribute(QgsField("hi_prec_fr", QVariant.Double, len=5, prec=2))
        basins.addAttribute(QgsField("hi_prec_du", QVariant.Double, len=4, prec=2))
        basins.addAttribute(QgsField("hi_prec_ti", QVariant.String, len=3))
        basins.addAttribute(QgsField("lo_prec_fr", QVariant.Double, len=6, prec=2))
        basins.addAttribute(QgsField("lo_prec_du", QVariant.Double, len=4, prec=2))
        basins.addAttribute(QgsField("lo_prec_ti", QVariant.String, len=3))
        
        # land class attributes (7)
        basins.addAttribute(QgsField("lc_dom", QVariant.Int, len=3, prec=0))
        basins.addAttribute(QgsField("agr_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("bare_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("forest_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("glac_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("lake_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("urban_fra", QVariant.Double, len=5, prec=3))
        
        # vegetation attributes (6)
        basins.addAttribute(QgsField("lai_max", QVariant.Double, len=4, prec=2))
        basins.addAttribute(QgsField("lai_diff", QVariant.Double, len=4, prec=2))
        basins.addAttribute(QgsField("ndvi_max", QVariant.Double, len=4, prec=2))
        basins.addAttribute(QgsField("ndvi_min", QVariant.Double, len=4, prec=2))
        basins.addAttribute(QgsField("gvf_max", QVariant.Double, len=4, prec=2))
        basins.addAttribute(QgsField("gvf_diff", QVariant.Double, len=4, prec=2))
        
        # soil attributes (10)
        basins.addAttribute(QgsField("bedrk_dep", QVariant.Double, len=5, prec=2)) # m
        basins.addAttribute(QgsField("root_dep", QVariant.Double, len=4, prec=2)) # m
        basins.addAttribute(QgsField("soil_poros", QVariant.Double, len=4, prec=2)) # -
        basins.addAttribute(QgsField("soil_condu", QVariant.Double, len=5, prec=3)) # /(100*24)
        basins.addAttribute(QgsField("soil_tawc", QVariant.Double, len=4, prec=2)) # m
        basins.addAttribute(QgsField("sand_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("silt_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("clay_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("grav_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("oc_fra", QVariant.Double, len=5, prec=3))
        
        # geological attributes from GLiM (16)
        basins.addAttribute(QgsField("gc_dom", QVariant.String, len=2))
        basins.addAttribute(QgsField("gc_ig_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_mt_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_pa_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_pb_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_pi_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_py_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_sc_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_sm_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_ss_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_su_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_va_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_vb_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gc_wb_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("geol_perme", QVariant.Double, len=5, prec=1))
        basins.addAttribute(QgsField("geol_poros", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("geo_perme2", QVariant.Double, len=5, prec=1))
        basins.addAttribute(QgsField("geo_poros2", QVariant.Double, len=5, prec=3))
        
        # more specific geological attributes from GLiM (layers 2 and 3)
        basins.addAttribute(QgsField("litho_dom", QVariant.String, len=2))
        basins.addAttribute(QgsField("litho_pavr", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("litho_pb", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("litho_vapy", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("litho_vb", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("litho_vbsr", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("litho_vbpy", QVariant.Double, len=5, prec=3))
        
        # geological attributes from Icelandic dataset(15)
        basins.addAttribute(QgsField("g_dom_NI", QVariant.String, len=2))
        basins.addAttribute(QgsField("g621_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("g701_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("g743_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gbinn_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gbnew_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gbold_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("ggnew_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("ggold_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("ghraun_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gmob_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gsgos_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gsinn_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gsn_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gsnew_fra", QVariant.Double, len=5, prec=3))
        basins.addAttribute(QgsField("gsold_fra", QVariant.Double, len=5, prec=3))
        
    ## Add results to attribute table

    # transfer calculated attribute values to fields
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("area_calc"), area/1000000) 
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("mvert_dist"), L/1000) 
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("mvert_ang"), angle) 
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("elon_ratio"), Re) 
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("strm_dens"), D*1000)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_dom"), glim_area_max_xx)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_ig_fra"), fra_gc_ig)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_mt_fra"), fra_gc_mt)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_pa_fra"), fra_gc_pa)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_pb_fra"), fra_gc_pb)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_pi_fra"), fra_gc_pi)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_py_fra"), fra_gc_py)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_sc_fra"), fra_gc_sc)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_sm_fra"), fra_gc_sm)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_ss_fra"), fra_gc_ss)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_su_fra"), fra_gc_su)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_va_fra"), fra_gc_va)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_vb_fra"), fra_gc_vb)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gc_wb_fra"), fra_gc_wb)
    # More specific geological attributes from GLiM layers 2 and 3
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("litho_dom"), glim_litho_area_max_xx)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("litho_pavr"), fra_litho_pavr)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("litho_pb"), fra_litho_pb)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("litho_vapy"), fra_litho_vapy)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("litho_vb"), fra_litho_vb)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("litho_vbsr"), fra_litho_vbsr)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("litho_vbpy"), fra_litho_vbpy)
    # Geological attributes from the Icelandic dataset:
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("g_dom_NI"), NIgeo_area_max_category)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("g621_fra"), fra_gc_621)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("g701_fra"), fra_gc_701)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("g743_fra"), fra_gc_743)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gbinn_fra"), fra_gc_binn)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gbnew_fra"), fra_gc_bnew)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gbold_fra"), fra_gc_bold)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("ggnew_fra"), fra_gc_gnew)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("ggold_fra"), fra_gc_gold)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("ghraun_fra"), fra_gc_hraun)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gmob_fra"), fra_gc_mob)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gsgos_fra"), fra_gc_sgos)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gsinn_fra"), fra_gc_sinn)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gsn_fra"), fra_gc_sn)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gsnew_fra"), fra_gc_snew)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("gsold_fra"), fra_gc_sold)
    # Land cover attributes
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("lc_dom"), clc_area_max_category)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("agr_fra"), agr_fra)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("bare_fra"), bare_fra)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("forest_fra"), forest_fra)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("glac_fra"), glac_fra)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("lake_fra"), lake_fra)
    basins.changeAttributeValue(basins.selectedFeatures()[0].id(), basins.fields().indexFromName("urban_fra"), urban_fra)  
    
    # end editing, save changes of layer
    basins.commitChanges()
        
    # unfilter layers
    basins.setSubsetString("")
    gauges.setSubsetString("")
    
    # remove selection of layers
    basins.removeSelection()
    gauges.removeSelection()

# delete unused temporary layers from registry
legend_layers = [i.layer() for i in QgsProject.instance().layerTreeRoot().children()]
registry_layers = QgsProject.instance().mapLayers().values()
for i in registry_layers:
    if not i in legend_layers:
        QgsProject.instance().removeMapLayer(i.id())

# print status end
print("Finished")