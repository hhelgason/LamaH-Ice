// Code obtained from Christoph Klingler Institute for Hydrology and Water Management, University of Natural Resources and Life Sciences, Vienna
// Adapted to Iceland by Hordur Helgason, University of Washington, 2022. Code accompanying the paper LamaH-Ice: LArge-SaMple DAta for Hydrology and Environmental Sciences for ICEland 
// This script runs on the Google Earth Engine, downloads MODIS LAI rasters

// Set point  for print purposes
var P1 = /* color: #d63000 */ee.Geometry.Point([-18.5, 64.8]);

// Import area of interest (aoi) as Feature Collection,
// aoi has to be uploaded first to GEE as shapefile (polygon which covers aoi) to Assests
var aoi = ee.FeatureCollection("users/helgason/ERA5L_grid");
            
// Set view to center of Raster
Map.centerObject(aoi, 7);

// Show aoi
Map.addLayer(aoi, {color: 'blue'},'AOI');

// Set needed Time Window and Time Step from Dataset
var start_dt = ee.Date("2012-01-01");
var end_dt = ee.Date("2012-12-31");

// Select dataset from GEE Data_Catalog
var dataset = ee.ImageCollection('MODIS/006/MCD12Q1')
  // Set timewindow, add 1 day to end date to make it inclusive
  .filterDate(start_dt, end_dt.advance(1, 'day'))
  // Selection of the needed bands
  .select('LC_Type1');
  
// Select first Image in ImageCollection
var image = dataset.first();
  
// Write to console
print(dataset, 'LC - IGBP');

// Check projektion before reprojection
print('CRS bev', image.projection());
// Check Scale before reprojection
print('Scale [m] bef', image.projection().nominalScale());

// Reproject from SR-ORG:6974 to EPSG:4326
var imager = image.reproject('EPSG:3057', null, 500);

// Check projektion after reprojection
print('CRS aft', imager.projection());
// Check Scale after reprojection
print('Scale [m] aft', imager.projection().nominalScale());

// Export the image LC, specifying scale and region.
Export.image.toDrive({
  image: imager,
  description: 'LC_IGBP_2012_epsg3057',
  scale: 500,
  region: aoi,
  fileFormat: 'GeoTIFF'
});