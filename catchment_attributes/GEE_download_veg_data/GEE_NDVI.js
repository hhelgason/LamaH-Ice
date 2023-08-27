// Code obtained from Christoph Klingler Institute for Hydrology and Water Management, University of Natural Resources and Life Sciences, Vienna
// Adapted to Iceland by Hordur Helgason, University of Washington, 2022. Code accompanying the paper LamaH-Ice: LArge-SaMple DAta for Hydrology and Environmental Sciences for ICEland
// This script runs on the Google Earth Engine, downloads NDVI rasters based on MODIS observations 

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
var start_dt = ee.Date("2000-04-01");
var end_dt = ee.Date("2022-03-31");

// Function to mask clouds based on the 'FparLai_QC' band of Modis data
var qualitymask = function(img) {
    var qa = img.select('State');
    // Check if images were produced at cloudy atm. conditions 
    // Bit 0 is true and Bit 1 is false --> = 1
    var qual = qa.bitwiseAnd(1 << 0)
               .and(qa.bitwiseAnd(0 << 1));
    // Return images, which were produced not at cloudy conditions
    return img.updateMask(qual.not());
};

// Select dataset from GEE Data_Catalog
var dataset = ee.ImageCollection('MODIS/061/MOD09Q1')  
  // Set timewindow, add 1 day to end date to make it inclusive
  .filterDate(start_dt, end_dt.advance(1, 'day'))
  // Selection of the needed bands
  .select(['sur_refl_b01', 'sur_refl_b02', 'State'])
  // Map function cloudmask over ImageCollection
  .map(qualitymask);
  
// Function to calculate and add an NDVI band to a single image
// NDVI = (NIR - Red) / (NIR + Red)
var addNDVI = function(image) {
  var NDVI = image.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']).rename('NDVI');
  return image.addBands(NDVI);
};

// Map the function addNDVI over the imagecollection
var datasetn = dataset.map(addNDVI);

// Create list with numbers 1 to 12
var months = ee.List.sequence(1, 12);

// Group by month, and then reduce within groups by mean(),
// the result is an ImageCollection with one image for each month;
var byMonth = ee.ImageCollection.fromImages(
    months.map(function (m) {
    return datasetn
      .select('NDVI')
      .filter(ee.Filter.calendarRange(m, m, 'month'))
      .mean()
      .set('month', m);
}));

// Calculate max monthly mean per pixel
var Max = byMonth.max();

// Calculate min monthly mean per pixel
var Min = byMonth.min();

// Calculate difference between max and monthly mean
var Diff = Max.subtract(Min);

// Write to console
print(datasetn, 'Modis with NDVI');
print(byMonth, 'NDVI monthly aggregated');

// Check projektion
print('CRS', Diff.projection());
// Check Scale
print('Scale [m]', Diff.projection().nominalScale());

// Print specific Image of ImageCollection
//Map.addLayer(byMonth.filterMetadata('month','equals',8), {}, 'Sel');

// Print Max monthly mean
Map.addLayer(Max.clip(P1.buffer(5000)), {}, 'Max');

// Print image with max pixel value over time period
//Map.addLayer(dataset.select('NDVI').max().clip(P1.buffer(5000)), {}, 'MPV');

// Reproject to EPSG:3057
var Maxr = Max.reproject('EPSG:3057', null, 250);
var Minr = Min.reproject('EPSG:3057', null, 250);

// Check projektion after reprojection
print('CRS aft', Maxr.projection());
print('CRS aft', Minr.projection());

// Check Scale after reprojection
print('Scale [m] aft', Maxr.projection().nominalScale());

// Export the image Max, specifying scale and region.
Export.image.toDrive({
  image: Max,
  description: 'NDVI_max_v61',
  scale: 250,
  region: aoi,
  fileFormat: 'GeoTIFF'
});

// Export the image Min, specifying scale and region.
Export.image.toDrive({
  image: Min,
  description: 'NDVI_min_v61',
  scale: 250,
  region: aoi,
  fileFormat: 'GeoTIFF'
});