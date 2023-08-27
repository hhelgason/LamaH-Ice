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
var start_dt = ee.Date("2002-08-01");
var end_dt = ee.Date("2021-07-31");

// Function to mask clouds based on the 'FparLai_QC' band of Modis data
var cloudmask = function(img) {
    var qa = img.select('FparExtra_QC');
    // Check if the Aerosol (Bit 3) or the Cirrus (Bit 4) or the
    // Internal cloud mask (Bit 5) or the Cloud shadow (Bit 6) is true
    var cloud = qa.bitwiseAnd(1 << 3)
                .or(qa.bitwiseAnd(1 << 4))
                .or(qa.bitwiseAnd(1 << 5))
                .or(qa.bitwiseAnd(1 << 6));
    // Return cloud-free images
    return img.updateMask(cloud.not());
};

// Select dataset from GEE Data_Catalog
var dataset = ee.ImageCollection('MODIS/061/MCD15A3H')
  // Set timewindow, add 1 day to end date to make it inclusive
  .filterDate(start_dt, end_dt.advance(1, 'day'))
  // Selection of the needed band
  .select(['Lai', 'FparExtra_QC'])
  // Map function cloudmask over ImageCollection
  .map(cloudmask);

// Create list with numbers 1 to 12
var months = ee.List.sequence(1, 12);

// Group by month, and then reduce within groups by mean(),
// the result is an ImageCollection with one image for each month;
var byMonth = ee.ImageCollection.fromImages(
    months.map(function (m) {
    return dataset
      .select('Lai')
      .filter(ee.Filter.calendarRange(m, m, 'month'))
      .mean()
      .set('month', m)
      // Divide by 10 (Scale = 0.1)
      .divide(10)
      // Set value to 0, where input mask is zero
      .unmask(0);
}));

// Calculate max monthly mean per pixel
var Max = byMonth.max();

// Calculate min monthly mean per pixel
var Min = byMonth.min();

// Calculate difference between max and monthly mean
var Diff = Max.subtract(Min);

// Write to console
print(dataset);
print(byMonth);

// Check projektion
print('CRS', Diff.projection());
// Check Scale
print('Scale [m]', Diff.projection().nominalScale());

// Print specific Image of ImageCollection
//Map.addLayer(byMonth.filterMetadata('month','equals',7), {}, 'Sel');

// Print Max monthly mean
Map.addLayer(Max.clip(P1.buffer(5000)), {}, 'Max');

// Print image with max pixel value over time period
//Map.addLayer(dataset.select('Lai').max().clip(P1.buffer(5000)), {}, 'MPV');

// Reproject to EPSG:3057
var Maxr = Max.reproject('EPSG:3057', null, 500);
var Minr = Min.reproject('EPSG:3057', null, 500);
var Diffr = Diff.reproject('EPSG:3057', null, 500);

// Check projektion after reprojection
print('CRS aft', Maxr.projection());
print('CRS aft', Minr.projection());
print('CRS aft', Diffr.projection());

// Check Scale after reprojection
print('Scale [m] aft', Maxr.projection().nominalScale());

// Export the image Max, specifying scale and region.
Export.image.toDrive({
  image: Max,
  description: 'LAI_max_v61',
  scale: 500,
  region: aoi,
  fileFormat: 'GeoTIFF'
});

// Export the image Min, specifying scale and region.
Export.image.toDrive({
  image: Min,
  description: 'LAI_min_v61',
  scale: 500,
  region: aoi,
  fileFormat: 'GeoTIFF'
});

// Export the image Diff, specifying scale and region.
Export.image.toDrive({
  image: Diff,
  description: 'LAI_diff_v61',
  scale: 500,
  region: aoi,
  fileFormat: 'GeoTIFF'
});