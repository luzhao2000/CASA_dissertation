// for 2014 Landsat 8 classification in GEE
var shanghai = ee.FeatureCollection('users/luluzhao2022/shanghai'); // read in Shanghai boundary shp as ROI

// Applies scaling factors.
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2); 
// Landsat Collection 2 surface reflectance has a scale factor of 0.0000275 and an additional offset of -0.2 per pixel.
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

var cloudMaskL457 = function(image) {
  var qa = image.select('QA_PIXEL');
//云层表示为第五位，云层置信度为6-7位，云阴影为第三位
//选择出有云并且云层置信度为中等，以及有云阴影覆盖的像元。
  var cloud = qa.bitwiseAnd(1 << 4).or (qa.bitwiseAnd(1 << 5))
          .and(qa.bitwiseAnd(1 << 7))
          .or(qa.bitwiseAnd(1 << 3));
  // 移除边界像元
  var mask2 = image.mask().reduce(ee.Reducer.min());
  //将检测有关云像元置为0，掩模保留位置不为0的数据。
  return image.updateMask(cloud.not()).updateMask(mask2);
};

var cloudMaskC2L7 = function(image) {
  var dilatedCloud = (1 << 1)
  var cloud = (1 << 3)
  var cloudShadow = (1 << 4)
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(dilatedCloud)
    .and(qa.bitwiseAnd(cloud))
    .or(qa.bitwiseAnd(cloudShadow))
  return image.updateMask(mask.not());
};


// surface reflectance Landsat 8 images
var Landsat_2014_raw = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2014-01-01', '2014-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 2))
    .map(applyScaleFactors);
// use waytwo to remove cloud cover in 2014 is better than wayone
var Landsat_2014_waytwo = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2014-01-01', '2014-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 2))
    .map(applyScaleFactors)
    .map(cloudMaskL457);
var Landsat_2014_wayone = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2014-01-01', '2014-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 2))
    .map(applyScaleFactors)
    .map(cloudMaskC2L7);    


//// to get the better images 
// apply the median reducer
// var Landsat_2014_1_median = Landsat_2014_wayone.reduce(ee.Reducer.median());
var Landsat_2014_2_median = Landsat_2014_waytwo.reduce(ee.Reducer.median());
var Landsat_2014_raw_median = Landsat_2014_raw.reduce(ee.Reducer.median());

var visualization = {
  bands: ['SR_B4_median', 'SR_B3_median', 'SR_B2_median'],
  min: 0.0,
  max: 0.3,
};

var vis1 = {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0.0,
  max: 0.3,
};
Map.setCenter(121.6720, 31.0030, 8); // focus on Shanghai
// Map.addLayer(Landsat_2014_wayone, vis1, 'Landsat_2014_wayone');
Map.addLayer(Landsat_2014_waytwo, vis1, 'Landsat_2014_waytwo');
// Map.addLayer(Landsat_2014_1_median, visualization, 'Landsat_2014_1_median');
Map.addLayer(Landsat_2014_2_median, visualization, 'Landsat_2014_2_median');
// Map.addLayer(Landsat_2014_raw, vis1, 'Landsat_2014_raw');
// Map.addLayer(Landsat_2014_raw_median, visualization, 'Landsat_2014_2_median');
Map.addLayer(shanghai);

/// mosaic two images together
// var Landsat_2015_scale_meanImage = Landsat_2015_raw_median.mean();
// Map.addLayer(Landsat_2015_scale_meanImage, vis1, 'Landsat_2015_scale_meanImage');

/// clip images
var Landsat_2014_scale_meanImage_clip = Landsat_2014_2_median.clip(shanghai)
  .select(['SR_B1_median', 'SR_B2_median', 'SR_B3_median', 'SR_B4_median', 'SR_B5_median', 'SR_B6_median', 'SR_B7_median']);


///填补去云后缺失缺失的部分
//由于2014年影像云覆盖去除效果较为显著，且出现漏洞的区域几乎为建筑，不影响绿地的提取

var vis_params3 = {
  bands: ['SR_B5_median', 'SR_B4_median', 'SR_B3_median'],
  min: 0,
  max: 0.3,
};
// map the layer
Map.addLayer(Landsat_2014_scale_meanImage_clip, vis_params3, 'Landsat_2014_scale_meanImage_clip');

//MNDVI, NDVI, and NDBI: 构建NDVI/NDBI/MNDWI作为光谱特征
//利用DEM数据的坡度和高程作为地形特征；然后将构建的这些特征作为影像的一个波段进行使用
var mndwi = Landsat_2014_scale_meanImage_clip.normalizedDifference(['SR_B3_median', 'SR_B6_median']).rename('MNDWI');//计算MNDWI
var ndbi = Landsat_2014_scale_meanImage_clip.normalizedDifference(['SR_B6_median', 'SR_B5_median']).rename('NDBI');//计算NDBI
var ndvi = Landsat_2014_scale_meanImage_clip.normalizedDifference(['SR_B5_median', 'SR_B4_median']).rename('NDVI');//计算NDVI
var strm=ee.Image("USGS/SRTMGL1_003");
var dem=ee.Algorithms.Terrain(strm);
var elevation=dem.select('elevation');
var slope=dem.select('slope');
Landsat_2014_scale_meanImage_clip = Landsat_2014_scale_meanImage_clip.addBands(ndvi).addBands(ndbi).addBands(mndwi).addBands(elevation.rename("ELEVATION")).addBands(slope.rename("SLOPE"));

// Then, select training data
// and make a FeatureCollection from the polygons
var polygons = ee.FeatureCollection([
  ee.Feature(green_space, {'class': 1}),
  ee.Feature(other_area, {'class': 2}), // Add a non-green space area with class '2'
  ee.Feature(water, {'class': 3}),
]);

// Use these bands for classification.
var bands = ['SR_B2_median', 'SR_B3_median', 'SR_B4_median', 'SR_B5_median', 'SR_B6_median','SR_B7_median','MNDWI','NDBI','NDVI','SLOPE', 'ELEVATION'];
// The name of the property on the points storing the class label.
var classProperty = 'class';

// we need to pull out the data from our training areas
var training = Landsat_2014_scale_meanImage_clip.select(bands).sampleRegions({
  collection: polygons,
  properties: [classProperty],
  scale: 30
});
print(training, "training");

/// machine learning classification method 1: CART
// Train a CART classifier
var classifier = ee.Classifier.smileCart().train({
  features: training,
  classProperty: classProperty,
});
// Print some info about the classifier (specific to CART)
print('CART, explained', classifier.explain());
// Classify the image
var CART_classified = Landsat_2014_scale_meanImage_clip.classify(classifier);
Map.addLayer(CART_classified, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_CART");

// Classify the reference data using the trained classifier
var classified_training = training.classify(classifier);
// Calculate the accuracy for each class
var accuracy = classified_training.errorMatrix('class', 'classification');
// Print the confusion matrix
print('CART Confusion Matrix:', accuracy);
// Calculate overall accuracy
var overallAccuracy = accuracy.accuracy();
print('CART Overall Accuracy:', overallAccuracy);


/// machine learning classification method 2: random forest
//Fist, add a column of random uniforms to the training dataset.
var withRandom = polygons.randomColumn('random');
print(withRandom, 'sample_polygon');
// We want to reserve some of the data for testing, to avoid overfitting the model (Roughly 70% training, 30% testing)
var split = 0.7;
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));
print(trainingPartition, "train");
print(testingPartition, "test");

// take samples from image for training and validation  
var training_rf = Landsat_2014_scale_meanImage_clip.select(bands).sampleRegions({
  collection: trainingPartition,
  properties: [classProperty],
  scale: 30,
});
var validation_rf = Landsat_2014_scale_meanImage_clip.select(bands).sampleRegions({
  collection: testingPartition,
  properties: [classProperty],
  scale: 30,
});
// Random Forest Classification
var rf1 = ee.Classifier.smileRandomForest(100)
  .train(training_rf, 'class', bands); //classifier
var rf2 = Landsat_2014_scale_meanImage_clip.classify(rf1);//validation
// Classify the test FeatureCollection.
var test = validation_rf.classify(rf1);

/// evaluate the extraction accuracy
var testAccuracy = test.errorMatrix('class', 'classification');
var consumers = testAccuracy.consumersAccuracy();
print('Validation error matrix: ', testAccuracy);
print('Validation overall accuracy: ', testAccuracy.accuracy());
print('Validation consumer accuracy: ', consumers);
// print('kappa accuracy', consumers.kappa());
Map.addLayer(rf2, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "RF");


/// random forest pixel approach
var pixel_number = 1000;
var green_space_points=ee.FeatureCollection.randomPoints(green_space, pixel_number).map(function(i){
  return i.set({'class': 1})})
var other_area_points=ee.FeatureCollection.randomPoints(other_area, pixel_number).map(function(i){
  return i.set({'class': 2})})
var water_points=ee.FeatureCollection.randomPoints(water, pixel_number).map(function(i){
  return i.set({'class': 3})})
var point_sample=ee.FeatureCollection([green_space_points,
                                  other_area_points,
                                  water_points])
                                  .flatten()
                                  .randomColumn();
// assign 70% of training points to validation 
var split=0.7
var training_sample = point_sample.filter(ee.Filter.lt('random', split));
var validation_sample = point_sample.filter(ee.Filter.gte('random', split));
// take samples from image for training and validation  
var training = Landsat_2014_scale_meanImage_clip.select(bands).sampleRegions({
  collection: training_sample,
  properties: ['class'],
  scale: 30,
});
var validation = Landsat_2014_scale_meanImage_clip.select(bands).sampleRegions({
  collection: validation_sample,
  properties: ['class'],
  scale: 30
});
// Random Forest Classification
var rf1_pixel = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
// Get information about the trained classifier.
print('Results of RF trained classifier', rf1_pixel.explain());

// conduct classification
var rf2_pixel = Landsat_2014_scale_meanImage_clip.classify(rf1_pixel);
Map.addLayer(rf2_pixel, {min: 1, max: 3, 
  palette: ['1c5f2c', 'b3ac9f','466b9f']},
  "RF_pixel");

/// assess accuracy
var rf_trainAccuracy = rf1_pixel.confusionMatrix();
print('Resubstitution error matrix: ', rf_trainAccuracy);
print('Training overall accuracy: ', rf_trainAccuracy.accuracy());
var validated = validation.classify(rf1_pixel);
var testAccuracy = validated.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('RF Validation error matrix: ', testAccuracy);
print('RF Validation overall accuracy: ', testAccuracy.accuracy())
print('RF Validation consumer accuracy: ', consumers);

// // Define the output file name and path
// var outputFileName = 'rf_pixel_classification_result_2014';
// var outputFileFolder = 'GEE_exports';

// // Set the scale and region of interest
// var scale = 30; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: rf2_pixel, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });

// // Define the output file name and path
// var outputFileName = 'CART_classification_result_2014';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 30; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: CART_classified, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });


//////////////////// for 2015 Landsat 8 classification in GEE////////////////////////////////
var shanghai = ee.FeatureCollection('users/luluzhao2022/shanghai'); // read in Shanghai boundary shp as ROI

// Applies scaling factors.
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2); 
// Landsat Collection 2 surface reflectance has a scale factor of 0.0000275 and an additional offset of -0.2 per pixel.
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

var cloudMaskL457 = function(image) {
  var qa = image.select('QA_PIXEL');
//云层表示为第五位，云层置信度为6-7位，云阴影为第三位
//选择出有云并且云层置信度为中等，以及有云阴影覆盖的像元。
  var cloud = qa.bitwiseAnd(1 << 4).or (qa.bitwiseAnd(1 << 5))
          .and(qa.bitwiseAnd(1 << 7))
          .or(qa.bitwiseAnd(1 << 3));
  // 移除边界像元
  var mask2 = image.mask().reduce(ee.Reducer.min());
  //将检测有关云像元置为0，掩模保留位置不为0的数据。
  return image.updateMask(cloud.not()).updateMask(mask2);
};

var cloudMaskC2L7 = function(image) {
  var dilatedCloud = (1 << 1)
  var cloud = (1 << 3)
  var cloudShadow = (1 << 4)
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(dilatedCloud)
    .and(qa.bitwiseAnd(cloud))
    .or(qa.bitwiseAnd(cloudShadow))
  return image.updateMask(mask.not());
};


// surface reflectance Landsat 8 images
var Landsat_2015_raw = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2015-01-01', '2015-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 0.5))
    .map(applyScaleFactors);
var Landsat_2015_waytwo = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2015-01-01', '2015-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 0.5))
    .map(applyScaleFactors)
    .map(cloudMaskL457);
var Landsat_2015_wayone = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2015-01-01', '2015-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 0.5))
    .map(applyScaleFactors)
    .map(cloudMaskC2L7);    


//// to get the better images 
// apply the median reducer
var Landsat_2015_1_median = Landsat_2015_wayone.reduce(ee.Reducer.median());
var Landsat_2015_2_median = Landsat_2015_waytwo.reduce(ee.Reducer.median());
var Landsat_2015_raw_median = Landsat_2015_raw.reduce(ee.Reducer.median());

var visualization = {
  bands: ['SR_B4_median', 'SR_B3_median', 'SR_B2_median'],
  min: 0.0,
  max: 0.3,
};

var vis1 = {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0.0,
  max: 0.3,
};
Map.setCenter(121.6720, 31.0030, 8); // focus on Shanghai
// Map.addLayer(Landsat_2015_wayone, vis1, 'Landsat_2015_wayone');
// Map.addLayer(Landsat_2015_waytwo, vis1, 'Landsat_2015_waytwo');
// Map.addLayer(Landsat_2015_1_median, visualization, 'Landsat_2015_1_median');
// Map.addLayer(Landsat_2015_2_median, visualization, 'Landsat_2015_2_median');
Map.addLayer(Landsat_2015_raw, vis1, 'Landsat_2015_raw');
Map.addLayer(Landsat_2015_raw_median, visualization, 'Landsat_2015_2_median');
Map.addLayer(shanghai);

/// mosaic two images together
// var Landsat_2015_scale_meanImage = Landsat_2015_raw_median.mean();
// Map.addLayer(Landsat_2015_scale_meanImage, vis1, 'Landsat_2015_scale_meanImage');

/// clip images
var Landsat_2015_scale_meanImage_clip = Landsat_2015_raw_median.clip(shanghai)
  .select(['SR_B1_median', 'SR_B2_median', 'SR_B3_median', 'SR_B4_median', 'SR_B5_median', 'SR_B6_median', 'SR_B7_median']);


///填补去云后缺失缺失的部分
//由于2015年影像云覆盖可忽略不计，所以省略这一步

var vis_params3 = {
  bands: ['SR_B5_median', 'SR_B4_median', 'SR_B3_median'],
  min: 0,
  max: 0.3,
};
// map the layer
Map.addLayer(Landsat_2015_scale_meanImage_clip, vis_params3, 'Landsat_2015_scale_meanImage_clip');

//MNDVI, NDVI, and NDBI: 构建NDVI/NDBI/MNDWI作为光谱特征
//利用DEM数据的坡度和高程作为地形特征；然后将构建的这些特征作为影像的一个波段进行使用
var mndwi = Landsat_2015_scale_meanImage_clip.normalizedDifference(['SR_B3_median', 'SR_B6_median']).rename('MNDWI');//计算MNDWI
var ndbi = Landsat_2015_scale_meanImage_clip.normalizedDifference(['SR_B6_median', 'SR_B5_median']).rename('NDBI');//计算NDBI
var ndvi = Landsat_2015_scale_meanImage_clip.normalizedDifference(['SR_B5_median', 'SR_B4_median']).rename('NDVI');//计算NDVI
var strm=ee.Image("USGS/SRTMGL1_003");
var dem=ee.Algorithms.Terrain(strm);
var elevation=dem.select('elevation');
var slope=dem.select('slope');
Landsat_2015_scale_meanImage_clip = Landsat_2015_scale_meanImage_clip.addBands(ndvi).addBands(ndbi).addBands(mndwi).addBands(elevation.rename("ELEVATION")).addBands(slope.rename("SLOPE"));

// Then, select training data
// and make a FeatureCollection from the polygons
var polygons = ee.FeatureCollection([
  ee.Feature(green_space, {'class': 1}),
  ee.Feature(other_area, {'class': 2}), // Add a non-green space area with class '2'
  ee.Feature(water, {'class': 3}),
]);

// Use these bands for classification.
var bands = ['SR_B2_median', 'SR_B3_median', 'SR_B4_median', 'SR_B5_median', 'SR_B6_median','SR_B7_median','MNDWI','NDBI','NDVI','SLOPE', 'ELEVATION'];
// The name of the property on the points storing the class label.
var classProperty = 'class';

// we need to pull out the data from our training areas
var training = Landsat_2015_scale_meanImage_clip.select(bands).sampleRegions({
  collection: polygons,
  properties: [classProperty],
  scale: 30
});
print(training, "training");

/// machine learning classification method 1: CART
// Train a CART classifier
var classifier = ee.Classifier.smileCart().train({
  features: training,
  classProperty: classProperty,
});
// Print some info about the classifier (specific to CART)
print('CART, explained', classifier.explain());
// Classify the image
var CART_classified = Landsat_2015_scale_meanImage_clip.classify(classifier);

Map.addLayer(CART_classified, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_CART");

// Classify the reference data using the trained classifier
var classified_training = training.classify(classifier);
// Calculate the accuracy for each class
var accuracy = classified_training.errorMatrix('class', 'classification');
// Print the confusion matrix
print('VART Confusion Matrix:', accuracy);
// Calculate overall accuracy
var overallAccuracy = accuracy.accuracy();
print('CART Overall Accuracy:', overallAccuracy);


/// machine learning classification method 2: random forest
//Fist, add a column of random uniforms to the training dataset.
var withRandom = polygons.randomColumn('random');
print(withRandom, 'sample_polygon');
// We want to reserve some of the data for testing, to avoid overfitting the model (Roughly 70% training, 30% testing)
var split = 0.7;
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));
print(trainingPartition, "train");
print(testingPartition, "test");

// take samples from image for training and validation  
var training_rf = Landsat_2015_scale_meanImage_clip.select(bands).sampleRegions({
  collection: trainingPartition,
  properties: [classProperty],
  scale: 30,
});
var validation_rf = Landsat_2015_scale_meanImage_clip.select(bands).sampleRegions({
  collection: testingPartition,
  properties: [classProperty],
  scale: 30,
});
// Random Forest Classification
var rf1 = ee.Classifier.smileRandomForest(100)
  .train(training_rf, 'class', bands); //classifier
var rf2 = Landsat_2015_scale_meanImage_clip.classify(rf1);//validation
// Classify the test FeatureCollection.
var test = validation_rf.classify(rf1);

/// evaluate the extraction accuracy
var testAccuracy = test.errorMatrix('class', 'classification');
var consumers = testAccuracy.consumersAccuracy();
print('Validation error matrix: ', testAccuracy);
print('Validation overall accuracy: ', testAccuracy.accuracy());
print('Validation consumer accuracy: ', consumers);
// print('kappa accuracy', consumers.kappa());
Map.addLayer(rf2, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "RF");


/// random forest pixel approach
var pixel_number = 1000;
var green_space_points=ee.FeatureCollection.randomPoints(green_space, pixel_number).map(function(i){
  return i.set({'class': 1})})
var other_area_points=ee.FeatureCollection.randomPoints(other_area, pixel_number).map(function(i){
  return i.set({'class': 2})})
var water_points=ee.FeatureCollection.randomPoints(water, pixel_number).map(function(i){
  return i.set({'class': 3})})
var point_sample=ee.FeatureCollection([green_space_points,
                                  other_area_points,
                                  water_points])
                                  .flatten()
                                  .randomColumn();
// assign 70% of training points to validation 
var split=0.7
var training_sample = point_sample.filter(ee.Filter.lt('random', split));
var validation_sample = point_sample.filter(ee.Filter.gte('random', split));
// take samples from image for training and validation  
var training = Landsat_2015_scale_meanImage_clip.select(bands).sampleRegions({
  collection: training_sample,
  properties: ['class'],
  scale: 30,
});
var validation = Landsat_2015_scale_meanImage_clip.select(bands).sampleRegions({
  collection: validation_sample,
  properties: ['class'],
  scale: 30
});
// Random Forest Classification
var rf1_pixel = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
// Get information about the trained classifier.
print('Results of RF trained classifier', rf1_pixel.explain());

// conduct classification
var rf2_pixel = Landsat_2015_scale_meanImage_clip.classify(rf1_pixel);
Map.addLayer(rf2_pixel, {min: 1, max: 3, 
  palette: ['1c5f2c', 'b3ac9f','466b9f']},
  "RF_pixel");

/// assess accuracy
var trainAccuracy = rf1_pixel.confusionMatrix();
print('Resubstitution error matrix: ', trainAccuracy);
print('Training overall accuracy: ', trainAccuracy.accuracy());
var validated = validation.classify(rf1_pixel);
var testAccuracy = validated.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('RF Validation error matrix: ', testAccuracy);
print('RF Validation overall accuracy: ', testAccuracy.accuracy())
print('RF Validation consumer accuracy: ', consumers);


// // Define the output file name and path
// var outputFileName = 'rf_pixel_classification_result_2015';
// var outputFileFolder = 'GEE_exports';

// // Set the scale and region of interest
// var scale = 30; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: rf2_pixel, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });

// // Define the output file name and path
// var outputFileName = 'CART_classification_result_2015';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 30; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: CART_classified, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });



//////////////////// for 2016 Landsat 8 classification in GEE////////////////////////////////
var shanghai = ee.FeatureCollection('users/luluzhao2022/shanghai'); // read in Shanghai boundary shp as ROI

// Applies scaling factors.
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2); 
// Landsat Collection 2 surface reflectance has a scale factor of 0.0000275 and an additional offset of -0.2 per pixel.
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

var cloudMaskL457 = function(image) {
  var qa = image.select('QA_PIXEL');
//云层表示为第五位，云层置信度为6-7位，云阴影为第三位
//选择出有云并且云层置信度为中等，以及有云阴影覆盖的像元。
  var cloud = qa.bitwiseAnd(1 << 4).or (qa.bitwiseAnd(1 << 5))
          .and(qa.bitwiseAnd(1 << 7))
          .or(qa.bitwiseAnd(1 << 3));
  // 移除边界像元
  var mask2 = image.mask().reduce(ee.Reducer.min());
  //将检测有关云像元置为0，掩模保留位置不为0的数据。
  return image.updateMask(cloud.not()).updateMask(mask2);
};

var cloudMaskC2L7 = function(image) {
  var dilatedCloud = (1 << 1)
  var cloud = (1 << 3)
  var cloudShadow = (1 << 4)
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(dilatedCloud)
    .and(qa.bitwiseAnd(cloud))
    .or(qa.bitwiseAnd(cloudShadow))
  return image.updateMask(mask.not());
};


// surface reflectance Landsat 8 images
var Landsat_2016_raw = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2016-01-01', '2016-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 2.5))
    .map(applyScaleFactors);
var Landsat_2016_waytwo = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2016-01-01', '2016-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 2.5))
    .map(applyScaleFactors)
    .map(cloudMaskL457);
var Landsat_2016_wayone = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2016-01-01', '2016-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 2.5))
    .map(applyScaleFactors)
    .map(cloudMaskC2L7);    


//// to get the better images 
// apply the median reducer
var Landsat_2016_1_median = Landsat_2016_wayone.reduce(ee.Reducer.median());
var Landsat_2016_2_median = Landsat_2016_waytwo.reduce(ee.Reducer.median());
var Landsat_2016_raw_median = Landsat_2016_raw.reduce(ee.Reducer.median());

var visualization = {
  bands: ['SR_B4_median', 'SR_B3_median', 'SR_B2_median'],
  min: 0.0,
  max: 0.3,
};

var vis1 = {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0.0,
  max: 0.3,
};
Map.setCenter(121.6720, 31.0030, 8); // focus on Shanghai
// Map.addLayer(Landsat_2016_wayone, vis1, 'Landsat_2016_wayone');
// Map.addLayer(Landsat_2016_waytwo, vis1, 'Landsat_2016_waytwo');
// Map.addLayer(Landsat_2016_1_median, visualization, 'Landsat_2016_1_median');
// Map.addLayer(Landsat_2016_2_median, visualization, 'Landsat_2016_2_median');
Map.addLayer(Landsat_2016_raw, vis1, 'Landsat_2016_raw');
Map.addLayer(Landsat_2016_raw_median, visualization, 'Landsat_2016_2_rawmedian');
Map.addLayer(shanghai);

/// clip images
var Landsat_2016_scale_meanImage_clip = Landsat_2016_raw_median.clip(shanghai)
  .select(['SR_B1_median', 'SR_B2_median', 'SR_B3_median', 'SR_B4_median', 'SR_B5_median', 'SR_B6_median', 'SR_B7_median']);


///填补去云后缺失缺失的部分
//由于2016年影像云覆盖可忽略不计，所以省略这一步

var vis_params3 = {
  bands: ['SR_B5_median', 'SR_B4_median', 'SR_B3_median'],
  min: 0,
  max: 0.3,
};
// map the layer
Map.addLayer(Landsat_2016_scale_meanImage_clip, vis_params3, 'Landsat_2016_scale_meanImage_clip');

//MNDVI, NDVI, and NDBI: 构建NDVI/NDBI/MNDWI作为光谱特征
//利用DEM数据的坡度和高程作为地形特征；然后将构建的这些特征作为影像的一个波段进行使用
var mndwi = Landsat_2016_scale_meanImage_clip.normalizedDifference(['SR_B3_median', 'SR_B6_median']).rename('MNDWI');//计算MNDWI
var ndbi = Landsat_2016_scale_meanImage_clip.normalizedDifference(['SR_B6_median', 'SR_B5_median']).rename('NDBI');//计算NDBI
var ndvi = Landsat_2016_scale_meanImage_clip.normalizedDifference(['SR_B5_median', 'SR_B4_median']).rename('NDVI');//计算NDVI
var strm=ee.Image("USGS/SRTMGL1_003");
var dem=ee.Algorithms.Terrain(strm);
var elevation=dem.select('elevation');
var slope=dem.select('slope');
Landsat_2016_scale_meanImage_clip = Landsat_2016_scale_meanImage_clip.addBands(ndvi).addBands(ndbi).addBands(mndwi).addBands(elevation.rename("ELEVATION")).addBands(slope.rename("SLOPE"));

// Then, select training data
// and make a FeatureCollection from the polygons
var polygons = ee.FeatureCollection([
  ee.Feature(green_space, {'class': 1}),
  ee.Feature(other_area, {'class': 2}), // Add a non-green space area with class '2'
  ee.Feature(water, {'class': 3}),
]);

// Use these bands for classification.
var bands = ['SR_B2_median', 'SR_B3_median', 'SR_B4_median', 'SR_B5_median', 'SR_B6_median','SR_B7_median','MNDWI','NDBI','NDVI','SLOPE', 'ELEVATION'];
// The name of the property on the points storing the class label.
var classProperty = 'class';

// we need to pull out the data from our training areas
var training = Landsat_2016_scale_meanImage_clip.select(bands).sampleRegions({
  collection: polygons,
  properties: [classProperty],
  scale: 30
});
print(training, "training");

/// machine learning classification method 1: CART
// Train a CART classifier
// var classifier = ee.Classifier.smileCart().train({
//   features: training,
//   classProperty: classProperty,
// });
// // Print some info about the classifier (specific to CART)
// print('CART, explained', classifier.explain());
// // Classify the image
// var CART_classified = Landsat_2016_scale_meanImage_clip.classify(classifier);
// Map.addLayer(CART_classified, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_CART");


// Split the data into training and validation sets
var split = 0.7; // Percentage for training set
var trainingPartition = training.filter(ee.Filter.lt('random', split));
var validationPartition = training.filter(ee.Filter.gte('random', split));
// Train a CART classifier on the training partition
var classifier = ee.Classifier.smileCart().train({
  features: trainingPartition,
  classProperty: classProperty,
});
// Classify the validation partition
var validationClassified = validationPartition.classify(classifier);
// Get the classified values and ground truth from the validation partition
var classifiedValues = validationClassified.select('classification');
var groundTruth = validationPartition.aggregate_array(classProperty);
// Compute the confusion matrix
var confusionMatrix = classifiedValues.errorMatrix(classProperty, groundTruth);
// Calculate overall accuracy
var overallAccuracy = confusionMatrix.accuracy();
// Print the confusion matrix and overall accuracy
print('Confusion Matrix:',confusionMatrix);
// print('Overall Accuracy:', overallAccuracy);


// // Classify the reference data using the trained classifier
// var classified_training = training.classify(classifier);
// // Calculate the accuracy for each class
// var accuracy = classified_training.errorMatrix('class', 'classification');
// // Print the confusion matrix
// print('CART Confusion Matrix:', accuracy);
// // Calculate overall accuracy
// var overallAccuracy = accuracy.accuracy();
// print('CART Overall Accuracy:', overallAccuracy);


// /// machine learning classification method 2: random forest
// //Fist, add a column of random uniforms to the training dataset.
// var withRandom = polygons.randomColumn('random');
// print(withRandom, 'sample_polygon');
// // We want to reserve some of the data for testing, to avoid overfitting the model (Roughly 70% training, 30% testing)
// var split = 0.7;
// var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
// var testingPartition = withRandom.filter(ee.Filter.gte('random', split));
// print(trainingPartition, "train");
// print(testingPartition, "test");

// // take samples from image for training and validation  
// var training_rf = Landsat_2016_scale_meanImage_clip.select(bands).sampleRegions({
//   collection: trainingPartition,
//   properties: [classProperty],
//   scale: 30,
// });
// var validation_rf = Landsat_2016_scale_meanImage_clip.select(bands).sampleRegions({
//   collection: testingPartition,
//   properties: [classProperty],
//   scale: 30,
// });
// // Random Forest Classification
// var rf1 = ee.Classifier.smileRandomForest(100)
//   .train(training_rf, 'class', bands); //classifier
// var rf2 = Landsat_2016_scale_meanImage_clip.classify(rf1);//validation
// // Classify the test FeatureCollection.
// var test = validation_rf.classify(rf1);

// /// evaluate the extraction accuracy
// var testAccuracy = test.errorMatrix('class', 'classification');
// var consumers = testAccuracy.consumersAccuracy();
// print('Validation error matrix: ', testAccuracy);
// print('Validation overall accuracy: ', testAccuracy.accuracy());
// print('Validation consumer accuracy: ', consumers);
// // print('kappa accuracy', consumers.kappa());
// Map.addLayer(rf2, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "RF");


// /// random forest pixel approach
// var pixel_number = 1000;
// var green_space_points=ee.FeatureCollection.randomPoints(green_space, pixel_number).map(function(i){
//   return i.set({'class': 1})})
// var other_area_points=ee.FeatureCollection.randomPoints(other_area, pixel_number).map(function(i){
//   return i.set({'class': 2})})
// var water_points=ee.FeatureCollection.randomPoints(water, pixel_number).map(function(i){
//   return i.set({'class': 3})})
// var point_sample=ee.FeatureCollection([green_space_points,
//                                   other_area_points,
//                                   water_points])
//                                   .flatten()
//                                   .randomColumn();
// // assign 70% of training points to validation 
// var split=0.7
// var training_sample = point_sample.filter(ee.Filter.lt('random', split));
// var validation_sample = point_sample.filter(ee.Filter.gte('random', split));
// // take samples from image for training and validation  
// var training = Landsat_2016_scale_meanImage_clip.select(bands).sampleRegions({
//   collection: training_sample,
//   properties: ['class'],
//   scale: 30,
// });
// var validation = Landsat_2016_scale_meanImage_clip.select(bands).sampleRegions({
//   collection: validation_sample,
//   properties: ['class'],
//   scale: 30
// });
// // Random Forest Classification
// var rf1_pixel = ee.Classifier.smileRandomForest(100)
//     .train(training, 'class');
// // Get information about the trained classifier.
// print('Results of RF trained classifier', rf1_pixel.explain());

// // conduct classification
// var rf2_pixel = Landsat_2016_scale_meanImage_clip.classify(rf1_pixel);
// Map.addLayer(rf2_pixel, {min: 1, max: 3, 
//   palette: ['1c5f2c', 'b3ac9f','466b9f']},
//   "RF_pixel");

// /// assess accuracy
// var trainAccuracy = rf1_pixel.confusionMatrix();
// print('Resubstitution error matrix: ', trainAccuracy);
// print('Training overall accuracy: ', trainAccuracy.accuracy());
// var validated = validation.classify(rf1_pixel);
// var testAccuracy = validated.errorMatrix('class', 'classification');
// var consumers=testAccuracy.consumersAccuracy()
// print('RF Validation error matrix: ', testAccuracy);
// print('RF Validation overall accuracy: ', testAccuracy.accuracy())
// print('RF Validation consumer accuracy: ', consumers);


// // Define the output file name and path
// var outputFileName = 'rf_pixel_classification_result_2016';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 30; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: rf2_pixel, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });
// // Define the output file name and path
// var outputFileName = 'CART_classification_result_2016';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 30; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: CART_classified, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });




//////////////////// for 2017 Landsat 8 classification in GEE////////////////////////////////
var shanghai = ee.FeatureCollection('users/luluzhao2022/shanghai'); // read in Shanghai boundary shp as ROI

// Applies scaling factors.
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2); 
// Landsat Collection 2 surface reflectance has a scale factor of 0.0000275 and an additional offset of -0.2 per pixel.
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

var cloudMaskL457 = function(image) {
  var qa = image.select('QA_PIXEL');
//云层表示为第五位，云层置信度为6-7位，云阴影为第三位
//选择出有云并且云层置信度为中等，以及有云阴影覆盖的像元。
  var cloud = qa.bitwiseAnd(1 << 4).or (qa.bitwiseAnd(1 << 5))
          .and(qa.bitwiseAnd(1 << 7))
          .or(qa.bitwiseAnd(1 << 3));
  // 移除边界像元
  var mask2 = image.mask().reduce(ee.Reducer.min());
  //将检测有关云像元置为0，掩模保留位置不为0的数据。
  return image.updateMask(cloud.not()).updateMask(mask2);
};

var cloudMaskC2L7 = function(image) {
  var dilatedCloud = (1 << 1)
  var cloud = (1 << 3)
  var cloudShadow = (1 << 4)
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(dilatedCloud)
    .and(qa.bitwiseAnd(cloud))
    .or(qa.bitwiseAnd(cloudShadow))
  return image.updateMask(mask.not());
};


// surface reflectance Landsat 8 images
var Landsat_2017_raw = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2017-01-01', '2017-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 0.5))
    .map(applyScaleFactors);
var Landsat_2017_waytwo = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2017-01-01', '2017-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 0.5))
    .map(applyScaleFactors)
    .map(cloudMaskL457);
var Landsat_2017_wayone = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2017-01-01', '2017-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 0.5))
    .map(applyScaleFactors)
    .map(cloudMaskC2L7);    


//// to get the better images 
// apply the median reducer
var Landsat_2017_1_median = Landsat_2017_wayone.reduce(ee.Reducer.median());
var Landsat_2017_2_median = Landsat_2017_waytwo.reduce(ee.Reducer.median());
var Landsat_2017_raw_median = Landsat_2017_raw.reduce(ee.Reducer.median());

var visualization = {
  bands: ['SR_B4_median', 'SR_B3_median', 'SR_B2_median'],
  min: 0.0,
  max: 0.3,
};

var vis1 = {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0.0,
  max: 0.3,
};
Map.setCenter(121.6720, 31.0030, 8); // focus on Shanghai
// Map.addLayer(Landsat_2017_wayone, vis1, 'Landsat_2017_wayone');
// Map.addLayer(Landsat_2017_waytwo, vis1, 'Landsat_2017_waytwo');
// Map.addLayer(Landsat_2017_1_median, visualization, 'Landsat_2017_1_median');
// Map.addLayer(Landsat_2017_2_median, visualization, 'Landsat_2017_2_median');
Map.addLayer(Landsat_2017_raw, vis1, 'Landsat_2017_raw');
Map.addLayer(Landsat_2017_raw_median, visualization, 'Landsat_2017_2_rawmedian');
Map.addLayer(shanghai);


/// clip images
var Landsat_2017_scale_meanImage_clip = Landsat_2017_raw_median.clip(shanghai)
  .select(['SR_B1_median', 'SR_B2_median', 'SR_B3_median', 'SR_B4_median', 'SR_B5_median', 'SR_B6_median', 'SR_B7_median']);


///填补去云后缺失缺失的部分
//由于2017年影像几乎无云覆盖，所以省略这一步

var vis_params3 = {
  bands: ['SR_B5_median', 'SR_B4_median', 'SR_B3_median'],
  min: 0,
  max: 0.3,
};
// map the layer
Map.addLayer(Landsat_2017_scale_meanImage_clip, vis_params3, 'Landsat_2017_scale_meanImage_clip');

//MNDVI, NDVI, and NDBI: 构建NDVI/NDBI/MNDWI作为光谱特征
//利用DEM数据的坡度和高程作为地形特征；然后将构建的这些特征作为影像的一个波段进行使用
var mndwi = Landsat_2017_scale_meanImage_clip.normalizedDifference(['SR_B3_median', 'SR_B6_median']).rename('MNDWI');//计算MNDWI
var ndbi = Landsat_2017_scale_meanImage_clip.normalizedDifference(['SR_B6_median', 'SR_B5_median']).rename('NDBI');//计算NDBI
var ndvi = Landsat_2017_scale_meanImage_clip.normalizedDifference(['SR_B5_median', 'SR_B4_median']).rename('NDVI');//计算NDVI
var strm=ee.Image("USGS/SRTMGL1_003");
var dem=ee.Algorithms.Terrain(strm);
var elevation=dem.select('elevation');
var slope=dem.select('slope');
Landsat_2017_scale_meanImage_clip = Landsat_2017_scale_meanImage_clip.addBands(ndvi).addBands(ndbi).addBands(mndwi).addBands(elevation.rename("ELEVATION")).addBands(slope.rename("SLOPE"));

// Then, select training data
// and make a FeatureCollection from the polygons
var polygons = ee.FeatureCollection([
  ee.Feature(green_space, {'class': 1}),
  ee.Feature(other_area, {'class': 2}), // Add a non-green space area with class '2'
  ee.Feature(water, {'class': 3}),
]);

// Use these bands for classification.
var bands = ['SR_B2_median', 'SR_B3_median', 'SR_B4_median', 'SR_B5_median', 'SR_B6_median','SR_B7_median','MNDWI','NDBI','NDVI','SLOPE', 'ELEVATION'];
// The name of the property on the points storing the class label.
var classProperty = 'class';

// we need to pull out the data from our training areas
var training = Landsat_2017_scale_meanImage_clip.select(bands).sampleRegions({
  collection: polygons,
  properties: [classProperty],
  scale: 30
});
print(training, "training");

/// machine learning classification method 1: CART
// Train a CART classifier
var classifier = ee.Classifier.smileCart().train({
  features: training,
  classProperty: classProperty,
});
// Print some info about the classifier (specific to CART)
print('CART, explained', classifier.explain());
// Classify the image
var CART_classified = Landsat_2017_scale_meanImage_clip.classify(classifier);

Map.addLayer(CART_classified, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_CART");

// Classify the reference data using the trained classifier
var classified_training = training.classify(classifier);
// Calculate the accuracy for each class
var accuracy = classified_training.errorMatrix('class', 'classification');
// Print the confusion matrix
print('CART Confusion Matrix:', accuracy);
// Calculate overall accuracy
var overallAccuracy = accuracy.accuracy();
print('CART Overall Accuracy:', overallAccuracy);


/// machine learning classification method 2: random forest
//Fist, add a column of random uniforms to the training dataset.
var withRandom = polygons.randomColumn('random');
print(withRandom, 'sample_polygon');
// We want to reserve some of the data for testing, to avoid overfitting the model (Roughly 70% training, 30% testing)
var split = 0.7;
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));
print(trainingPartition, "train");
print(testingPartition, "test");

// take samples from image for training and validation  
var training_rf = Landsat_2017_scale_meanImage_clip.select(bands).sampleRegions({
  collection: trainingPartition,
  properties: [classProperty],
  scale: 30,
});
var validation_rf = Landsat_2017_scale_meanImage_clip.select(bands).sampleRegions({
  collection: testingPartition,
  properties: [classProperty],
  scale: 30,
});
// Random Forest Classification
var rf1 = ee.Classifier.smileRandomForest(100)
  .train(training_rf, 'class', bands); //classifier
var rf2 = Landsat_2017_scale_meanImage_clip.classify(rf1);//validation
// Classify the test FeatureCollection.
var test = validation_rf.classify(rf1);

/// evaluate the extraction accuracy
var testAccuracy = test.errorMatrix('class', 'classification');
var consumers = testAccuracy.consumersAccuracy();
print('Validation error matrix: ', testAccuracy);
print('Validation overall accuracy: ', testAccuracy.accuracy());
print('Validation consumer accuracy: ', consumers);
// print('kappa accuracy', consumers.kappa());
Map.addLayer(rf2, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "RF");


/// random forest pixel approach
var pixel_number = 1000;
var green_space_points=ee.FeatureCollection.randomPoints(green_space, pixel_number).map(function(i){
  return i.set({'class': 1})})
var other_area_points=ee.FeatureCollection.randomPoints(other_area, pixel_number).map(function(i){
  return i.set({'class': 2})})
var water_points=ee.FeatureCollection.randomPoints(water, pixel_number).map(function(i){
  return i.set({'class': 3})})
var point_sample=ee.FeatureCollection([green_space_points,
                                  other_area_points,
                                  water_points])
                                  .flatten()
                                  .randomColumn();
// assign 70% of training points to validation 
var split=0.7
var training_sample = point_sample.filter(ee.Filter.lt('random', split));
var validation_sample = point_sample.filter(ee.Filter.gte('random', split));
// take samples from image for training and validation  
var training = Landsat_2017_scale_meanImage_clip.select(bands).sampleRegions({
  collection: training_sample,
  properties: ['class'],
  scale: 30,
});
var validation = Landsat_2017_scale_meanImage_clip.select(bands).sampleRegions({
  collection: validation_sample,
  properties: ['class'],
  scale: 30
});
// Random Forest Classification
var rf1_pixel = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
// Get information about the trained classifier.
print('Results of RF trained classifier', rf1_pixel.explain());

// conduct classification
var rf2_pixel = Landsat_2017_scale_meanImage_clip.classify(rf1_pixel);
Map.addLayer(rf2_pixel, {min: 1, max: 3, 
  palette: ['1c5f2c', 'b3ac9f','466b9f']},
  "RF_pixel");

/// assess accuracy
var trainAccuracy = rf1_pixel.confusionMatrix();
print('Resubstitution error matrix: ', trainAccuracy);
print('Training overall accuracy: ', trainAccuracy.accuracy());
var validated = validation.classify(rf1_pixel);
var testAccuracy = validated.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('RF Validation error matrix: ', testAccuracy);
print('RF Validation overall accuracy: ', testAccuracy.accuracy())
print('RF Validation consumer accuracy: ', consumers);


// // Define the output file name and path
// var outputFileName = 'rf_pixel_classification_result_2017';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 30; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: rf2_pixel, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });
// // Define the output file name and path
// var outputFileName = 'CART_classification_result_2017';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 30; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: CART_classified, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });





//////////////////// for 2018 Landsat 8 classification in GEE////////////////////////////////
var shanghai = ee.FeatureCollection('users/luluzhao2022/shanghai'); // read in Shanghai boundary shp as ROI

// Applies scaling factors.
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2); 
// Landsat Collection 2 surface reflectance has a scale factor of 0.0000275 and an additional offset of -0.2 per pixel.
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

var cloudMaskL457 = function(image) {
  var qa = image.select('QA_PIXEL');
//云层表示为第五位，云层置信度为6-7位，云阴影为第三位
//选择出有云并且云层置信度为中等，以及有云阴影覆盖的像元。
  var cloud = qa.bitwiseAnd(1 << 4).or (qa.bitwiseAnd(1 << 5))
          .and(qa.bitwiseAnd(1 << 7))
          .or(qa.bitwiseAnd(1 << 3));
  // 移除边界像元
  var mask2 = image.mask().reduce(ee.Reducer.min());
  //将检测有关云像元置为0，掩模保留位置不为0的数据。
  return image.updateMask(cloud.not()).updateMask(mask2);
};

var cloudMaskC2L7 = function(image) {
  var dilatedCloud = (1 << 1)
  var cloud = (1 << 3)
  var cloudShadow = (1 << 4)
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(dilatedCloud)
    .and(qa.bitwiseAnd(cloud))
    .or(qa.bitwiseAnd(cloudShadow))
  return image.updateMask(mask.not());
};


// surface reflectance Landsat 8 images
var Landsat_2018_raw = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2018-01-01', '2018-12-31')
    .filterBounds(shanghai)
    .filter(ee.Filter.lt("CLOUD_COVER", 1.5))
    .map(applyScaleFactors);
// var Landsat_2018_waytwo = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
//     .filterDate('2018-01-01', '2018-12-31')
//     .filterBounds(shanghai)
//     .filter(ee.Filter.lt("CLOUD_COVER", 1.5))
//     .map(applyScaleFactors)
//     .map(cloudMaskL457);
// var Landsat_2018_wayone = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
//     .filterDate('2018-01-01', '2018-12-31')
//     .filterBounds(shanghai)
//     .filter(ee.Filter.lt("CLOUD_COVER", 1.5))
//     .map(applyScaleFactors)
//     .map(cloudMaskC2L7);    


//// to get the better images 
// apply the median reducer
// var Landsat_2018_1_median = Landsat_2018_wayone.reduce(ee.Reducer.median());
// var Landsat_2018_2_median = Landsat_2018_waytwo.reduce(ee.Reducer.median());
var Landsat_2018_raw_median = Landsat_2018_raw.reduce(ee.Reducer.median());

var visualization = {
  bands: ['SR_B4_median', 'SR_B3_median', 'SR_B2_median'],
  min: 0.0,
  max: 0.3,
};

var vis1 = {
  bands: ['SR_B4', 'SR_B3', 'SR_B2'],
  min: 0.0,
  max: 0.3,
};
Map.setCenter(121.6720, 31.0030, 8); // focus on Shanghai
// Map.addLayer(Landsat_2018_wayone, vis1, 'Landsat_2018_wayone');
// Map.addLayer(Landsat_2018_waytwo, vis1, 'Landsat_2018_waytwo');
// Map.addLayer(Landsat_2018_1_median, visualization, 'Landsat_2018_1_median');
// Map.addLayer(Landsat_2018_2_median, visualization, 'Landsat_2018_2_median');
Map.addLayer(Landsat_2018_raw, vis1, 'Landsat_2018_raw');
Map.addLayer(Landsat_2018_raw_median, visualization, 'Landsat_2018_2_rawmedian');
Map.addLayer(shanghai);


/// clip images
var Landsat_2018_scale_meanImage_clip = Landsat_2018_raw_median.clip(shanghai)
  .select(['SR_B1_median', 'SR_B2_median', 'SR_B3_median', 'SR_B4_median', 'SR_B5_median', 'SR_B6_median', 'SR_B7_median']);


///填补去云后缺失缺失的部分
//由于2018年影像云覆盖可忽略不计，所以省略这一步

var vis_params3 = {
  bands: ['SR_B5_median', 'SR_B4_median', 'SR_B3_median'],
  min: 0,
  max: 0.3,
};
// map the layer
Map.addLayer(Landsat_2018_scale_meanImage_clip, vis_params3, 'Landsat_2018_scale_meanImage_clip');

//MNDVI, NDVI, and NDBI: 构建NDVI/NDBI/MNDWI作为光谱特征
//利用DEM数据的坡度和高程作为地形特征；然后将构建的这些特征作为影像的一个波段进行使用
var mndwi = Landsat_2018_scale_meanImage_clip.normalizedDifference(['SR_B3_median', 'SR_B6_median']).rename('MNDWI');//计算MNDWI
var ndbi = Landsat_2018_scale_meanImage_clip.normalizedDifference(['SR_B6_median', 'SR_B5_median']).rename('NDBI');//计算NDBI
var ndvi = Landsat_2018_scale_meanImage_clip.normalizedDifference(['SR_B5_median', 'SR_B4_median']).rename('NDVI');//计算NDVI
var strm=ee.Image("USGS/SRTMGL1_003");
var dem=ee.Algorithms.Terrain(strm);
var elevation=dem.select('elevation');
var slope=dem.select('slope');
Landsat_2018_scale_meanImage_clip = Landsat_2018_scale_meanImage_clip.addBands(ndvi).addBands(ndbi).addBands(mndwi).addBands(elevation.rename("ELEVATION")).addBands(slope.rename("SLOPE"));

// Then, select training data
// and make a FeatureCollection from the polygons
var polygons = ee.FeatureCollection([
  ee.Feature(green_space, {'class': 1}),
  ee.Feature(other_area, {'class': 2}), // Add a non-green space area with class '2'
  ee.Feature(water, {'class': 3}),
]);

// Use these bands for classification.
var bands = ['SR_B2_median', 'SR_B3_median', 'SR_B4_median', 'SR_B5_median', 'SR_B6_median','SR_B7_median','MNDWI','NDBI','NDVI','SLOPE', 'ELEVATION'];
// The name of the property on the points storing the class label.
var classProperty = 'class';

// we need to pull out the data from our training areas
var training = Landsat_2018_scale_meanImage_clip.select(bands).sampleRegions({
  collection: polygons,
  properties: [classProperty],
  scale: 30
});
print(training, "training");

/// machine learning classification method 1: CART
// Train a CART classifier
var classifier = ee.Classifier.smileCart().train({
  features: training,
  classProperty: classProperty,
});
// Print some info about the classifier (specific to CART)
print('CART, explained', classifier.explain());
// Classify the image
var CART_classified = Landsat_2018_scale_meanImage_clip.classify(classifier);

Map.addLayer(CART_classified, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_CART");

// Classify the reference data using the trained classifier
var classified_training = training.classify(classifier);
// Calculate the accuracy for each class
var accuracy = classified_training.errorMatrix('class', 'classification');
// Print the confusion matrix
print('CART Confusion Matrix:', accuracy);
// Calculate overall accuracy
var overallAccuracy = accuracy.accuracy();
print('CART Overall Accuracy:', overallAccuracy);

/// machine learning classification method 2: random forest
//Fist, add a column of random uniforms to the training dataset.
var withRandom = polygons.randomColumn('random');
print(withRandom, 'sample_polygon');
// We want to reserve some of the data for testing, to avoid overfitting the model (Roughly 70% training, 30% testing)
var split = 0.7;
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));
print(trainingPartition, "train");
print(testingPartition, "test");

// take samples from image for training and validation  
var training_rf = Landsat_2018_scale_meanImage_clip.select(bands).sampleRegions({
  collection: trainingPartition,
  properties: [classProperty],
  scale: 30,
});
var validation_rf = Landsat_2018_scale_meanImage_clip.select(bands).sampleRegions({
  collection: testingPartition,
  properties: [classProperty],
  scale: 30,
});
// Random Forest Classification
var rf1 = ee.Classifier.smileRandomForest(100)
  .train(training_rf, 'class', bands); //classifier
var rf2 = Landsat_2018_scale_meanImage_clip.classify(rf1);//validation
// Classify the test FeatureCollection.
var test = validation_rf.classify(rf1);

/// evaluate the extraction accuracy
var testAccuracy = test.errorMatrix('class', 'classification');
var consumers = testAccuracy.consumersAccuracy();
print('Validation error matrix: ', testAccuracy);
print('Validation overall accuracy: ', testAccuracy.accuracy());
print('Validation consumer accuracy: ', consumers);
// print('kappa accuracy', consumers.kappa());
Map.addLayer(rf2, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "RF");


/// random forest pixel approach
var pixel_number = 1000;
var green_space_points=ee.FeatureCollection.randomPoints(green_space, pixel_number).map(function(i){
  return i.set({'class': 1})})
var other_area_points=ee.FeatureCollection.randomPoints(other_area, pixel_number).map(function(i){
  return i.set({'class': 2})})
var water_points=ee.FeatureCollection.randomPoints(water, pixel_number).map(function(i){
  return i.set({'class': 3})})
var point_sample=ee.FeatureCollection([green_space_points,
                                  other_area_points,
                                  water_points])
                                  .flatten()
                                  .randomColumn();
// assign 70% of training points to validation 
var split=0.7
var training_sample = point_sample.filter(ee.Filter.lt('random', split));
var validation_sample = point_sample.filter(ee.Filter.gte('random', split));
// take samples from image for training and validation  
var training = Landsat_2018_scale_meanImage_clip.select(bands).sampleRegions({
  collection: training_sample,
  properties: ['class'],
  scale: 30,
});
var validation = Landsat_2018_scale_meanImage_clip.select(bands).sampleRegions({
  collection: validation_sample,
  properties: ['class'],
  scale: 30
});
// Random Forest Classification
var rf1_pixel = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
// Get information about the trained classifier.
print('Results of RF trained classifier', rf1_pixel.explain());

// conduct classification
var rf2_pixel = Landsat_2018_scale_meanImage_clip.classify(rf1_pixel);
Map.addLayer(rf2_pixel, {min: 1, max: 3, 
  palette: ['1c5f2c', 'b3ac9f','466b9f']},
  "RF_pixel");

/// assess accuracy
var trainAccuracy = rf1_pixel.confusionMatrix();
print('Resubstitution error matrix: ', trainAccuracy);
print('Training overall accuracy: ', trainAccuracy.accuracy());
var validated = validation.classify(rf1_pixel);
var testAccuracy = validated.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('RF Validation error matrix: ', testAccuracy);
print('RF Validation overall accuracy: ', testAccuracy.accuracy())
print('RF Validation consumer accuracy: ', consumers);



// // Define the output file name and path
// var outputFileName = 'rf_pixel_classification_result_2018';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 30; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: rf2_pixel, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });

// // Define the output file name and path
// var outputFileName = 'CART_classification_result_2018';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 30; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: CART_classified, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });



//////////////////// for 2019 Sentinel-2 classification in GEE////////////////////////////////
var shanghai = ee.FeatureCollection('users/luluzhao2022/shanghai'); // read in Shanghai boundary shp as ROI

//create the function, since the Sentinel data is scaled by 10000 so in order to display this between 0 and 1
var divide10000 = function(image) {
  return image.divide(10000);
};
var wayone = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterDate('2019-01-01', '2019-12-31')
                  .filterBounds(shanghai)  // Intersecting ROI
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 3)); // less than 1% cloud cover
var wayone_divide = wayone.map(divide10000);
var visualization = {
  min: 0.0,
  max: 0.3,
  bands: ['B4', 'B3', 'B2'],
};
Map.centerObject(shanghai, 10);
Map.addLayer(wayone_divide, visualization, 'wayoneRGB');

function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}
// var waytwo = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
//                   .filterDate('2019-01-01', '2019-12-31')
//                   .filterBounds(shanghai)  // Intersecting ROI
//                   // Pre-filter to get less cloudy granules.
//                   .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',5))
//                   .map(maskS2clouds);
// Map.addLayer(waytwo, visualization, 'waytwoRGB');
Map.addLayer(shanghai);

var way_one_median = wayone_divide.reduce(ee.Reducer.median());
// var way_two_median = waytwo.reduce(ee.Reducer.median());
var vis_params = {
  bands: ['B8_median', 'B4_median', 'B3_median'],
  min: 0.0, 
  max: 0.3,
};
var vis_params_T = {
  bands: ['B4_median', 'B3_median', 'B2_median'],
  min: 0.0, 
  max: 0.3,
};
// since the result of wayone is better, keep wayone to do further process
// Map.addLayer(way_one_median, vis_params, 'way_one_median_FalseColor');
Map.addLayer(way_one_median, vis_params_T, 'way_one_median_TrueColor');
var wayonemedian_clip = way_one_median.clip(shanghai)
Map.addLayer(wayonemedian_clip, vis_params, 'wayonemedian_clip');
// Map.addLayer(way_two_median, vis_params, 'way_two_median_TrueColor');
// var waytwomedian_clip = way_two_median.clip(shanghai)
// Map.addLayer(waytwomedian_clip, vis_params, 'waytwomedianRGB_clip');

// Make a FeatureCollection from the polygons
var polygons = ee.FeatureCollection([
  ee.Feature(green_space, {'class': 1}),
  ee.Feature(other_area, {'class': 2}),
  ee.Feature(water, {'class': 3}),
]);

var mndwi = wayonemedian_clip.normalizedDifference(['B3_median', 'B11_median']).rename('MNDWI');//计算MNDWI
var ndbi = wayonemedian_clip.normalizedDifference(['B11_median', 'B8A_median']).rename('NDBI');//计算NDBI
var ndvi = wayonemedian_clip.normalizedDifference(['B8A_median', 'B4_median']).rename('NDVI');//计算NDVI
var strm=ee.Image("USGS/SRTMGL1_003");
var dem=ee.Algorithms.Terrain(strm);
var elevation=dem.select('elevation');
var slope=dem.select('slope');
wayonemedian_clip = wayonemedian_clip.addBands(ndvi).addBands(ndbi).addBands(mndwi).addBands(elevation.rename("ELEVATION")).addBands(slope.rename("SLOPE"));


// Use these bands for classification.
var bands = ['B2_median', 'B3_median', 'B4_median', 'B5_median', 'B6_median', 'B7_median','B8A_median','MNDWI','NDBI','NDVI','SLOPE', 'ELEVATION'];
// The name of the property on the points storing the class label.
var classProperty = 'class';
// Sample the composite to generate training data.  Note that the
// class label is stored in the 'landcover' property.
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: polygons,
  properties: [classProperty],
  scale: 10
});
print(training, "training");

// Train a CART classifier.
var classifier = ee.Classifier.smileCart().train({
  features: training,
  classProperty: classProperty,
});
// Print some info about the classifier (specific to CART).
print('CART, explained', classifier.explain());
// Classify the image.
var CART_classified = wayonemedian_clip.classify(classifier);
// add output
Map.centerObject(shanghai);
Map.addLayer(CART_classified, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_CART");

// Classify the reference data using the trained classifier
var classified_training = training.classify(classifier);
// Calculate the accuracy for each class
var accuracy = classified_training.errorMatrix('class', 'classification');
// Print the confusion matrix
print('CART Confusion Matrix:', accuracy);
// Calculate overall accuracy
var overallAccuracy = accuracy.accuracy();
print('CART Overall Accuracy:', overallAccuracy);

/// random forest classifier:
// Optionally, do some accuracy assessment.  Fist, add a column of
// random uniforms to the training dataset.
var withRandom = polygons.randomColumn('random');
// We want to reserve some of the data for testing, to avoid overfitting the model.
var split = 0.5;  // Roughly 70% training, 30% testing.
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));
print(trainingPartition, "train")
print(testingPartition, "test")
// take samples from image for training and validation  
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: trainingPartition,
  properties: ['class'],
  scale: 10,
});
var validation = wayonemedian_clip.select(bands).sampleRegions({
  collection: testingPartition,
  properties: ['class'],
  scale: 10
});
// Classification
var rf1 = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
var rf2 = wayonemedian_clip.classify(rf1);
// Classify the test FeatureCollection.
var test = validation.classify(rf1);
var testAccuracy = test.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('Validation error matrix: ', testAccuracy); // accuracy doesn’t make sense
print('Validation overall accuracy: ', testAccuracy.accuracy())
print('Validation consumer accuracy: ', consumers);
Map.addLayer(rf2, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_RF");


/// random forest classification based on pixel scale
var pixel_number= 1000;
var green_space_points=ee.FeatureCollection.randomPoints(green_space, pixel_number).map(function(i){
  return i.set({'class': 1})})
var other_area_points=ee.FeatureCollection.randomPoints(other_area, pixel_number).map(function(i){
  return i.set({'class': 2})})
var water_points=ee.FeatureCollection.randomPoints(water, pixel_number).map(function(i){
  return i.set({'class': 3})})
var point_sample=ee.FeatureCollection([green_space_points,
                                  other_area_points,
                                  water_points])
                                  .flatten()
                                  .randomColumn();
// assign 70% of training points to validation 
var split=0.7
var training_sample = point_sample.filter(ee.Filter.lt('random', split));
var validation_sample = point_sample.filter(ee.Filter.gte('random', split));
// take samples from image for training and validation  
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: training_sample,
  properties: ['class'],
  scale: 10,
});
var validation = wayonemedian_clip.select(bands).sampleRegions({
  collection: validation_sample,
  properties: ['class'],
  scale: 10
});
// Random Forest Classification
var rf1_pixel = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
// Get information about the trained classifier.
print('Results of RF trained classifier', rf1_pixel.explain());
var rf2_pixel = wayonemedian_clip.classify(rf1_pixel);
Map.addLayer(rf2_pixel, {min: 1, max: 3, 
  palette: ['1c5f2c', 'b3ac9f','466b9f']},
  "Green_Space_RF_pixel");
var trainAccuracy = rf1_pixel.confusionMatrix();
print('Resubstitution error matrix: ', trainAccuracy);
print('Training overall accuracy: ', trainAccuracy.accuracy());
var validated = validation.classify(rf1_pixel);
var testAccuracy = validated.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('RF Validation error matrix: ', testAccuracy);
print('RF Validation overall accuracy: ', testAccuracy.accuracy())
print('RF Validation consumer accuracy: ', consumers);

// // Define the output file name and path
// var outputFileName = 'rf_pixel_classification_result_2019';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 15; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: rf2_pixel, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });
// // Define the output file name and path
// var outputFileName = 'CART_classification_result_2019';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 15; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: CART_classified, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });




//////////////////// for 2020 Sentinel-2 classification in GEE////////////////////////////////
var shanghai = ee.FeatureCollection('users/luluzhao2022/shanghai'); // read in Shanghai boundary shp as ROI

//create the function, since the Sentinel data is scaled by 10000 so in order to display this between 0 and 1
var divide10000 = function(image) {
  return image.divide(10000);
};
var wayone = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterDate('2020-01-01', '2020-12-31')
                  .filterBounds(shanghai)  // Intersecting ROI
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 1.5)); // less than 1% cloud cover
var wayone_divide = wayone.map(divide10000);
var visualization = {
  min: 0.0,
  max: 0.3,
  bands: ['B4', 'B3', 'B2'],
};
Map.centerObject(shanghai, 10);
Map.addLayer(wayone_divide, visualization, 'wayoneRGB');

function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}
var waytwo = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterDate('2020-01-01', '2020-12-31')
                  .filterBounds(shanghai)  // Intersecting ROI
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',5))
                  .map(maskS2clouds);
Map.addLayer(waytwo, visualization, 'waytwoRGB');
Map.addLayer(shanghai);

// var way_one_median = wayone_divide.reduce(ee.Reducer.median());
var way_two_median = waytwo.reduce(ee.Reducer.median());
var vis_params = {
  bands: ['B8_median', 'B4_median', 'B3_median'],
  min: 0.0, 
  max: 0.3,
};
var vis_params_T = {
  bands: ['B4_median', 'B3_median', 'B2_median'],
  min: 0.0, 
  max: 0.3,
};
// since the result of wayone is better, keep wayone to do further process
// Map.addLayer(way_one_median, vis_params, 'way_one_median_FalseColor');
// Map.addLayer(way_one_median, vis_params_T, 'way_one_median_TrueColor');
// var wayonemedian_clip = way_one_median.clip(shanghai)
// Map.addLayer(wayonemedian_clip, vis_params, 'wayonemedian_clip');
Map.addLayer(way_two_median, vis_params_T, 'way_two_median_TrueColor');
var waytwomedian_clip = way_two_median.clip(shanghai)
Map.addLayer(waytwomedian_clip, vis_params, 'waytwomedian_clip');

// Make a FeatureCollection from the polygons
var polygons = ee.FeatureCollection([
  ee.Feature(green_space, {'class': 1}),
  ee.Feature(other_area, {'class': 2}),
  ee.Feature(water, {'class': 3}),
]);

var mndwi = waytwomedian_clip.normalizedDifference(['B3_median', 'B11_median']).rename('MNDWI');//计算MNDWI
var ndbi = waytwomedian_clip.normalizedDifference(['B11_median', 'B8A_median']).rename('NDBI');//计算NDBI
var ndvi = waytwomedian_clip.normalizedDifference(['B8A_median', 'B4_median']).rename('NDVI');//计算NDVI
var strm=ee.Image("USGS/SRTMGL1_003");
var dem=ee.Algorithms.Terrain(strm);
var elevation=dem.select('elevation');
var slope=dem.select('slope');
waytwomedian_clip = waytwomedian_clip.addBands(ndvi).addBands(ndbi).addBands(mndwi).addBands(elevation.rename("ELEVATION")).addBands(slope.rename("SLOPE"));


// Use these bands for classification.
var bands = ['B2_median', 'B3_median', 'B4_median', 'B5_median', 'B6_median', 'B7_median','B8A_median','MNDWI','NDBI','NDVI','SLOPE', 'ELEVATION'];
// The name of the property on the points storing the class label.
var classProperty = 'class';
// Sample the composite to generate training data.  Note that the
// class label is stored in the 'landcover' property.
var training = waytwomedian_clip.select(bands).sampleRegions({
  collection: polygons,
  properties: [classProperty],
  scale: 10
});
print(training, "training");

// Train a CART classifier.
var classifier = ee.Classifier.smileCart().train({
  features: training,
  classProperty: classProperty,
});
// Print some info about the classifier (specific to CART).
print('CART, explained', classifier.explain());
// Classify the image.
var CART_classified = waytwomedian_clip.classify(classifier);
// add output
Map.centerObject(shanghai);
Map.addLayer(CART_classified, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_CART");

// Classify the reference data using the trained classifier
var classified_training = training.classify(classifier);
// Calculate the accuracy for each class
var accuracy = classified_training.errorMatrix('class', 'classification');
// Print the confusion matrix
print('CART Confusion Matrix:', accuracy);
// Calculate overall accuracy
var overallAccuracy = accuracy.accuracy();
print('CART Overall Accuracy:', overallAccuracy);

/// random forest classifier:
// Optionally, do some accuracy assessment.  Fist, add a column of
// random uniforms to the training dataset.
var withRandom = polygons.randomColumn('random');
// We want to reserve some of the data for testing, to avoid overfitting the model.
var split = 0.5;  // Roughly 70% training, 30% testing.
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));
print(trainingPartition, "train")
print(testingPartition, "test")
// take samples from image for training and validation  
var training = waytwomedian_clip.select(bands).sampleRegions({
  collection: trainingPartition,
  properties: ['class'],
  scale: 10,
});
var validation = waytwomedian_clip.select(bands).sampleRegions({
  collection: testingPartition,
  properties: ['class'],
  scale: 10
});
// Classification
var rf1 = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
var rf2 = waytwomedian_clip.classify(rf1);
// Classify the test FeatureCollection.
var test = validation.classify(rf1);
var testAccuracy = test.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('Validation error matrix: ', testAccuracy); // accuracy doesn’t make sense
print('Validation overall accuracy: ', testAccuracy.accuracy())
print('Validation consumer accuracy: ', consumers);
Map.addLayer(rf2, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_RF");


/// random forest classification based on pixel scale
var pixel_number= 1000;
var green_space_points=ee.FeatureCollection.randomPoints(green_space, pixel_number).map(function(i){
  return i.set({'class': 1})})
var other_area_points=ee.FeatureCollection.randomPoints(other_area, pixel_number).map(function(i){
  return i.set({'class': 2})})
var water_points=ee.FeatureCollection.randomPoints(water, pixel_number).map(function(i){
  return i.set({'class': 3})})
var point_sample=ee.FeatureCollection([green_space_points,
                                  other_area_points,
                                  water_points])
                                  .flatten()
                                  .randomColumn();
// assign 70% of training points to validation 
var split=0.7
var training_sample = point_sample.filter(ee.Filter.lt('random', split));
var validation_sample = point_sample.filter(ee.Filter.gte('random', split));
// take samples from image for training and validation  
var training = waytwomedian_clip.select(bands).sampleRegions({
  collection: training_sample,
  properties: ['class'],
  scale: 10,
});
var validation = waytwomedian_clip.select(bands).sampleRegions({
  collection: validation_sample,
  properties: ['class'],
  scale: 10
});
// Random Forest Classification
var rf1_pixel = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
// Get information about the trained classifier.
print('Results of RF trained classifier', rf1_pixel.explain());
var rf2_pixel = waytwomedian_clip.classify(rf1_pixel);
Map.addLayer(rf2_pixel, {min: 1, max: 3, 
  palette: ['1c5f2c', 'b3ac9f','466b9f']},
  "Green_Space_RF_pixel");
var trainAccuracy = rf1_pixel.confusionMatrix();
print('Resubstitution error matrix: ', trainAccuracy);
print('Training overall accuracy: ', trainAccuracy.accuracy());
var validated = validation.classify(rf1_pixel);
var testAccuracy = validated.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('RF Validation error matrix: ', testAccuracy);
print('RF Validation overall accuracy: ', testAccuracy.accuracy())
print('RF Validation consumer accuracy: ', consumers);

// // Define the output file name and path
// var outputFileName = 'rf_pixel_classification_result_2020';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 15; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: rf2_pixel, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });
// // Define the output file name and path
// var outputFileName = 'CART_classification_result_2020';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 15; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: CART_classified, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });




//////////////////// for 2021 Sentinel-2 classification in GEE////////////////////////////////
var shanghai = ee.FeatureCollection('users/luluzhao2022/shanghai'); // read in Shanghai boundary shp as ROI

//create the function, since the Sentinel data is scaled by 10000 so in order to display this between 0 and 1
var divide10000 = function(image) {
  return image.divide(10000);
};
var wayone = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterDate('2021-01-01', '2021-12-31')
                  .filterBounds(shanghai)  // Intersecting ROI
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 1.5)); // less than 1% cloud cover
var wayone_divide = wayone.map(divide10000);
var visualization = {
  min: 0.0,
  max: 0.3,
  bands: ['B4', 'B3', 'B2'],
};
Map.centerObject(shanghai, 10);
Map.addLayer(wayone_divide, visualization, 'wayoneRGB');

function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}
// var waytwo = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
//                   .filterDate('2021-01-01', '2021-12-31')
//                   .filterBounds(shanghai)  // Intersecting ROI
//                   // Pre-filter to get less cloudy granules.
//                   .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',5))
//                   .map(maskS2clouds);
// Map.addLayer(waytwo, visualization, 'waytwoRGB');
Map.addLayer(shanghai);

var way_one_median = wayone_divide.reduce(ee.Reducer.median());
// var way_two_median = waytwo.reduce(ee.Reducer.median());
var vis_params = {
  bands: ['B8_median', 'B4_median', 'B3_median'],
  min: 0.0, 
  max: 0.3,
};
var vis_params_T = {
  bands: ['B4_median', 'B3_median', 'B2_median'],
  min: 0.0, 
  max: 0.3,
};
// since the result of wayone is better, keep wayone to do further process
Map.addLayer(way_one_median, vis_params, 'way_one_median_FalseColor');
Map.addLayer(way_one_median, vis_params_T, 'way_one_median_TrueColor');
var wayonemedian_clip = way_one_median.clip(shanghai)
Map.addLayer(wayonemedian_clip, vis_params, 'wayonemedian_clip');
// Map.addLayer(way_two_median, vis_params_T, 'way_two_median_TrueColor');
// var waytwomedian_clip = way_two_median.clip(shanghai)
// Map.addLayer(waytwomedian_clip, vis_params, 'waytwomedian_clip');

// Make a FeatureCollection from the polygons
var polygons = ee.FeatureCollection([
  ee.Feature(green_space, {'class': 1}),
  ee.Feature(other_area, {'class': 2}),
  ee.Feature(water, {'class': 3}),
]);

var mndwi = wayonemedian_clip.normalizedDifference(['B3_median', 'B11_median']).rename('MNDWI');//计算MNDWI
var ndbi = wayonemedian_clip.normalizedDifference(['B11_median', 'B8A_median']).rename('NDBI');//计算NDBI
var ndvi = wayonemedian_clip.normalizedDifference(['B8A_median', 'B4_median']).rename('NDVI');//计算NDVI
var strm=ee.Image("USGS/SRTMGL1_003");
var dem=ee.Algorithms.Terrain(strm);
var elevation=dem.select('elevation');
var slope=dem.select('slope');
wayonemedian_clip = wayonemedian_clip.addBands(ndvi).addBands(ndbi).addBands(mndwi).addBands(elevation.rename("ELEVATION")).addBands(slope.rename("SLOPE"));


// Use these bands for classification.
var bands = ['B2_median', 'B3_median', 'B4_median', 'B5_median', 'B6_median', 'B7_median','B8A_median','MNDWI','NDBI','NDVI','SLOPE', 'ELEVATION'];
// The name of the property on the points storing the class label.
var classProperty = 'class';
// Sample the composite to generate training data.  Note that the
// class label is stored in the 'landcover' property.
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: polygons,
  properties: [classProperty],
  scale: 10
});
print(training, "training");

// Train a CART classifier.
var classifier = ee.Classifier.smileCart().train({
  features: training,
  classProperty: classProperty,
});
// Print some info about the classifier (specific to CART).
print('CART, explained', classifier.explain());
// Classify the image.
var CART_classified = wayonemedian_clip.classify(classifier);
// add output
Map.centerObject(shanghai);
Map.addLayer(CART_classified, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_CART");

// Classify the reference data using the trained classifier
var classified_training = training.classify(classifier);
// Calculate the accuracy for each class
var accuracy = classified_training.errorMatrix('class', 'classification');
// Print the confusion matrix
print('CART Confusion Matrix:', accuracy);
// Calculate overall accuracy
var overallAccuracy = accuracy.accuracy();
print('CART Overall Accuracy:', overallAccuracy);

/// random forest classifier:
// Optionally, do some accuracy assessment.  Fist, add a column of
// random uniforms to the training dataset.
var withRandom = polygons.randomColumn('random');
// We want to reserve some of the data for testing, to avoid overfitting the model.
var split = 0.5;  // Roughly 70% training, 30% testing.
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));
print(trainingPartition, "train")
print(testingPartition, "test")
// take samples from image for training and validation  
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: trainingPartition,
  properties: ['class'],
  scale: 10,
});
var validation = wayonemedian_clip.select(bands).sampleRegions({
  collection: testingPartition,
  properties: ['class'],
  scale: 10
});
// Classification
var rf1 = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
var rf2 = wayonemedian_clip.classify(rf1);
// Classify the test FeatureCollection.
var test = validation.classify(rf1);
var testAccuracy = test.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('Validation error matrix: ', testAccuracy); // accuracy doesn’t make sense
print('Validation overall accuracy: ', testAccuracy.accuracy())
print('Validation consumer accuracy: ', consumers);
Map.addLayer(rf2, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_RF");


/// random forest classification based on pixel scale
var pixel_number= 1000;
var green_space_points=ee.FeatureCollection.randomPoints(green_space, pixel_number).map(function(i){
  return i.set({'class': 1})})
var other_area_points=ee.FeatureCollection.randomPoints(other_area, pixel_number).map(function(i){
  return i.set({'class': 2})})
var water_points=ee.FeatureCollection.randomPoints(water, pixel_number).map(function(i){
  return i.set({'class': 3})})
var point_sample=ee.FeatureCollection([green_space_points,
                                  other_area_points,
                                  water_points])
                                  .flatten()
                                  .randomColumn();
// assign 70% of training points to validation 
var split=0.7
var training_sample = point_sample.filter(ee.Filter.lt('random', split));
var validation_sample = point_sample.filter(ee.Filter.gte('random', split));
// take samples from image for training and validation  
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: training_sample,
  properties: ['class'],
  scale: 10,
});
var validation = wayonemedian_clip.select(bands).sampleRegions({
  collection: validation_sample,
  properties: ['class'],
  scale: 10
});
// Random Forest Classification
var rf1_pixel = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
// Get information about the trained classifier.
print('Results of RF trained classifier', rf1_pixel.explain());
var rf2_pixel = wayonemedian_clip.classify(rf1_pixel);
Map.addLayer(rf2_pixel, {min: 1, max: 3, 
  palette: ['1c5f2c', 'b3ac9f','466b9f']},
  "Green_Space_RF_pixel");
var trainAccuracy = rf1_pixel.confusionMatrix();
print('Resubstitution error matrix: ', trainAccuracy);
print('Training overall accuracy: ', trainAccuracy.accuracy());
var validated = validation.classify(rf1_pixel);
var testAccuracy = validated.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('RF Validation error matrix: ', testAccuracy);
print('RF Validation overall accuracy: ', testAccuracy.accuracy())
print('RF Validation consumer accuracy: ', consumers);

// Define the output file name and path
var outputFileName = 'rf_pixel_classification_result_2021';
var outputFileFolder = 'GEE_exports';
// Set the scale and region of interest
var scale = 15; // Adjust the scale according to your analysis requirements
var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// Export the classified image as GeoTIFF
Export.image.toDrive({
  image: rf2_pixel, // Change to the appropriate classified image variable
  description: outputFileName,
  folder: outputFileFolder,
  fileNamePrefix: outputFileName,
  region: region,
  scale: scale,
  crs: 'EPSG:4326', // Adjust the CRS if needed
  fileFormat: 'GeoTIFF',
});
// Define the output file name and path
var outputFileName = 'CART_classification_result_2021';
var outputFileFolder = 'GEE_exports';
// Set the scale and region of interest
var scale = 15; // Adjust the scale according to your analysis requirements
var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// Export the classified image as GeoTIFF
Export.image.toDrive({
  image: CART_classified, // Change to the appropriate classified image variable
  description: outputFileName,
  folder: outputFileFolder,
  fileNamePrefix: outputFileName,
  region: region,
  scale: scale,
  crs: 'EPSG:4326', // Adjust the CRS if needed
  fileFormat: 'GeoTIFF',
});




//////////////////// for 2022 Sentinel-2 classification in GEE////////////////////////////////
var shanghai = ee.FeatureCollection('users/luluzhao2022/shanghai'); // read in Shanghai boundary shp as ROI

//create the function, since the Sentinel data is scaled by 10000 so in order to display this between 0 and 1
var divide10000 = function(image) {
  return image.divide(10000);
};
var wayone = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterDate('2022-01-01', '2022-12-31')
                  .filterBounds(shanghai)  // Intersecting ROI
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 1)); // less than 1% cloud cover
var wayone_divide = wayone.map(divide10000);
var visualization = {
  min: 0.0,
  max: 0.3,
  bands: ['B4', 'B3', 'B2'],
};
Map.centerObject(shanghai, 10);
Map.addLayer(wayone_divide, visualization, 'wayoneRGB');

function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}
// var waytwo = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
//                   .filterDate('2022-01-01', '2022-12-31')
//                   .filterBounds(shanghai)  // Intersecting ROI
//                   // Pre-filter to get less cloudy granules.
//                   .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',5))
//                   .map(maskS2clouds);
// Map.addLayer(waytwo, visualization, 'waytwoRGB');
Map.addLayer(shanghai);

var way_one_median = wayone_divide.reduce(ee.Reducer.median());
// var way_two_median = waytwo.reduce(ee.Reducer.median());
var vis_params = {
  bands: ['B8_median', 'B4_median', 'B3_median'],
  min: 0.0, 
  max: 0.3,
};
var vis_params_T = {
  bands: ['B4_median', 'B3_median', 'B2_median'],
  min: 0.0, 
  max: 0.3,
};
// since the result of wayone is better, keep wayone to do further process
Map.addLayer(way_one_median, vis_params, 'way_one_median_FalseColor');
Map.addLayer(way_one_median, vis_params_T, 'way_one_median_TrueColor');
var wayonemedian_clip = way_one_median.clip(shanghai)
Map.addLayer(wayonemedian_clip, vis_params, 'wayonemedian_clip');
// Map.addLayer(way_two_median, vis_params_T, 'way_two_median_TrueColor');
// var waytwomedian_clip = way_two_median.clip(shanghai)
// Map.addLayer(waytwomedian_clip, vis_params, 'waytwomedian_clip');

// Make a FeatureCollection from the polygons
var polygons = ee.FeatureCollection([
  ee.Feature(green_space, {'class': 1}),
  ee.Feature(other_area, {'class': 2}),
  ee.Feature(water, {'class': 3}),
]);

var mndwi = wayonemedian_clip.normalizedDifference(['B3_median', 'B11_median']).rename('MNDWI');//计算MNDWI
var ndbi = wayonemedian_clip.normalizedDifference(['B11_median', 'B8A_median']).rename('NDBI');//计算NDBI
var ndvi = wayonemedian_clip.normalizedDifference(['B8A_median', 'B4_median']).rename('NDVI');//计算NDVI
var strm=ee.Image("USGS/SRTMGL1_003");
var dem=ee.Algorithms.Terrain(strm);
var elevation=dem.select('elevation');
var slope=dem.select('slope');
wayonemedian_clip = wayonemedian_clip.addBands(ndvi).addBands(ndbi).addBands(mndwi).addBands(elevation.rename("ELEVATION")).addBands(slope.rename("SLOPE"));


// Use these bands for classification.
var bands = ['B2_median', 'B3_median', 'B4_median', 'B5_median', 'B6_median', 'B7_median','B8A_median','MNDWI','NDBI','NDVI','SLOPE', 'ELEVATION'];
// The name of the property on the points storing the class label.
var classProperty = 'class';
// Sample the composite to generate training data.  Note that the
// class label is stored in the 'landcover' property.
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: polygons,
  properties: [classProperty],
  scale: 10
});
print(training, "training");

// Train a CART classifier.
var classifier = ee.Classifier.smileCart().train({
  features: training,
  classProperty: classProperty,
});
// Print some info about the classifier (specific to CART).
print('CART, explained', classifier.explain());
// Classify the image.
var CART_classified = wayonemedian_clip.classify(classifier);
// add output
Map.centerObject(shanghai);
Map.addLayer(CART_classified, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_CART");

// Classify the reference data using the trained classifier
var classified_training = training.classify(classifier);
// Calculate the accuracy for each class
var accuracy = classified_training.errorMatrix('class', 'classification');
// Print the confusion matrix
print('CART Confusion Matrix:', accuracy);
// Calculate overall accuracy
var overallAccuracy = accuracy.accuracy();
print('CART Overall Accuracy:', overallAccuracy);

/// random forest classifier:
// Optionally, do some accuracy assessment.  Fist, add a column of
// random uniforms to the training dataset.
var withRandom = polygons.randomColumn('random');
// We want to reserve some of the data for testing, to avoid overfitting the model.
var split = 0.5;  // Roughly 70% training, 30% testing.
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));
print(trainingPartition, "train")
print(testingPartition, "test")
// take samples from image for training and validation  
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: trainingPartition,
  properties: ['class'],
  scale: 10,
});
var validation = wayonemedian_clip.select(bands).sampleRegions({
  collection: testingPartition,
  properties: ['class'],
  scale: 10
});
// Classification
var rf1 = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
var rf2 = wayonemedian_clip.classify(rf1);
// Classify the test FeatureCollection.
var test = validation.classify(rf1);
var testAccuracy = test.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('Validation error matrix: ', testAccuracy); // accuracy doesn’t make sense
print('Validation overall accuracy: ', testAccuracy.accuracy())
print('Validation consumer accuracy: ', consumers);
Map.addLayer(rf2, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_RF");


/// random forest classification based on pixel scale
var pixel_number= 1000;
var green_space_points=ee.FeatureCollection.randomPoints(green_space, pixel_number).map(function(i){
  return i.set({'class': 1})})
var other_area_points=ee.FeatureCollection.randomPoints(other_area, pixel_number).map(function(i){
  return i.set({'class': 2})})
var water_points=ee.FeatureCollection.randomPoints(water, pixel_number).map(function(i){
  return i.set({'class': 3})})
var point_sample=ee.FeatureCollection([green_space_points,
                                  other_area_points,
                                  water_points])
                                  .flatten()
                                  .randomColumn();
// assign 70% of training points to validation 
var split=0.7
var training_sample = point_sample.filter(ee.Filter.lt('random', split));
var validation_sample = point_sample.filter(ee.Filter.gte('random', split));
// take samples from image for training and validation  
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: training_sample,
  properties: ['class'],
  scale: 10,
});
var validation = wayonemedian_clip.select(bands).sampleRegions({
  collection: validation_sample,
  properties: ['class'],
  scale: 10
});
// Random Forest Classification
var rf1_pixel = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
// Get information about the trained classifier.
print('Results of RF trained classifier', rf1_pixel.explain());
var rf2_pixel = wayonemedian_clip.classify(rf1_pixel);
Map.addLayer(rf2_pixel, {min: 1, max: 3, 
  palette: ['1c5f2c', 'b3ac9f','466b9f']},
  "Green_Space_RF_pixel");
var trainAccuracy = rf1_pixel.confusionMatrix();
print('Resubstitution error matrix: ', trainAccuracy);
print('Training overall accuracy: ', trainAccuracy.accuracy());
var validated = validation.classify(rf1_pixel);
var testAccuracy = validated.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('RF Validation error matrix: ', testAccuracy);
print('RF Validation overall accuracy: ', testAccuracy.accuracy())
print('RF Validation consumer accuracy: ', consumers);

// // Define the output file name and path
// var outputFileName = 'rf_pixel_classification_result_2022';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 15; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: rf2_pixel, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });
// // Define the output file name and path
// var outputFileName = 'CART_classification_result_2022';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 15; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: CART_classified, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });




//////////////////// for 2023 Sentinel-2 classification in GEE////////////////////////////////
var shanghai = ee.FeatureCollection('users/luluzhao2022/shanghai'); // read in Shanghai boundary shp as ROI

//create the function, since the Sentinel data is scaled by 10000 so in order to display this between 0 and 1
var divide10000 = function(image) {
  return image.divide(10000);
};
var wayone = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterDate('2023-01-01', '2023-07-01')
                  .filterBounds(shanghai)  // Intersecting ROI
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 1)); // less than 1% cloud cover
var wayone_divide = wayone.map(divide10000);
var visualization = {
  min: 0.0,
  max: 0.3,
  bands: ['B4', 'B3', 'B2'],
};
Map.centerObject(shanghai, 10);
Map.addLayer(wayone_divide, visualization, 'wayoneRGB');

function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}
// var waytwo = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
//                   .filterDate('2023-01-01', '2023-07-01')
//                   .filterBounds(shanghai)  // Intersecting ROI
//                   // Pre-filter to get less cloudy granules.
//                   .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',5))
//                   .map(maskS2clouds);
// Map.addLayer(waytwo, visualization, 'waytwoRGB');
Map.addLayer(shanghai);

var way_one_median = wayone_divide.reduce(ee.Reducer.median());
// var way_two_median = waytwo.reduce(ee.Reducer.median());
var vis_params = {
  bands: ['B8_median', 'B4_median', 'B3_median'],
  min: 0.0, 
  max: 0.3,
};
var vis_params_T = {
  bands: ['B4_median', 'B3_median', 'B2_median'],
  min: 0.0, 
  max: 0.3,
};
// since the result of wayone is better, keep wayone to do further process
Map.addLayer(way_one_median, vis_params, 'way_one_median_FalseColor');
Map.addLayer(way_one_median, vis_params_T, 'way_one_median_TrueColor');
var wayonemedian_clip = way_one_median.clip(shanghai)
Map.addLayer(wayonemedian_clip, vis_params, 'wayonemedian_clip');
// Map.addLayer(way_two_median, vis_params_T, 'way_two_median_TrueColor');
// var waytwomedian_clip = way_two_median.clip(shanghai)
// Map.addLayer(waytwomedian_clip, vis_params, 'waytwomedian_clip');

// Make a FeatureCollection from the polygons
var polygons = ee.FeatureCollection([
  ee.Feature(green_space, {'class': 1}),
  ee.Feature(other_area, {'class': 2}),
  ee.Feature(water, {'class': 3}),
]);

var mndwi = wayonemedian_clip.normalizedDifference(['B3_median', 'B11_median']).rename('MNDWI');//计算MNDWI
var ndbi = wayonemedian_clip.normalizedDifference(['B11_median', 'B8A_median']).rename('NDBI');//计算NDBI
var ndvi = wayonemedian_clip.normalizedDifference(['B8A_median', 'B4_median']).rename('NDVI');//计算NDVI
var strm=ee.Image("USGS/SRTMGL1_003");
var dem=ee.Algorithms.Terrain(strm);
var elevation=dem.select('elevation');
var slope=dem.select('slope');
wayonemedian_clip = wayonemedian_clip.addBands(ndvi).addBands(ndbi).addBands(mndwi).addBands(elevation.rename("ELEVATION")).addBands(slope.rename("SLOPE"));


// Use these bands for classification.
var bands = ['B2_median', 'B3_median', 'B4_median', 'B5_median', 'B6_median', 'B7_median','B8A_median','MNDWI','NDBI','NDVI','SLOPE', 'ELEVATION'];
// The name of the property on the points storing the class label.
var classProperty = 'class';
// Sample the composite to generate training data.  Note that the
// class label is stored in the 'landcover' property.
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: polygons,
  properties: [classProperty],
  scale: 10
});
print(training, "training");

// Train a CART classifier.
var classifier = ee.Classifier.smileCart().train({
  features: training,
  classProperty: classProperty,
});
// Print some info about the classifier (specific to CART).
print('CART, explained', classifier.explain());
// Classify the image.
var CART_classified = wayonemedian_clip.classify(classifier);
// add output
Map.centerObject(shanghai);
Map.addLayer(CART_classified, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_CART");

// Classify the reference data using the trained classifier
var classified_training = training.classify(classifier);
// Calculate the accuracy for each class
var accuracy = classified_training.errorMatrix('class', 'classification');
// Print the confusion matrix
print('CART Confusion Matrix:', accuracy);
// Calculate overall accuracy
var overallAccuracy = accuracy.accuracy();
print('CART Overall Accuracy:', overallAccuracy);

/// random forest classifier:
// Optionally, do some accuracy assessment.  Fist, add a column of
// random uniforms to the training dataset.
var withRandom = polygons.randomColumn('random');
// We want to reserve some of the data for testing, to avoid overfitting the model.
var split = 0.5;  // Roughly 70% training, 30% testing.
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));
print(trainingPartition, "train")
print(testingPartition, "test")
// take samples from image for training and validation  
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: trainingPartition,
  properties: ['class'],
  scale: 10,
});
var validation = wayonemedian_clip.select(bands).sampleRegions({
  collection: testingPartition,
  properties: ['class'],
  scale: 10
});
// Classification
var rf1 = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
var rf2 = wayonemedian_clip.classify(rf1);
// Classify the test FeatureCollection.
var test = validation.classify(rf1);
var testAccuracy = test.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('Validation error matrix: ', testAccuracy); // accuracy doesn’t make sense
print('Validation overall accuracy: ', testAccuracy.accuracy())
print('Validation consumer accuracy: ', consumers);
Map.addLayer(rf2, {min: 1, max: 3, palette: ['1c5f2c', 'b3ac9f','466b9f']}, "Green_Space_RF");


/// random forest classification based on pixel scale
var pixel_number= 1000;
var green_space_points=ee.FeatureCollection.randomPoints(green_space, pixel_number).map(function(i){
  return i.set({'class': 1})})
var other_area_points=ee.FeatureCollection.randomPoints(other_area, pixel_number).map(function(i){
  return i.set({'class': 2})})
var water_points=ee.FeatureCollection.randomPoints(water, pixel_number).map(function(i){
  return i.set({'class': 3})})
var point_sample=ee.FeatureCollection([green_space_points,
                                  other_area_points,
                                  water_points])
                                  .flatten()
                                  .randomColumn();
// assign 70% of training points to validation 
var split=0.7
var training_sample = point_sample.filter(ee.Filter.lt('random', split));
var validation_sample = point_sample.filter(ee.Filter.gte('random', split));
// take samples from image for training and validation  
var training = wayonemedian_clip.select(bands).sampleRegions({
  collection: training_sample,
  properties: ['class'],
  scale: 10,
});
var validation = wayonemedian_clip.select(bands).sampleRegions({
  collection: validation_sample,
  properties: ['class'],
  scale: 10
});
// Random Forest Classification
var rf1_pixel = ee.Classifier.smileRandomForest(100)
    .train(training, 'class');
// Get information about the trained classifier.
print('Results of RF trained classifier', rf1_pixel.explain());
var rf2_pixel = wayonemedian_clip.classify(rf1_pixel);
Map.addLayer(rf2_pixel, {min: 1, max: 3, 
  palette: ['1c5f2c', 'b3ac9f','466b9f']},
  "Green_Space_RF_pixel");
var trainAccuracy = rf1_pixel.confusionMatrix();
print('Resubstitution error matrix: ', trainAccuracy);
print('Training overall accuracy: ', trainAccuracy.accuracy());
var validated = validation.classify(rf1_pixel);
var testAccuracy = validated.errorMatrix('class', 'classification');
var consumers=testAccuracy.consumersAccuracy()
print('RF Validation error matrix: ', testAccuracy);
print('RF Validation overall accuracy: ', testAccuracy.accuracy())
print('RF Validation consumer accuracy: ', consumers);

// // Define the output file name and path
// var outputFileName = 'rf_pixel_classification_result_2023';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 15; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: rf2_pixel, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });
// // Define the output file name and path
// var outputFileName = 'CART_classification_result_2023';
// var outputFileFolder = 'GEE_exports';
// // Set the scale and region of interest
// var scale = 15; // Adjust the scale according to your analysis requirements
// var region = shanghai.geometry(); // Assuming 'shanghai' is the region of interest feature collection
// // Export the classified image as GeoTIFF
// Export.image.toDrive({
//   image: CART_classified, // Change to the appropriate classified image variable
//   description: outputFileName,
//   folder: outputFileFolder,
//   fileNamePrefix: outputFileName,
//   region: region,
//   scale: scale,
//   crs: 'EPSG:4326', // Adjust the CRS if needed
//   fileFormat: 'GeoTIFF',
// });
