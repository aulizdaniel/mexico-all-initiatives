// AUX
var auxFuncs = require('users/ctroche/Mapbiomas:0.aux');

// USER
var name_x = 'cts';
var region_name = 3;

var sample_version_in = '2';
var stable_version_in = '2';
var classification_version_out = '2';

// ==============================
// INPUTS
// ==============================
var inputProperties = [
  'nir_median','nir_median_dry','nir_median_wet',
  'swir1_median','swir1_median_dry','swir1_median_wet',
  'ndvi_median','ndvi_median_dry','ndvi_median_wet',
  'ndvi_amp','evi2_median_dry','evi2_median_wet',
  'MNDWI_median_dry','MNDWI_median_wet',
  'LSWI_median_dry','LSWI_median_wet',
  'ndwi_median_dry','ndwi_median_wet',
  'MMRI_median_dry','MMRI_median_wet',
  'CMRI_median_dry','CMRI_median_wet',
  'gvs_median_dry','ndfi_median_dry','ndfi_median_wet',
  'wefi_median_wet','shade_median_dry'
];

// ==============================
// FUNCIÓN ÍNDICES (CLAVE)
// ==============================
function addMangrovIndices(img) {

  var green_dry  = img.select('green_median_dry');
  var green_wet  = img.select('green_median_wet');
  var nir_dry    = img.select('nir_median_dry');
  var nir_wet    = img.select('nir_median_wet');
  var swir1_dry  = img.select('swir1_median_dry');
  var swir1_wet  = img.select('swir1_median_wet');

  var MNDWI_dry = green_dry.subtract(swir1_dry).divide(green_dry.add(swir1_dry));
  var MNDWI_wet = green_wet.subtract(swir1_wet).divide(green_wet.add(swir1_wet));

  var LSWI_dry = nir_dry.subtract(swir1_dry).divide(nir_dry.add(swir1_dry));
  var LSWI_wet = nir_wet.subtract(swir1_wet).divide(nir_wet.add(swir1_wet));

  var ndvi_dry = img.select('ndvi_median_dry').toFloat().divide(65535).multiply(2).subtract(1);
  var ndvi_wet = img.select('ndvi_median_wet').toFloat().divide(65535).multiply(2).subtract(1);

  var MMRI_dry = MNDWI_dry.abs().subtract(ndvi_dry.abs())
                    .divide(MNDWI_dry.abs().add(ndvi_dry.abs()));

  var MMRI_wet = MNDWI_wet.abs().subtract(ndvi_wet.abs())
                    .divide(MNDWI_wet.abs().add(ndvi_wet.abs()));

  var CMRI_dry = ndvi_dry.subtract(MNDWI_dry);
  var CMRI_wet = ndvi_wet.subtract(MNDWI_wet);

  return img.addBands([
    MNDWI_dry.rename('MNDWI_median_dry'),
    MNDWI_wet.rename('MNDWI_median_wet'),
    LSWI_dry.rename('LSWI_median_dry'),
    LSWI_wet.rename('LSWI_median_wet'),
    MMRI_dry.rename('MMRI_median_dry'),
    MMRI_wet.rename('MMRI_median_wet'),
    CMRI_dry.rename('CMRI_median_dry'),
    CMRI_wet.rename('CMRI_median_wet')
  ], null, true);
}

// ==============================
// COLECCIÓN (AQUÍ ESTABA EL ERROR)
// ==============================
var Land_collection = ee.ImageCollection(
  'projects/mapbiomas-mosaics/assets/LANDSAT/LULC/MEXICO/mosaics-1'
).map(addMangrovIndices);   // 🔴 CLAVE

// ==============================
// REGIÓN
// ==============================
var regions = ee.FeatureCollection(
  'projects/mapbiomas-mexico/assets/LAND-COVER/COLLECTION-1/GENERAL/manglar_zonas'
);
var region = regions.filter(ee.Filter.eq('Zona', region_name));

Map.centerObject(region, 8);

// ==============================
// MUESTRAS
// ==============================
var sample_folder = 'projects/mapbiomas-mexico/assets/LAND-COVER/COLLECTION-1/GENERAL/SAMPLES/';

var trainedSamples = ee.FeatureCollection(
  sample_folder + 'trained_Samples_MEX_mangrove_z' +
  region_name + '_' + name_x + '_v' + sample_version_in
).filterBounds(region);

// ==============================
// STABLE
// ==============================
var stableName = sample_folder + 'MEX_Stable_Map_mangrove_z' +
                 region_name + '_' + name_x + '_v' + stable_version_in;

var stableMap = ee.Image(stableName).rename('class');

var stableSamples = stableMap.stratifiedSample({
  numPoints: 0,
  classBand: 'class',
  region: region,
  scale: 30,
  geometries: true
});

// ==============================
// AÑOS
// ==============================
var years_list = trainedSamples.aggregate_array('year')
  .distinct().sort().getInfo();

// ==============================
// LOOP
// ==============================
var classifiedList = years_list.map(function(year){

  var mosaic = Land_collection
    .filter(ee.Filter.eq('year', year))
    .mosaic();

  var samples_year = trainedSamples.filter(ee.Filter.eq('year', year));

  var stableSpectral = mosaic.sampleRegions({
    collection: stableSamples,
    scale: 30,
    geometries: true
  }).filter(ee.Filter.notNull(['red_median']));

  var totalSamples = samples_year.merge(stableSpectral);

  var clf = ee.Classifier.smileRandomForest({
    numberOfTrees: 100,
    seed: 15
  }).train({
    features: totalSamples,
    classProperty: 'class',
    inputProperties: inputProperties
  });

  var classified = mosaic.clip(region)
    .classify(clf)
    .rename('classification_' + year);

  var prob = mosaic.classify(
    clf.setOutputMode('MULTIPROBABILITY')
  );

var probMangrove = prob
  .arrayReduce(ee.Reducer.max(), [0])
  .arrayGet([0])  
  .rename('probability_' + year)
  .updateMask(classified.eq(5));

  return classified.addBands(probMangrove);
});

// ==============================
// STACK
// ==============================
var col = ee.ImageCollection.fromImages(classifiedList);

var classified_final = col.select('classification_.*').toBands()
var probability_raw = col.select('probability_.*').toBands();

// Obtener nombres actuales
var oldNames = probability_raw.bandNames();

// Crear nuevos nombres válidos (quitando índice inicial)
var newNames = oldNames.map(function(name){
  name = ee.String(name);
  // divide por "_" y elimina el primer elemento (índice)
  var parts = name.split('_');
  // reconstruir desde "probability_YYYY"
  return ee.String('probability_').cat(parts.get(-1));
});

// Renombrar
var probability_final = probability_raw.rename(newNames);
// ==============================
// VISUALIZACIÓN — FIX REAL
// ==============================

var visYear = '2025';

// --- Clasificación ---
var class2025 = classified_final
  .select('.*classification_' + visYear)  // 🔴 regex

Map.addLayer(
  class2025,
  {min: 1, max: 55, palette: ['#1a9641','#d7191c','#2c7bb6']},
  'Clasificación ' + visYear,
  true
);

// --- Máscara manglar ---
var mangroveMask = class2025.eq(5);

Map.addLayer(
  mangroveMask.selfMask(),
  {palette: ['red']},
  'Manglar mask ' + visYear,
  false
);

// --- Probabilidad (ya enmascarada en tu pipeline) ---
var prob2025 = probability_final
  .select('.*probability_' + visYear);

Map.addLayer(
  prob2025,
  {min: 0, max: 1, palette: ['white','yellow','red']},
  'Prob Manglar ' + visYear,
  true
);

// --- Debug ---
print('Bandas clasificación:', classified_final.bandNames());
print('Bandas probabilidad:', probability_final.bandNames());

// ==============================
// EXPORT
// ==============================
var class_folder  = 'projects/mapbiomas-mexico/assets/LAND-COVER/COLLECTION-1/GENERAL/classification/mangrove/';
var country_name = 'MEX';

var class_name = country_name + '_classification_mangrove_z' +
                 region_name + '_' + name_x + '_v' + classification_version_out;

var prob_name = country_name + '_probability_mangrove_z' +
                region_name + '_' + name_x + '_v' + classification_version_out;

Export.image.toAsset({
  image: classified_final,
  description: class_name,
  assetId: class_folder + class_name,
  scale: 30,
  pyramidingPolicy: {'.default': 'mode'},
  maxPixels: 1e13,
  region: region.geometry()
});

Export.image.toAsset({
  image: probability_final,
  description: prob_name,
  assetId: class_folder + prob_name,
  scale: 30,
  pyramidingPolicy: {'.default': 'mean'},
  maxPixels: 1e13,
  region: region.geometry()
});
