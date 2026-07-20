
// Modules
var mosaicProd = require("users/edimilsonrodriguessantos/mapbiomas:Col10/classificacao/mosaic_production.js");

// ============================================================================
// PARAMETERS
// ============================================================================

var version = 1;

// Output folder
var dirout = 'projects/mapbiomas-mexico/assets/Urban/COLLECTION-1/MOSAICS/CENTRO-SUR';


// Mosaic area


var gridMx = ee.FeatureCollection(
  'projects/mapbiomas-mexico/assets/Urban/COLLECTION-1/Samples/malla_geoestadistica_sel_id_zona_vecinos'
);

// Subzone within the grid
var SUBZONA = 'centro-sur';

// Short name for export
var regionName = 'centro_sur';

var gridCentroSur = gridMx.filter(ee.Filter.eq('subzonas', SUBZONA));

var region = gridCentroSur.geometry();


// Example: range 1985–2005
var activeYears = ee.List.sequence(2006, 2025).getInfo();


// ============================================================================
// OPTIONAL VISUALIZATION
// ============================================================================

/*
var test_year = 1985;
var testMosaic = mosaicProd.mosaicGen(test_year, region);

Map.centerObject(gridCentroSur, 6);
Map.addLayer(gridCentroSur, {color: 'red'}, 'Subzone');
Map.addLayer(
  testMosaic.clip(region),
  {bands: ['RED', 'GREEN', 'BLUE'], min: 1009, max: 1152},
  'Natural color ' + test_year
);
*/

// print('Active years:', activeYears);
// print('Subzone:', SUBZONA);
// print('Number of features:', gridCentroSur.size());



activeYears.forEach(function(year) {

  var mosaic = mosaicProd.mosaicGen(year, region).clip(region);

  var description = 'mosaic_mexico_' + regionName + '_urban_' + year + '_v' + version;
  var assetId = dirout + '/' + description;

  Export.image.toAsset({
    image: mosaic.set({
      'year': year,
      'version': version,
      'territory': 'MEXICO',
      'theme': 'Urban Area',
      'collection': 1,
      'subzona': SUBZONA
    }),
    description: description,
    assetId: assetId,
    region: region,
    scale: 30,
    maxPixels: 1e13,
    pyramidingPolicy: {'.default': 'mean'}
  });

  print('Tarea lanzada:', description);
});
