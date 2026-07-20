// ============================================================================
// PARAMETERS — EDIT HERE
// ============================================================================
/*
var celda_id        = 168;     // ID de la celda a procesar
var version         = 1;       // Versión (debe coincidir con clasificación de entrada)
var TEMPORAL_WINDOW = 2;       // Ventana temporal: ±N años para suavizado


// ============================================================================
// CONFIGURACIÓN TEMPORAL
// ============================================================================

var YEAR_START = 1985;
var YEAR_END   = 2025;
var year_viz   = [1986, 1990, 1995, 2000, 2010, 2020, 2025];

// --- Mosaicos pre-exportados (para visualización Natural Color) ---
var MOSAIC_DIR     = 'projects/mapbiomas-mosaics/assets/LANDSAT/LULC/MEXICO/Urban/COLLECTION-1/MOSAICS/CENTRO-SUR';
var MOSAIC_PREFIX  = 'mosaic_mexico_centro_sur_urban_';
var MOSAIC_VERSION = 1;

// ============================================================================
// FIN PARÁMETROS — NO MODIFICAR A PARTIR DE AQUÍ
// ============================================================================

// --- Rutas ---
var paths = {
  grid:            'projects/mapbiomas-mexico/assets/Urban/COLLECTION-1/Samples/malla_geoestadistica_sel_id_zona_vecinos',
  inputProba:      'projects/mapbiomas-mexico/assets/Urban/COLLECTION-1/Probabilities/',
  inputThresholds: 'projects/mapbiomas-mexico/assets/Urban/COLLECTION-1/Umbrales/',
  outputProba:     'projects/mapbiomas-mexico/assets/Urban/COLLECTION-1/Probabilities_Harmonized/',
  outputClass:     'projects/mapbiomas-mexico/assets/Urban/COLLECTION-1/Classification_Harmonized/'
};

// --- Cargar períodos y umbrales desde asset de clasificación ---
// (FeatureCollection pequeña: getInfo() es seguro y se ejecuta una sola vez)
var thresholdAssetId = paths.inputThresholds + 'umbrales_' + celda_id +
                       '_' + YEAR_START + '_' + YEAR_END + '_v' + version;

var thresholdsInfo = ee.FeatureCollection(thresholdAssetId)
  .sort('year_start')
  .toList(100)
  .getInfo();

if (thresholdsInfo.length === 0) {
  throw new Error('Asset de umbrales vacío o no encontrado: ' + thresholdAssetId);
}

var PERIODOS = thresholdsInfo.map(function(feat) {
  return {
    yearStart: feat.properties.year_start,
    yearEnd:   feat.properties.year_end,
    threshold: feat.properties.threshold,
    label:     feat.properties.periodo_label
  };
});

// --- Lista de años a procesar ---
var YEARS = [];
for (var y = YEAR_START; y <= YEAR_END; y++) {
  YEARS.push(y);
}

// --- Helper: encontrar umbral correspondiente a un año ---
function getThresholdForYear(year) {
  for (var i = 0; i < PERIODOS.length; i++) {
    if (year >= PERIODOS[i].yearStart && year <= PERIODOS[i].yearEnd) {
      return PERIODOS[i].threshold;
    }
  }
  throw new Error('No se encontró período para año ' + year);
}

// --- Helper: cargar mosaico Landsat de un año (solo bandas RGB para visualización) ---
function getMosaic(year) {
  var assetId = MOSAIC_DIR + '/' + MOSAIC_PREFIX + year + '_v' + MOSAIC_VERSION;
  return ee.Image(assetId).select(['RED', 'GREEN', 'BLUE']);
}

// ============================================================================
// 1. GEOMETRÍA DE CELDA
// ============================================================================

var gridMx    = ee.FeatureCollection(paths.grid);
var geomCelda = gridMx.filter(ee.Filter.eq('id', celda_id)).geometry();

print('═══════════════════════════════════════════════════════');
print('  ARMONIZACIÓN TEMPORAL — CELDA ' + celda_id);
print('  Ventana: ±' + TEMPORAL_WINDOW + ' años');
print('  Período: ' + YEAR_START + '–' + YEAR_END + ' (' + YEARS.length + ' años)');
print('  Versión entrada/salida: v' + version);
print('═══════════════════════════════════════════════════════');

print('──────────────────────────────────────────');
print('UMBRALES POR PERÍODO (heredados del clasificador):');
PERIODOS.forEach(function(p) {
  print('  ' + p.label + '  →  umbral = ' + p.threshold +
        ' (' + (p.threshold * 100) + '%)');
});

// ============================================================================
// 2. CARGAR PROBABILIDADES ORIGINALES (MULTIBANDA → COLECCIÓN)
// ============================================================================

var probaAssetId  = paths.inputProba + 'proba_' + celda_id +
                    '_' + YEAR_START + '_' + YEAR_END + '_v' + version;
var probaOriginal = ee.Image(probaAssetId);

print('──────────────────────────────────────────');
print('Cargando probabilidades de:');
print('  ' + probaAssetId);


// Convierte una imagen multibanda con bandas <prefix>YYYY a ImageCollection
// con una imagen por año y propiedad 'year' asignada (necesaria para join).

function multibandToCollection(img, prefix) {
  var images = YEARS.map(function(year) {
    return img.select(prefix + year)
      .rename('probability')
      .set('year', year)
      .set('celda_id', celda_id);
  });
  return ee.ImageCollection.fromImages(images);
}

var probaColOriginal = multibandToCollection(probaOriginal, 'probability_');

print('  ✓ Probabilidades originales:', probaColOriginal.size(), 'imágenes');

// ============================================================================
// 3. ARMONIZACIÓN TEMPORAL (MEDIA MÓVIL ±N AÑOS)
// ============================================================================


function temporalHarmonization(probCol) {
  print('Aplicando armonización con ventana ±' + TEMPORAL_WINDOW + ' años...');

  var join = ee.Join.saveAll({matchesKey: 'images'});

  var filter = ee.Filter.maxDifference({
    difference: TEMPORAL_WINDOW,
    leftField:  'year',
    rightField: 'year'
  });

  var joined = join.apply(probCol, probCol, filter);

  return ee.ImageCollection(joined.map(function(image) {
    var year      = image.get('year');
    var neighbors = ee.ImageCollection.fromImages(ee.List(image.get('images')));
    var meanProb  = neighbors.reduce(ee.Reducer.mean()).rename('probability');

    return meanProb.set({
      'year':            year,
      'celda_id':        celda_id,
      'processing':      'temporally_harmonized',
      'temporal_window': TEMPORAL_WINDOW,
      'n_images_used':   neighbors.size()
    });
  }));
}

var probaColHarmonized = temporalHarmonization(probaColOriginal);

// ============================================================================
// 4. RECOMPOSICIÓN A IMAGEN MULTIBANDA (PROBABILIDAD ARMONIZADA)
// ============================================================================

function collectionToMultiband(col, bandPrefix) {
  var bandNames = YEARS.map(function(year) {
    return bandPrefix + year;
  });

  var stack = ee.ImageCollection(YEARS.map(function(year) {
    return col.filter(ee.Filter.eq('year', year)).first()
      .rename(bandPrefix + year);
  })).toBands().rename(bandNames);

  return stack;
}

var probaHarmImage = collectionToMultiband(probaColHarmonized, 'probability_');

// ============================================================================
// 5. RECLASIFICACIÓN BINARIA CON UMBRALES POR PERÍODO
// ============================================================================

var classHarmBands = YEARS.map(function(year) {
  var threshold = getThresholdForYear(year);
  var thPct     = threshold * 100;

  return probaHarmImage.select('probability_' + year)
    .gte(thPct)
    .rename('classification_' + year);
});

var classHarmImage = ee.Image.cat(classHarmBands);

// ============================================================================
// 6. METADATOS Y CASTEO FINAL
// ============================================================================

probaHarmImage = probaHarmImage.toByte().set({
  'celda_id':        celda_id,
  'version':         version,
  'temporal_window': TEMPORAL_WINDOW,
  'processing':      'temporal_harmonization',
  'year_start':      YEAR_START,
  'year_end':        YEAR_END,
  'source_asset':    probaAssetId
});

classHarmImage = classHarmImage.toByte().set({
  'celda_id':                 celda_id,
  'version':                  version,
  'temporal_window':          TEMPORAL_WINDOW,
  'processing':               'temporal_harmonization_thresholded',
  'thresholds':               PERIODOS.map(function(p) { return p.threshold; }),
  'thresholds_source_asset':  thresholdAssetId,
  'year_start':               YEAR_START,
  'year_end':                 YEAR_END,
  'source_asset':             probaAssetId
});

// ============================================================================
// 7. DIAGNÓSTICO / CONTROL DE CALIDAD
// ============================================================================

print('──────────────────────────────────────────');
print('CONTROL DE CALIDAD:');
print('  Imágenes armonizadas:', probaColHarmonized.size());
print('  Años esperados:', YEARS.length);

print('──────────────────────────────────────────');
print('Vecinas usadas por año (efecto de borde de la ventana):');
[YEAR_START, YEAR_START + 1, YEAR_START + 2,
 YEAR_END - 2, YEAR_END - 1, YEAR_END].forEach(function(yr) {
  var img = probaColHarmonized.filter(ee.Filter.eq('year', yr)).first();
  print('  ' + yr + ' →', img.get('n_images_used'));
});

// ============================================================================
// 8. VISUALIZACIÓN COMPARATIVA
// ============================================================================

Map.centerObject(geomCelda, 10);

Map.addLayer(
  gridMx.filter(ee.Filter.eq('id', celda_id))
    .style({color: 'FFFFFF', fillColor: '00000000', width: 2.5}),
  {}, 'Celda ' + celda_id
);

var probaViz = {min: 0, max: 100, palette: ['68ff0a', 'fbff08', 'ff3406']};
var diffViz  = {min: -30, max: 30, palette: ['0000ff', 'ffffff', 'ff0000']};

year_viz.forEach(function(year) {
  var probOrig = probaColOriginal.filter(ee.Filter.eq('year', year)).first().clip(geomCelda);
  var probHarm = probaColHarmonized.filter(ee.Filter.eq('year', year)).first().clip(geomCelda);
  var diff     = probHarm.subtract(probOrig).rename('diff');

  var thPct   = getThresholdForYear(year) * 100;
  var binOrig = probOrig.gte(thPct);
  var binHarm = probHarm.gte(thPct);

  // Mosaico del año — Natural Color (RED-GREEN-BLUE)
  Map.addLayer(
    getMosaic(year).clip(geomCelda),
    {bands: ['RED', 'GREEN', 'BLUE'], min: 1004.18, max: 1204.82, gamma: 1.3},
    'Mosaico ' + year + ' (RED-GREEN-BLUE)',
    false
  );

  Map.addLayer(probOrig, probaViz, 'Prob. original ' + year, false);
  Map.addLayer(probHarm, probaViz, 'Prob. armonizada ' + year, false);
  Map.addLayer(diff,     diffViz,  'Diff (armon − orig) ' + year, false);
  Map.addLayer(binOrig.selfMask(), {palette: ['FF0000']},
               'Bin. original ' + year + ' (≥' + thPct + '%)', false);
  Map.addLayer(binHarm.selfMask(), {palette: ['FFAA00']},
               'Bin. armonizada ' + year + ' (≥' + thPct + '%)', false);
});

// ============================================================================
// 9. LEYENDA
// ============================================================================

var legend = ui.Panel({
  style: {
    position:        'bottom-left',
    padding:         '8px 12px',
    backgroundColor: 'white'
  }
});

legend.add(ui.Label({
  value: 'Armonización temporal — Celda ' + celda_id,
  style: {fontWeight: 'bold', fontSize: '13px', margin: '0 0 6px 0'}
}));

legend.add(ui.Label({
  value: 'Ventana: ±' + TEMPORAL_WINDOW + ' años',
  style: {fontSize: '11px', color: '666', margin: '0 0 6px 0', fontStyle: 'italic'}
}));

// Barra de probabilidad
var gradient = ui.Thumbnail({
  image: ee.Image.pixelLonLat().select('longitude')
    .multiply(100 / 2)
    .visualize({min: 0, max: 100, palette: ['68ff0a', 'fbff08', 'ff3406']}),
  params: {bbox: [0, 0, 2, 0.15], dimensions: '200x15'},
  style:  {stretch: 'horizontal', margin: '0 0 4px 0'}
});
legend.add(ui.Label({value: 'Probabilidad (0–100%)',
                     style: {fontSize: '11px', margin: '0'}}));
legend.add(gradient);

// Barra de diferencia
var diffBar = ui.Thumbnail({
  image: ee.Image.pixelLonLat().select('longitude')
    .multiply(60 / 2).subtract(30)
    .visualize({min: -30, max: 30, palette: ['0000ff', 'ffffff', 'ff0000']}),
  params: {bbox: [0, 0, 2, 0.15], dimensions: '200x15'},
  style:  {stretch: 'horizontal', margin: '4px 0 4px 0'}
});
legend.add(ui.Label({value: 'Diff armon − orig (−30 a +30)',
                     style: {fontSize: '11px', margin: '0'}}));
legend.add(diffBar);

function legendEntry(color, label) {
  var box = ui.Label({style: {
    backgroundColor: '#' + color, padding: '8px',
    margin: '4px 6px 4px 0', border: '1px solid #ccc'
  }});
  var text = ui.Label({value: label, style: {margin: '0 0 4px 0', fontSize: '12px'}});
  return ui.Panel({widgets: [box, text], layout: ui.Panel.Layout.Flow('horizontal')});
}

legend.add(legendEntry('FF0000', 'Bin. original'));
legend.add(legendEntry('FFAA00', 'Bin. armonizada'));

Map.add(legend);

// ============================================================================
// 10. EXPORTS
// ============================================================================
*/
/*
var baseName = celda_id + '_' + YEAR_START + '_' + YEAR_END + '_v' + version;

Export.image.toAsset({
  image:            probaHarmImage,
  description:      'proba_harm_' + baseName,
  assetId:          paths.outputProba + 'proba_harm_' + baseName,
  region:           geomCelda,
  scale:            30,
  maxPixels:        1e13,
  pyramidingPolicy: {'.default': 'mean'}
});

Export.image.toAsset({
  image:            classHarmImage,
  description:      'class_harm_' + baseName,
  assetId:          paths.outputClass + 'class_harm_' + baseName,
  region:           geomCelda,
  scale:            30,
  maxPixels:        1e13,
  pyramidingPolicy: {'.default': 'mode'}
});

print('──────────────────────────────────────────');
print('  Exportaciones programadas:');
print('    → proba_harm_' + baseName + '  (' + YEARS.length + ' bandas, byte)');
print('    → class_harm_' + baseName + '  (' + YEARS.length + ' bandas, byte)');
print('═══════════════════════════════════════════════════════');
print('  ARMONIZACIÓN COMPLETA');
print('═══════════════════════════════════════════════════════');*/


// ============================================================================
// MAPBIOMAS MEXICO — TEMPORAL HARMONIZATION IN BATCH (CENTRO-SUR)
// Collection 1 · v0.01
// Processes all cells with exported probabilities and queues a harmonized
// export per cell to the Harmonized_probabilities folder
// ============================================================================

// ============================================================================
// PARAMETERS — EDIT HERE
// ============================================================================

var version           = 1;     // Version of input and output assets
var TEMPORAL_WINDOW   = 2;     // ±N-year window for the moving average
var SKIP_EXISTING     = true;  // Skip already-harmonized cells (idempotence)
var MAX_CELLS_PER_RUN = 0;     // Defensive cap (0 = no limit)

// ============================================================================
// TEMPORAL CONFIGURATION
// ============================================================================

var YEAR_START = 1985;
var YEAR_END   = 2025;

var YEARS = [];
for (var y = YEAR_START; y <= YEAR_END; y++) {
  YEARS.push(y);
}

// ============================================================================
// PATHS
// ============================================================================

var paths = {
  grid:        'projects/mapbiomas-mexico/assets/Urban/COLLECTION-1/Samples/malla_geoestadistica_sel_id_zona_vecinos',
  inputProba:  'projects/mapbiomas-mexico/assets/Urban/COLLECTION-1/Probabilities/',
  outputProba: 'projects/mapbiomas-mexico/assets/Urban/COLLECTION-1/Harmonized_probabilities/'
};

// ============================================================================
// DISCOVERY HELPERS
// ============================================================================

/** Lists all assets in a folder, handling pagination (1000/call). */
function listAllAssets(folder) {
  var allAssets = [];
  var pageToken = null;
  var iterations = 0;
  var MAX_ITER = 10;
  do {
    var result = pageToken
      ? ee.data.listAssets(folder, {pageToken: pageToken})
      : ee.data.listAssets(folder);
    allAssets = allAssets.concat(result.assets || []);
    pageToken = result.nextPageToken;
    iterations++;
  } while (pageToken && iterations < MAX_ITER);
  return allAssets;
}

/** Filters assets whose id ends exactly in _v<version>. */
function filterByVersion(assets, versionStr) {
  var suffix = '_v' + versionStr;
  return assets.filter(function(a) {
    return a.id.slice(-suffix.length) === suffix;
  });
}

/** Extracts celda_id from an asset id like '.../proba_<id>_1985_2025_v1'. */
function extractCellIdFromProba(assetId) {
  var match = assetId.match(/proba_(\d+)_\d{4}_\d{4}_v/);
  return match ? parseInt(match[1], 10) : null;
}

/** Extracts celda_id from an asset id like '.../proba_harm_<id>_1985_2025_v1'. */
function extractCellIdFromHarm(assetId) {
  var match = assetId.match(/proba_harm_(\d+)_\d{4}_\d{4}_v/);
  return match ? parseInt(match[1], 10) : null;
}

// ============================================================================
// DISCOVER CELLS TO PROCESS
// ============================================================================

print('═══════════════════════════════════════════════════════');
print('  ARMONIZACIÓN TEMPORAL — BATCH');
print('  Ventana: ±' + TEMPORAL_WINDOW + ' años');
print('  Versión: v' + version);
print('  Output: ' + paths.outputProba);
print('═══════════════════════════════════════════════════════');

// Input: cells with exported probabilities
var probaAssets = filterByVersion(listAllAssets(paths.inputProba), version);
var allCellIds  = probaAssets
  .map(function(a) { return extractCellIdFromProba(a.id); })
  .filter(function(id) { return id !== null; });

print('Probabilidades de entrada disponibles:', allCellIds.length);

// Output: already-harmonized cells (to skip)
var existingCellIds = [];
if (SKIP_EXISTING) {
  try {
    var existingAssets = filterByVersion(listAllAssets(paths.outputProba), version);
    existingCellIds = existingAssets
      .map(function(a) { return extractCellIdFromHarm(a.id); })
      .filter(function(id) { return id !== null; });
    print('Ya armonizadas (se saltan):', existingCellIds.length);
  } catch (e) {
    print('Output folder vacío o no existe aún — se crearán todas.');
  }
}

// Final list to process
var cellsToProcess = allCellIds.filter(function(id) {
  return existingCellIds.indexOf(id) === -1;
});

if (MAX_CELLS_PER_RUN > 0 && cellsToProcess.length > MAX_CELLS_PER_RUN) {
  print('⚠ Limitando a primeras ' + MAX_CELLS_PER_RUN + ' de ' +
        cellsToProcess.length + ' celdas (MAX_CELLS_PER_RUN)');
  cellsToProcess = cellsToProcess.slice(0, MAX_CELLS_PER_RUN);
}

print('──────────────────────────────────────────');
print('Celdas a procesar (' + cellsToProcess.length + '):', cellsToProcess);

// ============================================================================
// HARMONIZATION FUNCTIONS
// ============================================================================

var gridMx = ee.FeatureCollection(paths.grid);

/** Converts a multiband image <prefix>YYYY to an ImageCollection (1 image/year). */
function multibandToCollection(img, prefix, celdaId) {
  var images = YEARS.map(function(year) {
    return img.select(prefix + year)
      .rename('probability')
      .set('year', year)
      .set('celda_id', celdaId);
  });
  return ee.ImageCollection.fromImages(images);
}

/** Moving average ±TEMPORAL_WINDOW years via a temporal join. */
function temporalHarmonization(probCol, celdaId) {
  var join = ee.Join.saveAll({matchesKey: 'images'});
  var filter = ee.Filter.maxDifference({
    difference: TEMPORAL_WINDOW,
    leftField:  'year',
    rightField: 'year'
  });
  var joined = join.apply(probCol, probCol, filter);

  return ee.ImageCollection(joined.map(function(image) {
    var year      = image.get('year');
    var neighbors = ee.ImageCollection.fromImages(ee.List(image.get('images')));
    var meanProb  = neighbors.reduce(ee.Reducer.mean()).rename('probability');

    return meanProb.set({
      'year':            year,
      'celda_id':        celdaId,
      'processing':      'temporally_harmonized',
      'temporal_window': TEMPORAL_WINDOW,
      'n_images_used':   neighbors.size()
    });
  }));
}

/** Recomposes an ImageCollection into a multiband image with bands <prefix>YYYY. */
function collectionToMultiband(col, bandPrefix) {
  var bandNames = YEARS.map(function(year) { return bandPrefix + year; });
  var stack = ee.ImageCollection(YEARS.map(function(year) {
    return col.filter(ee.Filter.eq('year', year)).first().rename(bandPrefix + year);
  })).toBands().rename(bandNames);
  return stack;
}

// ============================================================================
// PER-CELL PROCESSING
// ============================================================================

function processCell(celdaId) {
  var probaAssetId = paths.inputProba + 'proba_' + celdaId +
                     '_' + YEAR_START + '_' + YEAR_END + '_v' + version;

  var probaOriginal      = ee.Image(probaAssetId);
  var probaColOriginal   = multibandToCollection(probaOriginal, 'probability_', celdaId);
  var probaColHarmonized = temporalHarmonization(probaColOriginal, celdaId);

  var probaHarmImage = collectionToMultiband(probaColHarmonized, 'probability_')
    .toByte()
    .set({
      'celda_id':        celdaId,
      'version':         version,
      'temporal_window': TEMPORAL_WINDOW,
      'processing':      'temporal_harmonization',
      'year_start':      YEAR_START,
      'year_end':        YEAR_END,
      'source_asset':    probaAssetId
    });

  var geomCelda = gridMx.filter(ee.Filter.eq('id', celdaId)).geometry();
  var baseName  = celdaId + '_' + YEAR_START + '_' + YEAR_END + '_v' + version;

  Export.image.toAsset({
    image:            probaHarmImage,
    description:      'proba_harm_' + baseName,
    assetId:          paths.outputProba + 'proba_harm_' + baseName,
    region:           geomCelda,
    scale:            30,
    maxPixels:        1e13,
    pyramidingPolicy: {'.default': 'mean'}
  });
}

// ============================================================================
// EXECUTION
// ============================================================================

if (cellsToProcess.length === 0) {
  print('──────────────────────────────────────────');
  print('  Nada que procesar.');
  print('  Si esperabas tareas: revisar SKIP_EXISTING o el folder de entrada.');
} else {
  cellsToProcess.forEach(function(celdaId, idx) {
    print('  [' + (idx + 1) + '/' + cellsToProcess.length + '] → celda ' + celdaId);
    processCell(celdaId);
  });

  print('──────────────────────────────────────────');
  print('  ✓ Tareas encoladas: ' + cellsToProcess.length);
  print('  → Pestaña "Tasks" → Run para iniciar cada una');
  print('  (o usa el botón "Run all" si aparece)');
}

print('═══════════════════════════════════════════════════════');
