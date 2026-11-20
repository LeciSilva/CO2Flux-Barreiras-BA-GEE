/********************************************************************
/TCC LECI - ÍNDICE ESPECTRAL ASSOCIADO À ASSIMILAÇÃO DE CO2 NA ÁREA URBANA DE BARREIRAS
 //Período: 2016 a 2025
 //Sensor: Sentinel-2 Harmonized TOA
  */
// ================================================================
// 1. CONFIGURAÇÕES GERAIS
// ================================================================

var anoInicial = 2016;
var anoFinal = 2025;
var escala = 10;
var limiarNuvem = 40;
var nAnos = 10;

var roi = ee.FeatureCollection(
  'projects/elvisbergmm/assets/Limites_Bairros'
).geometry();

// ================================================================
// 2. COLEÇÕES SENTINEL-2 E PROBABILIDADE DE NUVENS
// ================================================================

var sentinel = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
  .filterBounds(roi)
  .filterDate(anoInicial + '-01-01', (anoFinal + 1) + '-01-01')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80));

var probNuvem = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
  .filterBounds(roi)
  .filterDate(anoInicial + '-01-01', (anoFinal + 1) + '-01-01');

// ================================================================
// 3. VINCULAR IMAGENS À PROBABILIDADE DE NUVENS (JUNÇÃO)
// ================================================================

var juncoes = ee.ImageCollection(
  ee.Join.saveFirst('cloud_mask').apply({
    primary: sentinel,
    secondary: probNuvem,
    condition: ee.Filter.equals({
      leftField: 'system:index',
      rightField: 'system:index'
    })
  })
).filter(ee.Filter.notNull(['cloud_mask']));

// ================================================================
// 4. MÁSCARA DE NUVENS E BORDAS
// ================================================================

function mascararNuvens(image) {
  var prob = ee.Image(image.get('cloud_mask')).select('probability');
  var mascara = prob.lt(limiarNuvem);

  mascara = mascara
    .and(image.select('B8A').mask())
    .and(image.select('B9').mask());

  return image
    .updateMask(mascara)
    .select(['B2', 'B3', 'B4', 'B8'])
    .divide(10000)
    .copyProperties(image, ['system:time_start']);
}

var sentinelLimpo = juncoes.map(mascararNuvens);

// ================================================================
// 5. COMPOSTOS ANUAIS E ÍNDICES ESPECTRAIS
// ================================================================

var anos = ee.List.sequence(anoInicial, anoFinal);

var colecaoAnual = ee.ImageCollection.fromImages(
  anos.map(function(ano) {
    ano = ee.Number(ano);
    var inicio = ee.Date.fromYMD(ano, 1, 1);
    var fim = inicio.advance(1, 'year');

    var imagensAno = sentinelLimpo.filterDate(inicio, fim);
    var composto = imagensAno.median().clip(roi);

    var ndvi = composto.normalizedDifference(['B8', 'B4']).rename('NDVI');
    var pri = composto.normalizedDifference(['B3', 'B2']).rename('PRI_proxy');
    var spri = pri.add(1).divide(2).rename('sPRI');
    var co2 = ndvi.multiply(spri).rename('CO2_Index');

    return co2
      .addBands(ndvi)
      .addBands(pri)
      .addBands(spri)
      .set('ano', ano)
      .set('num_cenas', imagensAno.size())
      .set('system:time_start', inicio.millis());
  })
).sort('ano');

// =========================================================================================
// 6. ESTATÍSTICAS ANUAIS (OTIMIZADAS PARA EXPORTAÇÃO) E QUANTIFICAÇÃO DE CENAS UTILIZADAS
// =========================================================================================

var tabelaAnual = ee.FeatureCollection(
  anos.map(function(ano) {
    ano = ee.Number(ano);
    var imagemAno = ee.Image(colecaoAnual.filter(ee.Filter.eq('ano', ano)).first()).select('CO2_Index');

    // tileScale aumentado para 16 evita estouro de memória no cálculo estatístico
    var stats = imagemAno.reduceRegion({
      reducer: ee.Reducer.mean().combine({
        reducer2: ee.Reducer.stdDev(), sharedInputs: true
      }).combine({
        reducer2: ee.Reducer.min(), sharedInputs: true
      }).combine({
        reducer2: ee.Reducer.max(), sharedInputs: true
      }),
      geometry: roi,
      scale: escala,
      maxPixels: 1e13,
      bestEffort: true,
      tileScale: 16
    });

    return ee.Feature(null, {
      'ano': ano,
      'num_cenas': imagemAno.get('num_cenas'),
      'media': stats.get('CO2_Index_mean'),
      'desvio_padrao': stats.get('CO2_Index_stdDev'),
      'minimo': stats.get('CO2_Index_min'),
      'maximo': stats.get('CO2_Index_max')
    });
  })
);

// Exibe no Console o total de cenas válidas no decênio e a distribuição anual
print('Total de cenas Sentinel-2 utilizadas (2016-2025):', sentinelLimpo.size());
print('Cenas por ano (2016-2025):', colecaoAnual.aggregate_array('ano'), colecaoAnual.aggregate_array('num_cenas'));
// ================================================================
// 7. ANÁLISE DE TENDÊNCIA: MANN-KENDALL E SEN'S SLOPE
// ================================================================

var listaImagens = colecaoAnual.select('CO2_Index').toList(nAnos);
var sInit = ee.Image.constant(0).rename('S');

for (var i = 0; i < nAnos - 1; i++) {
  var imgI = ee.Image(listaImagens.get(i));
  for (var j = i + 1; j < nAnos; j++) {
    var imgJ = ee.Image(listaImagens.get(j));
    var diff = imgJ.subtract(imgI);
    var sinal = diff.where(diff.gt(0), 1).where(diff.lt(0), -1).where(diff.eq(0), 0);
    sInit = sInit.add(sinal);
  }
}
var S = sInit.clip(roi);

var n = ee.Number(nAnos);
var varianciaS = n.multiply(n.subtract(1)).multiply(n.multiply(2).add(5)).divide(18);
var desvioPadraoS = ee.Image.constant(varianciaS).sqrt();

var Z = ee.Image.constant(0)
  .where(S.gt(0), S.subtract(1).divide(desvioPadraoS))
  .where(S.lt(0), S.add(1).divide(desvioPadraoS))
  .rename('Z_score');

var listaSlopes = ee.List([]);
for (var i = 0; i < nAnos - 1; i++) {
  var imgI = ee.Image(listaImagens.get(i));
  var anoI = ee.Number(ee.Image(listaImagens.get(i)).get('ano'));
  for (var j = i + 1; j < nAnos; j++) {
    var imgJ = ee.Image(listaImagens.get(j));
    var anoJ = ee.Number(ee.Image(listaImagens.get(j)).get('ano'));
    var dt = anoJ.subtract(anoI);
    var declividadePar = imgJ.subtract(imgI).divide(dt);
    listaSlopes = listaSlopes.add(declividadePar);
  }
}

var sensSlope = ee.ImageCollection.fromImages(listaSlopes)
  .reduce(ee.Reducer.median())
  .rename('Sens_Slope')
  .clip(roi);

var tendenciaCategorizada = ee.Image.constant(0)
  .where(Z.gt(1.96), 2)
  .where(Z.gt(0).and(Z.lte(1.96)), 1)
  .where(Z.lt(-1.96), -2)
  .where(Z.lt(0).and(Z.gte(-1.96)), -1)
  .rename('Tendencia_CO2')
  .clip(roi);

// ================================================================
// 8. GRÁFICO DA SÉRIE TEMPORAL DOS VALORES MÉDIOS
// ================================================================

var graficoSerie = ui.Chart.feature.byFeature({
  features: tabelaAnual,
  xProperty: 'ano',
  yProperties: ['media']
}).setOptions({
  title: 'Série Temporal do CO2_Index (Valores Médios - 2016 a 2025)',
  hAxis: {title: 'Ano', format: '####', gridlines: {count: 10}},
  vAxis: {title: 'Média do CO2_Index'},
  lineWidth: 2,
  pointSize: 5,
  colors: ['#1a9850'],
  trendlines: {
    0: {color: 'red', lineWidth: 1, visibleInLegend: true, labelInLegend: 'Tendência Linear'}
  },
  legend: {position: 'bottom'}
});

print(graficoSerie);

// ================================================================
// 9. GRÁFICO DE HISTOGRAMA DA DISTRIBUIÇÃO ESPACIAL (MÉDIA HISTÓRICA 2016-2025)
// ================================================================

// Calcula a média de longo prazo de cada pixel para integrar toda a série temporal
var imgMediaSerie = colecaoAnual.select('CO2_Index').mean();

var histograma = ui.Chart.image.histogram({
  image: imgMediaSerie,
  region: roi,
  scale: 150, // Escala ampliada para 150m: evita estouro de memória mantendo o formato do gráfico
  maxPixels: 1e13
}).setOptions({
  title: 'Distribuição Espacial e Frequência do CO2_Index (Média Histórica 2016 - 2025)',
  hAxis: {title: 'Intervalos do CO2_Index (Média Temporal)', viewWindow: {min: -0.2, max: 0.8}},
  vAxis: {title: 'Número de pixels (Amostrados)'},
  colors: ['#006837'],
  legend: {position: 'none'}
});

print(histograma);


// =====================================================================
/// =====================================================================
// 10. EXPORTAÇÕES DIRETAS PARA O GOOGLE DRIVE (SEM MAPA INTERATIVO)
// =====================================================================

// 1. Exportar a tabela CSV de estatísticas anuais
Export.table.toDrive({
  collection: tabelaAnual,
  description: 'SerieTemporal_CO2_2016_2025',
  fileFormat: 'CSV',
  selectors: ['ano', 'media', 'desvio_padrao', 'minimo', 'maximo']
});

// 2. Exportar o Mapa de Tendência e Significância espacial (Resolução: 10m)
Export.image.toDrive({
  image: sensSlope.addBands(Z).addBands(tendenciaCategorizada).toFloat(),
  description: 'Tendencia_Significancia_CO2_2016_2025',
  scale: escala, 
  region: roi,
  maxPixels: 1e13
});

// =====================================================================
// EXPORTAR A IMAGEM DA MÉDIA HISTÓRICA (2016-2025) PARA O DRIVE
// =====================================================================

Export.image.toDrive({
  image: imgMediaSerie, // Utiliza a variável que já calcula a média no seu script
  description: 'Media_Historica_CO2_Index_2016_2025',
  scale: escala, // Mantém os 10 metros de resolução do Sentinel-2
  region: roi,
  maxPixels: 1e13
});


// ================================================================
// 11. QUANTIFICAÇÃO ESPACIAL DAS ÁREAS CRÍTICAS (HECTARES E %)
// ================================================================

// 1. Calcula a área de cada pixel individual em metros quadrados (m²)
var imagemAreaPixel = ee.Image.pixelArea();

// 2. Multiplica a imagem de tendências pela área do pixel para obter a métrica espacial
// Adicionamos bandas mascaradas para cada classe crítica de interesse
var imagemAreasClasses = ee.Image.cat([
  imagemAreaPixel.updateMask(tendenciaCategorizada.eq(2)).rename('Alta_Significativa'),
  imagemAreaPixel.updateMask(tendenciaCategorizada.eq(-2)).rename('Baixa_Significativa'),
  imagemAreaPixel.updateMask(tendenciaCategorizada.eq(1)).rename('Alta_Nao_Significativa'),
  imagemAreaPixel.updateMask(tendenciaCategorizada.eq(-1)).rename('Baixa_Nao_Significativa'),
  imagemAreaPixel.updateMask(tendenciaCategorizada.eq(0)).rename('Sem_Tendencia')
]);

// 3. Soma a área de todos os pixels dentro da sua Região de Interesse (ROI)
var somaAreasM2 = imagemAreasClasses.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: roi,
  scale: escala, // Mantém os 10 metros do Sentinel-2
  maxPixels: 1e13,
  tileScale: 16
});

// 4. Converte os valores de metros quadrados (m²) para Hectares (ha) dividindo por 10.000
var areaAltaSigHa = ee.Number(somaAreasM2.get('Alta_Significativa')).divide(10000);
var areaBaixaSigHa = ee.Number(somaAreasM2.get('Baixa_Significativa')).divide(10000);
var areaAltaNaoSigHa = ee.Number(somaAreasM2.get('Alta_Nao_Significativa')).divide(10000);
var areaBaixaNaoSigHa = ee.Number(somaAreasM2.get('Baixa_Nao_Significativa')).divide(10000);
var areaSemTendHa = ee.Number(somaAreasM2.get('Sem_Tendencia')).divide(10000);

// 5. Calcula a área total da ROI somando todas as classes anteriores
var areaTotalRoiHa = areaAltaSigHa
  .add(areaBaixaSigHa)
  .add(areaAltaNaoSigHa)
  .add(areaBaixaNaoSigHa)
  .add(areaSemTendHa);

// 6. Calcula as porcentagens de cada classe em relação ao total do território
var pctAltaSig = areaAltaSigHa.divide(areaTotalRoiHa).multiply(100);
var pctBaixaSig = areaBaixaSigHa.divide(areaTotalRoiHa).multiply(100);
var pctAltaNaoSig = areaAltaNaoSigHa.divide(areaTotalRoiHa).multiply(100);
var pctBaixaNaoSig = areaBaixaNaoSigHa.divide(areaTotalRoiHa).multiply(100);
var pctSemTend = areaSemTendHa.divide(areaTotalRoiHa).multiply(100);

// 7. Imprime os resultados formatados prontos para copiar para a redação do TCC
print('--- QUANTIFICAÇÃO ESPACIAL DAS TENDÊNCIAS (RESULTADOS TCC) ---');
print('Área Total Analisada (ha):', areaTotalRoiHa);
print('------------------------------------------------------------');
print('Classe [2] - Alta Significativa (ha):', areaAltaSigHa, 'Porcentagem (%):', pctAltaSig);
print('Classe [1] - Alta Não Significativa (ha):', areaAltaNaoSigHa, 'Porcentagem (%):', pctAltaNaoSig);
print('Classe [0] - Sem Tendência Clara (ha):', areaSemTendHa, 'Porcentagem (%):', pctSemTend);
print('Classe [-1] - Baixa Não Significativa (ha):', areaBaixaNaoSigHa, 'Porcentagem (%):', pctBaixaNaoSig);
print('Classe [-2] - Baixa Significativa (ha):', areaBaixaSigHa, 'Porcentagem (%):', pctBaixaSig);

// =====================================================================
// =====================================================================
// 11. CAMADAS DE VISUALIZAÇÃO NO MAPA (OTIMIZADAS)
// =====================================================================

// Centraliza o mapa na sua área de estudo (ROI)
Map.centerObject(roi, 12);

// Parâmetros de visualização
var visCO2 = {
  min: 0.1,
  max: 0.6,
  palette: ['#d73027', '#f46d43', '#fdae61', '#fee08b', '#d9ef8b', '#a6d96a', '#66bd63', '#1a9850']
};

var visTendencia = {
  min: -2,
  max: 2,
  palette: ['#d73027', '#f46d43', '#e0e0e0', '#a6d96a', '#1a9850']
};

// OTIMIZAÇÃO: Reprojetar nominalmente para a escala de exibição ou usar pyramiding
// Reduz a carga de memória ao carregar as camadas no mapa interativo
var imgMediaSerieExibicao = imgMediaSerie.reproject({
  crs: 'EPSG:4326',
  scale: 30 // Suaviza levemente a exibição para 30m na tela (evita o erro)
});

var tendenciaExibicao = tendenciaCategorizada.reproject({
  crs: 'EPSG:4326',
  scale: 30
});

// =====================================================================
// EXPORTAÇÃO COMPLEMENTAR DAS IMAGENS DE 2016 E 2025
// =====================================================================

// 1. Filtrar a imagem correspondente ao ano de 2016 na sua coleção anual
var img2016 = ee.Image(colecaoAnual.filter(ee.Filter.eq('ano', 2016)).first()).select('CO2_Index');

// Exportar a imagem de 2016 para o Google Drive
Export.image.toDrive({
  image: img2016,
  description: 'CO2_Index_Ano_2016',
  scale: escala,
  region: roi,
  maxPixels: 1e13
});

// 2. Filtrar a imagem correspondente ao ano de 2025 na sua coleção anual
var img2025 = ee.Image(colecaoAnual.filter(ee.Filter.eq('ano', 2025)).first()).select('CO2_Index');

// Exportar a imagem de 2025 para o Google Drive
Export.image.toDrive({
  image: img2025,
  description: 'CO2_Index_Ano_2025',
  scale: escala,
  region: roi,
  maxPixels: 1e13
});
