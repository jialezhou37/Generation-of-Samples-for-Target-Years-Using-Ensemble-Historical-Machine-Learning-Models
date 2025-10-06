var studyArea = ee.FeatureCollection("projects/ee-zjl/assets/USA_quhua/studyArea"),
    S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED");


Map.centerObject(studyArea,7)
Map.addLayer(studyArea,{},'studyArea')

/********************************************************************************
step1:generate 2023 othercrop map and CDL_maize_soybean in studyArea
*********************************************************************************/

var ESA10v200 = ee.ImageCollection("ESA/WorldCover/v100").first()
var Cropland = ee.Image.constant(1)
  .updateMask(ESA10v200.eq(40)).clip(studyArea)


var CDL_2023=ee.ImageCollection("USDA/NASS/CDL").filterDate('2023-01-01','2023-12-31').first().select('cropland').clip(studyArea)


//maize
var maize_2023=CDL_2023.eq(1).selfMask()
Map.addLayer(maize_2023,{min:1,max:1,palette:['ffd300']},'maize_2023')
//soybean
var soybean_2023=CDL_2023.eq(5).selfMask()
Map.addLayer(soybean_2023,{min:1,max:1,palette:['267000']},'soybean_2023')
//othercrop
var maize_soybean_2023=ee.ImageCollection.fromImages([maize_2023,soybean_2023]).mosaic().rename('label_1').unmask(0).updateMask(Cropland)
var othercrop_2023=maize_soybean_2023.eq(0).selfMask()
Map.addLayer(othercrop_2023,{min:1,max:1,palette:['00a8e2']},'othercrop_2023')


/********************************************************************************
step3:prepare 2023 featureImg to classiication
*********************************************************************************/

        var interpolFunction = require('users/cxh1216/GEEcourse:function/Linear_Interpol');
        var s2tbx = require('users/BNU_LXY/pkgs:s2_toolbox_ukr')
        var year = 2023;
        var startDate = ee.Date('2023-05-01'),endDate = ee.Date('2023-10-01'),SY = ee.Date('2023-01-01');
        var interval = 10;
        var band = ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'];
        
        var S2_Filtered = S2.filterDate(startDate,endDate).filterBounds(studyArea);
        print(S2_Filtered,'S2_Filtered')
        var S2_noCloud = S2_Filtered.map(function(image){
          var qa = image.select('QA60'); 
          var cloudBitMask = 1 << 10;
          var cirrusBitMask = 1 << 11;
          var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(qa.bitwiseAnd(cirrusBitMask).eq(0));
          return image.updateMask(mask);
        }).select(band);
        
        var S2_Composition = interpolFunction.TemporalComposite(S2_noCloud,startDate,endDate,10)
        
        var S2_Interpol = interpolFunction.LinearInterpol(S2_Composition);
        
            S2_Interpol = S2_Interpol.map(function(image){
            var NDVI = image.normalizedDifference(['B8','B4']).rename('NDVI');
            var LSWI=image.normalizedDifference(['B8','B11']).rename('LSWI');
            var NDSVI=image.normalizedDifference(['B11','B4']).rename('NDSVI');
            var NDTI=image.normalizedDifference(['B11','B12']).rename('NDTI')
            var RENDVI=image.normalizedDifference(['B8','B6']).rename('RENDVI')
            var B = image.float().select('B2').divide(10000)
            var R = image.float().select('B4').divide(10000)
            var N = image.float().select('B8').divide(10000)
            var EVI = image.expression(
              '2.5*(N-R)/(N+6*R-7.5*B+1)',
              {
                'B':B, 'R':R, 'N':N
              }
            ).rename('EVI')
            var RE1 = image.float().select('B5').divide(10000)
            var RE2 = image.float().select('B6').divide(10000)
            var RE3 = image.float().select('B7').divide(10000)
            var REP = image.expression(
              '(705+35*(0.5*(RE3+R)-RE1)/(RE2-RE1))/1000',
              {
                'R':R, 'RE1':RE1, 'RE2':RE2, 'RE3':RE3
              }
            ).rename('REP')
            var B61=image.float().select('B6').divide(10000).rename('B61')
            var B111=image.float().select('B11').divide(10000).rename('B111')
            var B121=image.float().select('B12').divide(10000).rename('B121')
          return image.addBands(B61).addBands(B111).addBands(B121).addBands(NDVI).addBands(LSWI).addBands(NDSVI).addBands(NDTI).addBands(RENDVI).addBands(EVI).addBands(REP)
        })
        
        
        var featureImg=S2_Interpol.select(['B61','B111','B121','NDVI','LSWI','NDSVI','NDTI','RENDVI','EVI','REP']).toBands().clip(studyArea)
        var featureBand=featureImg.bandNames()
        print(featureBand,'featureBand')
        var featureImg=featureImg.select(featureBand,featureBand.map(function(str){
          var str=ee.String(str)
          var str_1='_'
          var part=str.split(str_1)
          var slicestr=part.get(1)
          var Index=part.get(0)
          var slicestr=ee.String(slicestr)
          var finalstr=slicestr.cat('_').cat(Index)
          return finalstr
        }))
        
        



        var featureImg_2023=featureImg.updateMask(Cropland)
        var featureBand=featureImg_2023.bandNames()


 



/********************************************************************************
step4:classifer2
*********************************************************************************/

var maizesample=featureImg_2023.sampleRegions({
    collection:maizesample,
    properties:['Label'],
    scale:10,
    tileScale:16,
    geometries:true,
})
var maizesample_expand=featureImg_2023.sampleRegions({
    collection:maizesample_expand,
    properties:['Label'],
    scale:10,
    tileScale:16,
    geometries:true,
})
var soybeansample=featureImg_2023.sampleRegions({
    collection:soybeansample,
    properties:['Label'],
    scale:10,
    tileScale:16,
    geometries:true,
}).limit(500)
var trainingSample_2023=maizesample.merge(soybeansample).merge(othercropsample).merge(maizesample_expand)
var trainingFeature=featureBand
        
        
var classifier_2023 = ee.Classifier.smileRandomForest(100)
                                .train({
                                  features:trainingSample_2023,
                                  classProperty:'Label',
                                  inputProperties:trainingFeature
                                })
                                
                                
var classification_2023=featureImg_2023.classify(classifier_2023,'Label').rename('label')

Map.addLayer(classification_2023,{min:0,max:2,palette:['00a8e2','ffd300','267000']},'classification_2023_2')