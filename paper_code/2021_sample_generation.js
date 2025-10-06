var USA = ee.FeatureCollection("projects/ee-zjl/assets/USA_quhua/USA"),
    S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED"),
    studyArea = ee.FeatureCollection("projects/ee-zjl/assets/USA_quhua/studyArea");
Map.centerObject(studyArea,7)
Map.addLayer(studyArea,{},'studyArea')

/********************************************************************************
step1:generate othercrop map and CDL_maize_soybean in Indiana
*********************************************************************************/
var ESA10v200 = ee.ImageCollection("ESA/WorldCover/v100").first()
var Cropland = ee.Image.constant(1)
  .updateMask(ESA10v200.eq(40)).clip(studyArea)


var CDL_2021=ee.ImageCollection("USDA/NASS/CDL").filterDate('2021-01-01','2021-12-31').first().select('cropland').clip(studyArea)


//maize
var maize_2021=CDL_2021.eq(1).selfMask()
Map.addLayer(maize_2021,{min:1,max:1,palette:['ffd300']},'maize_2021')
//soybean
var soybean_2021=CDL_2021.eq(5).selfMask()
Map.addLayer(soybean_2021,{min:1,max:1,palette:['267000']},'soybean_2021')
//othercrop
var maize_soybean_2021=ee.ImageCollection.fromImages([maize_2021,soybean_2021]).mosaic().rename('label_1').unmask(0).updateMask(Cropland)
var othercrop_2021=maize_soybean_2021.eq(0).selfMask()
Map.addLayer(othercrop_2021,{min:1,max:1,palette:['00a8e2']},'othercrop_2021')
/**************************************************************************
step2:find optimal focal radii
1、generate 10000 random points in StudyArea for each crop
2、select points which locate in the area of each crop(crop points)
3、caculate the mean and var of NDVI in July for each crop points in 16m 32m 48m 64m 80m 96m......
4、choose the stabilize use the function findoptimalradii 


result:maize:48m,soybean:96m,othercrop:0m
***************************************************************************/

      var maizesample_r=ee.FeatureCollection.randomPoints({
        region:studyArea,
        points:10000,
        seed:1,
        maxError:1
      })
      var soybeansample_r=ee.FeatureCollection.randomPoints({
        region:studyArea,
        points:10000,
        seed:5,
        maxError:1
      })
      var othercropsample_r=ee.FeatureCollection.randomPoints({
        region:studyArea,
        points:10000,
        seed:0,
        maxError:1
      })
      
      //define focal radii model
      
      function neighFun(img,kernalRadius,roi){
          var kernel = ee.Kernel.square(kernalRadius,'pixels',false);
          var kernelArea = (ee.Number(kernalRadius).multiply(2).add(1)).pow(2);
          var imgNeibor = ee.Image(img).convolve(kernel)
                            .eq(kernelArea)
                            .set("system:footprint",roi.geometry());
          return img.updateMask(imgNeibor);
      }
      
      
      
      //define July NDVI mean
      
      var NDVI_July=S2.filterDate('2021-07-01','2021-08-01')
                      .filterBounds(studyArea)
                      .map(function(image){
                        var NDVI=image.normalizedDifference(['B8','B4']).rename('NDVI')
                        return image.addBands(NDVI)
                      })
                      .select(['NDVI'])
                      .mean()
                      .clip(studyArea)
      Map.addLayer(NDVI_July,{},'NDVI_July')
      
      //maize_process
      
      var maize_feature_list=[]
      var maize_mean_list=[]
      
      for(var i=0;i<=10;i++){
        
        var maize_mask=neighFun(maize_2021,i,studyArea)
        var inward_maize=maize_2021.updateMask(maize_mask).int16().selfMask()
        
        var maizesample=inward_maize.reduceRegions({
          collection:maizesample_r,
          reducer:ee.Reducer.mean(),
          scale:16
        })
        
        var maizesample=maizesample.filter(ee.Filter.notNull(['mean']))
        
        var maizesample_statistics=NDVI_July.reduceRegions({
          collection:maizesample,
          reducer:ee.Reducer.mean(),
          scale:10
        })
        
        var reducer = ee.Reducer.mean().combine({
          reducer2: ee.Reducer.stdDev(),
          sharedInputs: true
        })
        
        var maizesample_statistics_dic=maizesample_statistics.reduceColumns({
          reducer:reducer,
          selectors:['mean']
        })
        var mean=maizesample_statistics_dic.get('mean')
        var stdDev=maizesample_statistics_dic.get('stdDev')
        
        var feature=ee.Feature(null,{
          'scale':i,
          'NDVI_mean':maizesample_statistics_dic.get('mean'),
          'NDVI_stdDev':maizesample_statistics_dic.get('stdDev')
        })
        
        maize_feature_list.push(feature)
        maize_mean_list.push(mean)
      }
      var maize_mean_list=ee.List(maize_mean_list)
     
      
      
      var feature_list=ee.List(maize_feature_list)
      var feature=ee.FeatureCollection(feature_list)
      var Maize_chart=ui.Chart.feature.byFeature(feature,'scale',['NDVI_mean','NDVI_stdDev']).setOptions({
        title:'Maize optimal focal radii',
        hAxis: {title: 'Scale'},
        vAxis:{title: 'Index Value'},
        lineWidth: 2,
        pointSize: 4,
        series: {
          0: {color: 'red'},
          1: {color: 'green'},
        }
      })
      print(Maize_chart)
      
      
      
      
      
      //soybean_process

      var soybean_feature_list=[]
      var soybean_mean_list=[]
      for(var i=0;i<=10;i++){
        
        var soybean_mask=neighFun(soybean_2021,i,studyArea)
        var inward_soybean=soybean_2021.updateMask(soybean_mask).int16().selfMask()
        
        var soybeansample=inward_soybean.reduceRegions({
          collection:soybeansample_r,
          reducer:ee.Reducer.mean(),
          scale:16
        })
        
        var soybeansample=soybeansample.filter(ee.Filter.notNull(['mean']))
        
        var soybeansample_statistics=NDVI_July.reduceRegions({
          collection:soybeansample,
          reducer:ee.Reducer.mean(),
          scale:10
        })
        
        var reducer = ee.Reducer.mean().combine({
          reducer2: ee.Reducer.stdDev(),
          sharedInputs: true
        })
        
        var soybeansample_statistics_dic=soybeansample_statistics.reduceColumns({
          reducer:reducer,
          selectors:['mean']
        })
        var mean=soybeansample_statistics_dic.get('mean')
        var stdDev=soybeansample_statistics_dic.get('stdDev')
        
        var feature=ee.Feature(null,{
          'scale':i,
          'NDVI_mean':soybeansample_statistics_dic.get('mean'),
          'NDVI_stdDev':soybeansample_statistics_dic.get('stdDev')
        })
        
        soybean_feature_list.push(feature)
        soybean_mean_list.push(mean)
      }
      var soybean_mean_list=ee.List(soybean_mean_list)
      
      
      
      var feature_list=ee.List(soybean_feature_list)
      var feature=ee.FeatureCollection(feature_list)
      var Soybean_chart=ui.Chart.feature.byFeature(feature,'scale',['NDVI_mean','NDVI_stdDev']).setOptions({
        title:'Soybean optimal focal radii',
        hAxis: {title: 'Scale'},
        vAxis:{title: 'Index Value'},
        lineWidth: 2,
        pointSize: 4,
        series: {
          0: {color: 'red'},
          1: {color: 'green'},
        }
      })
      print(Soybean_chart)
      
      
      
      
/********************************************************************************
step3:use focal radii to maize,soybean and othercrop
*********************************************************************************/

var maize_mask=neighFun(maize_2021,1,studyArea)
var inward_maize=maize_2021.updateMask(maize_mask).int16().selfMask().rename('Label')

var soybean_mask=neighFun(soybean_2021,1,studyArea)
var inward_soybean=soybean_2021.updateMask(soybean_mask).int16().selfMask().multiply(2).rename('Label')

var othercrop_mask=neighFun(othercrop_2021,0,studyArea)
var inward_othercrop=othercrop_2021.updateMask(othercrop_mask).int16().selfMask().multiply(0).rename('Label')

Map.addLayer(inward_maize,{min:1,max:1,palette:['FF0000']},'inward_maize')
Map.addLayer(inward_soybean,{min:1,max:1,palette:['22ff00']},'inward_soybean')
Map.addLayer(inward_othercrop,{min:1,max:1,palette:['1500ff']},'inward_othercrop')


var numberofsample=5000


var maizesample=inward_maize.stratifiedSample({
  numPoints:numberofsample,
  classBand:'Label',
  region:studyArea,
  scale:10,
  tileScale:16,
  geometries:true
})

var soybeansample=inward_soybean.stratifiedSample({
  numPoints:numberofsample,
  classBand:'Label',
  region:studyArea,
  scale:10,
  tileScale:16,
  geometries:true
})

var othercropsample=inward_othercrop.stratifiedSample({
  numPoints:numberofsample,
  classBand:'Label',
  region:studyArea,
  scale:10,
  tileScale:16,
  geometries:true
})





/********************************************************************************
step4:sample purify by NDVI_mean±NDVI_std
*********************************************************************************/


      
      var interpolFunction = require('users/cxh1216/GEEcourse:function/Linear_Interpol');
      var s2tbx = require('users/BNU_LXY/pkgs:s2_toolbox_ukr')
      var year = 2021;
      var startDate = ee.Date('2021-05-01'),endDate = ee.Date('2021-10-01'),SY = ee.Date('2021-01-01');
      var interval = 30;
      var band = ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'];
      
      var S2_Filtered = S2.filterDate(startDate,endDate).filterBounds(studyArea);
      
      var S2_noCloud = S2_Filtered.map(function(image){
        var qa = image.select('QA60'); 
        var cloudBitMask = 1 << 10;
        var cirrusBitMask = 1 << 11;
        var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(qa.bitwiseAnd(cirrusBitMask).eq(0));
        return image.updateMask(mask);
      }).select(band);
      
      var S2_Composition = interpolFunction.TemporalComposite(S2_noCloud,startDate,endDate,30)
      
      var S2_Interpol = interpolFunction.LinearInterpol(S2_Composition);
      
          S2_Interpol = S2_Interpol.map(function(image){
          var NDVI = image.normalizedDifference(['B8','B4']).rename('NDVI');
          var LSWI=image.normalizedDifference(['B8','B11']).rename('LSWI');
          var NDSVI=image.normalizedDifference(['B11','B4']).rename('NDSVI');
          var NDTI=image.normalizedDifference(['B11','B12']).rename('NDTI')
          var RENDVI=image.normalizedDifference(['B8','B6']).rename('RENDVI')
        return image.addBands(NDVI).addBands(LSWI).addBands(NDSVI).addBands(NDTI).addBands(RENDVI)
      }).map(s2tbx.addEVI)
      
      
      var featureImg=S2_Interpol.select(['NDVI']).toBands().clip(studyArea)
      var featureBand=featureImg.bandNames()
      
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
      
      var ESA10v200 = ee.ImageCollection("ESA/WorldCover/v200").first()
      var Cropland = ee.Image.constant(1)
        .clip(studyArea)
        .updateMask(ESA10v200.eq(40))


  
var featureImg=featureImg.updateMask(Cropland)
var featureBand=featureImg.bandNames()

var sample_list=ee.List([maizesample,soybeansample,othercropsample])




var purify_sample_list=sample_list.map(function(sample){
  var sample2=featureImg.sampleRegions({
    collection:sample,
    properties: ['Label'],
    scale:10,
    tileScale:16,
    geometries:true,
  })
  
  var Mean=sample2.reduceColumns({
    reducer:ee.Reducer.mean().repeat(6),
    selectors:featureBand
  })
  var Mean=ee.List(ee.Dictionary(Mean).get('mean'))
  var StdDev=sample2.reduceColumns({
    reducer:ee.Reducer.stdDev().repeat(6),
    selectors:featureBand
  })
  var StdDev=ee.List(ee.Dictionary(StdDev).get('stdDev'))
  var MeanStdDevPlus=ee.List.sequence(0,5).map(function(item){
    var MeanStdDevPlus = ee.Number(Mean.get(item)).add(StdDev.get(item));
    return MeanStdDevPlus;
  })
  
  var MeanStdDevMinus = ee.List.sequence(0,5).map(function(item){
    var MeanStdDevMinus = ee.Number(Mean.get(item)).subtract(StdDev.get(item));
    return MeanStdDevMinus;
  });
  
  var sample3=sample2.map(function(fea){
    var feaValue = ee.List(ee.Feature(fea).toDictionary(featureBand).values());
    var flag = ee.List.sequence(0,5).map(function(item){
      var a_b = ee.Number(MeanStdDevMinus.get(item)).subtract(feaValue.get(item));
      var c_b = ee.Number(MeanStdDevPlus.get(item)).subtract(feaValue.get(item));
      return ee.Number(a_b).multiply(c_b);
    });
  
    return fea.set('flag0',flag.get(0)).set('flag1',flag.get(1)).set('flag2',flag.get(2))
              .set('flag3',flag.get(3)).set('flag4',flag.get(4)).set('flag5',flag.get(5))
            
  })
  
  var purify_sample=sample3.filter(ee.Filter.lt("flag0",0))
                            .filter(ee.Filter.lt("flag1",0))
                            .filter(ee.Filter.lt("flag2",0))
                            .filter(ee.Filter.lt("flag3",0))
                            .filter(ee.Filter.lt("flag4",0))
                            .filter(ee.Filter.lt("flag5",0))

  return purify_sample
})

var purify_maize_sample=purify_sample_list.get(0)
var purify_soybean_sample=purify_sample_list.get(1)
var purify_othercrop_sample=purify_sample_list.get(2)
/********************************************************************************
step5:extract feature to sample
*********************************************************************************/
        var interpolFunction = require('users/cxh1216/GEEcourse:function/Linear_Interpol');
        var s2tbx = require('users/BNU_LXY/pkgs:s2_toolbox_ukr')
        var year = 2021;
        var startDate = ee.Date('2021-05-01'),endDate = ee.Date('2021-10-01'),SY = ee.Date('2021-01-01');
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
        
        



        var featureImg=featureImg.updateMask(Cropland)
        var featureBand=featureImg.bandNames()
        
      var maizesample_purify_2021_withfeature=featureImg.sampleRegions({
        collection:purify_maize_sample,
        properties:['Label'],
        scale: 10,
        tileScale:16,
        geometries:true
      })
      var soybeansample_purify_2021_withfeature=featureImg.sampleRegions({
        collection:purify_soybean_sample,
        properties:['Label'],
        scale: 10,
        tileScale:16,
        geometries:true
      })
      var othercropsample_purify_2021_withfeature=featureImg.sampleRegions({
        collection:purify_othercrop_sample,
        properties:['Label'],
        scale: 10,
        tileScale:16,
        geometries:true
      })        