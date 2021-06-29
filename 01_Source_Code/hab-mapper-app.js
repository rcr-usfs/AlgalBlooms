var app = {};
/////////////////////////////////////////////////////////
//Initial params
app.availableStates = ['OR','WA','WY'];
app.selectedState = 'WA';//Available states are OR, WA, and WY
app.defaultStartYear = 2017;
app.defaultEndYear = 2020;
app.defaultStartMonth = 7;
app.defaultEndMonth = 10;

app.availableStartYear = 2010;
app.availableEndYear = 2020;
app.availableStartMonth = 3;
app.availableEndMonth = 10;

app.minZThresh = 0.5;
app.maxZThresh = 4;
app.defaultZThresh = 1;

app.defaultPctHABThresh = 0.5;
app.minPctHABThresh = 0.1;
app.maxPctHABThresh = 25;

var summary_areas_dict = {};
summary_areas_dict['WA'] = ee.FeatureCollection('projects/gtac-algal-blooms/assets/ancillary/WA_FS_Named_Recreation_Lakes_v2');
summary_areas_dict['OR'] = ee.FeatureCollection('projects/gtac-algal-blooms/assets/ancillary/OR_FS_Named_Recreation_Lakes_v2');
summary_areas_dict['WY'] = ee.FeatureCollection('projects/gtac-algal-blooms/assets/ancillary/WY_FS_Named_Recreation_Lakes_v2');

//-------- Map Panel-----------------
// app.center = {lon:-118.1803, lat:36.8341,zoom:12};//Zoomed in center
// app.center = {lon:-119, lat:38.2486,zoom:6};//Zoomed to all CA
app.map = ui.Map();
app.map.style().set({cursor:'crosshair'});
app.map.setOptions('HYBRID');


app.blackline = function(){return ui.Label('__________________________________________')};
app.wideBlackline = function(){return ui.Label('____________________________________________________________________________')};
// -------------Legend Panel----------------------------------------
// Creates a color bar thumbnail image for use in legend from the given color palette.
function makeColorBarParams(palette, min, max) {
  return {
    bbox: [0, 0, 1, 0.1],
    dimensions: '150x10',
    format: 'png',
    min: min,
    max: max,
    palette: palette,
  };
}

// Create a panel with three numbers for the legend.
app.addToLegend = function(palette, min, middle, max, title){
  // Create the color bar for the legend.
  var colorBar = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: makeColorBarParams(palette),
    style: {stretch: 'horizontal', margin: '0px 0px', maxHeight: '20px', position: 'bottom-right'},
  });
  
  // Create a panel with three numbers for the legend.
  var legendLabels = ui.Panel({
    widgets: [
      ui.Label(min, {margin: '2px 2px',fontSize:'8pt'}),
      ui.Label(
          (middle),
          {margin: '2px 2px', textAlign: 'center', stretch: 'horizontal',fontSize:'8pt'}),
      ui.Label(max, {margin: '2px 2px',fontSize:'8pt'})
    ],
    layout: ui.Panel.Layout.flow('horizontal')
  });
  
  var legendTitle = ui.Label({
    value: title,
    style: {fontWeight: 'bold',margin:'1px 2px',fontSize:'8pt'}
  });
  
  var legendPanel = ui.Panel({widgets: [legendTitle, colorBar, legendLabels], style: {position: 'bottom-right',width:'250px',margin:'2px'}});
  // return legendPanel;
  app.legendPanel.add(legendPanel);
};
////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//Sync year range sliders
app.yearSliderSyncer = function(){
  var startYear = app.startYearSlider.getValue();
  var endYear = app.endYearSlider.getValue();
  
  if(startYear > endYear){app.endYearSlider.setValue(startYear,false)};
  if(endYear < startYear){app.startYearSlider.setValue(endYear,false)};
};
app.monthSliderSyncer = function(){
  var startMonth = app.startMonthSlider.getValue();
  var endMonth = app.endMonthSlider.getValue();
  // print(startMonth,endMonth)
  if(startMonth > endMonth){
    app.endMonthSlider.setValue(startMonth,false)};
  if(endMonth < startMonth){app.startMonthSlider.setValue(endMonth,false)};
};
app.updateState = function(){
  app.selectedState = app.stateSelect.getValue();
  app.summary_areas = summary_areas_dict[app.selectedState];
  app.map.centerObject(app.summary_areas.geometry().bounds());
  app.run();
  // app.map.centerObject(app.tables.sort('Pct_HAB',false).limit(1).geometry().bounds());
}

app.run = function(){
  
  app.map.clear();
  app.map.setOptions('HYBRID');
  app.legendPanel.clear();
  var startYear = app.startYearSlider.getValue();
  var endYear = app.endYearSlider.getValue();
  var startMonth = app.startMonthSlider.getValue();
  var endMonth = app.endMonthSlider.getValue();
  var state = app.selectedState;
  var nameEnd = state + ' yrs '+startYear.toString()+'-'+endYear.toString() + ' mths '+startMonth.toString() + '-' + endMonth.toString();
  
  var state = app.stateSelect.getValue();
  var years = ee.List.sequence(startYear,endYear).getInfo();
  var months = ee.List.sequence(startMonth,endMonth).getInfo();
  var zThresh = app.zThreshSlider.getValue();
  var pctHABThresh = app.pctHABThreshSlider.getValue();
  //Bring in tables and filter to study area and year
  app.tables = ee.data.listAssets('projects/gtac-algal-blooms/assets/outputs/HAB-Summary-Tables').assets;
  app.tables = app.tables.map(function(t){return t.id});
  app.tables = app.tables.filter(function(t){return  t.indexOf(state + '_') >-1});
  app.tables = app.tables.filter(function(t){return years.indexOf(parseInt(t.split('_yr')[t.split('_yr').length-1]))>-1});
  app.tables = app.tables.filter(function(t){return months.indexOf(parseInt(t.split('_m')[t.split('_m').length-1]))>-1});
  // print(app.tables)
  //Read in tables
  app.tables = app.tables.map(function(t){
    var f = ee.FeatureCollection(t);
    var m = t.split('_m')[t.split('_m').length-1];
    var y = t.split('_yr')[t.split('_yr').length-1].split('_')[0];
    f = f.map(function(f){return f.set({'year':parseInt(y),'month':parseInt(m)})});
    return f;
  });
  app.tables = ee.FeatureCollection(app.tables).flatten();

  // print(app.tables.size());
  // print(tables.filter(ee.Filter.eq('Pct_HAB',0)).size());
  // Map.addLayer(tables.filter(ee.Filter.eq('Pct_HAB',0)).reduceToImage(['Pct_HAB'],ee.Reducer.mean()),{min:0,max:0,palette:'00F'},'Zero')

  //Bring in z layers
  var z = ee.ImageCollection('projects/gtac-algal-blooms/assets/outputs/HAB-Z-Images')
  .filter(ee.Filter.calendarRange(startYear,endYear,'year'))
  .filter(ee.Filter.calendarRange(startMonth,endMonth,'month'))
  .filter(ee.Filter.eq('studyAreaName',state)).select([0],['HAB Z-Score'])
  .map(function(img){return img.divide(1000).copyProperties(img,['system:time_start'])});

  var summerWaterMask = z.filter(ee.Filter.calendarRange(7,9,'month')).map(function(img){return img.mask().selfMask()}).count().gt(2).selfMask();
  z = z.map(function(img){return img.updateMask(summerWaterMask)});
  
  app.zTimeSeries = z;

  var hab = z.map(function(img){return img.gte(zThresh)});
  // Map.addLayer(summerWaterMask,{'min':1,'max':1,'palette':'00F'},'summer water mask')
  // app.map.addLayer(z,{'min':0,'max':2000,'palette':'00F,0F0,F00'},'Z Time Series '+nameEnd,false);


  app.map.addLayer(z.max(),{'min':0,'max':2,'palette':'00F,0F0,F00'},'Max Z '+nameEnd,false);
  app.map.addLayer(hab.sum(),{'min':0,'max':months.length*years.length,'palette':'00F,0F0,F00'},'HAB Count '+nameEnd,false);

// // var groupByFieldOptions = ['REGION','FORESTNAME','DISTRICTNA','GNIS_Name'];
// // var groupByField = groupByFieldOptions[1];
  var pctHAB = app.tables.reduceToImage(['Pct_HAB'],ee.Reducer.max());
  var habSites = pctHAB.updateMask(pctHAB.gte(pctHABThresh));
  var cleanSites = pctHAB.updateMask(pctHAB.lt(pctHABThresh));
  app.map.addLayer(habSites,{min:pctHABThresh,max:25,palette:'FF0,F00'},'HAB Rec Sites '+nameEnd);
  app.map.addLayer(cleanSites,{min:0,max:pctHABThresh,palette:'00F,0FF'},'Clean Rec Sites '+nameEnd);
  
  app.addToLegend('00F,0FF', '0 % HAB', null, pctHABThresh.toString() + ' % HAB', 'Clean Rec Sites '+nameEnd);
  app.addToLegend('FF0,F00', pctHABThresh.toString() + ' % HAB', null,  '25 % HAB', 'HAB Rec Sites '+nameEnd);
  app.addToLegend('00F,0F0,F00', '0 HAB Obs', null, (months.length*years.length) + ' HAB Obs', 'HAB Count '+nameEnd);
  app.addToLegend('00F,0F0,F00', '0 (Z-Score)', null, '2 (Z-Score)', 'Max Z '+nameEnd);
  
  app.map.onClick(app.tableQuery)
};
/////////////////////////////////////////////////////////
app.tableQuery = function(event){
  var pt = ee.Geometry.Point([event.lon,event.lat]);
  try{app.map.remove(app.selectedPt);}
  catch(err){print(err)};
  app.selectedPt = ui.Map.Layer(pt,{},'Click Location');
  app.map.add(app.selectedPt);
  var sites = app.tables.filterBounds(pt);
  app.chartPanel.clear();
  // print(sites.size())
  if(sites.size().getInfo()>0){
    sites = sites.map(function(site){return site.set('date',ee.Date.fromYMD(site.get('year'),site.get('month'),1).millis())});
    sites = sites.sort('date');
    
    try{app.map.remove(app.selectedOutline);}
    catch(err){print(err)};
    app.selectedOutline = ui.Map.Layer(ee.Image().paint(sites.limit(1),null,1),{min:1,max:1,palette:'0DD'},'Selected Rec Site');
    app.map.add(app.selectedOutline);
    
    var siteName = sites.first().get('GNIS_Name').getInfo();
    
    var pcts = sites.toList(1000,0).map(function(f){
      f = ee.Feature(f);
      var date = ee.Date(f.get('date')).format();
      return f.get('Pct_HAB')});
    var dates = sites.toList(1000,0).map(function(f){
      f = ee.Feature(f);
      var date = ee.Date.fromYMD(f.get('year'),f.get('month'),1).format('YYYY-MM');
      return date});
    // print(pcts);
    // print(dates);
    var chart = ui.Chart.array.values(pcts,0,dates);
    chart.setOptions({
      lineWidth: 2,
      color:'080',
      title:siteName+' Rec Site Percent Mapped HAB',
      legend: {position: 'none'},
      height:325,
      chartArea: {'width': '85%', 'height': '60%'},
      hAxis: {title: "Date" , direction:1, slantedText:true, slantedTextAngle:45 },
      vAxis: {title: "% HAB" }
    });
    
  
    
    
    app.chartPanel.add(app.wideBlackline());
    app.chartPanel.add(chart);
    var zChart = ui.Chart.image.series(app.zTimeSeries, sites, ee.Reducer.mean(),30);
    zChart = zChart.setOptions({
    height:325,
    title: siteName+' Mean HAB Z-Score',
    'legend': {'position': 'bottom'},
    'chartArea': {'width': '85%', 'height': '70%'},
    vAxis : { 
        textStyle : {
            fontSize: 10 
        },title:'HAB Z-Score'
    },
  });
  app.chartPanel.add(app.wideBlackline());
  app.chartPanel.add(zChart);
  }else{
    app.chartPanel.clear();
    app.chartPanel.add(app.wideBlackline());
    app.chartPanel.add(ui.Label('No rec sites where you clicked. Please click on a rec site',{fontSize:'10pt'}));
  
    var zChart = ui.Chart.image.series(app.zTimeSeries, pt, ee.Reducer.first(),30);
    zChart = zChart.setOptions({
    height:325,
    title: 'HAB Z-Score Click Point',
    'legend': {'position': 'bottom'},
    'chartArea': {'width': '85%', 'height': '70%'},
    vAxis : { 
        textStyle : {
            fontSize: 10 
        },title:'HAB Z-Score'
    },
  });
  app.chartPanel.add(app.wideBlackline());
  app.chartPanel.add(zChart);
  }
  
  
}
/////////////////////////////////////////////////////////
//----- Options & Filters Panel---------------------------
// Left side of window, use to select filters and view options
app.optionsFiltersPanel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {width: '300px', position: 'top-left'}
});


app.optionsFiltersPanel.add(ui.Label('HAB Mapper 5k',{'fontSize':'20pt','fontWeight':'bold','color':'008800'}));
app.optionsFiltersPanel.add(app.blackline());
app.optionsFiltersPanel.add(ui.Label('What is this app?',{'fontWeight':'bold'}));
app.optionsFiltersPanel.add(ui.Label('This is a beta tool intended to depict potential presence of harmful algal blooms.'));
app.optionsFiltersPanel.add(app.blackline());
app.stateSelect = ui.Select(app.availableStates, 'Select a state: ',app.selectedState, app.updateState); 
app.optionsFiltersPanel.add(app.stateSelect);
app.optionsFiltersPanel.add(app.blackline());
app.optionsFiltersPanel.add(ui.Label('Choose years to monitor:'));
app.startYearSlider =ui.Slider(app.availableStartYear, app.availableEndYear, app.defaultStartYear,1,null,'horizontal',null,{width: '250px'});
app.endYearSlider =ui.Slider(app.availableStartYear, app.availableEndYear, app.defaultEndYear,1,null,'horizontal',null,{width: '250px'});
app.optionsFiltersPanel.add(app.startYearSlider);
app.optionsFiltersPanel.add(app.endYearSlider);
app.startYearSlider.onSlide(app.yearSliderSyncer);
app.endYearSlider.onSlide(app.yearSliderSyncer);

app.optionsFiltersPanel.add(app.blackline());
app.optionsFiltersPanel.add(ui.Label('Choose months to monitor:'));
app.startMonthSlider =ui.Slider(app.availableStartMonth, app.availableEndMonth, app.defaultStartMonth,1,null,'horizontal',null,{width: '250px'});
app.endMonthSlider =ui.Slider(app.availableStartMonth, app.availableEndMonth, app.defaultEndMonth,1,null,'horizontal',null,{width: '250px'});
app.optionsFiltersPanel.add(app.startMonthSlider);
app.optionsFiltersPanel.add(app.endMonthSlider);
app.startMonthSlider.onSlide(app.monthSliderSyncer);
app.endMonthSlider.onSlide(app.monthSliderSyncer);

app.optionsFiltersPanel.add(app.blackline());
app.optionsFiltersPanel.add(ui.Label('Choose Z-Score threshold for spatially explicit HAB locations:'));
app.zThreshSlider =ui.Slider(app.minZThresh, app.maxZThresh, app.defaultZThresh,0.5,null,'horizontal',null,{width: '250px'});
app.optionsFiltersPanel.add(app.zThreshSlider);
app.optionsFiltersPanel.add(app.blackline());
app.optionsFiltersPanel.add(ui.Label('Choose % of rec area that had to have HAB to be shown:'));
app.pctHABThreshSlider =ui.Slider(app.minPctHABThresh, app.maxPctHABThresh, app.defaultPctHABThresh,0.1,null,'horizontal',null,{width: '250px'});
app.optionsFiltersPanel.add(app.pctHABThreshSlider);
app.optionsFiltersPanel.add(app.blackline());
app.runButton = ui.Button('Update Map', app.run)
app.optionsFiltersPanel.add(app.runButton);
//////////////////////////////////////////////////////////////////////////////////
//----- Plot Panel---------------------------
// Right side of window, shows charts for clicked polygons
app.plotPanel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {width: '400'}
});
app.legendPanel = ui.Panel(null, null, {position: 'bottom-right',margin:'5px'});
app.chartPanel = ui.Panel(null, null, {stretch: 'horizontal'});
app.plotPanel.add(ui.Label('Click on map to chart HAB time series'));

app.plotPanel.add(app.chartPanel);
app.plotPanel.add(app.wideBlackline());
app.plotPanel.add(ui.Label('Map Key',{fontWeight:'bold',fontSize:'12pt'}));
app.plotPanel.add(app.legendPanel);
app.plotPanel.add(app.wideBlackline());
/////////////////////////////////////////////////////////////
//Put app together
ui.root.clear();
ui.root.add(app.optionsFiltersPanel);

ui.root.add(app.map);
ui.root.add(app.plotPanel);
// ui.root.add(app.plotPanel);
app.updateState();

