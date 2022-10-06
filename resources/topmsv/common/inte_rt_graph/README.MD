# intensity-retentiontime graph library

inte_rt_graph is a library for drawing intensity-retentiontime svg graph on html page.

## Dependencies

D3

## Usage

Include library js file in your html page.
```html
<script src="your_dir/inte_rt_graph.js"></script>
```


InteRtGraph(svgID, data, rt_ID, inte_ID, height, width, padding)

Create a inte_rt_graph object, only svgID and data are required, other parameters will be default val if you don't set manually.
svgID is a string which is the id for svg element in html
data is an array which contains rt and inteSum, for example, [{rt:539.830854,inteSum:37038836.45996094},{rt:540.060222,inteSum:38193423.51464844}].
rt_ID is a string which is the id for retentiontime element in html, default value: 'rt-hover'
inte_ID is a string which is the id for intensity element in html, default value: 'intensity-hover'
height is a number which is svg graph height, default value: 120
width is a number which is svg graph width, default value: 1100
padding is a object which contains top, bottom, left, right property, default value: {top: 10, right: 10, bottom: 50, left: 80}


drawGraph()

Draw graph by given paramenters.


moveLine(rt)

Move the line on svg graph to given retention time location.
rt should be retentiontime in hour.


```javascript
let graphObj = new InteRtGraph("rt-sum",data);
graphObj.drawGraph();

graphObj.moveLine(rt);
```


## Example

There is a example.html in the same directory for reference.