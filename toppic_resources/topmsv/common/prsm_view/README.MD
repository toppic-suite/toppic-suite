# PrsmGraph Library

PrsmGraph is a library for drawing Amino Acid Sequence table on html page.

## Dependencies

D3, jQuery
Topview library: proteoform

## Usage

Include all files from prsm_graph folder to the head of html file:
```html
<script src="your_dir/prsm_graph/draw_prsm.js"></script>
<script src="your_dir/prsm_graph/prsm_data.js"></script>
<script src="your_dir/prsm_graph/prsm_graph.js"></script>
<script src="your_dir/prsm_graph/prsm_para.js"></script>
```

Give svg_id to svg tag of html:
```html
<svg id = "your_svg_id" ></svg>
```

Javascript
```javascript
let prsmGraph = new PrsmGraph("svg_id", prsm);
prsmGraph.redraw();
```

## Example

Inside example directory, you can find a simple example with prsm data included.