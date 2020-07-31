# svgAAS

svgAAS is a library for drawing Amino Acid Sequence table on svg.

## Dependencies

D3 and jQuery
Topview library: prsm_data_util

## Usage

Include all files from svgAAS folder to head of html file:
```html
<script src="your_dir/svgAAS/coordinates_util.js"></script>
<script src="your_dir/svgAAS/draw_svg.js"></script>
<script src="your_dir/svgAAS/svgAAS.js"></script>
```

Give svg_id to svg tag of html:
```html
<svg id = "your_svg_id" ></svg>
```

Call function svgAAS(prsm_data, svg_id) to draw svg with given svg_id:
```javascript
svgAAS(your_prsm_data, "your_svg_id");
```

## Example

Inside example directory, you can find a simple example with prsm data included.