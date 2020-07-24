# svgAAS

svgAAS is a library for drawing Amino Acid Sequence table on svg.

## Dependencies

D3 and jQuery

## Usage

Include all files from svgAAS folder to head of html file:
```html
<script src="your_dir/svgAAS/coordinates_util.js"></script>
<script src="your_dir/svgAAS/drawSvg.js"></script>
<script src="your_dir/svgAAS/prsmDataUtil.js"></script>
<script src="your_dir/svgAAS/SvgAAS.js"></script>
```

Give svg_id to svg tag of html:
```html
<svg id = "your_svg_id" ></svg>
```

Call function svgAAS(prsm_data, svg_id) to draw svg with given svg_id:
```javascript
svgAAS(your_prsm_data, "your_svg_id");
```