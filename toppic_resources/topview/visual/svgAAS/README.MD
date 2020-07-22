# svgAAS

svgAAS is a library for drawing Amino Acid Sequence table on svg.

## Dependencies

D3 and jQuery

## Usage

Include all files from svgAAS folder to head of html file:
```
<script src="./svgAAS/coordinates_util.js"></script>
<script src="./svgAAS/drawSvg.js"></script>
<script src="./svgAAS/prsmDataUtil.js"></script>
<script src="./svgAAS/SvgAAS.js"></script>
```

Give svg_id to svg tag of html:
```
<svg id = "your_svg_id" ></svg>
```

Call function svgAAS(prsm_data, svg_id) to draw svg with given svg_id:
```
svgAAS(your_prsm_data, "your_svg_id");
```