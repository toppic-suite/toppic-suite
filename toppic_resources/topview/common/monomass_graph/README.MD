# monomassGraph

monomassGraph is a library for drawing mono mass graph on html page.

## Dependencies

D3, jQuery, spectrum_graph library in topview

## Usage

Include javascript files
```html
<script src="your_dir/multiscan.js"></script>
<script src="your_dir/calculate_prefix_and_suffix_mass.js"></script>
<script src="your_dir/generate_spectrum_graphs.js"></script>
<script src="your_dir/prsm_data_util.js"></script>
<script src="your_dir/draw_monomass_spectrum.js"></script>
<script src="your_dir/amino_acid_distribution.js"></script>
<script src="your_dir/monomass_graph.js"></script>
<script src="your_dir/peakData.js"></script>
```

Build monomass_div and monoMass_nav inside
```html
<div id="monomass_div">
    <ul class="nav nav-tabs graph_nav" id="monoMass_nav"></ul>
</div>
```

Call monomassGraph with your prsm data (there is an example in example folder called prsm133.js)
```javascript
monomassGraph(your_prsm_data);
```

## Example

Please open testMonomass.html in example folder to check how it works.