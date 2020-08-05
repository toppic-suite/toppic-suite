# monomassGraph

monomassGraph is a library for drawing mono mass graph on html page.

## Dependencies

D3, jQuery
Topview library: spectrum_graph, calculate_prefix_and_suffix_mass, amino_acid_distribution

## Usage

Include javascript files
```html
<script src="your_dir/monomass_graph.js"></script>
```

Build monomass_div and monoMass_nav in html page, and give div_id and nav_id to monomass_graph object.
```html
<div id="monomass_div">
    <ul class="nav nav-tabs graph_nav" id="monoMass_nav"></ul>
</div>
```

Call monomassGraph with your prsm data(required), nav_id(optional, default value: "monoMass_nav"), div_id(optional, default value: "monomass_div")

```javascript
let monomassObj = new MonomassGraph(prsm_data, "monoMass_nav", "monomass_div");
monomassObj.drawGraph();
```

## Example

Please open testMonomass.html in example folder to check how it works.