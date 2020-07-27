# monomassGraph

monomassGraph is a library for drawing mono mass graph on html page.

## Dependencies

D3, jQuery, spectrum graph library in topview

## Usage

Include javascript files
```html
<script src="your_dir/monomassGraph/multiscan.js"></script>
<script src="your_dir/monomassGraph/calculateprefixandsuffixmass.js"></script>
<script src="your_dir/monomassGraph/generateSpectrumGraphs.js"></script>
<script src="your_dir/monomassGraph/prsmDataUtil.js"></script>
<script src="your_dir/monomassGraph/drawMonoMassSpectrum.js"></script>
<script src="your_dir/monomassGraph/aminoaciddistribution.js"></script>
<script src="your_dir/monomassGraph/monomassGraph.js"></script>
```

Build monomass_div and monoMass_nav inside
```html
<div id="monomass_div">
    <ul class="nav nav-tabs graph_nav" id="monoMass_nav"></ul>
</div>
```

Call monomassGraph with your prsm data (there is an example in test folder called prsm133.js)
```javascript
monomassGraph(your_prsm_data);
```

## Example

Please open testMonomass.html in test folder to check how it works.