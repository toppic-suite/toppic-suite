function getSVGString( svgNode ) {
	svgNode.setAttribute('xlink', 'http://www.w3.org/1999/xlink');
	var cssStyleText = getCSSStyles( svgNode );
	appendCSS( cssStyleText, svgNode );
	var serializer = new XMLSerializer();
	var svgString = serializer.serializeToString(svgNode);
	svgString = svgString.replace(/(\w+)?:?xlink=/g, 'xmlns:xlink='); // Fix root xlink without namespace
	svgString = svgString.replace(/NS\d+:href/g, 'xlink:href'); // Safari NS namespace fix
	return svgString;
	function getCSSStyles( parentElement ) {
		var selectorTextArr = [];
		// Add Parent element Id and Classes to the list
		selectorTextArr.push( '#'+parentElement.id );
		for (var c = 0; c < parentElement.classList.length; c++)
				if ( !contains('.'+parentElement.classList[c], selectorTextArr) )
					selectorTextArr.push( '.'+parentElement.classList[c] );
		// Add Children element Ids and Classes to the list
		var nodes = parentElement.getElementsByTagName("*");
		for (var i = 0; i < nodes.length; i++) {
			var id = nodes[i].id;
			if ( !contains('#'+id, selectorTextArr) )
				selectorTextArr.push( '#'+id );
			var classes = nodes[i].classList;
			for (var c = 0; c < classes.length; c++)
				if ( !contains('.'+classes[c], selectorTextArr) )
					selectorTextArr.push( '.'+classes[c] );
		}
		// Extract CSS Rules
		var extractedCSSText = "";
		for (var i = 0; i < document.styleSheets.length; i++) {
			var s = document.styleSheets[i];
			
			try {
			    if(!s.cssRules) continue;
			} catch( e ) {
		    		if(e.name !== 'SecurityError') throw e; // for Firefox
		    		continue;
		    	}
			var cssRules = s.cssRules;
			for (var r = 0; r < cssRules.length; r++) {
				if ( contains( cssRules[r].selectorText, selectorTextArr ) )
					extractedCSSText += cssRules[r].cssText;
			}
		}
		
		return extractedCSSText;
		
		function contains(str,arr) {
			return arr.indexOf( str ) === -1 ? false : true;
		}
	}
	function appendCSS( cssText, element ) {
		var styleElement = document.createElement("style");
		styleElement.setAttribute("type","text/css"); 
		styleElement.innerHTML = cssText;
		var refNode = element.hasChildNodes() ? element.children[0] : null;
		element.insertBefore( styleElement, refNode );
	}
}

function svgString2Image( svgString, width, height, format, callback ) {
	var format = format ? format : 'png';
	var imgsrc = 'data:image/svg+xml;base64,'+ btoa( unescape( encodeURIComponent( svgString ) ) ); // Convert SVG string to data URL
	var canvas = document.createElement("canvas");
	var context = canvas.getContext("2d");
	canvas.width = width;
	canvas.height = height;
	var image = new Image();
	image.onload = function() {
		context.clearRect ( 0, 0, width, height );
		context.drawImage(image, 0, 0, width, height);
		
		canvas.toBlob( function(blob) {
			var filesize = Math.round( blob.length/1024 ) + ' KB';
			if ( callback ) callback( blob, filesize );
		});
		
	};
	image.src = imgsrc;
}

function svg2svg(svgNode,name)
{
	var serializer = new XMLSerializer();
	var svgData  = serializer.serializeToString(svgNode);
	var svgBlob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"});
	var svgUrl = URL.createObjectURL(svgBlob);
	var downloadLink = document.createElement("a");
	downloadLink.href = svgUrl;
	downloadLink.download = name;
	document.body.appendChild(downloadLink);
	downloadLink.click();
	document.body.removeChild(downloadLink);
}
function popupnamewindow(type,svgtype,id,x,y){
	console.log("in popup : ", id);
	d3.selectAll("#tooltip_imagename").remove() ;
	var div = d3.select("body").append("div")
	.attr("class", "tooltip")
	.attr("id","tooltip_imagename")
	.style("opacity", 1);

	div.transition()
	.duration(200)
	.style("opacity", .9);
	div.html(
			'<input type="text" placeholder="Image Name" id="imagename" />'+
			'<button id="saveimage" style = "none" type="button">save</button>'
			)
	.style("left", (x - 30) + "px")             
	.style("top", (y - 60) + "px")
	// .style("transform","translateX(-35%)!important")
	.attr("box-sizing","border")
	.attr("display","inline-block")
	.attr("min-width","1.5em")
	.attr("padding","2px")
	.attr("margin-left","0px")
	.attr("text-align","center")
	.attr("text-decoration","none")
	.attr("border","1px solid #111111")
	.attr("background-color","white");

	$("#saveimage").click(function(){
		let imagename = $("#imagename").val();
		if(svgtype == "seq" && (imagename == null || imagename == ""))
		{
			imagename = "seq";
		}else{
			imagename = "spectrum";
		}
		if(type == "svg"){
			d3.selectAll("#tooltip_imagename").remove() ;
			let svg_element = d3.selectAll("#"+id).node();
			svg2svg(svg_element,imagename);
		}
		if(type == "png"){
			d3.selectAll("#tooltip_imagename").remove() ;
			let l_svgContainer = d3.select("#"+id);
			let svgString = getSVGString(l_svgContainer.node());
			let svg_element = document.getElementById(id);
			let bBox = svg_element.getBBox();
			let width = bBox.width;
			let height = bBox.height ;
			svgString2Image( svgString, 2*width, 2*height, 'png', save ); 
			function save( dataBlob, filesize ){
				saveAs( dataBlob, imagename ); 
			}
		}
	})
}
