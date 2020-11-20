/**
 * Function to produce a pop up window to provide name and set name to 
 * the image while downloading the image of Graph SVG and Sequence SVG
 * @param {String} type - Provides if the image is downloaded as svg or png
 * @param {String} id - Provides the id of the svg to be downloaded
 * @param {Float} x - Provides coordinate on where to show a tooltip block to enter name of the image to be downloaded
 * @param {Float} y - Provides coordinate on where to show a tooltip block to enter name of the image to be downloaded
 */
function popupNameWindow(type,id,x,y){
	d3.selectAll("#tooltip_imagename").remove() ;
	var div = d3.select("body").append("div")
	.attr("class", "tooltip")
	.attr("id","tooltip_imagename")
	.style("opacity", 1);

	// Provides a tooltip to enter a name for the image to be downloaded
	div.transition()
	.duration(200)
	.style("opacity", .9);
	div.html( 
			'<input type="text" placeholder="Image Name" id="imagename" />'+
			'<button id="saveimage" style = "none" type="button">save</button>'
			)
	.style("left", (x - 30) + "px")  // x Coordinate of the position of the tooltip           
	.style("top", (y - 60) + "px")	// y Coordinate of the position of the tooltip 
	.attr("box-sizing","border")
	.attr("display","inline-block")
	.attr("min-width","1.5em")
	.attr("padding","2px")
	.attr("margin-left","0px")
	.attr("text-align","center")
	.attr("text-decoration","none")
	.attr("border","1px solid #111111")
	.attr("background-color","white");
	
	// On click action to save the image on click of download button
	$("#saveimage").click(function(){
		let imagename = $("#imagename").val();
		if( imagename == null || imagename == "")
		{
			imagename = "spectrum";
		}
		// Check if the image needs to be downloaded as svg
		if(type == "svg"){
			d3.selectAll("#tooltip_imagename").remove() ;
      console.log(id);
			let svgContainer = d3.select("#"+id);
			let svgElement = svgContainer.node();
			svg2svg(svgElement,imagename);
		}
		// Check if the image needs to be downloaded as png
		if(type == "png"){
			d3.selectAll("#tooltip_imagename").remove() ;
			let l_svgContainer = d3.select("#"+id);
			let svgString = getSVGString(l_svgContainer.node());
			let specParams =  new SpectrumParameters();
			let width = specParams.svgWidth;
			let height = specParams.svgHeight ;
			svgString2Image( svgString, 2*width, 2*height, 'png', save ); 
			function save( dataBlob, filesize ){
				saveAs( dataBlob, imagename ); 
			}
		}
	})
}