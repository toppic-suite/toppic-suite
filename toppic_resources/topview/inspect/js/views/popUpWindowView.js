function popupnamewindow(type,svgtype,id,x,y){
	// console.log("in popup : ", id);
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
		if(svgtype === "seq" && (imagename === null || imagename === ""))
		{
			imagename = "seq";
		}else{
			imagename = "spectrum";
		}
		if(type === "svg"){
			d3.selectAll("#tooltip_imagename").remove() ;
			let svg_element = d3.selectAll("#"+id).node();
			svg2svg(svg_element,imagename);
		}
		if(type === "png"){
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
