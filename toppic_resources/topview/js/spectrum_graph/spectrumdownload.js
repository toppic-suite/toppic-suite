/**
 * Class to execute download operations for the spectrum graph
 */
class SpectrumDownload {
    //Unique Id values for graph and download buttons 
    static SPECDOWNLOADPNG = "spectrum_download_png";
    static SPECDOWNLOADSVG = "spectrum_download_svg";
    static downloadbuttonX = 41;
    static spanWidth = 50;
    static downloadButtonWidth = 30;
    static downloadButtonHeight = 30;
    static buttonWidth = 50;
    static buttonHeight = 28;
    static buttonOne_Y = 33;
    static buttonTwo_Y = 60;
    //default costructor
    constructor(){
    }
    /**
     * Draw a rectangular block to keep download image
     */
    static addDownloadRect = function(svgId,spectrumParameters){
        d3.select("#g_spectrum_button").remove();
        let svg = d3.select("body").select(svgId);
        let group = svg.append("g").attr("id","g_spectrum_button");
    
        group.append("svg:image").attr("id","spectrum_download")
             .attr('x', (spectrumParameters.svgWidth-SpectrumDownload.downloadbuttonX))
             .attr('y', 0)//top most point in svg
             .attr('width', SpectrumDownload.downloadButtonWidth)
             .attr('height', SpectrumDownload.downloadButtonHeight)
             .attr("xlink:href", "../js/spectrum_graph/images/download.png")	
    
        d3.selectAll("#g_spectrum_button").on('mouseover',function(){
                d3.select(this).style("cursor", "pointer"); 
                })
                .on('mouseout',function(){
                    d3.select(this).style("cursor", "default"); 
                });
        $("#spectrum_download").on("click",function(){
            d3.selectAll(".downloadgraph_button").remove();
            SpectrumDownload.addDownloadbuttons(svgId,spectrumParameters);
        });
    }
    /**
     * Setting images as download buttons
     */
    static addDownloadbuttons = function(svgId,spectrumParameters){
        let group = d3.select("#g_spectrum_button");
        group.append("svg:image").attr("id",SpectrumDownload.SPECDOWNLOADSVG)
            .attr("class","downloadgraph_button")
            .attr('x', (spectrumParameters.svgWidth-SpectrumDownload.spanWidth))
            .attr('y', SpectrumDownload.buttonOne_Y)
            .attr('width', SpectrumDownload.buttonWidth)
            .attr('height', SpectrumDownload.buttonHeight)
            .attr("xlink:href", "../js/spectrum_graph/images/svg.png")
            .on("click",function(){
                SpectrumDownload.popup("svg",spectrumParameters,svgId);
            });
        
        group.append("svg:image").attr("id",SpectrumDownload.SPECDOWNLOADPNG)
            .attr("class","downloadgraph_button")
            .attr('x', (spectrumParameters.svgWidth-SpectrumDownload.spanWidth))
            .attr('y', SpectrumDownload.buttonTwo_Y)
            .attr('width', SpectrumDownload.buttonWidth)
            .attr('height', SpectrumDownload.buttonHeight)
            .attr("xlink:href", "../js/spectrum_graph/images/png.png")
            .on("click",function(){
                SpectrumDownload.popup("png",spectrumParameters,svgId);
            });
        //SpectrumDownload.download(svgId,spectrumParameters);
    }
    /**
     * On click action to download spectrum graph as SVG/PNG
     */
    static download = function(type,spectrumParameters,svgId,imagename){
       if(type == "svg")
       {
            $("#g_spectrum_button").remove();
            let  name = imagename+".svg";
            let svg_element = d3.selectAll(svgId).node();
            svg2svg(svg_element,name);
            
       }else{
            $("#g_spectrum_button").remove();
            let l_svgContainer = d3.select(svgId);
            let svgString = getSVGString(l_svgContainer.node());
            //let svg_element = document.getElementById(SpectrumDownload.SPECTRUMGRAPHID);
            //let bBox = svg_element.getBBox();
            let name = imagename+".png";
            let width = spectrumParameters.svgWidth;
            let height = spectrumParameters.svgHeight ;
            svgString2Image( svgString, 2*width, 2*height, 'png', save ); 
            function save( dataBlob, filesize ){
                saveAs( dataBlob, name ); 
            }
       } 
       SpectrumDownload.addDownloadRect(svgId,spectrumParameters);
       $('#tooltip_imagename').remove();   
            
    }
    static popup =  function(type,spectrumParameters,svgId){
        console.log("in popup");
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
        .style("left", (d3.event.pageX - 30) + "px")             
        .style("top", (d3.event.pageY - 45) + "px")
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
            if( imagename == null || imagename == "")
            {
                imagename = "spectrum";
            }
            SpectrumDownload.download(type,spectrumParameters,svgId, imagename);
        })
        
    }
}