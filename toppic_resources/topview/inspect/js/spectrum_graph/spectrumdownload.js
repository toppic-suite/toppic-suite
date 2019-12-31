/**
 * Class to execute download operations for the spectrum graph
 */
class SpectrumDownload {
    //Unique Id values for graph and download buttons 
    static SPECDOWNLOADPNG = "spectrum_download_png";
    static SPECDOWNLOADSVG = "spectrum_download_svg";
    static SPECTRUMGRAPHID = "spectrum";
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
             .attr("xlink:href", "resources/js/spectrum_graph/images/download.png")	
    
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
            .attr("xlink:href", "resources/js/spectrum_graph/images/svg.png");
        
        group.append("svg:image").attr("id",SpectrumDownload.SPECDOWNLOADPNG)
            .attr("class","downloadgraph_button")
            .attr('x', (spectrumParameters.svgWidth-SpectrumDownload.spanWidth))
            .attr('y', SpectrumDownload.buttonTwo_Y)
            .attr('width', SpectrumDownload.buttonWidth)
            .attr('height', SpectrumDownload.buttonHeight)
            .attr("xlink:href", "resources/js/spectrum_graph/images/png.png");
        SpectrumDownload.download(svgId,spectrumParameters);
    }
    /**
     * On click action to download spectrum graph as SVG/PNG
     */
    static download = function(svgId,spectrumParameters){
        //On click action to download spectrum graph SVG in .svg format
        $("#"+SpectrumDownload.SPECDOWNLOADSVG).click(function(){
            
            $("#g_spectrum_button").remove();
            let name = "spectrum.svg"
            let svg_element = d3.selectAll("#"+SpectrumDownload.SPECTRUMGRAPHID).node();
            svg2svg(svg_element,name);
            SpectrumDownload.addDownloadRect(svgId,spectrumParameters);
        })
        //On click action to download spectrum graph PNG in .png format
        $("#"+SpectrumDownload.SPECDOWNLOADPNG).click(function(){
            $("#g_spectrum_button").remove();
            let l_svgContainer = d3.select("#"+SpectrumDownload.SPECTRUMGRAPHID);
            let svgString = getSVGString(l_svgContainer.node());
            //let svg_element = document.getElementById(SpectrumDownload.SPECTRUMGRAPHID);
            //let bBox = svg_element.getBBox();
            let width = spectrumParameters.svgWidth;
            let height = spectrumParameters.svgHeight ;
            
            svgString2Image( svgString, 2*width, 2*height, 'png', save ); 
            function save( dataBlob, filesize ){
                saveAs( dataBlob, 'spectrum.png' ); 
            }
            SpectrumDownload.addDownloadRect(svgId,spectrumParameters);
        })
    }
}