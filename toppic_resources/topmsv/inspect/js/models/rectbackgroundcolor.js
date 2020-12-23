/**
 * Class to add background color to the sequence when mass shift added by clicking on the acid
 */
class rectBGColor{
    constructor(bgColorList,massShiftList){
        this.bgColorList = bgColorList;
    }
    /**
     * Function to get all the positions and color to be add background color
     * @param {Integer} position - Contains position at which mass shif to be added  
     * @param {String} color - Contians background color(comes from dropdown list) to be added to the mass shift
     */
    getMassListToColorBG(position, color){
        let bgColorList = this.bgColorList;
        let len = bgColorList.length;
        if(len > 0)
        {
            let match = false;
            for(let i=0; i<len ;i++)
            {
                match = false;
                if(bgColorList[i].position == position)
                {
                    match = true;
                    bgColorList[i].color = color;
                    break;
                }
            }
            if(!match)
            {
                bgColorList.push({position:position,color:color})
            }
        }else{
            bgColorList.push({position:position,color:color})
        }

        return bgColorList;
    }
    /**
     * Function Adds selected background color to the amino acid
     * @param {Array} bgColorList - Set all the aminoacids with desired background colors selected from drop down
     */
    setBackGroundColorOnMassShift(bgColorList){
        let len = bgColorList.length;
        let para = parameters();
        for(let i=0;i<len;i++)
        {
            let x,y;
            x = getX(para,bgColorList[i].position);
            y = getY(para,bgColorList[i].position);
            this.rect_Backgroundcolor(x,y,"seqsvg",bgColorList[i].color)
        }
    }
    
    /**
     * Function to add background color at occurence acids
     * @param {Integer} x - Starting x position to add background color
     * @param {Integer} y - Strating y position to add background color
     * @param {String} id - Contains id of svg tag from html
     * @param {String} color - Background color to be added
     */
    rect_Backgroundcolor(x,y,id,color){
        /*	font-size 16px is equal to 12pt	*/
        let font_width = 15 ;
        /*	to draw the rect color uniformly */
        let font_height = 15 ;
        let svgContainer = d3.select("#"+id);
        x = x-3;
        svgContainer.append("rect")
                        .attr("class","rectBGColor")
                        .attr("x", x)
                        .attr("y", y-font_height)
                        .attr("width", font_width)
                        .attr("height", 20)
                        .attr("dy", "0em")
                        .style("fill", color)
                        .style("fill-opacity", ".6")
    }
}