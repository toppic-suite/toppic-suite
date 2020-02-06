class rectBGColor{
    constructor(bgColorList,massShiftList){
        this.bgColorList = bgColorList;
    }

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
    
    /* Code to color the background of a occurence acids */
    rect_Backgroundcolor(x,y,id,color){
        /*	font-size 16px is equal to 12pt	*/
        let font_width = 12 ;
        /*	to draw the rect color uniformly */
        let font_height = 15 ;
        let svgContainer = d3.select("#"+id);
        
        //Temporary Fix
        console.log("svgContainer : ", svgContainer);
    
        svgContainer.append("rect")
                        .attr("class","rectBGColor")
                        .attr("x", x)
                        .attr("y", y-font_height)
                        .attr("width", font_width)
                        .attr("height", 20)
                        .attr("dy", "0em")
                        .style("fill", color)
                        .style("fill-opacity", ".6")
                        // .style("stroke-width", "1.5px");
    }
}