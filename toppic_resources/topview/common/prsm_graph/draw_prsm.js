function drawPrsm(svgId, para, prsm) {
  let svg = d3.select("body").select("#"+svgId);
  // Removes all the elements under SVG group 'svgGroup' everytime there this function is called
  $("#" + svgId).children().remove();
  initSvg(svg, para, prsm);
  if (prsm.showStartSkipped) {
    addStartSkippedLine(svg,para,prsm.startSkippedInfo);
  }
  if (prsm.showEndSkipped) {
    addEndSkippedLine(svg,para,prsm.displayFirstPos, prsm.displayLastPos + 1, prsm.endSkippedInfo, prsm.showStartSkipped);
  }
  addAnnos(svg, para, prsm.displayFirstPos, prsm.annotations, prsm.showStartSkipped);
  addAnnoBackground(svg, para, prsm.displayFirstPos, prsm.annotations, prsm.showStartSkipped);
  addAminoAcids(svg, para, prsm.residues, prsm.displayFirstPos, prsm.displayLastPos, prsm.showStartSkipped);
  if (para.showNum) {
    addPosNums(svg, para, prsm.displayFirstPos, prsm.displayLastPos, prsm.rowNum, prsm.showStartSkipped);
  }
  drawStartSymbol(svg, para, prsm.displayFirstPos, prsm.formFirstPos, prsm.showStartSkipped);
  drawEndSymbol(svg, para, prsm.displayFirstPos, prsm.displayLastPos, prsm.formLastPos, prsm.showStartSkipped);
  drawBreakPoints(svg, para, prsm.displayFirstPos, prsm.formFirstPos, prsm.breakPoints, prsm.showStartSkipped);
}

function initSvg(svg, para, prsm) {
  let [width, height] = para.getSvgSize(prsm.rowNum, prsm.showStartSkipped, prsm.showEndSkipped);
  //console.log("width", width, "height", height);
  svg.attr("width", width)
    .attr("height", height)
    .attr("font-family","FreeMono,Miltonian,monospace")
    .attr("font-size","16px")
    .style("fill", para.svgBackgroundColor);
}

/**
 * Get the terminated/skipped acid information on to the svg
 * @param {String} svg - Contians id of the SVG tag from html.
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} info - Contains the skipped information. 
 */
function addStartSkippedLine(svg, para, info) {
  let [x,y] = para.getSkipStartCoordinates();
  //	Get the coordinates to write the skip information at the start of acid
  svg.append("text")
    .attr("x", x)
    .attr("y", y)
    .attr("class","prsm_svg_group")
    .style("fill","black")
    .text(info);
}

/**
 * Get the terminated/skipped acid information on to the svg
 * @param {String} svg - Contians id of the SVG tag from html.
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} info - Contains the skipped information. 
 */
function addEndSkippedLine(svg, para, startPos, linePos, info, isStartSkipped) {
  let [x,y] = para.getSkipEndCoordinates(linePos, startPos);
  //console.log("x", x, "y", y, "info", info);

  if (isStartSkipped){
    y = y + para.middleMargin;
  }

  //	Get the coordinates to write the skip information at the start of acid
  svg.append("text")
    .attr("x", x)
    .attr("y", y)
    .attr("class","prsm_svg_group")
    .style("fill","black")
    .text(info);
}

/**
 * Draw the sequence on to svg
 * @param {String} svg - Svg element.
 * @param {String} seq - Contains the protein sequence
 */
function addAminoAcids (svg, para, residues, start, end, isStartSkipped) {
  let aaGroup = svg.append("g")
    .attr("id", "amino_acid")
    .attr("class","prsm_svg_group")

  for (let i = 0; i < residues.length; i++) {
    let pos = residues[i].position;
    if (pos >= start && pos <= end) {
      let letter = residues[i].acid;
      let xPos = para.getX(pos, start); 
      let yPos = para.getY(pos, start);
      let color = residues[i].color;

      if (isStartSkipped){
        yPos = yPos + para.middleMargin;
      }
      
      //console.log(letter, xPos, yPos, color);
      aaGroup.append("text")
        .attr("x", xPos)
        .attr("y", yPos) 
        .style("fill", color)
        .text(letter);
    }
  }
}

/**
 * Put the numerical positions at the start and end of each row of the sequence
 * @param {*} para Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} prsm - Contains the complete information of prsm. 
 * @param {String} id - Contians id of the SVG tag from html.
 */
function addPosNums (svg, para, startPos, lastPos, rowNum, isStartSkipped) {
  let numGroup = svg.append("g")
    .attr("id", "pos_num")
    .attr("class","prsm_svg_group");

  for (let i = 0; i < rowNum; i++) {
    let leftPos = startPos + i * para.rowLength;
    let [lx,ly] = para.getLeftNumCoordinates (leftPos, startPos);
    let leftText = leftPos + 1;
   
    if (isStartSkipped){
      ly = ly + para.middleMargin;
    }

    numGroup.append("text")
      .attr("x",lx)
      .attr("y",ly)
      .text(leftText)
      .style("text-anchor", "end")
      .style("fill", "black");

    let rightPos = leftPos + para.rowLength - 1;
    if (rightPos > lastPos) {
      rightPos = lastPos;
    }
    let [rx,ry] = para.getRightNumCoordinates (rightPos, startPos);
    let rightText = rightPos + 1;

    if (isStartSkipped){
      ry = ry + para.middleMargin;
    }
    
    numGroup.append("text")
      .attr("x",rx)
      .attr("y",ry)
      .text(rightText)
      .style("text-anchor", "start")
      .style("fill", "black");
  }
}

/**
 * Draw the annotations to show the start and end position of the sequence
 * @param {String} svg - Contians id of the SVG tag from html.
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 */
function drawStartSymbol(svg, para, dispFirstPos, formFirstPos, isStartSkipped) {
  if (dispFirstPos !== formFirstPos) {
		let [x,y] = para.getAACoordinates(formFirstPos, dispFirstPos);
		x = x - (para.letterWidth/2) ;

    if (isStartSkipped){
      y = y + para.middleMargin;
    }
    
		let coordinates = (x)+","+(y+2)+ " " +(x+5)+","+ (y+2)+" "+(x+5)+","+(y-12)+ " "+(x) + ","+(y-12);
    svg.append("polyline")
      .attr("class","prsm_svg_group")
      .attr("points", coordinates)
      .style("fill", "none")
      .style("stroke", "red")
      .style("stroke-width", "1.3");
  }
}

function drawEndSymbol(svg, para, dispFirstPos, dispLastPos, formLastPos, isStartSkipped) {
  if (dispLastPos !== formLastPos) {
		let [x,y] = para.getAACoordinates(formLastPos, dispFirstPos);
		x = x + (para.letterWidth/2) ;
    
    if (isStartSkipped){
      y = y + para.middleMargin;
    }
    
    let coordinates = (x+7)+","+(y-12)+ " " +(x+2)+","+ (y-12)+" "+(x+2)+","+(y+2)+ " "+(x+7) + ","+(y+2);

    svg.append("polyline")
          .attr("class","prsm_svg_group")
					.attr("points", coordinates)
					.style("fill", "none")
					.style("stroke", "red")
					.style("stroke-width", "1.3") ;
	}
}


/**
 * Invoke drawBreakPoints method 
 * @param {String} svg - Contians id of the SVG tag from html.
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} breakPoints - Json object with information of the position
 */
function drawBreakPoints(svg, para, dispFirstPos, formFirstPos, breakPoints, isStartSkipped) {
  let bpGroup = svg.append("g")
    .attr("id", "break_point")
    .attr("class","prsm_svg_group");
  for (let i = 0; i < breakPoints.length; i++) {
    let bp = breakPoints[i];
    //console.log(bp.position, dispFirstPos);
    let [x,y] = para.getBpCoordinates(bp.position, dispFirstPos);

    if (isStartSkipped){
      y = y + para.middleMargin;
    }
    
    // Setting polyline coordinates
    let coordinates; 
    if (bp.existNIon && !bp.existCIon) {
      coordinates = (x-2)+","+(y-13)+ " " +(x+4)+","+ (y-11)+" "+(x+4)+","+(y+2);
    }
    else if (!bp.existNIon && bp.existCIon) {
      coordinates = (x+4)+","+ (y-11)+" "+(x+4)+","+(y+2)+ " "+(x+10) + ","+(y+5);
    }
    else {
	    coordinates =  (x-2)+","+(y-13)+ " " + (x+4)+","+ (y-11)+" "+(x+4)+","+(y+2)+ " "+(x+10) + ","+(y+5);
    }
    bpGroup.append("polyline")
      .attr("points", coordinates)
      .style("fill", "none")
      .style("stroke", "#1e90ff")
      .style("stroke-width", "1");	
    //	Rectangle to have flexible on click and on mouse actions	
    let anno = bp.anno;
    let bpIonPos = bp.position - formFirstPos;
    bpGroup.append("rect")
      .attr("class","break_point")
      .attr("x", x)
      .attr("y", y-14)
      .attr("width", 13)
      .attr("height", 23)
      .attr("ion_pos", bpIonPos)
      .style("opacity", 0)
      .on("mouseover", function(){
        appendBpAnno(anno);
      })
      .on("mouseout", function(d){
        removeBpAnno();	
      }); 
  }
}

/**
 * Function to add annotation to the polylines on mouseOver
 * @param {String} anno - Contains consolidated annotation to display
 */
function appendBpAnno(anno) {
  // tooltip is a bootstrap class
  let div = d3.select("body").append("div")	
    .attr("id", "break_point_anno")
    .attr("class", "tooltip")				
    .style("opacity", 0); 
  div.transition()		
    .duration(10)		
    .style("opacity", .9);
  div.html(anno)	
    .style("left", (d3.event.pageX)  + "px")		
    .style("top", (d3.event.pageY - 28)+ "px") ;
}

/**
 * Function to remove annotation to the polylines on mouseOver
 */
function removeBpAnno() {
	d3.selectAll("#break_point_anno").remove();
}


/**
 * Color the background of the occurence
 * @param {String} svg - Contians id of the SVG tag from html.
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 * @param {Object} shifts - Contains the complete information of mass shifts. 
 */
function addAnnos(svg, para, firstPos, annos, isStartSkipped) {
  let annoGroup = svg.append("g")
    .attr("id", "mass_shift_anno")
    .attr("class","prsm_svg_group");
  let yShiftIdx = 0;
  for (let i = 0; i < annos.length; i++) {
    let leftPos = annos[i].leftPos;
    let [x1,y1] = para.getAACoordinates(leftPos, firstPos);
    //console.log(leftPos, x1, y1);
    let overlap = false;
    if (i > 0) {
      let prevPos = parseInt(annos[i-1].leftPos);
      let [x2,y2] = para.getAACoordinates(prevPos, firstPos);
			// subtract -2 for CSS and alignment purpose
      let annoLen = annos[i-1].annoText.length * (para.fontWidth-2);
      //console.log(shifts[i-1].anno, annoLen, x2, x1);
			if(y1 == y2 && (x1-x2 < annoLen)) {
        overlap = true;
			}
    }
    // change y shift based on overlapping
    if (overlap) {
      yShiftIdx = (yShiftIdx + 1) % 2;
    }
    else {
      yShiftIdx = 0;
    }
    y1 = y1 + para.modAnnoYShifts[yShiftIdx];
    let annoText = annos[i].annoText;
    if (isStartSkipped){
      y1 = y1 + para.middleMargin;
    }

    annoGroup.append("text")
      .attr("x", x1)
      .attr("y", y1)
      .text(annoText)
      .attr("fill","black")
      .attr("font-size","15px");
  }
}

/**
 * Code to color the background of a occurence acids
 * @param {Integer} x - Contains x coordinate of the start postion to add background color
 * @param {Integer} y - Contains y coordinate of the end position to add background color
 * @param {String} id - Contains id of the SVG tag from html.
 * @param {Integer} width - Contains width to whuich backgroud color has to be added
 * @param {Object} para - Contains parameters of width, letter space etc., to draw SVG
 */
function addAnnoBackground(svg, para, firstPos, annos, isStartSkipped) {
  let bgGroup = svg.append("g")
    .attr("id", "background")
    .attr("class","prsm_svg_group");
  for (let i = 0; i < annos.length; i++) {
    let leftPos = annos[i].leftPos;
    let rightPos = annos[i].rightPos;
    //console.log("left", leftPos, "right", rightPos);
    let startRow = Math.floor((leftPos -firstPos)/para.rowLength);
    let endRow = Math.floor((rightPos - 1 - firstPos)/para.rowLength);
    //console.log("start row", startRow, "end row", endRow);
    for (let j = startRow; j <= endRow; j++) {
      let rowLeft = firstPos + j * para.rowLength;
      let rowRight = rowLeft + para.rowLength - 1;
      if (j == startRow) {
        rowLeft = leftPos;
      }
      if (j == endRow) {
        rowRight = rightPos-1;
      }
      //console.log("row left", rowLeft, "row right", rowRight);
      let [x1,y1] = para.getAACoordinates(rowLeft, firstPos);
      let [x2,y2] = para.getAACoordinates(rowRight, firstPos);

      if (isStartSkipped){
        y1 = y1 + para.middleMargin;
      }

      if (isStartSkipped){
        y2 = y2 + para.middleMargin;
      }
    
      bgGroup.append("rect")
        .attr("x", x1-para.fontWidth*0.1)
        .attr("y", y1-para.fontHeight*0.7)
        .attr("width", x2-x1 + para.fontWidth)
        .attr("height", para.fontHeight)
        .style("fill", para.massShiftColor)
        .style("fill-opacity", ".4")
        .style("stroke-width", "1.5px");
    }
  }
}

