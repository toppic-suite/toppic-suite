class InteRtGraph {
    xScale_g;
    fixedLine_g;

    inteRtArray;
    svg_ID;
    rt_ID;
    inte_ID;

    padding;
    width;
    height;

    constructor(svg_ID, inteRtArray, rt_ID = 'rt-hover', inte_ID = 'intensity-hover', height = 120, width = 1100, padding = {top: 10, right: 10, bottom: 50, left: 80}) {
        this.inteRtArray = inteRtArray;

        this.svg_ID = "#"+svg_ID;
        this.rt_ID = rt_ID;
        this.inte_ID = inte_ID;

        this.width = width;
        this.height = height;
        this.padding = padding;
    }

    set padding(obj) {
        if (!obj.top || !obj.bottom || !obj.left || !obj.right) {
            console.log("padding should be an object that contanins top, bottom, left and right");
            return;
        }
        this.svg_padding = obj;
    }

    get padding() {
        this.svg_padding;
    }

    set width(value) {
        if (value > this.svg_padding.left + this.svg_padding.right) {
            this.svg_width = value;
        } else {
            console.log("width should be larger than padding left plus padding right!");
        }
    }

    get width() {
        return this.svg_width;
    }

    set height(value) {
        if (value > this.svg_padding.top + this.svg_padding.bottom) {
            this.svg_height = value;
        } else {
            console.log("height should be larger than padding top plus padding bottom");
        }
    }

    get height() {
        return this.svg_height;
    }

    drawGraph() {

        let inteRtArray = this.inteRtArray;
        let padding = this.padding;
        let rt_ID = this.rt_ID;
        let inte_ID = this.inte_ID;
        let width = this.width;
        let height = this.height;

        let maxInte = d3.max(this.inteRtArray, function(d) {
            return d.inteSum;
        });
    
        let formatPercent = d3.format(".0%");
    
        this.inteRtArray.forEach(function (element) {
            element.rt = element.rt/60;
            element.intePercentage = element.inteSum/maxInte;
        });
        this.inteRtArray.sort((a,b) => {a.rt > b.rt ? 1:-1});
    
        let min = d3.min(this.inteRtArray, function(d) {
            return d.intePercentage;
        });
        let max = d3.max(this.inteRtArray, function(d) {
            return d.intePercentage;
        });
    
        let minRT = d3.min(this.inteRtArray, function(d) {
            return d.rt;
        });
        let maxRT = d3.max(this.inteRtArray, function(d) {
            return d.rt;
        });
        let xScale = d3.scaleLinear()
            .domain([0, maxRT+5])
            .range([0, this.width - this.padding.left - this.padding.right]);
        this.xScale_g = xScale;

        let yScale = d3.scaleLinear()
            .domain([0, max])
            .range([this.height - this.padding.top - this.padding.bottom, 0]);
    
        let svg = d3.select(this.svg_ID)
            .append('svg')
            .attr('viewBox', "0 0 "+ this.width + " "+this.height)
            .attr('preserveAspectRatio', 'xMidYMid meet')
            .attr('width', '100%')
            .attr('height', '100%');
    
        let xAxis = d3.axisBottom()
            .scale(xScale)
            .ticks(20);
    
        let yAxis = d3.axisLeft()
            .scale(yScale)
            .tickFormat(formatPercent)
            .ticks(5);
    
        svg.append('g')
            .attr('class', 'axis')
            .attr('transform', 'translate(' + this.padding.left + ',' + (this.height - this.padding.bottom) + ')')
            .call(xAxis);
        // text label for the x axis
        svg.append("text")
            // .attr("fill", "black")//set the fill here
            .attr("transform",
                "translate(" + ((this.width+this.padding.left-this.padding.right)/2) + " ," +
                (this.height - this.padding.bottom + 35) + ")")
            .style("text-anchor", "middle")
            .text("Retention Time (mins)");
    
        svg.append('g')
            .attr('class', 'axis')
            .attr('transform', 'translate(' + this.padding.left + ',' + this.padding.top + ')')
            .call(yAxis);
        // text label for the y axis
        svg.append("text")
            .attr("transform", "rotate(-90)")
            .attr("y", 20)
            .attr("x",0 - (this.height / 2) + 20)
            .attr("dy", "1em")
            .style("text-anchor", "middle")
            .text("Intensity");
    
        let linePath = d3.line()
            .x(function(d){ return xScale(d.rt) })
            .y(function(d){ return yScale(d.intePercentage) }).curve(d3.curveBasis);
    
        svg.append('g')
            .append('path')
            .attr('class', 'line-path')
            .attr('transform', 'translate(' + this.padding.left + ',' + this.padding.top + ')')
            .attr('d', linePath(this.inteRtArray))
            .attr('fill', 'none')
            .attr('stroke-width', 1)
            .attr('stroke', 'black');

    
        //Line chart mouse over
        let hoverLineGroup = svg.append("g")
            .attr("class", "hover-line");
        let hoverLine = hoverLineGroup
            .append("line")
            .attr("stroke", "#ff0000")
            .attr("x1", this.padding.left).attr("x2", this.padding.left)
            .attr("y1", this.padding.top).attr("y2", this.height-this.padding.bottom);
        let fixedLine = hoverLineGroup
            .append("line")
            .attr("stroke", "#ff8000")
            .attr("x1", this.padding.left).attr("x2", this.padding.left)
            .attr("y1", this.padding.top).attr("y2", this.height-this.padding.bottom);
        this.fixedLine_g = fixedLine;
    
        hoverLine.style("opacity", 1e-6);

        svg
            .on("mouseout", hoverMouseOff)
            .on("mouseover mousemove touchmove", hoverMouseOn)
            .on("click", mouseClick);
    
    
        let bisectRT = d3.bisector(function(d) { return d.rt; }).right;
    
        function mouseClick() {
            let mouse_x = d3.mouse(this)[0];
            let mouse_y = d3.mouse(this)[1];
            let maxMouse = xScale(maxRT);
            let mouseRT = xScale.invert(mouse_x-padding.left);
            let i = bisectRT(inteRtArray, mouseRT); // returns the index to the current data item
    
            if(i>0 && i < inteRtArray.length && mouse_y < height-padding.bottom && mouse_y > padding.top) {
                let d0 = inteRtArray[i - 1];
                let d1 = inteRtArray[i];
                // work out which date value is closest to the mouse
                let d = mouseRT - d0.rt > d1.rt - mouseRT ? d1 : d0;
                fixedLine.attr("x1", mouse_x).attr("x2", mouse_x);
                fixedLine.style("opacity", 1);
                // findNextLevelOneScan(d.scanNum);
            } else if (i === inteRtArray.length && mouse_x -padding.left<= maxMouse+1 && mouse_y < height-padding.bottom && mouse_y > padding.top)
            {
                let d = inteRtArray[i-1];
                fixedLine.attr("x1", mouse_x).attr("x2", mouse_x);
                fixedLine.style("opacity", 1);
                // findNextLevelOneScan(d.scanNum);
            } else {
                //fixedLine.style("opacity", 1e-6);
            }
        }

        function hoverMouseOn() {
            let mouse_x = d3.mouse(this)[0];
            let mouse_y = d3.mouse(this)[1];
            let maxMouse = xScale(maxRT);
            hoverLine.attr("x1", mouse_x).attr("x2", mouse_x);
            hoverLine.style("opacity", 1);
            let graph_y = yScale.invert(mouse_y);
            let graph_x = xScale.invert(mouse_x-padding.left);
    
            let mouseRT = xScale.invert(mouse_x-padding.left);
            let i = bisectRT(inteRtArray, mouseRT); // returns the index to the current data item
            if(i>0 && i < inteRtArray.length && mouse_y < height-padding.bottom && mouse_y > padding.top) {
                let d0 = inteRtArray[i - 1];
                let d1 = inteRtArray[i];
                // work out which date value is closest to the mouse
                let d = mouseRT - d0.rt > d1.rt - mouseRT ? d1 : d0;
                if(document.getElementById(rt_ID)) {
                    document.getElementById(rt_ID).innerHTML = Math.round(d.rt * 100)/100;
                }
                if (document.getElementById(inte_ID)) {
                    document.getElementById(inte_ID).innerHTML = d.inteSum.toExponential(2);
                }
    
                hoverLine.style("opacity", 1);
            } else if (i === inteRtArray.length&& mouse_x-padding.left <= maxMouse+1 && mouse_y < height-padding.bottom && mouse_y > padding.top)
            {
                let d = inteRtArray[i-1];
                if(document.getElementById(rt_ID)) {
                    document.getElementById(rt_ID).innerHTML = Math.round(d.rt * 100)/100;
                }
                if (document.getElementById(inte_ID)) {
                    document.getElementById(inte_ID).innerHTML = d.inteSum.toExponential(2);
                }
    
                hoverLine.style("opacity", 1);
            } else {
                if(document.getElementById(rt_ID)) {
                    document.getElementById(rt_ID).innerHTML = 0;
                }
                if (document.getElementById(inte_ID)) {
                    document.getElementById(inte_ID).innerHTML = 0;
                }
                hoverLine.style("opacity", 0);
            }
        }
        function hoverMouseOff() {
            hoverLine.style("opacity", 1e-6);
        }
    }

    moveLine(rt) {
        let newX = this.xScale_g(rt) + this.padding.left;
        this.fixedLine_g.attr("x1", newX).attr("x2", newX);
        this.fixedLine_g.style("opacity", 1);
    }
}