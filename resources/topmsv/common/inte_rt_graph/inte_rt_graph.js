"use strict";
class InteRtGraph {
    constructor(svg_ID, inteRtArray, onClickFunc = () => { }, scanNum_ID = 'scan-hover', rt_ID = 'rt-hover', inte_ID = 'intensity-hover', height = 120, width = 1100, padding = { left: 80, right: 10, head: 10, bottom: 50 }) {
        this.rawInteRtArray_ = [];
        this.inteRtArray_ = [];
        this.rawInteRtArray_ = inteRtArray;
        this.svg_ID_ = "#" + svg_ID;
        this.rt_ID_ = rt_ID;
        this.inte_ID_ = inte_ID;
        this.scanNum_ID_ = scanNum_ID;
        this.width_ = width;
        this.height_ = height;
        this.padding_ = padding;
        this.onClickFunc = onClickFunc;
    }
    /*
        set padding(obj: Padding) {
            if (!obj.head || !obj.bottom || !obj.left || !obj.right) {
                console.log("padding should be an object that contanins top, bottom, left and right");
                return;
            }
            this.svg_padding_ = obj;
        }
    
        get padding(): Padding {
            return this.svg_padding_;
        }
    
        set width(value) {
            if (value > this.svg_padding_.left + this.svg_padding_.right) {
                this.svg_width_ = value;
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
    */
    drawGraph() {
        let inteRtArray = this.inteRtArray_;
        let padding = this.padding_;
        let rt_ID = this.rt_ID_;
        let inte_ID = this.inte_ID_;
        let scanNum_ID = this.scanNum_ID_;
        let width = this.width_;
        let height = this.height_;
        if (!this.rawInteRtArray_.length) {
            return;
        }
        this.rawInteRtArray_.forEach(function (element) {
            let inteRtSingle = {
                inteSum: parseFloat(element.inteSum),
                rt: parseFloat(element.rt),
                intePercentage: -1,
                scanNum: element.scanNum
            };
            inteRtArray.push(inteRtSingle);
        });
        let maxInte = d3.max(this.inteRtArray_, function (d) {
            return d.inteSum;
        });
        let formatPercent = d3.format(".0%");
        this.inteRtArray_.forEach(function (element) {
            if (maxInte == undefined) {
                console.error("ERROR: invalid intensity in inte-rt graph");
                return;
            }
            element.rt = element.rt / 60;
            element.intePercentage = element.inteSum / maxInte;
        });
        //this.inteRtArray_.sort((a,b) => {a.rt > b.rt ? 1:-1});
        let min = d3.min(this.inteRtArray_, function (d) {
            return d.intePercentage;
        });
        let max = d3.max(this.inteRtArray_, function (d) {
            return d.intePercentage;
        });
        let minRT = d3.min(this.inteRtArray_, function (d) {
            return d.rt;
        });
        let maxRT = d3.max(this.inteRtArray_, function (d) {
            return d.rt;
        });
        if (!min || !max || !minRT || !maxRT) {
            console.error("ERROR: invalid intensity or rt in inte-rt graph");
            return;
        }
        let xScale = d3.scaleLinear()
            .domain([0, maxRT + 5])
            .range([0, this.width_ - this.padding_.left - this.padding_.right]);
        this.xScale_g_ = xScale;
        let yScale = d3.scaleLinear()
            .domain([0, max])
            .range([this.height_ - this.padding_.head - this.padding_.bottom, 0]);
        let svg = d3.select(this.svg_ID_)
            .append('svg')
            .attr('viewBox', "0 0 " + this.width_ + " " + this.height_)
            .attr('preserveAspectRatio', 'xMidYMid meet')
            .attr('width', '100%')
            .attr('height', '100%');
        //@ts-ignore   
        let xAxis = d3.axisBottom()
            //@ts-ignore   
            .scale(xScale)
            .ticks(20);
        //@ts-ignore   
        let yAxis = d3.axisLeft()
            //@ts-ignore   
            .scale(yScale)
            //@ts-ignore   
            .tickFormat(formatPercent)
            .ticks(5);
        svg.append('g')
            .attr('class', 'axis')
            .attr('transform', 'translate(' + this.padding_.left + ',' + (this.height_ - this.padding_.bottom) + ')')
            .call(xAxis);
        // text label for the x axis
        svg.append("text")
            // .attr("fill", "black")//set the fill here
            .attr("transform", "translate(" + ((this.width_ + this.padding_.left - this.padding_.right) / 2) + " ," +
            (this.height_ - this.padding_.bottom + 35) + ")")
            .style("text-anchor", "middle")
            .text("Retention Time (mins)");
        svg.append('g')
            .attr('class', 'axis')
            .attr('transform', 'translate(' + this.padding_.left + ',' + this.padding_.head + ')')
            .call(yAxis);
        // text label for the y axis
        svg.append("text")
            .attr("transform", "rotate(-90)")
            .attr("y", 20)
            .attr("x", 0 - (this.height_ / 2) + 20)
            .attr("dy", "1em")
            .style("text-anchor", "middle")
            .text("Intensity");
        let linePath = d3.line()
            //@ts-ignore   
            .x(function (d) { return xScale(d.rt); })
            //@ts-ignore   
            .y(function (d) { return yScale(d.intePercentage); }).curve(d3.curveBasis);
        svg.append('g')
            .append('path')
            .attr('class', 'line-path')
            .attr('transform', 'translate(' + this.padding_.left + ',' + this.padding_.head + ')')
            //@ts-ignore   
            .attr('d', linePath(this.inteRtArray_))
            .attr('fill', 'none')
            .attr('stroke-width', 1)
            .attr('stroke', 'black');
        //Line chart mouse over
        let hoverLineGroup = svg.append("g")
            .attr("class", "hover-line");
        let hoverLine = hoverLineGroup
            .append("line")
            .attr("stroke", "#ff0000")
            .attr("x1", this.padding_.left).attr("x2", this.padding_.left)
            .attr("y1", this.padding_.head).attr("y2", this.height_ - this.padding_.bottom);
        let fixedLine = hoverLineGroup
            .append("line")
            .attr("stroke", "#ff8000")
            .attr("x1", this.padding_.left).attr("x2", this.padding_.left)
            .attr("y1", this.padding_.head).attr("y2", this.height_ - this.padding_.bottom);
        this.fixedLine_g_ = fixedLine;
        hoverLine.style("opacity", 1e-6);
        let self = this;
        svg
            .on("mouseout", hoverMouseOff)
            .on("mouseover mousemove touchmove", hoverMouseOn)
            .on("click", mouseClick);
        let bisectRT = d3.bisector(function (d) { return d.rt; }).right;
        function mouseClick() {
            let mouse_x = d3.mouse(d3.event.currentTarget)[0];
            let mouse_y = d3.mouse(d3.event.currentTarget)[1];
            //@ts-ignore 
            let maxMouse = xScale(maxRT); // maxRT is already checked if it is undefined
            let mouseRT = xScale.invert(mouse_x - padding.left);
            let i = bisectRT(inteRtArray, mouseRT); // returns the index to the current data item
            if (!maxMouse || !mouseRT) {
                console.error("ERROR: invalid rt value in intensity-rt graph");
                return;
            }
            if (i > 0 && i < inteRtArray.length && mouse_y < height - padding.bottom && mouse_y > padding.head) {
                let d0 = inteRtArray[i - 1];
                let d1 = inteRtArray[i];
                // work out which date value is closest to the mouse
                let d = mouseRT - d0.rt > d1.rt - mouseRT ? d1 : d0;
                fixedLine.attr("x1", mouse_x).attr("x2", mouse_x);
                fixedLine.style("opacity", 1);
                //@ts-ignore 
                self.onClickFunc(d.scanNum);
            }
            else if (i === inteRtArray.length && mouse_x - padding.left <= maxMouse + 1 && mouse_y < height - padding.bottom && mouse_y > padding.head) {
                let d = inteRtArray[i - 1];
                fixedLine.attr("x1", mouse_x).attr("x2", mouse_x);
                fixedLine.style("opacity", 1);
                //@ts-ignore 
                self.onClickFunc(d.scanNum);
            }
            else {
                //fixedLine.style("opacity", 1e-6);
            }
        }
        function hoverMouseOn() {
            let mouse_x = d3.mouse(d3.event.currentTarget)[0];
            let mouse_y = d3.mouse(d3.event.currentTarget)[1];
            //@ts-ignore 
            let maxMouse = xScale(maxRT); // maxRT is already checked if it is undefined
            hoverLine.attr("x1", mouse_x).attr("x2", mouse_x);
            hoverLine.style("opacity", 1);
            let graph_y = yScale.invert(mouse_y);
            let graph_x = xScale.invert(mouse_x - padding.left);
            let mouseRT = xScale.invert(mouse_x - padding.left);
            if (!maxMouse || !mouseRT) {
                console.error("ERROR: invalid rt value in intensity-rt graph");
                return;
            }
            let i = bisectRT(inteRtArray, mouseRT); // returns the index to the current data item
            if (i > 0 && i < inteRtArray.length && mouse_y < height - padding.bottom && mouse_y > padding.head) {
                let d0 = inteRtArray[i - 1];
                let d1 = inteRtArray[i];
                // work out which date value is closest to the mouse
                let d = mouseRT - d0.rt > d1.rt - mouseRT ? d1 : d0;
                if (document.getElementById(rt_ID)) {
                    document.getElementById(rt_ID).innerHTML = (Math.round(d.rt * 100) / 100).toString();
                }
                if (document.getElementById(inte_ID)) {
                    document.getElementById(inte_ID).innerHTML = d.inteSum.toExponential(2);
                }
                if (document.getElementById(inte_ID)) {
                    document.getElementById(scanNum_ID).innerHTML = d.scanNum;
                }
                hoverLine.style("opacity", 1);
            }
            else if (i === inteRtArray.length && mouse_x - padding.left <= maxMouse + 1 && mouse_y < height - padding.bottom && mouse_y > padding.head) {
                let d = inteRtArray[i - 1];
                if (document.getElementById(rt_ID)) {
                    document.getElementById(rt_ID).innerHTML = (Math.round(d.rt * 100) / 100).toString();
                }
                if (document.getElementById(inte_ID)) {
                    document.getElementById(inte_ID).innerHTML = d.inteSum.toExponential(2);
                }
                if (document.getElementById(inte_ID)) {
                    document.getElementById(scanNum_ID).innerHTML = d.scanNum;
                }
                hoverLine.style("opacity", 1);
            }
            else {
                if (document.getElementById(rt_ID)) {
                    document.getElementById(rt_ID).innerHTML = "0";
                }
                if (document.getElementById(inte_ID)) {
                    document.getElementById(inte_ID).innerHTML = "0";
                }
                if (document.getElementById(inte_ID)) {
                    document.getElementById(scanNum_ID).innerHTML = "0";
                }
                hoverLine.style("opacity", 0);
            }
        }
        function hoverMouseOff() {
            hoverLine.style("opacity", 1e-6);
        }
    }
    moveLine(rt) {
        let newX = this.xScale_g_(rt) + this.padding_.left;
        this.fixedLine_g_.attr("x1", newX).attr("x2", newX);
        this.fixedLine_g_.style("opacity", 1);
    }
}
