class PrsmGraph {

  constructor(svgId, prsm){
    this.id = svgId;
    this.para = new PrsmPara();
    this.data = new PrsmData(); 
    this.data.initData(prsm, this.para); 
  }

  redraw = function(){
    //console.log(this.envList);
    this.data.updatePara(this.para);
    drawPrsm(this.id, this.para, this.data); 
  }

}
