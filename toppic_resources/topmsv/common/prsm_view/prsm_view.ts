class PrsmView {
  private id_: string;
  private para_: PrsmPara;
  private data_: PrsmViewData | null;

  constructor(svgId: string, prsmObj: Prsm | null, graphData: PrsmViewData | null = null, allowMod: boolean = false) {
    this.id_ = svgId;
    this.para_ = new PrsmPara(allowMod);
    this.data_ = graphData;

    if (!this.data_) {
      this.data_ = new PrsmViewData();
      if (!prsmObj) {
        console.error("ERROR: PrsmObj is empty!");
        return;
      }
      this.data_.initData(prsmObj, this.para_);   
    }
  }
  getData(): PrsmViewData | null {
    return this.data_;
  }
  getPara(): PrsmPara {
    return this.para_;
  }
  redraw(): void {
    if (!this.data_) {
      console.error("PrsmView data is empty");
      return;
    }
    this.data_.updatePara(this.para_);
    drawPrsm(this.id_, this.para_, this.data_); 
  }
}
