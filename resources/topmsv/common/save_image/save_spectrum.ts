class SaveSpectrum{
  private massGraphList_: SpectrumView[];
  private spectrumGraphList_: SpectrumView[];

  private graphModalObj_: SpectrumView | undefined;

  constructor(spectrumGraphList: SpectrumView[], massGraphList: SpectrumView[]){
    this.spectrumGraphList_ = spectrumGraphList;
    this.massGraphList_ = massGraphList;
  }

  addSpectrumModal = (): void => {
    let spectrumModal = 
    `<div class="modal" id="ms2_graph_popup_window" role="dialog">
      <div class="modal-dialog modal-sm" role="document">
        <div class="modal-content">
          <div class="modal-header ">
            <h3>Save MS/MS Spectrum</h3>
            <button type="button" class="close" data-dismiss="modal">&times;</button>
          </div>
          <div class="modal-body">
            <div> 
              <table class="table table-sm table_modal">
                <tbody>
                  <tr>
                    <div>
                      <td class="td_popup_variables" id="popup-env-btns" style="text-align:left">Show envelopes: &nbsp;&nbsp;&nbsp;
                        <input type="radio" name="show_envelopes" id="show_envelopes_yes"
                                            class = "show_envelopes" id="show_envelopes_no"
                                            checked>Yes</input>
                        <input type="radio"
                              name="show_envelopes" 
                              class = "show_envelopes">No</input>
                      </td>
                      <td class="td_popup_variables" id="popup-ion-btns">Show ions: &nbsp;&nbsp;&nbsp;
                        <input type="radio" name="show_ions" class = "show_ions" id="show_ions_yes" checked>Yes
                        <input type="radio" name="show_ions" class = "show_ions" id="show_ions_yes">No
                      </td>
                      <td class="td_popup_variables" id="show-error-plot-btns">Show error plot: &nbsp;&nbsp;&nbsp;
                        <input type="radio" name="show_error" class = "show_error" id="show_error_yes" checked>Yes
                        <input type="radio" name="show_error" class = "show_error" id="show_error_yes">No
                      </td>
                      <td class="td_popup_button" style="text-align:right">
                        <button type = "button" class="btn btn-primary btn-sm "  id ="ms2_popup_redraw_btn" >Redraw</button>
                    </td>
                    </div>
                  </tr>
                </tbody>
              </table>
            </div>
          </div>
          <div class="modal-body" id="ms2_graph_popup_svg_div">
            <svg id="popup_ms2_svg" style="background-color:#F8F8F8;"></svg>
          </div>
          <div class="modal-footer">
            <button type="button" class="btn btn-primary btn-sm custom " id =
            "download_ms2_graph_png_btn" >
              <i class="fa fa-download"></i><span>&nbsp;&nbsp;PNG</span>
            </button>
            <button type="button" class="btn btn-primary btn-sm custom " id =
            "download_ms2_graph_svg_btn" >
              <i class="fa fa-download"></i><span>&nbsp;&nbsp;SVG</span>
            </button>
            <button type="button" class="btn btn-primary btn-sm custom " data-dismiss="modal">Close</button>
          </div>
        </div>
      </div>
    </div>`
  
    //append to body
    $("body").append(spectrumModal);
  }
  createModalGraph = (): void => {
      // Get the Scan number of the Current spectrum graph showing on the screen
      let currentSpectrum: SpectrumView;
      let ms2Id: string = $('.ms2_graph_list.active')[0].id;
      let ms2Split: string[] = ms2Id.split("_");
      let type: string = ms2Split[ms2Split.length-2];
      let ms2Index: number = parseInt(ms2Split[ms2Split.length-1]);
      let svgId: string = "popup_ms2_svg";

      if (type == "graphlist") {
        currentSpectrum = this.spectrumGraphList_[ms2Index];
        this.graphModalObj_ = new SpectrumView(svgId,currentSpectrum.getPeakList());
        this.graphModalObj_.addRawSpectrumAnno(currentSpectrum.getEnvList(), currentSpectrum.getIonList());//ionList of currentSpectrum cannot be null
        this.graphModalObj_.getPara().updateScale(currentSpectrum.getPara().getWinMinMz(), currentSpectrum.getPara().getWinMaxMz(), currentSpectrum.getPara().getWinMaxInte());
        drawBasicSpectrum(this.graphModalObj_.getSvgId(), this.graphModalObj_.getPara(), this.graphModalObj_.getPeakList(), this.graphModalObj_.getIonList());
        $("#popup-env-btns").show();
        $("#show-error-plot-btns").hide();

        //set the buttons to default value
        $("#show_envelopes_yes").prop("checked", true);
        $("#show_envelopes_no").prop("checked", false);
 
        $("#show_ions_yes").prop("checked", true);
        $("#show_ions_no").prop("checked", false);
      }else{
        currentSpectrum = this.massGraphList_[ms2Index];
        this.graphModalObj_ = new SpectrumView(svgId,currentSpectrum.getPeakList());
		    this.graphModalObj_.addMonoMassSpectrumAnno(currentSpectrum.getIonList(), currentSpectrum.getProteoform()!, currentSpectrum.getNIon(), currentSpectrum.getCIon());
        this.graphModalObj_.getPara().setMonoMassGraph(true);
        this.graphModalObj_.getPara().updateScale(currentSpectrum.getPara().getWinMinMz(), currentSpectrum.getPara().getWinMaxMz(), currentSpectrum.getPara().getWinMaxInte());
        drawMonoMassSpectrum(this.graphModalObj_.getSvgId(), this.graphModalObj_.getPara(), this.graphModalObj_.getProteoform(), this.graphModalObj_.getNMassList(), 
        this.graphModalObj_.getCMassList(), (this.graphModalObj_.getIonList()) ? this.graphModalObj_.getIonList()! : []) //this.graphModalObj_.getIonList() returns matchedIon[] or null
        $("#popup-env-btns").hide();
        $("#show-error-plot-btns").show();

        //set the buttons to default value
        $("#show_error_yes").prop("checked", true);
        $("#show_error_no").prop("checked", false);

        $("#show_ions_yes").prop("checked", true);
        $("#show_ions_no").prop("checked", false);

        $("#show_envelopes_yes").prop("checked", true);
        $("#show_envelopes_no").prop("checked", false);
      }
  }
  initSpectrumModalEventHandler = () => {
    // MS2 graph save button 
    $("#ms2_graph_save_btn").click(() => {
      this.createModalGraph();
      this.drawModalGraph();
      //@ts-ignore
      $("#ms2_graph_popup_window").draggable({
        appendTo: "body"
      });
    })

    // Click of redraw invoke this and generate the popup spectrum
    $("#ms2_popup_redraw_btn").click(() => {
      if (!this.graphModalObj_) {
        console.error("ERROR: spectrum graph object is undefined");
        return;
      }
      let para = this.graphModalObj_.getPara(); 
      para.setShowEnvelopes_($("#show_envelopes_yes").prop("checked"));
      para.setShowIons($("#show_ions_yes").prop("checked"));
      para.setShowError($("#show_error_yes").prop("checked"));

      let svg: HTMLElement | null = document.getElementById("popup_ms2_svg");
      if (svg) {
        if (!para.getShowError() && para.getIsMonoMassGraph()) {
          svg.setAttribute("viewBox", '0 0 ' + para.getSVGWidth().toString() + ' ' + (para.getSVGHeight() - 40).toString());
        }
        else {
          svg.setAttribute("viewBox", '0 0 ' + para.getSVGWidth().toString() + ' ' + para.getSVGHeight().toString());
        }
      }
      this.drawModalGraph();
    })

    // Coordinates at which a pop window to be launched to give name for the image to be downloaded
    // All the actions for different image buttons on different spectrums
    d3.select("#download_ms2_graph_png_btn").on("click", () => {
      let x = d3.event.pageX;
      let y = d3.event.pageY;
      // function inside prsmtohtml.js to pop up a window
      popupNameWindow("png", "popup_ms2_svg",x,y)
    })
    d3.select("#download_ms2_graph_svg_btn").on("click", () => {
      let x = d3.event.pageX;
      let y = d3.event.pageY;
      popupNameWindow("svg", "popup_ms2_svg",x,y)
    })
  }
  drawModalGraph = () => {
		this.graphModalObj_!.redraw();
  }
  main = () => {
    this.addSpectrumModal();
    this.initSpectrumModalEventHandler();
  }
}

