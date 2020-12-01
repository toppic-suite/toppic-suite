class SavePrsm{
  prsmGraphObj;
  prsmModalGraphObj;

  constructor(prsmGraphObj){
    this.prsmGraphObj = prsmGraphObj;
  }

  addPrsmModal = () => {
    let prsmModal = `<div class="modal" id="save_prsm_popup_window" role="dialog">
      <div class="modal-dialog modal-lg" role="document">
        <div class="modal-content" id ="save-image-modal">
          <div class="modal-header ">
            <!-- Your first column here -->
            <h5 class="modal-title ml-auto">Save PrSM Image</h5>
            <button type="button" class="close" data-dismiss="modal">&times;</button>
          </div>
          <div classs="modal-header modal2-header pb-0">
            <table class="table table-sm table_modal">
              <tr>
                <td>Amino acids per row:</td>
                <td><input id="row-size" type="text" size="7" value=30></td>
                <td>Letter width:</td>
                <td><input id="letter-width" type="text" size="7" value=28></td>
                <td>Row height:</td>
                <td><input id="row-height" type="text" size="7" value=40></td>
              </tr>
              <tr>
                <td>Gap between blocks:</td>
                <td><input id="block-width" type="text" size="7" value=20></td>
                <td>Gap between num and seq:</td>
                <td><input id="num-width" type="text" size="7" value="20"></td>
                <td>Show numbers:</td>
                <td>
                  <input type="radio" name="show-num" id = "show-num" checked>Yes
                  <input type="radio" name="show-num" id = "show-num" >No
                </td>
              </tr>
              <tr>
                <td>Show skipping information:</td>
                <td>
                  <input type="radio" name="show-skipped-lines" id=
                  "show-skipped-lines" checked>Yes
                  <input type="radio" name="show-skipped-lines" id= "show-skipped-lines">No
                </td>
                <td></td>
                <td>
                  <button type = "button" class="btn btn-primary btn-sm "  id =
                  "prsm_graph_redraw_btn" style="width:77%" >Redraw</button>
                </td>
              </tr>
            </table>
          </div>
          <div class="modal-body">
            <svg id = "prsm_popup_svg" class="prsm_popup_svg" style="background-color:white"></svg>
          </div>
          <div class="modal-footer">
            <button class="btn btn-primary btn-sm custom "  
                    id="prsm_popup_help_btn" data-toggle="modal"
                                            data-target="#prsm_help_popup_window">Help</button>
            <button type="button" class="btn btn-primary btn-sm custom " id = "prsm_popup_png_btn" >
            <i class="fa fa-download"></i><span>&nbsp;&nbsp;PNG</span>
            </button>
            <button type="button" class="btn btn-primary btn-sm custom " id = "prsm_popup_svg_btn" >
            <i class="fa fa-download"></i><span>&nbsp;&nbsp;SVG</span>
            </button>
            <button type="button" class="btn btn-primary btn-sm custom " data-dismiss="modal">Close</button>
          </div>
        </div>
      </div>
    </div>
    <div class="modal" id="prsm_help_popup_window" role="dialog">
    <div class="modal-dialog modal-sm" role="document">
      <div class="modal-content help-window">
        <div class="modal-header ">
          <h5 class="modal-title ml-auto">Help</h5>
          <button type="button" class="close" id = "prsm-help-window-close-btn" data-dismiss="modal">&times;</button>
        </div>
        <div class="modal-body" >
          <ul>
            <li>Show skipping information: When the protein sequence is not a whole
              sequence, the option will add additional lines with information about
              the number of amino acids not included in the image.</li>
            <!--<li>White b/g: it will add a white background to the image.</li> -->
          </ul> 
        </div> 
      </div>
    </div>
  </div>
  `
  
    //append to body
    $("body").append(prsmModal);
  }
  createPrsmModalGraph = () => {
    this.prsmModalGraphObj = new PrsmGraph("prsm_popup_svg",null,this.prsmGraphObj.data);

  }
  initPrsmModalEventHandler = () => {
    // Save PrSM popup window
    d3.select('#save_prsm_btn').on("click", () => {
      this.prsmModalGraphObj.redraw();
      //	set the dimensions of popup svg to default values
      let para = this.prsmModalGraphObj.para;
      document.getElementById("row-size").value = para.rowLength ;
      document.getElementById("letter-width").value = para.letterWidth ;
      document.getElementById("row-height").value = para.rowHeight ;
      document.getElementById("block-width").value = para.gapWidth ;
      document.getElementById("num-width").value = para.numericalWidth ;
      document.getElementsByName("show-num")[0].checked = para.showNum ;
      document.getElementsByName("show-skipped-lines")[0].checked = para.showSkippedLines ;

      // Allows to drag the pop up windows
      $("#save_prsm_popup_window").draggable({
        appendTo: "body"
      });
    });

    d3.select('#prsm_graph_redraw_btn').on("click", () => {
      let para = this.prsmModalGraphObj.para;
      para.rowLength = parseInt(document.getElementById("row-size").value);
      para.letterWidth = parseInt(document.getElementById("letter-width").value) ;
      para.rowHeight = parseInt(document.getElementById("row-height").value) ;
      para.gapWidth = parseInt(document.getElementById("block-width").value) ;
      para.numericalWidth = parseInt(document.getElementById("num-width").value) ;

      //	Check whether show numbers is checked
      if(document.getElementsByName("show-num")[0].checked) {
        para.setShowNum(true);
      }
      else {
        para.setShowNum(false);
      }
      //	Check to show skipped lines 
      if(document.getElementsByName("show-skipped-lines")[0].checked) {
        para.showSkippedLines = true ;
      }
      else {
        para.showSkippedLines = false ;
      }
      this.prsmModalGraphObj.redraw();
    });	

    $("#prsm_popup_help_btn").on("click",() => {
      $("#prsm_help_popup_window").draggable({
        appendTo: "body"
      });
      //grey out other parts of screen
      $(".modal-backdrop").css('z-index',3000);      
      $("#prsm_help_popup_window").css('z-index',3001);     
    });

    $("#prsm-help-window-close-btn").on("click",() => {
      $(".modal-backdrop").css('z-index',1040);  
    })
      //	Download the svg as ".svg" image
    d3.select('#prsm_popup_svg_btn').on("click",() => {
      let svgId = "prsm_popup_svg" ;
      let x = d3.event.pageX;
      let y = d3.event.pageY;
      popupNameWindow("svg",svgId,x,y);
    });
    //	Download svg as PNG Image
    d3.select('#prsm_popup_png_btn').on("click", () => {
      let svgId = "prsm_popup_svg" ;
      let x = d3.event.pageX;
      let y = d3.event.pageY;
      popupNameWindow("png",svgId,x,y);
    });
  }
  drawPrsmModalGraph = () => {
		this.prsmModalGraphObj.para.rowLength = 40;
		this.prsmModalGraphObj.para.letterWidth = 28;
		this.prsmModalGraphObj.redraw();
  }
  main = () => {
    this.addPrsmModal();
    this.createPrsmModalGraph();
    this.initPrsmModalEventHandler();
    this.drawPrsmModalGraph();
  }
}

