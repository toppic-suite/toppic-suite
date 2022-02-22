/**
 * Starting point of building spectrums page.
 * Gets the data of all spectrum and shows on the html
 * @param {String} folderName - Provides the path to the data folder.
 * Provides path which helps to navigate to protein and prsm pages.
 */
/**
 * generate a search box for ms1 and ms2 spectra
 */
function searchBoxEventListner(): void {
  //find scan using ID in the page
  let scanNum: HTMLInputElement = <HTMLInputElement> document.getElementById("ms2_scanid");
  if (scanNum) {
    let searchResult: HTMLElement | null = document.getElementById("ms2_"+scanNum.value);
    if (searchResult) {
      let navbar: HTMLElement | null = document.getElementById("nav-div");  
      if (navbar) {
        let topOfElement = searchResult.offsetTop - navbar.clientHeight;
        window.scroll({top: topOfElement});
      }
    }
  }
}
function addSearchBoxEvent(): void {
  let submitBtn: HTMLElement | null = document.getElementById("ms2_submit");
  if (submitBtn) {
      submitBtn.addEventListener("click", searchBoxEventListner);
  }
}

/**
 * convert the json prsm data into HTML and create links for each spectra to navigate
 * @param {object} prsm - Contains data of a single prsm
 * @param {String} folderName - Provides path to build navigation links
 */
 function prsmToHtml(prsm: any, folderName: string): HTMLLIElement {
  let protein = prsm.annotated_protein;
  let div: HTMLLIElement = document.createElement('li');
  let id: string = "p" + prsm.prsm_id;
  div.setAttribute("id", id);

  let childDiv: HTMLDivElement = document.createElement('div');
  childDiv.style.display = "flex";
  let ms1Spec: HTMLHeadingElement = document.createElement('p');
  ms1Spec.innerHTML = "MS1 Scan #" + prsm.ms.ms_header.ms1_scans;
  ms1Spec.classList.add("scan-title");
  ms1Spec.id = "ms1_" + prsm.ms.ms_header.ms1_scans;

  let ms2Spec: HTMLHeadingElement = document.createElement('p');
  ms2Spec.innerHTML = "MS2 Scan #" + prsm.ms.ms_header.scans;
  ms2Spec.classList.add("scan-title");
  ms2Spec.id = "ms2_" + prsm.ms.ms_header.scans;

  let protName: HTMLHeadingElement = document.createElement('p');
  protName.innerHTML = protein.sequence_name;
  protName.classList.add("scan-title");

  childDiv.appendChild(ms2Spec);
  childDiv.appendChild(ms1Spec);
  childDiv.appendChild(protName);
  div.appendChild(childDiv);

 // let residues: {"position": string, "acid": string}[] = protein.annotation.residue;
  
  /*for (let i = protein.annotation.first_residue_position; i < protein.annotation.last_residue_position; i++) {
    sequence = sequence.concat(residues[i].acid);
  }*/

  let sequence: string = protein.annotation.annotated_seq;

  let p: HTMLParagraphElement = document.createElement('p');
  let a: HTMLAnchorElement = document.createElement('a');
  a.href = "prsm" + ".html" + "?folder=" + folderName + "&protein=" + prsm.prsm_id;
  //a.innerHTML = protein.sequence_name + " " + "first residue = " + protein.annotation.first_residue_position + " " + "last residue = " + protein.annotation.last_residue_position;
  a.innerHTML = sequence;
  p.appendChild(a);
  div.appendChild(p);
  return div;
}
 function allSpectrum(folderName: string) {
  //@ts-ignore
  let l_prsms: any = prsm_data.prsms.prsm;
  let count: number = 1;
  // get the div container 
  let div: Element = document.getElementsByClassName("container")[0];

  //switch between spectrum identification and protein identification
  let p: HTMLHeadingElement = document.createElement('p');//need to assign class 
  let x: string = location.href;
  let l_split: string = x.split(/[?#]+/)[0];
  let idx = l_split.lastIndexOf('\\');
  if (idx < 0) {
    idx = l_split.lastIndexOf('/');
  }
  let newAddress: string = l_split.slice(0,idx + 1) + "proteins.html?folder=" + folderName;
  p.innerHTML = '<a href=' + newAddress + '>Switch to Protein Identification</a>';
  //div.prepend(p);

  let titleDiv: HTMLElement | null = document.getElementById('title-div');
  let h2: HTMLElement | null = document.getElementById('id-title');

  if (!h2) {
      console.error("ERROR: Identification title header doesn't exist");
      return;
  }
  if (!titleDiv) {
    console.error("ERROR: Title div doesn't exist");
    return;
  }
  // Check to see if protein variable inside l_proteins is an array.
  // Checks for multiple proteins
  if (Array.isArray(l_prsms)) {
    count = l_prsms.length;
    document.title = count + " spectra are identified";
    h2.innerHTML = count + " spectra are identified.";
  }
  else {
    document.title = count + " spectrum is identified";
    h2.innerHTML = count + " spectrum is identified.";
  }
  let br: HTMLBRElement = document.createElement('br');
  // create header with protein count 
  div.appendChild(br);
  // Creating ordered list
  let ol: HTMLOListElement = document.createElement('ol');
  // get the best prsm for each protein and form unique links for all the proteins
  // Check to see if protein variable inside l_proteins is an array.
  if (Array.isArray(l_prsms)) {
    l_prsms.forEach(function (prsm: any, index: number) {
      let div_temp: HTMLLIElement = prsmToHtml(prsm, folderName);
      let br1: HTMLBRElement = document.createElement('br');
      div_temp.appendChild(br1);
      ol.appendChild(div_temp);
    });
  }
  else {
    let prsm: any = l_prsms;
    let div_temp: HTMLLIElement = proteinToHtml(prsm, folderName);
    let br1: HTMLBRElement = document.createElement('br');
    div_temp.appendChild(p);
    div_temp.appendChild(br1);
    ol.appendChild(div_temp);
  }
  div.appendChild(ol);

  addSearchBoxEvent();
}
