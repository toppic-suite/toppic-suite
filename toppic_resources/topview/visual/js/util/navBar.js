//draw navigation bar while avoiding cross origin issue

var drawNav = function(){
	//depending on where it is being called... is it index html?
	let x = location.href;
	let n = x.lastIndexOf("/");
	let htmlName = x.substring(n + 1, x.length)

	let relPath;//the path to the other pages are different for index.html only
	
	if (htmlName == "index.html"){
		relPath = " ";
	}
	else{
		relPath = "../";
	}

	let navCode = '<nav class="navbar navbar-expand-lg fixed-top navbar-dark bg-dark"><div id="for-flex" style="display:flex; width: 100%"><div class=" navcontainer"><button class="navbar-toggler" type="button" data-toggle="collapse" data-target=".navbar-collapse" aria-controls="navbar-collapse" aria-expanded="false" aria-label="Toggle navigation"><span class="navbar-toggler-icon"></span></button><div class="collapse navbar-collapse"><a class="navbar-brand logo" href="' + relPath + 'index.html"><span class="headBlock"><h3><strong id="toppic_icon">T</strong>opView</h3></span></a><ul class="navbar-nav mr-auto mt-2 mt-lg-0"><li class="navtab">|</li><li class="nav-item"><a class="nav-link" href="' + relPath + 'visual/proteins.html"><mark>Identifications</mark></a></li><li class="navtab">|</li><li class="nav-item"><a class="nav-link" href="' + relPath + 'inspect/spectrum.html">Visual Inspection</a></li><li class="navtab">|</li><li class="nav-item"><a class="nav-link" href="' + relPath + 'spectrum/ms.html">Spectra</a></li></ul><ul class="navbar-nav navbar-right"><li class="nav-item"><a class="nav-link" href="https://soic.iupui.edu/people/xiaowen-liu/" target="#"">Contact</a></li><li class="navtab">|</li><li><a class="nav-link" href="' + relPath + 'inspect/help.html" target = "#" >Help</a></li></ul></div></div></div></nav>'
	
	let aboutBtn = '<li class="navtab">|</li><li class="nav-item"><a class="nav-link" href="/Home/About">About</a></li>'//this seems to be broken. Where is this file??

	 $(function(){
	 $("#nav-bar").html(navCode)
	 });
}();
