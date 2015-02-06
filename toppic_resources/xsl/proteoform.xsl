<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

<xsl:template match="compatible_proteoform">
<html>
  <title>Proteoform #<xsl:value-of select="proteoform_id"/> from <xsl:value-of select="sequence_name"/><xsl:text> </xsl:text>
		<xsl:value-of select="sequence_description"/>
  
  
  </title>
<link rel="stylesheet" type="text/css" href="../resources/media/css/jquery.dataTables.css"></link>
<link rel="stylesheet" type="text/css" href="../resources/bootstrap.min.css"></link>
<link rel="stylesheet" type="text/css" href="../resources/prsm.css"></link>
<script type="text/javascript" language="javascript" src="../resources/media/js/jquery.js"></script>
<script type="text/javascript" language="javascript" src="../resources/media/js/jquery.dataTables.js"></script>
<style type="text/css">
	@font-face {
		font-family: "FreeMono";
		src: url("../resources/FreeMono.ttf")
	}
</style>
<body>
<div class="container">
<xsl:call-template name="navigation"/>
<h2>Proteoform #<xsl:value-of select="proteoform_id"/></h2>
<p style="font-size:16px;"><xsl:value-of select="prsm_number"/> PrSM<xsl:if test="prsm_number > 1">s</xsl:if> for this proteoform </p>
<table class="table table-hover table-bordered table-striped table-nonfluid">
<tr>
<th>Scan</th>
<th>E-value</th>
<th># all peaks</th>
<th># matched peaks</th>
<th># matched fragment ions</th>
<th>Link</th>
</tr>
<xsl:apply-templates select="prsm">
<xsl:sort select="e_value" order="ascending" data-type="number"/>
</xsl:apply-templates>
</table>
<xsl:call-template name="navigation"/>
</div>
</body>
</html>

</xsl:template>

<xsl:template match="prsm">
<tr>
<td align="center"> <xsl:value-of select="ms/ms_header/scans"/> </td>
<td align="center"> <xsl:value-of select="e_value"/> </td>
<td align="center"><xsl:value-of select="count(ms/peaks/peak)"/></td>
<td align="center"><xsl:value-of select="matched_peak_number"/></td>
<td align="center"><xsl:value-of select="matched_fragment_number"/></td>
<td align="center"><a onfocus="cur = {position()-1};selectRow();" href="../prsms/prsm{prsm_id}.html" id="link{position()}">See PrSM>></a></td>
</tr>
</xsl:template>


<xsl:template name="navigation">
<p style="font-size:15px;">
	<a href="../proteins.html">All proteins</a> /
  <a href="../proteins/protein{sequence_id}.html"><xsl:value-of select="sequence_name"/><xsl:text> </xsl:text>
		<xsl:value-of select="sequence_description"/>
  </a>
</p>
</xsl:template>

</xsl:stylesheet>
