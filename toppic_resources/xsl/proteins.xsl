<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

<xsl:template match="proteins">
<html>
<title><xsl:value-of select="count(protein)"/> proteins were identified</title>
<link rel="stylesheet" type="text/css" href="./resources/media/css/jquery.dataTables.css"></link>
<link rel="stylesheet" type="text/css" href="./resources/bootstrap.min.css"></link>
<link rel="stylesheet" type="text/css" href="./resources/prsm.css"></link>
<script type="text/javascript" language="javascript" src="./resources/media/js/jquery.js"></script>
<script type="text/javascript" language="javascript" src="./resources/media/js/jquery.dataTables.js"></script>
<style type="text/css">
	@font-face {
		font-family: "FreeMono";
		src: url("./resources/FreeMono.ttf")
	}
</style>
<body>
<div class="container">
<h2><xsl:value-of select="count(protein)"/> proteins were identified.</h2>
<br/>
<xsl:apply-templates select="protein"></xsl:apply-templates>
</div>
</body>
</html>
</xsl:template>

<xsl:template match="protein">
<div id="p{sequence_id}">
  <p style="font-size:16px;">
    <a href="proteins/protein{sequence_id}.html"><xsl:value-of select="sequence_name"/><xsl:text> </xsl:text>
		<xsl:value-of select="sequence_description"/>
  </a></p>
<xsl:apply-templates select="compatible_proteoform/prsm"></xsl:apply-templates>
<br/>
</div>
</xsl:template>



<xsl:template match="compatible_proteoform/prsm">
<xsl:if test="position()=1">
<xsl:choose>
	<xsl:when test="count(../prsm) = 1">
		<p style="font-size:16px;">There is only <a href="prsms/prsm{prsm_id}.html">1 PrSM</a> with an E-value <xsl:value-of select="e_value"/>.</p>
	</xsl:when>
    <xsl:otherwise>
		<p style="font-size:16px;">The <a href="prsms/prsm{prsm_id}.html">best PrSM</a> has an E-value <xsl:value-of select="e_value"/>. There
        <xsl:choose>
			<xsl:when test="count(../../compatible_proteoform) = 1">is <a href="species/species{annotated_protein/species_id}.html">1 proteoform</a>.  </xsl:when>
            <xsl:otherwise>are <xsl:value-of select="count(../../compatible_proteoform)"/> proteoforms.
            </xsl:otherwise>
        </xsl:choose></p>
	</xsl:otherwise>
</xsl:choose>
</xsl:if>
</xsl:template>
</xsl:stylesheet>
