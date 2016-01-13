<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
<xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

<xsl:include href="prsm_anno.xsl"/>

<xsl:template match="prsm">
<html>
<head>
<title>Protein-Spectrum-Match for Spectrum #<xsl:value-of select="ms/ms_header/ids"/></title>
<link rel="stylesheet" type="text/css" href="../resources/media/css/jquery.dataTables.css"></link>
<link rel="stylesheet" type="text/css" href="../resources/bootstrap.min.css"></link>
<link rel="stylesheet" type="text/css" href="../resources/prsm.css"></link>
<script type="text/javascript" language="javascript" src="../resources/media/js/jquery.js"></script>
<script type="text/javascript" language="javascript" src="../resources/media/js/jquery.dataTables.js"></script>
<script type="text/javascript">
$(document).ready(function() {
	$('#spectrum').dataTable( {
		"scrollY":        "400px",
		"scrollCollapse": true,
		"paging":         false,
		"order": [[ 0, "asc" ]],
		"bSortClasses": false
	} );
} );
		
var peaksCount = <xsl:value-of select="count(ms/peaks/peak)"/>;

function showMatchedPeaks() {
  var elems = document.getElementsByClassName("matched_peak");
  for(var i = 0; elems.length > i; i++) {
    elems[i].style.display = "";
  }
  elems = document.getElementsByClassName("unmatched_peak");
  for(var i = 0; elems.length > i; i++) {
    elems[i].style.display = "none";
  }
  $('div.dataTables_scrollBody').height(400);
}

function showNotMatchedPeaks() {
  var elems = document.getElementsByClassName("matched_peak");
  for(var i = 0; elems.length > i; i++) {
    elems[i].style.display = "none";
  }
  elems = document.getElementsByClassName("unmatched_peak");
  for(var i = 0; elems.length > i; i++) {
    elems[i].style.display = "";
  }
  $('div.dataTables_scrollBody').height(400);
}

function showAllPeaks() {
  var elems = document.getElementsByClassName('matched_peak');
  for(var i = 0; elems.length > i; i++) {
    elems[i].style.display = '';
  }
  elems = document.getElementsByClassName('unmatched_peak');
  for(var i = 0; elems.length > i; i++) {
    elems[i].style.display = '';
  }
  $('div.dataTables_scrollBody').height(400);
}

function showIonPeaks(ids) {
  var elems = document.getElementsByClassName('matched_peak');
  for(var i = 0; elems.length > i; i++) {
    elems[i].style.display = 'none';
  }
  elems = document.getElementsByClassName('unmatched_peak');
  for(var i = 0; elems.length > i; i++) {
    elems[i].style.display = 'none';
  }
			 
	var s = ids.split(",");
		   
	for (i = 0; s.length>i; i++) {
    elems = document.getElementsByName(s[i]);
    for(var j = 0; elems.length > j; j++) {
		  elems[j].style.display  =  "";
      elems[j].style.background  =  "#BEECFF";
    }
	}
	$('div.dataTables_scrollBody').height(400);
}
-->
		
</script>

<style type="text/css">
	@font-face {
		font-family: "FreeMono";
		src: url("../resources/FreeMono.ttf")
	}
</style>
</head>
<body>
<div class="container">
	<xsl:call-template name="navigation"/>
	<h2>Protein-Spectrum-Match #<xsl:value-of select="prsm_id"/> for Spectrum #<xsl:value-of select="ms/ms_header/ids"/></h2>
	<table class="table table-borderless" style="font-size:15px">
	<tr>
		<td>PrSM ID: </td><td> <xsl:value-of select="prsm_id"/> </td>
		<td>Scan(s): </td><td> <xsl:value-of select="ms/ms_header/scans"/> </td>
		<td>Precursor charge:</td><td> <xsl:value-of select="ms/ms_header/precursor_charge"/> </td>
	</tr>
	<tr>
		<td>Precursor m/z: </td><td> <xsl:value-of select="ms/ms_header/precursor_mz"/> </td>
		<td>Precursor mass:</td><td> <xsl:value-of select="ms/ms_header/precursor_mono_mass"/> </td>
		<td>Proteoform mass:</td><td> <xsl:value-of select="annotated_protein/proteoform_mass"/> </td>
	</tr>
	<tr>
		<td># matched peaks:</td><td> <xsl:value-of select="matched_peak_number"/> </td>
		<td># matched fragment ions:</td><td> <xsl:value-of select="matched_fragment_number"/> </td>
		<td># unexpected PTMs:</td><td> <xsl:value-of select="annotated_protein/unexpected_change_number"/> </td>
	</tr>
	<tr>
		<td>E-value:</td><td> <xsl:value-of select="e_value"/> </td>
		<td>P-value:</td><td> <xsl:value-of select="p_value"/> </td>
		<td>Spectral FDR:</td><td> <xsl:value-of select="fdr"/> </td>
	</tr>
	</table>
	<br/>
	<xsl:apply-templates select="annotated_protein/annotation" mode="prsm"/> 
	<br/>
	<div style="font-size:16px;">
		<xsl:text>&#160;&#160;&#160;&#160;&#160;</xsl:text>
		<a href="#" onclick="showAllPeaks();"><xsl:text>All peaks (</xsl:text><xsl:value-of select="count(ms/peaks/peak)"/><xsl:text>)</xsl:text></a>
		<xsl:text>&#160;&#160;</xsl:text>
		<a href="#" onclick="showMatchedPeaks();"><xsl:text>Matched peaks (</xsl:text><xsl:value-of select="count(ms/peaks/peak[count(matched_ions/matched_ion)>0])"/><xsl:text>)</xsl:text></a>
		<xsl:text>&#160;&#160;</xsl:text>
		<a href="#" onclick="showNotMatchedPeaks();"><xsl:text>Not matched peaks (</xsl:text><xsl:value-of select="count(ms/peaks/peak[count(matched_ions/matched_ion)=0])"/><xsl:text>)</xsl:text></a>
	</div>
	<br/>
	<table id="spectrum" class="display" cellspacing="0" width="100%">
	<thead>
	<xsl:call-template name="peaks-header"/>
	</thead>
	<tbody>
		<xsl:apply-templates select="ms/peaks/peak"/>
	</tbody>
	</table>
	<br/>
	<xsl:call-template name="navigation"/>
</div>	
</body>
</html>

</xsl:template>

<xsl:template match="ms/peaks/peak/matched_ions/matched_ion" mode="peakMatch">
	if (peakId == <xsl:value-of select="concat(../../id, ion_type)"/>)
		return true;
</xsl:template>

<xsl:template name="navigation">
	<p style="font-size:15px;">
		<a href="../proteins.html">All proteins</a> /
		<a href="../proteins/protein{annotated_protein/sequence_id}.html">
		<xsl:value-of select="annotated_protein/sequence_name"/><xsl:text> </xsl:text>
		<xsl:value-of select="annotated_protein/sequence_description"/></a> /
		<a href="../proteoforms/proteoform{annotated_protein/proteoform_id}.html">
        Proteoform #<xsl:value-of select="annotated_protein/proteoform_id"/></a>
	</p>
</xsl:template>

<xsl:template match="annotated_protein/annotation" mode="prsm">
	<div id="alignment" style="font-family: 'FreeMono'; font-size:16px;line-height:2.5;background-color:#FFF;">
		<xsl:variable name="row_residue_num" select="30"/>
			<xsl:variable name="residue_num" select="protein_length"/>
			<xsl:variable name="table_first_row">
			<xsl:variable name="tmp" select="floor(first_residue_position div $row_residue_num)"/>
			<xsl:choose>
			<xsl:when test="$tmp &gt; 0">
				<xsl:value-of select="$tmp - 1"/>
			</xsl:when>
			<xsl:otherwise> 
				<xsl:value-of select="$tmp"/>
			</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>

		<xsl:variable name="table_total_row_num" select="ceiling($residue_num div $row_residue_num)"/>
			<xsl:variable name="table_last_row">
			<xsl:variable name="tmp" select="floor(last_residue_position div $row_residue_num)"/>
			<xsl:choose>
				<xsl:when test="$table_total_row_num &gt; $tmp + 1">
					<xsl:value-of select="$tmp + 1"/>
				</xsl:when>
				<xsl:otherwise> 
					<xsl:value-of select="$tmp"/>
				</xsl:otherwise>
			</xsl:choose>
		</xsl:variable>

		<table  class="table table-borderless table-condensed nopadding table-bold" >
        <xsl:if test="$table_first_row &gt; 0">
			<tr>
				<td colspan="2"></td>
				<td colspan="62" align="center"><xsl:text>... </xsl:text><xsl:value-of select="$table_first_row * 30"/><xsl:text> amino acid residues are skipped at the N-terminus ...</xsl:text></td>
				<td colspan="3"></td>
			</tr>
        </xsl:if>
        <xsl:call-template name="add_table_row">
			<xsl:with-param name="i" select="$table_first_row"/>
			<xsl:with-param name="table_last_row" select="$table_last_row"/>
        </xsl:call-template>

        <xsl:if test="$table_last_row * 30 + 30 &lt; $residue_num">
			<tr>
				<td>&#160; </td>
			</tr>
			<tr>
				<td colspan="2"></td>
				<td colspan="62" align="center"><xsl:text>... </xsl:text><xsl:value-of select="$residue_num - $table_last_row * 30 - 30"/><xsl:text> amino acid residues are skipped at the C-terminus ...</xsl:text></td>
				<td colspan="3"></td>
			</tr>
        </xsl:if>
		</table>         
	</div>

	<div style="font-size:16px;">
		<xsl:variable name="fix_ptm_num" select="count(expected_change)"/>
		<xsl:if test="$fix_ptm_num &gt; 0">
			<br/>
			<xsl:text>&#160;&#160;&#160;&#160;&#160;Fixed PTMs: </xsl:text> 
			<xsl:apply-templates select="expected_change"/>
			<br/>
		</xsl:if>
		<xsl:variable name="variable_ptm_num" select="count(unexpected_change)"/>
		<xsl:if test="$variable_ptm_num &gt; 0">
			<br/>
			<xsl:text>&#160;&#160;&#160;&#160;&#160;Variable PTMs: </xsl:text> 
			<xsl:apply-templates select="unexpected_change"/>
			<br/>
		</xsl:if>
	</div>
</xsl:template>

<xsl:template match="cleavage">
	<xsl:text disable-output-escaping="yes"><![CDATA[ style="]]></xsl:text>
	<xsl:if test="cleavage_type = 'seq_start'">
		<xsl:text>font-weight:bold;color:red;</xsl:text>
    </xsl:if>
    <xsl:if test="cleavage_type = 'seq_end'">
		<xsl:text>font-weight:bold;color:red;</xsl:text>
    </xsl:if>
    <xsl:if test="is_unexpected_change = '1'">
		<xsl:if test="unexpected_change_color = 0">
			<xsl:text>background:#DFFFFF;</xsl:text>
        </xsl:if>
        <xsl:if test="unexpected_change_color = 1">
			<xsl:text>background:#CECEF6;</xsl:text>
        </xsl:if>
	</xsl:if>
    <xsl:text disable-output-escaping="yes"><![CDATA[">]]></xsl:text>
	<xsl:choose>
		<xsl:when test="exist_n_ion = '0' and exist_c_ion = '0'">
			<xsl:choose>
				<xsl:when test="cleavage_type = 'seq_start'">
					<xsl:text>]</xsl:text>
				</xsl:when>
				<xsl:when test="cleavage_type = 'seq_end'">
					<xsl:text>[</xsl:text>
				</xsl:when>
            	<xsl:otherwise>
					<xsl:text>&#160;</xsl:text>
				</xsl:otherwise>
			</xsl:choose>
        </xsl:when>
        <xsl:otherwise>
			<a style="text-decoration:none" href="#">
            <xsl:attribute name="title"> 
				<xsl:apply-templates select="matched_peaks" mode="title"/>
            </xsl:attribute>
            <xsl:attribute name="onclick"> 
				<xsl:text>showIonPeaks('</xsl:text>
				<xsl:apply-templates select="matched_peaks" mode="prsm"/>
				<xsl:text>')</xsl:text>
            </xsl:attribute>
            <xsl:choose>
				<xsl:when test="exist_n_ion = '1' and exist_c_ion = '0'">
					<xsl:text disable-output-escaping="yes">&amp;#x23AB;</xsl:text>
				</xsl:when>
				<xsl:when test="exist_n_ion = '0' and exist_c_ion = '1'">
					<xsl:text disable-output-escaping="yes">&amp;#x23A9;</xsl:text>
				</xsl:when>
				<xsl:when test="exist_n_ion = '1' and exist_c_ion = '1'">
					<xsl:text disable-output-escaping="yes">&amp;#x23B1;</xsl:text>
				</xsl:when>
            </xsl:choose>
			</a>
		</xsl:otherwise>
	</xsl:choose>
	<xsl:text disable-output-escaping="yes"><![CDATA[]]></xsl:text>
</xsl:template>

<xsl:template name="peaks-header">
<tr>
<th width="25">Spec</th>
<th width="25">Peak</th>
<th width="90">Mono mass</th>
<th width="90">Mono m/z</th>
<th width="80" style="vertical-align:middle">Intensity</th>
<th width="75" style="vertical-align:middle">Charge</th>
<th width="103">Theoretical mass</th>
<th width="50" style="vertical-align:middle">Ion</th>
<th width="70" style="vertical-align:middle">Pos</th>
<th width="95">Mass error</th>
<th width="80">PPM error</th>
</tr>
</xsl:template>

<xsl:template match="peak">
<xsl:variable name="peakID" select="peak_id" as="xs:integer"/>
<xsl:if test="count(matched_ions/matched_ion) = 0">
<tr id="spec{spec_id}peak{peak_id}" class="unmatched_peak">
<td align="center"><xsl:value-of select="spec_id"/></td>
<td align="center"><xsl:value-of select="$peakID + 1"/></td>
<td align="center"><xsl:value-of select="monoisotopic_mass"/></td>
<td align="center"><xsl:value-of select="monoisotopic_mz"/></td>
<td align="center"><xsl:value-of select="intensity"/></td>
<td align="center"><xsl:value-of select="charge"/></td>
<td></td><td></td><td></td><td></td><td></td></tr>
</xsl:if>
<xsl:if test="not(count(matched_ions/matched_ion) = 0)">
<xsl:apply-templates select="matched_ions/matched_ion"/>
</xsl:if>
</xsl:template>

<xsl:template match="matched_ion">
<xsl:variable name="peakID" select="../../peak_id" as="xs:integer"/>
<tr id="spec{../../spec_id}peak{../../peak_id}{ion_type}" class="matched_peak" name="{ion_position}">
<td align="center"><xsl:value-of select="../../spec_id"/></td>
<td align="center"><xsl:value-of select="$peakID + 1"/></td>
<td align="center"><xsl:value-of select="../../monoisotopic_mass"/></td>
<td align="center"><xsl:value-of select="../../monoisotopic_mz"/></td>
<td align="center"><xsl:value-of select="../../intensity"/></td>
<td align="center"><xsl:value-of select="../../charge"/></td>	
<td align="center"><xsl:value-of select="theoretical_mass"/></td>
<td align="center"><xsl:value-of select="ion_type"/><xsl:value-of select="ion_display_position"/></td>
<td align="center"><xsl:value-of select="ion_left_position"/></td>
<td align="center"><xsl:value-of select="mass_error"/></td>
<td align="center"><xsl:value-of select="ppm"/></td>
</tr>
</xsl:template>

</xsl:stylesheet>

