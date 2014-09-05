<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

  <xsl:variable name="alignWidth">40</xsl:variable>

  <xsl:template match="prsm" mode="basic">
    <p style="font-family:monospace;font-size:17;line-height: 20pt; ">
      <nobr>
        <xsl:apply-templates select="align/align-start"/>
        <xsl:apply-templates select="align/wrong-shift-start"/>
        <xsl:apply-templates select="align/residue"/>
        <xsl:apply-templates select="align/wrong-shift-finish"/>
        <xsl:apply-templates select="align/align-finish"/><br/>
      </nobr>
    </p>
  </xsl:template>

  <xsl:template match="wrong-shift-start">
    <xsl:call-template name="wrong-shift"/>
  </xsl:template>

  <xsl:template match="wrong-shift-finish">
    <xsl:call-template name="wrong-shift"/>
  </xsl:template>

  <xsl:template name="wrong-shift">
    <span  style="position: relative;"><span style="position: absolute; top:-10pt; font-size: 8pt; color:red; text-decoration:none;">
        <xsl:value-of select="format-number(text(), '0.00')"/>
    </span><a href="#"><font color="red">&#160;</font></a></span>
  </xsl:template>

  <xsl:template match="residue">
    <xsl:apply-templates select="shift-start"/>

    <xsl:value-of select="acid"/>

    <xsl:apply-templates select="shift-finish"/>
    <xsl:if test="residue-id mod 10 = 9">&#160;</xsl:if>
    <xsl:if test="residue-id mod $alignWidth = ($alignWidth -1)"><br/></xsl:if>
  </xsl:template>

  <xsl:template match="residue" mode="alignStart">
    <xsl:value-of select="acid"/>
    <xsl:if test="residue-id mod 10 = 9">&#160;</xsl:if>
    <xsl:if test="residue-id mod $alignWidth = ($alignWidth -1)"><br/></xsl:if>
  </xsl:template>

  <xsl:template match="residue" mode="alignFinish">
    <xsl:value-of select="acid"/>
    <xsl:if test="residue-id mod 10 = 9">&#160;</xsl:if>
    <xsl:if test="residue-id mod $alignWidth = ($alignWidth -1)"><br/></xsl:if>
  </xsl:template>


  <xsl:template match="align-start">
    <font color="gray">
      <xsl:apply-templates select="residue" mode="alignStart"/>
    </font>
  </xsl:template>

  <xsl:template match="align-finish">
    <font color="gray">
      <xsl:apply-templates select="residue" mode="alignFinish"/>
    </font>
  </xsl:template>

  <xsl:template match="shift-start">
    <xsl:text disable-output-escaping="yes"><![CDATA[<span  style="position: relative;"><span style="position: absolute; top:-10pt; font-size: 8pt; color:red; text-decoration:none;"></xsl:text>
        ]]></xsl:text><xsl:value-of select="format-number(text(), '0.00')"/><xsl:text disable-output-escaping="yes"><![CDATA[</span><a href="#"><font color="red">]]></xsl:text>
      </xsl:template>

      <xsl:template match="shift-finish"><xsl:text disable-output-escaping="yes"><![CDATA[</font></a></span>]]></xsl:text></xsl:template>

<xsl:template name="tag-view">
  <xsl:param name="value"/>
  <xsl:choose>
    <xsl:when test="$value='nme1'">No N-terminal modification(canonical)</xsl:when>
    <xsl:when test="$value='nme2'">No N-terminal modification(non-canonical)</xsl:when>
    <xsl:when test="$value='nme3'">N-Terminal Methionine Excision(canonical)</xsl:when>
    <xsl:when test="$value='nme4'">N-Terminal Methionine Excision(non-canonical)</xsl:when>
    <xsl:when test="$value='nme5'">N-Terminal Methionine Excision(canonical) and N-Terminal Acetylation</xsl:when>
    <xsl:when test="$value='nme6'">N-Terminal Methionine Excision(non-canonical) and N-Terminal Acetylation</xsl:when>
    <xsl:when test="$value='nme7'">N-Terminal Acetylation</xsl:when>
    <xsl:when test="$value='spt1'">Signal Peptide Removal</xsl:when>
    <xsl:when test="$value='spt2'">Removal of protein prefix that does not follow the canonical AXA signal peptide motif</xsl:when>
    <xsl:when test="$value='unusual-nmod'">Unusual N-terminal modification</xsl:when>
    <xsl:when test="$value='c-truncation'">Removal of protein suffix</xsl:when>
    <xsl:when test="$value='internal'">A modification (or multiple modifications) on internal residue</xsl:when>
    <xsl:when test="$value='alternative'">Potential alternative Start Codon as indicated by a
      truncated prefix that either ends in Met or is followed by Met. In the
      former case, Met is followed by one of the seven amino acids
      describing the NME specificity rule: G,A,P,S,T,V,C.
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="$value"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template name="tag-script">
  <script type="text/javascript">
    function defaultCheck(id, tags) {
    if (!(tags instanceof Array)) {
    tags = [tags];
    }
    for (var i = 0; tags.length > i; i++) {
    if (!(tags[i] instanceof Array)) {
    tags[i] = [tags[i]];
    }
    if (isForTags(id, tags[i]) == tags[i].length) {
    return true;
    }
    }
    return false;
    }

    function showTags(tags, link, f) {
    var i;
    if (f == null)
    f = defaultCheck;
    for (i = 0; ids.length>i; i++) {
    var id = ids[i];
    var st = defaultCheck(id, tags) ? "" : "none";
    document.getElementById('p' + id).style.display = st;
    var align = document.getElementById('pa' + id);
    if (align) {
    align.style.display = st;
    }
    }
    markLink(link);
    }

    var links = new Array();

    function markLink(link) {
    if (link == null)
    return;
    if (0 > links.indexOf(link)) {
    links.push(link);
    }
    for (var i = 0; links.length > i;i++) {
    var curLink = links[i];
    curLink.style.fontWeight = curLink == link ? 'bold' : 'normal';
    }
    }


    function showAll(link) {
    var i;
    for (i = 0; ids.length>i; i++) {
    var id = ids[i];
    document.getElementById('p' + id).style.display  = "";
    }
    markLink(link);
    }

  </script>
</xsl:template>

<xsl:template match="tag" mode="check">
  if (tags[i] == '<xsl:value-of select="value"/>') score++;
</xsl:template>

<xsl:template match="value">
  <xsl:if test="position() > 1"><br/></xsl:if>
  <xsl:call-template name="tag-view">
    <xsl:with-param name="value" select="."/>
  </xsl:call-template>
</xsl:template>

</xsl:stylesheet>

