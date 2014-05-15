<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:variable name="alignWidth">40</xsl:variable>

    <xsl:template match="prsm" mode="basic">
        <p style="font-family:monospace;font-size:16;line-height:2.5; ">
        <xsl:apply-templates select="annotated_protein/annotation/character" mode="basic"/>
        <xsl:if test="annotated_protein/db_acid_number > annotated_protein/last_residue_position">
          <br/>
          <!--ignore acids first:<xsl:value-of select="floor(annotated_protein/first_residue_position div $alignWidth)*$alignWidth"/>; end:<xsl:value-of select="annotated_protein/db_acid_number - annotated_protein/last_residue_position"/>-->
          display seq start:<xsl:value-of select="floor(annotated_protein/first_residue_position div $alignWidth)*$alignWidth+1"/>; end:<xsl:value-of select="floor(annotated_protein/last_residue_position div $alignWidth)*$alignWidth + $alignWidth"/>
        </xsl:if>
        </p>
    </xsl:template>

    <xsl:template match="prsm" mode="species">
        <xsl:if test="position()=1">
            <p>
            <xsl:choose>
                <xsl:when test="count(../prsm) > 1">
                    <a href="../prsms/prsm{prsm_id}.html">Best PrSM</a> has E-value <xsl:value-of select="e_value"/>
                    and precursor mass <xsl:value-of select="ms/ms_header/precursor_mass"/>,
                    there are <a href="../species/species{../species_id}.html"><xsl:value-of select="count(../prsm)"/> PrSMs</a>.</xsl:when>
                <xsl:otherwise>
                   There is only <a href="../prsms/prsm{prsm_id}.html">1 PrSM</a>
                    with E-value <xsl:value-of select="e_value"/> and precursor mass <xsl:value-of select="precursor_mass"/>.
                </xsl:otherwise>
            </xsl:choose>
            </p>
            <table border="0"  cellspacing="0px" cellpadding="0px">
            <xsl:text disable-output-escaping="yes"><![CDATA[<tr><td height='16px' colspan='44'></td></tr><tr><td height='16px' colspan='44'></td></tr><tr>]]></xsl:text><!--  -->
            <xsl:apply-templates select="." mode="basic"/>
            <xsl:text disable-output-escaping="yes"><![CDATA[</tr>]]></xsl:text></table>

        </xsl:if>
    </xsl:template>

    <xsl:template match="character" mode="basic">
        
        <xsl:if test="type = 'cleavage'">
                   <xsl:if test="cleavage_type = 'unexpected_shift' and shift_no_letter != 0">
                        <xsl:choose>
                            <xsl:when  test="display_position = '0' ">
                                <span  style="position: relative;">
                                    <span style="position: absolute; top:-24pt; font-size: 8pt; color:red; text-decoration:none;">
                                        <xsl:value-of select="shift_no_letter"/>
                                    </span>
                                </span>
                            </xsl:when>
                            <xsl:when  test="display_position = '1' ">
                                <span  style="position: relative;">
                                    <span style="position: absolute; top:-32pt; font-size: 8pt; color:blue; text-decoration:none;">
                                        <xsl:value-of select="shift_no_letter"/>
                                    </span>
                                </span>
                            </xsl:when>
                        </xsl:choose>
                    </xsl:if>
        </xsl:if>
        <xsl:if test="type = 'residue' and position >= floor(../../first_residue_position div $alignWidth)*$alignWidth and (floor(../../last_residue_position div $alignWidth)*$alignWidth+$alignWidth)>position">
            <xsl:if test="position  mod $alignWidth = 0 and position != 0"><!--  -->
                <xsl:text disable-output-escaping="yes"><![CDATA[</tr><tr><td height='16px' colspan='44'></td></tr><tr><td height='16px' colspan='44'></td></tr><tr>]]></xsl:text>
                
            </xsl:if>
            <xsl:if test="position  mod 10 = 0 "><!-- and position != 0 -->
                <xsl:if test="residue_type = 'unexpected_shift' and display_background = '0'">
                <td width='8px' height='16px' bgcolor="#F6CECE">
                    <span style ="color:gray; background:#F6CECE">
                    <xsl:text> </xsl:text>
                    </span>
                    </td>
                </xsl:if>
                <xsl:if test="residue_type = 'unexpected_shift' and display_background = '1'">
                <td width='8px' height='16px' bgcolor="#CECEF6">
                    <span style ="color:gray; background:#CECEF6">
                    <xsl:text> </xsl:text>
                    </span>
                    </td>
                </xsl:if>
                <xsl:if test="residue_type != 'unexpected_shift'">
                <td width='8px' height='16px' style=''>
                    <span style ="color:gray;">
                    <xsl:text> </xsl:text>
                    </span>
                    </td>
                </xsl:if>
<!--
                <xsl:if test="../character[postion()-1]/residue_type = 'unexpected_shift'">
                    <span style ="color:gray; background:#F6CECE">
                    <xsl:text>+</xsl:text>
                    </span>
                </xsl:if>
                <xsl:if test="../character[postion()-1]/residue_type != 'unexpected_shift'">
                    <span style ="color:gray;">
                    <xsl:text>-</xsl:text>
                    </span>
                </xsl:if>
-->
            </xsl:if>
            <xsl:choose>
                <xsl:when test="residue_type = 'n_trunc'">
                <td width='8px' height='16px' style=''>
                    <span style ="color:gray">
                        <xsl:value-of select="acid"/>
                    </span>
                </td>
                </xsl:when>
                <xsl:when test="residue_type = 'c_trunc'">
                <td width='8px' height='16px' style=''>
                    <span style ="color:gray">
                        <xsl:value-of select="acid"/>
                    </span>
                </td>
                </xsl:when>
                <xsl:when test="residue_type = 'unexpected_shift'">
                    <!--<xsl:if test="is_modification = '1'">
                        <xsl:choose>
                            <xsl:when  test="display_position = '0' and display_background = '0'">
                                <div style="position: relative;">
                                    <div style="position: absolute; top:-36px;left:-8px; width:45px; font-size: 8pt; color:red; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                            <xsl:when  test="display_position = '0' and display_background = '1'">
                                <div style="position: relative;">
                                    <div style="position: absolute; top:-36px;left:-8px; width:45px; font-size: 8pt; color:blue; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                            <xsl:when  test="display_position = '1' and display_background = '0'">
                                <div style="position: relative;">
                                    <div id="{floor(position div 30)}" shift="{display_position}" style="position: absolute; top:-56px; width:45px; font-size: 8pt; color:blue; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                            <xsl:when  test="display_position = '1' and display_background = '1'">
                                <div style="position: relative;">
                                    <div id="{floor(position div 30)}" shift="{display_position}" style="position: absolute; top:-56px; width:45px; font-size: 8pt; color:red; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                        </xsl:choose>
                    </xsl:if>-->
                    <xsl:choose>
                        <xsl:when  test="display_background = '0' and is_expected = '0'">
                        <td width='8px' height='16px' bgcolor="#F6CECE">
                        <span style ="color:black; background:#F6CECE">
                            <xsl:value-of select="acid"/>
                        </span>
                        <xsl:choose>
                            <xsl:when  test="is_modification = '1' and display_position = '0' and display_background = '0'">
                                <div style="position: relative;">
                                    <div style="position: absolute; top:-36px; width:100px; font-size: 8pt; color:red; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                            <xsl:when  test="is_modification = '1' and display_position = '1' and display_background = '0'">
                                <div style="position: relative;">
                                    <div id="{floor(position div 30)}" shift="{display_position}" style="position: absolute; top:-56px; width:100px; font-size: 8pt; color:red; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                        </xsl:choose>
                        </td>
                        </xsl:when>
                        <xsl:when  test="display_background = '1' and is_expected = '0'">
                        <td width='8px' height='16px' bgcolor="#CECEF6">
                        <span style ="color:black; background:#CECEF6">
                            <xsl:value-of select="acid"/>
                        </span>
                        <xsl:choose>
                            <xsl:when  test="is_modification = '1' and display_position = '0' and display_background = '1'">
                                <div style="position: relative;">
                                    <div style="position: absolute; top:-36px; width:100px; font-size: 8pt; color:blue; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                            <xsl:when  test="is_modification = '1' and display_position = '1' and display_background = '1'">
                                <div style="position: relative;">
                                    <div id="{floor(position div 30)}" shift="{display_position}" style="position: absolute; top:-56px; width:100px; font-size: 8pt; color:blue; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                        </xsl:choose>
                        </td>
                        </xsl:when>
                        <xsl:when  test="display_background = '0' and is_expected = '1'">
                        <td width='8px' height='16px' bgcolor="#F6CECE">
                        <span style ="color:{shift_style}; background:#F6CECE">
                            <xsl:value-of select="acid"/>
                        </span>
                        <xsl:choose>
                            <xsl:when  test="is_modification = '1' and display_position = '0' and display_background = '0'">
                                <div style="position: relative;">
                                    <div style="position: absolute; top:-36px; width:100px; font-size: 8pt; color:red; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                            <xsl:when  test="is_modification = '1' and display_position = '1' and display_background = '0'">
                                <div style="position: relative;">
                                    <div id="{floor(position div 30)}" shift="{display_position}" style="position: absolute; top:-56px; width:100px; font-size: 8pt; color:red; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                        </xsl:choose>
                        </td>
                        </xsl:when>
                        <xsl:when  test="display_background = '1' and is_expected = '1'">
                        <td width='8px' height='16px' bgcolor="#CECEF6">
                        <span style ="color:{shift_style}; background:#CECEF6">
                            <xsl:value-of select="acid"/>
                        </span>
                        <xsl:choose>
                            <xsl:when  test="is_modification = '1' and display_position = '0' and display_background = '1'">
                                <div style="position: relative;">
                                    <div style="position: absolute; top:-36px; width:100px; font-size: 8pt; color:blue; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                            <xsl:when  test="is_modification = '1' and display_position = '1' and display_background = '1'">
                                <div style="position: relative;">
                                    <div id="{floor(position div 30)}" shift="{display_position}" style="position: absolute; top:-56px; width:100px; font-size: 8pt; color:blue; text-decoration:none;">
                                        <xsl:value-of select="shift"/>
                                    </div>
                                </div>
                            </xsl:when>
                        </xsl:choose>
                        </td>
                        </xsl:when>
                    </xsl:choose>
                </xsl:when>
                <xsl:when test="is_expected = '1'">
                    <td width='8px' height='16px' bgcolor="#DFFFFF">
                    <span style ="color:{shift_style}; background:#DFFFFF">
                        <xsl:value-of select="acid"/>
                    </span>
                    </td>
<!--
                    <xsl:if test="is_modification = '2'">
                        <xsl:choose>
                            <xsl:when  test="display_position = '0'">
                                <span  style="position: relative;">
                                    <span style="position: absolute; top:-24pt; font-size: 8pt; color:blue; text-decoration:none;">
                                        <xsl:value-of select="format-number(shift, '0.00')"/>
                                    </span>
                                </span>
                            </xsl:when>
                            <xsl:when  test="display_position = '1'">
                                <span  style="position: relative;">
                                    <span style="position: absolute; top:-32pt; font-size: 8pt; color:blue; text-decoration:none;">
                                        <xsl:value-of select="format-number(shift, '0.00')"/>
                                    </span>
                                </span>
                            </xsl:when>
                        </xsl:choose>
                    </xsl:if>
                    <span style ="color:black; background:#F6CECE">
                        <xsl:value-of select="acid"/>
                    </span>
-->
                </xsl:when>
                <xsl:when test="residue_type = 'normal'">
                    <td width='8px' height='16px' bgcolor="#DFFFFF">
                    <span style ="color:black; background:#DFFFFF">
                        <xsl:value-of select="acid"/>
                    </span>
                    </td>
                </xsl:when>

            </xsl:choose>
        </xsl:if>
    </xsl:template>


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

