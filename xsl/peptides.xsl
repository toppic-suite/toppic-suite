<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:include href="basic-prsm.xsl"/>
    <xsl:include href="tags.xsl"/>

    <xsl:template match="proteins">
        <html>
            <title><xsl:value-of select="count(protein/peptide)"/> protein species discovered</title>
            <script type="text/javascript">
                var ids  = new Array(<xsl:apply-templates select="protein/peptide" mode="id"/>);

                var parents = new Array(<xsl:apply-templates select="protein" mode="id"/>);

                function isForTags(id, tags) {
                    var score = 0;
                    if (!(tags instanceof Array)) {
                        tags = [tags];
                    }
                    var i;
                    <xsl:apply-templates select="protein/peptide" mode="tag"/>
                    return score;
                }

                function isProteinVisible(id) {
                    <xsl:apply-templates select="protein/peptide" mode="protein"/>
                    return false;
                }

                function showHideProteins() {
                    for (i = 0; parents.length>i; i++) {
                        var id = parents[i];
                        document.getElementById('pp' + id).style.display  = isProteinVisible(id) ?"" : "none";
                    }
                }

                function showTagsWrapper(tags, link) {
                    showTags(tags, link);
                    showHideProteins();
                }
            </script>
            <xsl:call-template name="tag-script"/>
            <body>
                <h2><xsl:value-of select="count(protein/peptide)"/> protein species discovered</h2>
                <a href="#" onclick="showAll(this);showHideProteins();">Show all</a>

                <xsl:call-template name="tags"/>

                <br/><br/>
                <table border="1">
                    <tr><th width="70">Protein species #</th><th width="70">Best PrSM #</th>
                        <th  width="70">E-value</th><th  width="70"># of PrSMs</th>
                        <th  width="105">Scan number</th>
                        <th  width="70">Precursor mass</th>
                        <th  width="70"># all peaks</th>
                        <th  width="70">Duplicate score</th>
                        <th  width="70">Unique score</th>
                        <th  width="70">Alignment type</th>
                        <th  width="70"># modifications</th>
                        <td>&#160;</td></tr>
                <xsl:apply-templates select="protein">
                    <xsl:sort select="min(peptide/prsm/e-value)" order="ascending" data-type="number"/>
                </xsl:apply-templates>
                </table>
            </body>
        </html>
    </xsl:template>

    <xsl:template name="info">
        <xsl:param name="description"/>
        <xsl:param name="canonical"/>
        <xsl:param name="nonCanonical"></xsl:param>
        <xsl:call-template name="info-core">
            <xsl:with-param name="description" select="$description"/>
            <xsl:with-param name="canonical" select="$canonical"/>
            <xsl:with-param name="nonCanonical" select="$nonCanonical"/>
            <xsl:with-param name="objects" select="protein/peptide"/>
        </xsl:call-template>
    </xsl:template>

    <xsl:template match="protein">
        <tr id="pp{protein-id}">
            <td colspan="12">
                <a href="proteins/protein{protein-id}.html"><xsl:value-of select="protein-name"/></a>
            </td>
         </tr>
        <xsl:apply-templates select="peptide">
            <xsl:sort select="min(prsm/e-value)" order="ascending" data-type="number"/>
        </xsl:apply-templates>

    </xsl:template>

    <xsl:template match="peptide">
            <xsl:apply-templates select="prsm">
                  <xsl:sort select="e-value" order="ascending" data-type="number"/>
            </xsl:apply-templates>
    </xsl:template>



    <xsl:template match="peptide" mode="id">
        <xsl:value-of select="peptide-id"/><xsl:if test="last()>position()">,</xsl:if>
    </xsl:template>

    <xsl:template match="protein" mode="id">
        <xsl:value-of select="protein-id"/><xsl:if test="last()>position()">,</xsl:if>
    </xsl:template>

    <xsl:template match="peptide" mode="tag">
        if (<xsl:value-of select="peptide-id"/> == id) {
            for (i=0; tags.length > i; i++) {
                <xsl:apply-templates select="tag" mode="check"/>
            }
        }
    </xsl:template>

    <xsl:template match="peptide" mode="protein">
        if (<xsl:value-of select="../protein-id"/> == id) {
            if (document.getElementById('p<xsl:value-of select="peptide-id"/>').style.display == "") {
                return true;
            }
        }
    </xsl:template>

    <xsl:template match="prsm">
        <xsl:if test="position()=1">

            <tr id="p{../peptide-id}">
                <td align="center">
                    <a href="peptides/peptide{peptide-id}.html">
                        <xsl:value-of select="peptide-id"/>
                    </a>
                </td>

                <td align="center">
                    <a href="prsms/prsm{spectrum-id}.html">
                        <xsl:value-of select="spectrum-id"/>
                    </a>
                </td>
                <td align="right">
                    <xsl:value-of select="e-value"/>
                </td>
                <td align="center">
                    <xsl:value-of select="count(../prsm)"/>
                </td>
                <td align="center">
                    <xsl:value-of select="scans"/>
                </td>
                <td align="center">
                    <xsl:value-of select="precursor-mass"/>
                </td>
                <td align="center">
                    <xsl:value-of select="count(peaks/peak)"/>
                </td>
                <td align="center">
                    <xsl:value-of select="duplicate-score"/>
                </td>
                <td align="center">
                    <xsl:value-of select="unique-score"/>
                </td>
                <td align="center">
                    <xsl:value-of select="alignment-type"/>
                </td>
                <td align="center">
                    <xsl:value-of select="count(align/residue/shift-start)"/>
                </td>
                <td>&#160;</td>

            </tr>
            <tr id="pa{../peptide-id}">
                <td colspan="12">
                    <xsl:apply-templates select="." mode="basic"/>
                </td>
            </tr>
        </xsl:if>
    </xsl:template>


</xsl:stylesheet>
