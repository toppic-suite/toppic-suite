<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:include href="basic_prsm.xsl"/>

    <xsl:template match="proteins">
        <html>
            <title><xsl:value-of select="count(protein)"/> proteins discovered</title>
            <body>
                <h2><xsl:value-of select="count(protein)"/> proteins discovered</h2>

                <br/>
                <xsl:apply-templates select="protein">
<!--
                    <xsl:sort select="min(species/prsm/e_value)" order="ascending" data-type="number"/>
-->
                </xsl:apply-templates>
            </body>
        </html>
    </xsl:template>


    <xsl:template match="protein">
        <div id="p{sequence_id}">

        <a href="proteins/protein{sequence_id}.html"><xsl:value-of select="sequence_name"/></a>    <br/>

            <xsl:apply-templates select="species/prsm">
<!--
                  <xsl:sort select="e_value" order="ascending" data-type="number"/>
-->
            </xsl:apply-templates>
            <br/>&#160;
        </div>
    </xsl:template>



    <xsl:template match="prsm">
        <xsl:if test="position()=1">
            <xsl:choose>
                <xsl:when test="count(../../species/prsm) = 1">There is only <a href="prsms/prsm{prsm_id}.html">1 PrSM</a> with E-value <xsl:value-of select="e_value"/>.</xsl:when>
                <xsl:otherwise>
                    <a href="prsms/prsm{prsm_id}.html">Best PrSM</a> has E-value <xsl:value-of select="e_value"/>, there
                    <xsl:choose>
                        <xsl:when test="count(../../species) = 1">is <a href="species/species{annotated_protein/species_id}.html">1 protein species</a>.</xsl:when>
                        <xsl:otherwise>are <xsl:value-of select="count(../../species)"/> protein species.</xsl:otherwise>
                    </xsl:choose>
                </xsl:otherwise>
            </xsl:choose>
        </xsl:if>
    </xsl:template>


</xsl:stylesheet>
