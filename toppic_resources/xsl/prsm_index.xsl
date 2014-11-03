<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>
    <xsl:template match="peptides">
        <html>
            <title>Peptides list</title>

            <body>
                <xsl:call-template name="navigation"/>

        <xsl:apply-templates select="peptide">
            <xsl:sort select="match-count" order="descending" data-type="number"/>
        </xsl:apply-templates>
                <xsl:call-template name="navigation"/>

            </body>
        </html>
    </xsl:template>

    <xsl:template match="peptide">
        <xsl:value-of select="match-count"/> PrSM<xsl:if test="match-count > 1">s</xsl:if>
        for
        <a href="peptides/peptide{peptide-id}.html">
            <xsl:value-of select="sequence-name"/>
        </a>
        <br/>
    </xsl:template>

    <xsl:template name="navigation">
        <p>
            <a href="index.html">Protein-Spectrum-Matches</a>&#160;Protein Species
        </p>
    </xsl:template>

</xsl:stylesheet>
