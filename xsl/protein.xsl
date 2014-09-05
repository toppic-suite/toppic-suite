<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:include href="basic_prsm.xsl"/>

    <xsl:template match="protein">
        <html>
            <title>Protein species for protein <xsl:value-of select="sequence_name"/>
            </title>
            <body>
                <xsl:call-template name="navigation"/>
                <p>
                    <xsl:value-of select="species_number"/> proteins species for the
                    <strong><xsl:value-of select="sequence_name"/></strong>
                </p>

<!--
                <xsl:apply-templates select="species">
                    <xsl:sort select="min(prsm/e_value)" order="ascending" data-type="number"/>
                </xsl:apply-templates>
-->

                <xsl:call-template name="navigation"/>
            </body>
        </html>
    </xsl:template>


    <xsl:template match="species">
       <div id="p{species_id}">

            <h4><xsl:value-of select="position()"/> Protein species #<xsl:value-of select="species_id"/></h4>
            <xsl:apply-templates select="prsm" mode="species">
<!--
                <xsl:sort select="min(e-value)" order="ascending" data-type="number"/>
-->
            </xsl:apply-templates>
        </div>

    </xsl:template>

    <xsl:template name="navigation">
        <p>
            <a href="../proteins.html">All proteins</a>&#160;
        </p>
    </xsl:template>

</xsl:stylesheet>
