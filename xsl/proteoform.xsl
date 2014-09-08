<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:include href="basic_prsm.xsl"/>

    <xsl:template match="species">
        <html>
            <title>Protein species
                #<xsl:value-of select="species_id"/> from <xsl:value-of select="sequence-name"/>
            </title>

            <body>

                <xsl:call-template name="navigation"/>
                <h2>Protein Species #<xsl:value-of select="species_id"/></h2>


                <h3><xsl:value-of select="prsm_number"/> PrSM<xsl:if test="prsm_number > 1">s</xsl:if> 
                    for this protein species </h3>

                <xsl:apply-templates select="prsm[position()=1]/tag/value"/>
                <xsl:apply-templates select="prsm[position()=1]" mode="species"/>
                <table border="1">
                    <tr>
                        <th>E-value</th>
                        <th># all peaks</th>
                        <th># matched peaks</th>
                        <th># matched fragments</th>
                        <th>Link</th>
                    </tr>
                    <xsl:apply-templates select="prsm" mode="id">
                        <xsl:sort select="e_value" order="ascending" data-type="number"/>
                    </xsl:apply-templates>
                </table>

                <xsl:apply-templates select="species"/>


                <xsl:call-template name="navigation"/>

            </body>
        </html>
    </xsl:template>

    <xsl:template match="prsm" mode="id">
        <tr>
            <td align="right">
                <xsl:value-of select="e_value"/>
            </td>
            <td align="center"><xsl:value-of select="count(ms/peaks/peak)"/></td>
            <td align="center"><xsl:value-of select="matched_peak_number"/></td>
            <td align="center"><xsl:value-of select="matched_fragment_number"/></td>
            <td align="center"><a onfocus="cur = {position()-1};selectRow();" href="../prsms/prsm{prsm_id}.html" id="link{position()}">See PrSM>></a></td>
        </tr>
    </xsl:template>


    <xsl:template name="navigation">
        <p>
             <a href="../proteins.html">All proteins</a> /
             <a href="../proteins/protein{sequence_id}.html"><xsl:value-of select="sequence_name"/></a>
        </p>
    </xsl:template>
</xsl:stylesheet>
