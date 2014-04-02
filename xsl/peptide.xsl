<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:include href="basic-prsm.xsl"/>

    <xsl:template match="peptide-page">
        <html>
            <title>Protein species
                #<xsl:value-of select="peptide/peptide-id"/> from <xsl:value-of select="peptide/sequence-name"/>
            </title>

            <script>
                var rowCount = <xsl:value-of select="count(peptide/prsm)"/>;
                var cur = 0;

                function selectRow() {
                    for (var i = 0; rowCount > i; i++) {
                        var row = document.getElementById("row" + (i+1));
                        row.style.backgroundColor = cur == i ? "yellow" : "white";
                    }
                    document.getElementById("link" + (cur+1)).focus();
                }

                function keyCheck(event) {
                    var KeyID = event.keyCode;
                    if (event.ctrlKey) {
                        switch(KeyID) {
                            case 37:
                            case 38: cur--; break;
                            case 39:
                            case 40: cur++; break;
                        }
                        cur = (cur + rowCount)%rowCount;
                        selectRow();
                    }
                 }
            </script>

            <body>

                <xsl:call-template name="navigation"/>

                        <xsl:apply-templates select="peptide"/>

                <p>Use <strong>Ctrl+Arrows</strong> for navigation between PrSMs, <strong>Enter</strong> for going to the selected one.</p>

                <xsl:call-template name="navigation"/>

                <script type="text/javascript">
                    selectRow();
                    document.onkeyup = keyCheck;
                </script>
            </body>
        </html>
    </xsl:template>


    <xsl:template match="peptide">

        <h2>Protein Species #<xsl:value-of select="peptide-id"/></h2>

        <xsl:apply-templates select="prsm[position()=1]" mode="basic"/>

        <h3><xsl:value-of select="count(prsm)"/> PrSM<xsl:if test="count(prsm)>1">s</xsl:if> for this protein species </h3>

        <div>
        <xsl:apply-templates select="prsm[position()=1]/tag/value"/>
        </div>
        <br/>

        <table border="1">
            <tr>
                <th>E-value</th>
                <th>Duplicate score</th>
                <th>Unique score</th>
                <th># all peaks</th>
                <th># matched peaks</th>
                <th># not matched peaks</th>
                <th>Link</th>
            </tr>
            <xsl:apply-templates select="prsm" mode="id">
                <xsl:sort select="e-value" order="ascending" data-type="number"/>
            </xsl:apply-templates>
        </table>

    </xsl:template>

    <xsl:template match="prsm" mode="id">
        <tr onclick="document.location='../prsms/prsm{spectrum-id}.html'" id="row{position()}">
            <td align="right">
                <xsl:value-of select="e-value"/>
            </td>
            <td align="center"><xsl:value-of select="duplicate-score"/></td>
            <td align="center"><xsl:value-of select="unique-score"/></td>
            <td align="center"><xsl:value-of select="total"/></td>
            <td align="center"><xsl:value-of select="matched"/></td>
            <td align="center"><xsl:value-of select="unmatched"/></td>
            <td align="center"><a onfocus="cur = {position()-1};selectRow();" href="../prsms/prsm{spectrum-id}.html" id="link{position()}">See PrSM>></a></td>
        </tr>
    </xsl:template>


    <xsl:template name="navigation">
        <p>
             <a href="../proteins.html">All proteins</a> /
             <a href="../proteins/protein{peptide/sequence-id}.html"><xsl:value-of select="peptide/sequence-name"/></a>
        </p>
    </xsl:template>
</xsl:stylesheet>
