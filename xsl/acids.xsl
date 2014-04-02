<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>


    <xsl:template match="amino_acid_list">
        var acids = [
            <xsl:apply-templates select="amino_acid"/>
        ]
    </xsl:template>

    <xsl:template match="amino_acid">
        {
            name:'<xsl:value-of select="name"/>',
            letter:'<xsl:value-of select="one_letter"/>',
            code:'<xsl:value-of select="three_letter"/>',
            composition:'<xsl:value-of select="composition"/>',
            monoMass:<xsl:value-of select="mono_mass"/>,
            averageMass:<xsl:value-of select="average_mass"/>,
            frequence:<xsl:value-of select="frequency"/>
        },

    </xsl:template>
</xsl:stylesheet>
