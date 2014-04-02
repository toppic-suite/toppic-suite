<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:template name="tags">

        <xsl:call-template name="info">
            <xsl:with-param name="description">without N-terminal modification</xsl:with-param>
            <xsl:with-param name="canonical" select="'nme1'"/>
            <xsl:with-param name="nonCanonical" select="'nme2'"/>
        </xsl:call-template>

        <xsl:call-template name="info">
            <xsl:with-param name="description">with N-Terminal Methionine Excision</xsl:with-param>
            <xsl:with-param name="canonical" select="'nme3'"/>
            <xsl:with-param name="nonCanonical" select="'nme4'"/>
        </xsl:call-template>

        <xsl:call-template name="info">
            <xsl:with-param name="description">with N-Terminal Methionine Excision and N-Terminal Acetylation
            </xsl:with-param>
            <xsl:with-param name="canonical" select="'nme5'"/>
            <xsl:with-param name="nonCanonical" select="'nme6'"/>
        </xsl:call-template>

        <xsl:call-template name="info">
            <xsl:with-param name="description">with N-Terminal Acetylation</xsl:with-param>
            <xsl:with-param name="canonical" select="'nme7'"/>
        </xsl:call-template>

        <xsl:call-template name="info">
            <xsl:with-param name="description">with Signal Peptide Removal</xsl:with-param>
            <xsl:with-param name="canonical" select="'spt1'"/>
        </xsl:call-template>

        <xsl:call-template name="info">
            <xsl:with-param name="description">with Removal of protein prefix that does not follow the canonical AXA
                signal peptide motif
            </xsl:with-param>
            <xsl:with-param name="canonical" select="'spt2'"/>
        </xsl:call-template>

        <xsl:call-template name="info">
            <xsl:with-param name="description">with Unusual N-terminal modification</xsl:with-param>
            <xsl:with-param name="canonical" select="'unusual-nmod'"/>
        </xsl:call-template>

        <xsl:call-template name="info">
            <xsl:with-param name="description">with potential alternative Start Codon as indicated by a truncated prefix
                that either ends in Met or is followed by Met. In the former case, Met is followed by one of the seven
                amino acids describing the NME specificity rule: G, A, P, S, T, V, C
            </xsl:with-param>
            <xsl:with-param name="canonical" select="'alternative'"/>
        </xsl:call-template>


    </xsl:template>


    <xsl:template name="info-core">
        <xsl:param name="description"/>
        <xsl:param name="canonical"/>
        <xsl:param name="nonCanonical"></xsl:param>
        <xsl:param name="objects"/>
        <xsl:variable name="n" select="count($objects/tag[value=$canonical or value=$nonCanonical])"/>
        <xsl:if test="$n > 0">
            <br/>
            <a href="#">
                <xsl:attribute name="onclick">showTagsWrapper([['<xsl:value-of select="$canonical"/>'], ['<xsl:value-of select="$nonCanonical"/>']], this);</xsl:attribute>
                <xsl:value-of select="$n"/>
                <xsl:text> </xsl:text>
                <xsl:value-of select="$description"/>
            </a>
            <xsl:variable name="nc" select="count($objects/tag[value=$nonCanonical])"/>
            <xsl:variable name="ct"
                          select="count($objects[count(tag[value='c-truncation']) > 0 and (count(tag[value=$canonical or value=$nonCanonical]) > 0)])"/>
            <xsl:variable name="internal"
                          select="count($objects[count(tag[value='internal']) > 0 and (count(tag[value=$canonical or value=$nonCanonical]) > 0)])"/>
            <xsl:if test="$n = $nc">, non-canonical</xsl:if>
            <xsl:if test="$n = $ct">, with removal of protein suffix</xsl:if>
            <xsl:if test="$n = $internal">, with internal modifications</xsl:if>

            <xsl:variable name="includeNc" select="($nc > 0 and $n > $nc)"/>
            <xsl:variable name="includeCt" select="($ct > 0 and $n > $ct)"/>
            <xsl:variable name="includeInternal" select="($internal > 0 and $n > $internal)"/>
            <xsl:if test="$includeNc or $includeCt or $includeInternal">
                <i>, including </i>

                <xsl:if test="$includeNc">
                    <a href="#">
                        <xsl:attribute name="onclick">showTagsWrapper(['<xsl:value-of select="$nonCanonical"/>'], this);</xsl:attribute>
                        <xsl:value-of select="$nc"/>
                        non-canonical
                    </a>
                </xsl:if>

                <xsl:if test="$includeCt">
                    <xsl:if test="$includeNc">, </xsl:if>
                    <a href="#">
                        <xsl:attribute name="onclick">showTagsWrapper([['<xsl:value-of select="$canonical"/>', 'c-truncation'], ['<xsl:value-of select="$nonCanonical"/>', 'c-truncation']], this);
                        </xsl:attribute>
                        <xsl:value-of select="$ct"/>
                        with removal of protein suffix
                    </a>
                </xsl:if>

                <xsl:if test="$includeInternal">
                    <xsl:if test="$includeNc or $includeCt">, </xsl:if>
                    <a href="#">
                        <xsl:attribute name="onclick">showTagsWrapper([['<xsl:value-of select="$canonical"/>', 'internal'], ['<xsl:value-of select="$nonCanonical"/>', 'internal']], this);
                        </xsl:attribute>
                        <xsl:value-of select="$internal"/>
                        with internal modifications
                    </a>
                </xsl:if>
            </xsl:if>
        </xsl:if>
    </xsl:template>
</xsl:stylesheet>

