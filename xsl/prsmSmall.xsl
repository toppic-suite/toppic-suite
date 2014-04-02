<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>
    <xsl:template match="prsm">
        <html>
            <body>
                <xsl:apply-templates select="align"/>
            </body>
        </html>
    </xsl:template>


    <xsl:template match="align">
        <p style="font-family:Monospace;font-size:17;line-height: 18pt">
            <xsl:apply-templates select="align-start"/>
            <xsl:apply-templates select="residue"/>
            <xsl:apply-templates select="align-finish"/>
        </p>
    </xsl:template>



    <xsl:template match="residue">
        <xsl:apply-templates select="shift-start"/>

        <xsl:value-of select="acid"/>

        <xsl:apply-templates select="shift-finish"/>
        <xsl:choose>
            <xsl:when test="not(in-shift)">
            <a style="text-decoration:none" href="#">
                <xsl:if test="last()>position()">
                    <xsl:attribute name="title">b<xsl:value-of select="position()"/> y<xsl:value-of
                            select="last()-position()"/>
                    </xsl:attribute>
                </xsl:if>
                <xsl:choose>
                    <xsl:when test="count(ion)>0">
                        <span id="break{residue-id}">
                            <xsl:choose>
                            <xsl:when test="(count(ion/b-ion) > 0) and (count(ion/y-ion) > 0)">
                                <xsl:text disable-output-escaping="yes">&amp;#x23B1;</xsl:text>
                            </xsl:when>
                            <xsl:when test="count(ion/b-ion) > 0">
                                <xsl:text disable-output-escaping="yes">&amp;#x23AB;</xsl:text>
                            </xsl:when>
                            <xsl:when test="count(ion/y-ion) > 0">
                                <xsl:text disable-output-escaping="yes">&amp;#x23A9;</xsl:text>
                            </xsl:when>
                        </xsl:choose>
                        </span>

                    </xsl:when>
                    <xsl:when test="position() = last() and ../align-finish">
                        <font color="red">[</font>
                    </xsl:when>
                    <xsl:otherwise>&#160;</xsl:otherwise>
                </xsl:choose>
            </a>
        </xsl:when>
            <xsl:otherwise>&#160;</xsl:otherwise>
        </xsl:choose>


        <xsl:if test="residue-id mod 60 = 59"><br/></xsl:if>
    </xsl:template>

    <xsl:template match="residue" mode="alignStart">
        <xsl:value-of select="acid"/>
        <xsl:if test="last()>position()">&#160;</xsl:if>
        <xsl:if test="residue-id mod 60 = 59"><br/></xsl:if>
    </xsl:template>

    <xsl:template match="residue" mode="alignFinish">
        <xsl:value-of select="acid"/>&#160;<xsl:if test="residue-id mod 60 = 59"><br/></xsl:if>
    </xsl:template>


    <xsl:template match="align-start">
        <font color="gray">
            <xsl:apply-templates select="residue" mode="alignStart"/>
        </font>
        <font color="red">]</font>
    </xsl:template>

    <xsl:template match="align-finish">
        <font color="gray">
            <xsl:apply-templates select="residue" mode="alignFinish"/>
        </font>
    </xsl:template>

    <xsl:template match="shift-start">
        <xsl:text disable-output-escaping="yes"><![CDATA[<span  style="position: relative;"><span style="position: absolute; top:-10pt; font-size: 8pt; color:red; text-decoration:none;"></xsl:text>
        ]]></xsl:text><xsl:value-of select="."/><xsl:text disable-output-escaping="yes"><![CDATA[</span><a href="#"><font color="red">]]></xsl:text>
    </xsl:template>
    <xsl:template match="shift-finish">
        <xsl:text disable-output-escaping="yes"><![CDATA[</font></a></span>]]></xsl:text>
    </xsl:template>

</xsl:stylesheet>
