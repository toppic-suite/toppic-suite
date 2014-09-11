<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>
   
  <xsl:include href="prsm_anno.xsl"/>

  <xsl:template match="protein">
    <html>
      <title>Proteoforms for protein <xsl:value-of select="sequence_name"/> </title>
      <body>
        <xsl:call-template name="navigation"/>

        <p>
          <xsl:value-of select="compatible_proteoform_number"/> proteoforms for protein <xsl:value-of select="sequence_name"/>
        </p>
      </body>

      <xsl:apply-templates select="compatible_proteoform">
        <!-- sort not working  <xsl:sort select="min(proteoform_id)" order="ascending" data-type="number"/> -->
      </xsl:apply-templates>
      <xsl:call-template name="navigation"/>
    </html>
  </xsl:template>

  <xsl:template name="navigation">
    <p>
      <a href="../proteins.html">All proteins</a>&#160;
    </p>
  </xsl:template>

  <xsl:template match="compatible_proteoform">
    <div id="p{proteoform_id}">
      <h4><xsl:value-of select="position()"/> Proteoform #<xsl:value-of select="proteoform_id"/></h4>
      <xsl:apply-templates select="prsm" mode="protein">
      </xsl:apply-templates>
    </div>
  </xsl:template>


  <xsl:template match="prsm" mode="protein">
    <xsl:if test="position()=1">
      <p>
        <xsl:choose>
          <xsl:when test="count(../prsm) > 1">
            The <a href="../prsms/prsm{prsm_id}.html">best PrSM</a> has an E-value <xsl:value-of select="e_value"/>
            and a precursor mass <xsl:value-of select="ms/ms_header/precursor_mass"/>.
            There are <a href="../proteoforms/proteoform{../proteoform_id}.html"><xsl:value-of select="count(../prsm)"/> PrSMs</a> in total.
          </xsl:when>
          <xsl:otherwise>
            There is only <a href="../prsms/prsm{prsm_id}.html">1 PrSM</a>
            with an E-value <xsl:value-of select="e_value"/> and a precursor mass <xsl:value-of select="ms/ms_header/precursor_mass"/>.
          </xsl:otherwise>
        </xsl:choose>
      </p>
      <xsl:apply-templates select="annotated_protein/annotation" mode="protein">
        <!--<xsl:sort select="min(e-value)" order="ascending" data-type="number"/> -->
      </xsl:apply-templates>
    </xsl:if>
  </xsl:template>

  <xsl:template match="annotated_protein/annotation" mode="protein">
    <div id="alignment" style="font-family: 'FreeMono', Miltonian, monospace; font-size:16;line-height:2.5;background-color:#FFF">
      <xsl:variable name="row_residue_num" select="30"/>
      <xsl:variable name="residue_num" select="protein_length"/>
      <xsl:variable name="table_first_row">
        <xsl:variable name="tmp" select="floor(first_residue_position div $row_residue_num)"/>
        <xsl:choose>
          <xsl:when test="$tmp &gt; 0">
            <xsl:value-of select="$tmp - 1"/>
          </xsl:when>
          <xsl:otherwise> 
            <xsl:value-of select="$tmp"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>

      <xsl:variable name="table_total_row_num" select="ceiling($residue_num div $row_residue_num)"/>
      <xsl:variable name="table_last_row">
        <xsl:variable name="tmp" select="floor(last_residue_position div $row_residue_num)"/>
        <xsl:choose>
          <xsl:when test="$table_total_row_num &gt; $tmp + 1">
            <xsl:value-of select="$tmp + 1"/>
          </xsl:when>
          <xsl:otherwise> 
            <xsl:value-of select="$tmp"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>

      <table border="0"  cellspacing="0px" cellpadding="0px">
        <xsl:if test="$table_first_row &gt; 0">
          <tr>
            <td colspan="2"></td>
            <td colspan="62" align="center">
              <xsl:text>... </xsl:text> 
              <xsl:value-of select="$table_first_row * 30"/>
              <xsl:text> amino acid residues are skipped at the N-terminus ...</xsl:text> 
            </td>
            <td colspan="3"></td>
          </tr>
        </xsl:if>
        <xsl:call-template name="add_table_row">
          <xsl:with-param name="i" select="$table_first_row"/>
          <xsl:with-param name="table_last_row" select="$table_last_row"/>
        </xsl:call-template>

        <xsl:if test="$table_last_row * 30 + 30 &lt; $residue_num">
          <tr> <td>&#160; </td></tr>
          <tr>
            <td colspan="2"></td>
            <td colspan="62" align="center">
              <xsl:text>... </xsl:text> 
              <xsl:value-of select="$residue_num - $table_last_row * 30 - 30"/>
              <xsl:text> amino acid residues are skipped at the C-terminus ...</xsl:text> 
            </td>
            <td colspan="3"></td>
          </tr>
        </xsl:if>
      </table>         
    </div>
  </xsl:template>

  <xsl:template match="cleavage">
    <xsl:text disable-output-escaping="yes"><![CDATA[<span style="]]></xsl:text>
      <xsl:if test="cleavage_type = 'seq_start'">
        <xsl:text>font-weight:bold;color:red;</xsl:text>
      </xsl:if>
      <xsl:if test="cleavage_type = 'seq_end'">
        <xsl:text>font-weight:bold;color:red;</xsl:text>
      </xsl:if>

      <xsl:if test="is_unexpected_change = '1'">
        <xsl:if test="unexpected_change_color = 0">
          <xsl:text>background:#DFFFFF;</xsl:text>
        </xsl:if>
        <xsl:if test="unexpected_change_color = 1">
          <xsl:text>background:#CECEF6;</xsl:text>
        </xsl:if>
      </xsl:if>
      <xsl:text disable-output-escaping="yes"><![CDATA[">]]></xsl:text>

      <xsl:choose>
        <xsl:when test="cleavage_type = 'seq_start'">
          <xsl:text>]</xsl:text>
        </xsl:when>
        <xsl:when test="cleavage_type = 'seq_end'">
          <xsl:text>[</xsl:text>
        </xsl:when>
        <xsl:otherwise>
          <xsl:text>&#160;</xsl:text>
        </xsl:otherwise>
      </xsl:choose>
      <xsl:text disable-output-escaping="yes"><![CDATA[</span>]]></xsl:text>
  </xsl:template>

</xsl:stylesheet>

