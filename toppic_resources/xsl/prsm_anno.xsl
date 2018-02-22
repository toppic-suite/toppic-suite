<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
  <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

  <xsl:template match="residue">
    <xsl:text disable-output-escaping="yes"><![CDATA[ style="]]></xsl:text>
    <xsl:if test="residue_type = 'known_change'">
      <xsl:text>font-weight:bold;color:red;</xsl:text>
    </xsl:if>
    <xsl:if test="residue_type = 'n_truncation'">
      <xsl:text>color:grey;</xsl:text>
    </xsl:if>
    <xsl:if test="residue_type = 'c_truncation'">
      <xsl:text>color:grey;</xsl:text>
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
    <xsl:if test="is_unexpected_change = '1'">
      <a style="text-decoration:none" href="#">
        <xsl:attribute name="title">
          <xsl:value-of select="anno" />
        </xsl:attribute>
        <xsl:if test="not(residue_type = 'known_change')">
          <xsl:if test="possible_pos_color = 1">
            <font color="red">
              <xsl:value-of select="acid" />
            </font>
          </xsl:if>
          <xsl:if test="possible_pos_color = 0">
            <font color="black">
              <xsl:value-of select="acid" />
            </font>
          </xsl:if>					
        </xsl:if>
        <xsl:if test="residue_type = 'known_change'">
          <font color="red">
            <xsl:value-of select="acid" />
          </font>					
        </xsl:if>
      </a>
    </xsl:if>
    <xsl:if test="not(is_unexpected_change = '1')">
      <xsl:value-of select="acid" />
    </xsl:if>
    <xsl:text disable-output-escaping="yes"><![CDATA[]]></xsl:text>
  </xsl:template>

  <xsl:template match="matched_peaks" mode="prsm">
    <xsl:for-each select="matched_peak">
      <xsl:value-of select="ion_position"/><xsl:text>,</xsl:text>
    </xsl:for-each>
  </xsl:template>

  <xsl:template match="matched_peaks" mode="title">
    <xsl:for-each select="matched_peak">
      <xsl:value-of select="ion_type"/>
      <xsl:value-of select="ion_display_position"/>
      <xsl:text>&#160;</xsl:text>
      <xsl:value-of select="peak_charge"/>
      <xsl:text>+&#160;</xsl:text>
    </xsl:for-each>
  </xsl:template>

  <xsl:template name="add_one_letter">
    <xsl:param name="i" />
    <xsl:param name="j" />
    <xsl:param name="k" />
    <xsl:if test="$k &lt; 10">
      <xsl:variable name="pos" select="$i * 30 + $j * 10 + $k"/>
      <xsl:text disable-output-escaping="yes"><![CDATA[<td width="14px" ]]></xsl:text>
        <xsl:if test="$pos &lt;= protein_length + 1">
          <xsl:apply-templates select="cleavage[position = $pos]"/>
        </xsl:if>
        <xsl:text disable-output-escaping="yes"><![CDATA[</td>]]></xsl:text>
      <xsl:text disable-output-escaping="yes"><![CDATA[<td width="14px" ]]></xsl:text>
        <xsl:if test="$pos &lt;= protein_length">
          <xsl:apply-templates select="residue[position = $pos]"/>
        </xsl:if>
        <xsl:text disable-output-escaping="yes"><![CDATA[</td>]]></xsl:text>
      <xsl:call-template name="add_one_letter">
        <xsl:with-param name="i" select="$i"/>
        <xsl:with-param name="j" select="$j"/>
        <xsl:with-param name="k" select="$k+1"/>
      </xsl:call-template>
    </xsl:if>
  </xsl:template>

  <xsl:template name="add_ten_letters">
    <xsl:param name="i" />
    <xsl:param name="j" />

    <xsl:if test="$j &lt; 3">
      <xsl:call-template name="add_one_letter">
        <xsl:with-param name="i" select="$i"/>
        <xsl:with-param name="j" select="$j"/>
        <xsl:with-param name="k" select="0"/>
      </xsl:call-template>

      <xsl:if test="$j &lt; 2">
        <td width="16px"></td>
      </xsl:if>

      <xsl:call-template name="add_ten_letters">
        <xsl:with-param name="i" select="$i"/>
        <xsl:with-param name="j" select="$j+1"/>
      </xsl:call-template>
    </xsl:if>
  </xsl:template>

  <xsl:template name="add_begin_position">
    <xsl:param name="position" />

    <xsl:text disable-output-escaping="yes"><![CDATA[<td width="40px" align="right">]]></xsl:text>
      <xsl:value-of select="$position + 1"/>
      <xsl:text disable-output-escaping="yes"><![CDATA[</td>]]></xsl:text>

    <xsl:text disable-output-escaping="yes"><![CDATA[<td width="16px"></td>]]></xsl:text>
  </xsl:template>

  <xsl:template name="add_end_position">
    <xsl:param name="position" />

    <xsl:text disable-output-escaping="yes"><![CDATA[<td width="16px"></td>]]></xsl:text>
    <xsl:text disable-output-escaping="yes"><![CDATA[<td width="40px" align="left">]]></xsl:text>
      <xsl:choose>
        <xsl:when test="$position &lt; protein_length">
          <xsl:value-of select="$position + 1"/>
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="protein_length"/>
        </xsl:otherwise>
      </xsl:choose>
      <xsl:text disable-output-escaping="yes"><![CDATA[</td>]]></xsl:text>
  </xsl:template>


  <xsl:template match="unexpected_change" mode="first_column">
    <xsl:param name="pos" />

    <xsl:text disable-output-escaping="yes"><![CDATA[<td align="left"; colspan="]]></xsl:text>

      <xsl:variable name="end_pos">
        <xsl:choose>
          <xsl:when test="right_position  &gt;= $pos + 60">
            <!-- 59 + 2 blanks + 2 right end num + 1 cell for large numbers -->
            <xsl:value-of select="$pos + 64"/>
          </xsl:when>
          <xsl:otherwise> 
            <xsl:value-of select="right_position"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>

      <xsl:variable name="tmp" select="$end_pos - $pos + 1"/>
      <xsl:variable name="blank_num" select="floor(($tmp - 1) div 20)"/>
      <xsl:variable name="span" select="$tmp + $blank_num"/>
      <xsl:value-of select="$span"/>

      <xsl:text disable-output-escaping="yes"><![CDATA[">]]></xsl:text>

      <xsl:if test="$pos = left_position and segment_type ='SHIFT'">

        <xsl:choose>
          <xsl:when test="unexpected_change_color = 0">
            <xsl:text disable-output-escaping="yes"><![CDATA[<font color="#2D3333">]]></xsl:text>
              <xsl:value-of select="mass_shift"/>
              <xsl:text disable-output-escaping="yes"><![CDATA[</font>]]></xsl:text>
          </xsl:when>
          <xsl:when test="unexpected_change_color = 1">
            <xsl:text disable-output-escaping="yes"><![CDATA[<font color="#3E3E4A">]]></xsl:text>
              <xsl:value-of select="mass_shift"/>
              <xsl:text disable-output-escaping="yes"><![CDATA[</font>]]></xsl:text>
          </xsl:when>
        </xsl:choose>
      </xsl:if>

      <xsl:text disable-output-escaping="yes"><![CDATA[</td>]]></xsl:text>

  </xsl:template>

  <xsl:template match="unexpected_change" mode="non_first_column">
    <xsl:param name="pos" />
    <xsl:param name="row_start_pos" />

    <xsl:variable name="end_pos">
      <xsl:choose>
        <xsl:when test="right_position  &gt;= $row_start_pos + 60">
          <!-- 59 + 2 blanks + 2 right end num + 1 cell for large numbers -->
          <xsl:value-of select="$row_start_pos + 64"/>
        </xsl:when>
        <xsl:otherwise> 
          <xsl:value-of select="right_position"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <xsl:text disable-output-escaping="yes"><![CDATA[<td align="left"  colspan="]]></xsl:text>

      <xsl:variable name="tmp" select="$end_pos - $pos + 1"/>
      <xsl:variable name="blank_num" select="floor(($end_pos mod 60) div 20) - floor(($pos mod 60)  div 20)"/>
      <xsl:variable name="span" select="$tmp + $blank_num"/>
      <xsl:value-of select="$span"/>
      <xsl:text disable-output-escaping="yes"><![CDATA[" ]]></xsl:text>

      <xsl:text disable-output-escaping="yes"><![CDATA[>]]></xsl:text>

      <xsl:if test="$pos = left_position and segment_type ='SHIFT'">

        <xsl:choose>
          <xsl:when test="unexpected_change_color = 0">
            <xsl:text disable-output-escaping="yes"><![CDATA[<font color="#2D3333">]]></xsl:text>
              <xsl:value-of select="mass_shift"/>
              <xsl:text disable-output-escaping="yes"><![CDATA[</font>]]></xsl:text>
          </xsl:when>
          <xsl:when test="unexpected_change_color = 1">
            <xsl:text disable-output-escaping="yes"><![CDATA[<font color="#3E3E4A">]]></xsl:text>
              <xsl:value-of select="mass_shift"/>
              <xsl:text disable-output-escaping="yes"><![CDATA[</font>]]></xsl:text>
          </xsl:when>
        </xsl:choose>
      </xsl:if>

      <xsl:text disable-output-escaping="yes"><![CDATA[</td>]]></xsl:text>

  </xsl:template>

  <xsl:template name="add_shift">
    <xsl:param name="i" />
    <xsl:param name="j" />

    <xsl:if test="$j &lt; 60">
      <xsl:variable name="pos" select="$i * 60 + $j"/>
      <xsl:choose> 
        <xsl:when test="$j = 0">
          <xsl:apply-templates select="unexpected_change[left_position &lt;= $pos and right_position &gt;= $pos ]" mode="first_column">
            <xsl:with-param name="pos" select="$pos"/>
          </xsl:apply-templates>
        </xsl:when>

        <xsl:when test="$j &gt; 0">
          <xsl:apply-templates select="unexpected_change[left_position = $pos]" mode="non_first_column">
            <xsl:with-param name="pos" select="$pos"/>
            <xsl:with-param name="row_start_pos" select="$i * 60"/>
          </xsl:apply-templates>
        </xsl:when>
      </xsl:choose>
      <xsl:call-template name="add_shift">
        <xsl:with-param name="i" select="$i"/>
        <xsl:with-param name="j" select="$j+1"/>
      </xsl:call-template>
    </xsl:if>

  </xsl:template>


  <xsl:template name="add_table_row">
    <xsl:param name="i" />
    <xsl:param name="table_last_row" />
    <xsl:if test="$table_last_row &gt;= $i">

      <!-- shifts -->
      <tr>
        <xsl:text>&#10;</xsl:text>
        <xsl:text disable-output-escaping="yes"><![CDATA[<td colspan="2">&#160;</td>]]></xsl:text>

        <xsl:call-template name="add_shift">
          <xsl:with-param name="i" select="$i"/>
          <xsl:with-param name="j" select="0"/>
        </xsl:call-template>

        <xsl:text>&#10;</xsl:text>

      </tr>
      <xsl:text>&#10;</xsl:text>

      <tr>
        <xsl:text>&#10;</xsl:text>
        <xsl:variable name="row_residue_num" select="30"/>
        <xsl:variable name="begin_pos" select="$i * $row_residue_num"/>

        <xsl:call-template name="add_begin_position">
          <xsl:with-param name="position" select="$begin_pos"/>
        </xsl:call-template>


        <xsl:call-template name="add_ten_letters">
          <xsl:with-param name="i" select="$i"/>
          <xsl:with-param name="j" select="0"/>
        </xsl:call-template>

        <xsl:variable name="end_pos" select="$i * $row_residue_num + $row_residue_num - 1"/>

        <xsl:call-template name="add_end_position">
          <xsl:with-param name="position" select="$end_pos"/>
        </xsl:call-template>
        <!-- an extra cell for large shift numbers -->
        <td></td>

        <xsl:text>&#10;</xsl:text>
      </tr>
      <xsl:text>&#10;</xsl:text>

      <xsl:call-template name="add_table_row">
        <xsl:with-param name="i" select="$i + 1"/>
        <xsl:with-param name="table_last_row" select="$table_last_row"/>
      </xsl:call-template>
    </xsl:if>
  </xsl:template>

  <xsl:template match="occurence">
    <xsl:value-of select="acid_letter"/>
    <xsl:value-of select="position + 1"/>
    <xsl:text>&#160;</xsl:text>
  </xsl:template>

  <xsl:template match="expected_change">
    <xsl:element name="a">
      <xsl:attribute name="href">
        <xsl:text>http://www.unimod.org/modifications_view.php?editid1=</xsl:text>
        <xsl:value-of select="ptm/unimod" />
      </xsl:attribute>
      <xsl:attribute name="target">
        <xsl:text>_blank</xsl:text>
      </xsl:attribute>
      <font color="red">
        <xsl:value-of select="ptm/abbreviation"/>
        <xsl:text>&#160;[</xsl:text>
        <xsl:apply-templates select="occurence"/>
        <xsl:text>]</xsl:text>
      </font>
    </xsl:element>
    <xsl:text>&#160;</xsl:text>
  </xsl:template>

  <xsl:template match="unexpected_change">
    <xsl:text>&#160;</xsl:text>
    <xsl:variable name="s_type" select="segment_type" />
    <xsl:if test="contains($s_type, 'SHIFT')">
      <xsl:variable name="known_ptm" select="count(ptm)" />
      <xsl:if test="$known_ptm &gt; 0">
        <xsl:element name="a">
          <xsl:attribute name="href">
            <xsl:text>http://www.unimod.org/modifications_view.php?editid1=</xsl:text>
            <xsl:value-of select="ptm/unimod" />
          </xsl:attribute>
          <xsl:attribute name="target">
            <xsl:text>_blank</xsl:text>
          </xsl:attribute>
          <font color="red">
            <xsl:value-of select="ptm/abbreviation" />
            <xsl:text>&#160;[</xsl:text>
            <xsl:value-of select="occurence" />
            <xsl:text>]&#160;&#160;</xsl:text>
          </font>
        </xsl:element>
      </xsl:if>
      <xsl:if test="$known_ptm = 0">
        <xsl:choose>
          <xsl:when test="occurence != ''">
            <font color="red">
              <xsl:text>Unknown [</xsl:text>
              <xsl:value-of select="occurence" />
              <xsl:text>]</xsl:text>
            </font>
          </xsl:when>
          <xsl:otherwise>
            <font color="red">
              <xsl:text>Unknown [</xsl:text>
              <xsl:value-of select="mass_shift" />
              <xsl:text>]</xsl:text>  
            </font>  
          </xsl:otherwise>
        </xsl:choose>
      </xsl:if>
      <xsl:text>&#160;</xsl:text>			
    </xsl:if>
  </xsl:template>

</xsl:stylesheet>

