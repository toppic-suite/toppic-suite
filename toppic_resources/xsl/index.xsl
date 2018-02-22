<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>
  <xsl:template match="index">
    <html>
      <title>All Protein-Spectrum-Matches (PrSMs) list</title>
      <script>

        function autoResize(id){
          frame = document.getElementById(id);
          frame.height = frame.contentWindow.document.body.scrollHeight;
        }

        function addRemoveRow(prsmId) {
          tbl = document.getElementById('prsmTable');
          lastRow = tbl.rows.length;

          imageRow = document.getElementById('image' + prsmId);
          if (imageRow != null) {
            for (i =0; lastRow > i; i++) {
              if (tbl.rows[i].id=='image' + prsmId) {
  
              }
            }
            tbl.deleteRow(i);
            return;
          }

          for (i =0; lastRow > i; i++) {
            if (tbl.rows[i].id=='row' + prsmId) {
              break;
            }
          }

          row = tbl.insertRow(i + 1);

          row.id = 'image' + prsmId;

          cell = row.insertCell(0);
          cell.colSpan = 6;

          align = document.createElement('iframe');
          id = 'align' + prsmId;
          align.id = id;
          align.src = 's/prsm' + prsmId + '.html';
          align.frameBorder = 0;
          align.scrolling = 'no';
          align.width = '1320';
          align.height = '10';
          align.onload = function(){autoResize(id);};
          cell.appendChild(align);
        }
      </script>
      <body>
        <xsl:call-template name="navigation" />

        <br/>
        <table border="0">
          <tr>
            <td align="right">Identified spectras:</td>
            <td align="left"><xsl:value-of select="total-prsm"/></td>
          </tr>
          <tr>
            <td align="right">Identified proteins:</td>
            <td align="left"><xsl:value-of select="total-sequence"/></td>
          </tr>
          <tr>
            <td align="right">Identified protein species:</td>
            <td align="left"><xsl:value-of select="total-peptide"/></td>
          </tr>
        </table>
        <br/>
        <table border="0">
          <xsl:apply-templates select="modification-stat"/>
        </table>
        <br/>

        <xsl:apply-templates select="spectrums"/>

        <xsl:call-template name="navigation" />
      </body>
    </html>
  </xsl:template>

  <xsl:template match="modification-stat">
    <tr>
      <td>Identified protein species with <xsl:value-of select="number"/> modifications: </td>
      <td><xsl:value-of select="count"/></td>
    </tr>
  </xsl:template>

  <xsl:template match="spectrums">
    <table border="1" id="prsmTable">
      <thead><tr><th>Spectrum ID</th><th>Protein name</th><th>Unique score</th><th># of modifications</th>
          <th>E-value</th><th>Protein mass</th><th>Scans</th></tr></thead>
      <xsl:apply-templates select="spectrum/prsm"/>
    </table>
  </xsl:template>

  <xsl:template match="prsm">
    <tr id="row{spectrum-id}">
      <td>
        <a href="prsms/prsm{spectrum-id}.html"><xsl:value-of  select="spectrum-id"/></a>
        &#160;
        <a href="#" onclick="addRemoveRow({spectrum-id}); return false;">+/-</a>
      </td>
      <td><xsl:value-of select="substring-before(substring-after(sequence-name, '|'),'|')"/></td>
      <td><xsl:value-of select="unique-score"/></td>
      <td><xsl:value-of select="count(align/residue/shift-start)"/></td>
      <td><xsl:value-of select="e-value"/></td>
      <td align="right"><xsl:value-of select="sequence-mass"/></td>
      <td><xsl:value-of select="scans"/></td>

    </tr>

  </xsl:template>

  <xsl:template name="navigation">
    <p>
      Protein-Spectrum-Matches <a href="index-peptide.html">Protein Species</a>
    </p>
  </xsl:template>

</xsl:stylesheet>
