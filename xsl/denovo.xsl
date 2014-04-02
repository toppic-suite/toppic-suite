<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:variable name="countBY"
                  select="count(prsm/peaks/peak/breaks/break[type='B']) + count(prsm/peaks/peak/breaks/break[type='Y'])"/>
    <xsl:variable name="beginIon">
        <xsl:choose>
            <xsl:when test="$countBY > 0">B</xsl:when>
            <xsl:otherwise>C</xsl:otherwise>
        </xsl:choose>
    </xsl:variable>
    <xsl:variable name="endIon">
        <xsl:choose>
            <xsl:when test="$countBY > 0">Y</xsl:when>
            <xsl:otherwise>Z</xsl:otherwise>
        </xsl:choose>
    </xsl:variable>

    <xsl:template match="prsm">
        <html>
<head>
            <title>De Novo Approach for Spectrum #<xsl:value-of select="spectrum-id"/> -
                Peptide #<xsl:value-of select="peptide-id"/> for <xsl:value-of select="sequence-name"/>
            </title>
    <script type="text/javascript" src="../acids.js"/>
    <style>
        span.match,td.break{font-weight:bold;};
        {font-weight:bold;};
    </style>

</head>
            <body>

                <xsl:call-template name="navigation"/>

                <h2>De Novo Approach for  Spectrum #<xsl:value-of select="spectrum-id"/>
                </h2>

                <table cellpadding="3">
                    <tr>
                        <td>PrSM ID: </td>
                        <td>
                            <xsl:value-of select="prsm-id"/>
                        </td>
                        <td>
                             Protein mass:
                        </td>
                        <td>
                            <xsl:value-of select="sequence-mass"/>
                        </td>

                        <td>
                             Precursor m/z:
                        </td>
                        <td>
                            <xsl:value-of select="precursor-mz"/>
                        </td>
                        <td>Precursor charge:</td>
                        <td>
                            <xsl:value-of select="precursor-charge"/>
                        </td>
                        <td>Precursor mass:</td>
                        <td>
                            <xsl:value-of select="precursor-mass"/>
                        </td>
                    </tr>
                    <tr>
                        <td colspan="8">
                            Scan(s): <xsl:value-of select="scans"/>
                        </td>
                    </tr>
                </table>
                <br/>

                <script>
                    var total = <xsl:value-of select="precursor-mass"/>;
                    var peaks = [
                        <xsl:apply-templates select="peaks/peak"/>
                    ];

                    var all = [];

                    for (var i = 0; peaks.length > i; i++) {
                        all[peaks[i].id] = peaks[i];
                    }

                    var protein = '<xsl:apply-templates select="align/align-start/residue/acid"/><xsl:apply-templates select="align/residue/acid"/><xsl:apply-templates select="align/align-finish/residue/acid"/>';

                    function getMassLocal(p) {
                        var m = p.mass;
                        if (p.mod == "ammonia") {
                            m += ammonia;
                        }
                        if (p.mod == "water") {
                            m += water;
                        }
                        if (p.mod == "plus") {
                            m += 1;
                        }
                        if (p.mod == "minus") {
                            m -= 1;
                        }
                        return  p.reverse == 1 ? total - m : m;
                    }

                    function getMass(p) {
                        var ans = getMassLocal(p)
                        var count = 1;
                        var a = p.attach;
                        for (var i = 0; a.length > i; i++) {
                            ans += getMass(a[i]);
                            count++;
                        }
                        return ans / count;
                    }

                    function getMassStr(peak) {
                        var ans = "";
                        ans += peak.mass;

                        if (peak.matched) {
                            ans+='*';
                        }
                        if (peak.selected) {
                            ans += '!';
                        }
                        ans += '&lt;br&gt;';
                        return ans;
                    }

                    function getMassMod(peak) {
                        var ans = "";
                        if (peak.mod == null) {
                            ans+= '&lt;a href="javascript:plus(' + peak.id + ');"&gt;+&lt;/a&gt; ';
                            ans+= '&lt;a href="javascript:addAmmonia(' + peak.id + ');"&gt;NH3&lt;/a&gt; ';
                            ans+= '&lt;a href="javascript:addWater(' + peak.id + ');"&gt;H20&lt;/a&gt; ';
                            ans+= '&lt;a href="javascript:minus(' + peak.id + ');"&gt;-&lt;/a&gt; ';
                        }  else {
                            if (peak.mod == "ammonia") {
                                ans += "NH3"
                            }
                            if (peak.mod == "water") {
                                ans += "H20";
                            }
                            if (peak.mod == "plus") {
                                ans += "+";
                            }
                            if (peak.mod == "minus") {
                                ans += "-";
                            }

                            ans+= ' &lt;a href="javascript:clear(' + peak.id + ');"&gt;clear&lt;/a&gt; ';
                        }
                        ans += '&lt;br&gt;';
                        return ans;
                    }

                    function getAllMass(peak) {
                        var ans = getMassStr(peak);
                        var a = peak.attach;
                        for (var i = 0; a.length > i; i++) {
                            ans += " " + getMassStr(a[i]);
                        }
                        return ans;
                    }

                    function getAllMassMod(peak) {
                        var ans = getMassMod(peak);
                        var a = peak.attach;
                        for (var i = 0; a.length > i; i++) {
                            ans += " " + getMassMod(a[i]);
                        }
                        return ans;
                    }

                    function getIntensity(peak) {
                        var ans = peak.intensity;
                        var a = peak.attach;

                        for (var i = 0; a.length > i; i++) {
                            ans += " " + a[i].intensity;
                        }
                        return ans;
                    }


                    function getMassRepresentation (p, stop) {
                        var ans = getMassLocal(p).toFixed(3);
                        if (p.reverse == 1) {
                            ans =  "(" + ans + ")";
                        }
                        if (stop)
                            return ans;
                        var a = p.attach;
                        for (var i = 0; a.length > i; i++) {
                            ans += " " + getMassRepresentation (a[i], true);
                        }
                        return ans;
                    }
                    function comparePeaks(p1, p2) {
                        return getMass(p1) - getMass(p2);
                    }

                    var proton = 1.007276;
                    var isotope = 1.00235;
                    var water = 18.010565;
                    var ammonia = 17.0265;
                    var oxygen = 15.9949;

                    function appendInfo(cell, p1, p2) {
                        var delta = getMass(p1) - getMass(p2);
                        var text = delta.toFixed(3) + " ";
                        if (delta>20 &amp;&amp; 260 > delta) {
                            var pos = -1;
                            var best = 1000000;
                            for (var i = 0; acids.length > i; i++) {
                                var acid = acids[i];
                                var nbest = Math.abs(acid.monoMass - delta);
                                nbest = Math.min(nbest, Math.abs(acid.monoMass - delta - 1));
                                nbest = Math.min(nbest, Math.abs(acid.monoMass - delta + 1));
                                if (best > nbest) {
                                    best = nbest;
                                    pos = i;
                                }
                                if (0.1 > nbest) {
                                    text += acids[i].letter  +  (delta - acids[i].monoMass).toFixed(3) + " ";
                                }
                            }
                            if (0.1 > best) {
                                cell.style.fontWeight='bold';
                            }
                        }
                        cell.appendChild(document.createTextNode(text));
                    }

                    function appendDiff(cell, p1, p2) {
                        var t = 0;
                        if (p2 == null) {
                            return;
                        }

                        var m =   getMass(p1) - getMass(p2);
                        var cur = p2.residueId  + 1;

                        var a  = "";
                        var text = m.toFixed(3) + "&lt;br&gt;";
                        while (m > t) {
                            var ch = protein.charAt(cur);
                            t += acidMass(ch);
                            a += ch;
                            cur++;

                            var s = a;

                            var delta =  m - t;
                            if (300 > Math.abs(delta)) {
                                s = a + " " + delta.toFixed(3);
                                if (0.1 > Math.abs(delta)) {
                                    s = "&lt;b&gt;" + s + "&lt;/b&gt;";
                                }
                                text += s + "&lt;br&gt;";
                            }
                        }

                        cell.innerHTML = text;
                    }

                    function populate() {
                        peaks.sort(comparePeaks);
                        var table = document.getElementById("peaks");
                        while (table.rows.length > 1) {
                            table.deleteRow(1);
                        }

                        var cur = 0;
                        var log = new Object();
                        log.changes = [];
                        for (var i = 0;  all.length > i; i++) {
                            var change = new Object();
                            var peak = all[i];
                            if (peak != null) {
                                change.mass = peak.mass;
                                change.reverse = peak.reverse;
                                var mod =  peak.mod;
                                if (mod != null) {
                                    change.mod = mod;
                                }
                                log.changes.push(change);
                            }
                        }

                        document.getElementById("log").value = JSON.stringify(log);

                        var prev = null;

                        for (var i = 0;  peaks.length > i; i++) {
                            var peak = peaks[i];

                            var row = table.insertRow(i+1);

                            var c = [];
                            for (var j = 0; 14 > j; j++) {
                                c[j] = row.insertCell(j);
                            }

                            addLink(c[0], 'reverse',  "reverse("+ i +")");
                            if (i>0) {
                                addLink(c[0], 'join up',  "unite("+ i +")");
                            }
                            if (peak.attach.length > 0) {
                                addLink(c[0], 'separate',  "separate("+ i +")");
                            }

                            if (peak.matched == 1) {
                                addLink(c[0], 'process',  "process("+ i +")");
                            }

                            c[1].innerHTML = getAllMass(peak);
                            c[2].innerHTML = getAllMassMod(peak);
                            c[3].appendChild(document.createTextNode(getIntensity(peak)));
                            c[4].appendChild(document.createTextNode(getMassRepresentation(peak)));

                            if (peak.context != null) {
                                var context = peak.context + " " + peak.contextFull;
                                c[5].appendChild(document.createTextNode(context));
                                c[6].appendChild(document.createTextNode(
                                        acidMass(context.charAt(0)).toFixed(3) + " " + acidMass(context.charAt(1)).toFixed(3)
                                    )
                                );
                            }


                            appendDiff(c[7], peak, prev);

                            for (var k = 1; 7 >k; k++) {
                                if (i >= k) {
                                    appendInfo(c[k+7], peak, peaks[i-k]);
                                }
                            }
                            if (peak.matched == 1) {
                                prev = peak;
                            } else {
                                var a = peak.attach;
                                for (var j = 0; a.length > j; j++) {
                                    if (a[j].matched == 1) {
                                        prev = a[j];
                                    }
                                }
                            }
                        }
                    }

                    function addLink(cell, text, action) {
                        var link = document.createElement('a');
                        link.href= "javascript:" + action + ";";
                        link.appendChild(document.createTextNode(text));
                        cell.appendChild(document.createTextNode(' '));
                        cell.appendChild(link);
                    }

                    function acidMass(a) {
                        for (var i = 0; acids.length > i; i++) {
                                var acid = acids[i];
                                if (acid.letter == a) {
                                    return acid.monoMass;
                                }
                        }
                    }

                    function reverse(id) {
                        var peak = peaks[id];
                        peak.reverse = 1 - peak.reverse;
                        var a = peak.attach;
                        for (var i = 0; a.length > i; i++) {
                            a[i].reverse = 1 - a[i].reverse;
                        }
                        populate();
                    }

                    function unite(id) {
                        var t = peaks[id-1];
                        var s = peaks[id];
                        doUnite(t, s);
                        peaks.splice(id,1);
                        populate();
                    }

                    function doUnite(t, s) {
                        if (s.context != null) {
                            t.context = s.context;
                            t.contextFull = s.contextFull;
                        }
                        t.attach.push(s);
                        var a = s.attach;
                        for (var i = 0; a.length > i; i++) {
                            t.attach.push(a[i]);
                        }
                        s.attach = new Array();
                    }

                    function addAmmonia(id) {
                       all[id].mod = "ammonia";
                       populate();
                    }
                    function addWater(id) {
                       all[id].mod = "water";
                       populate();
                    }

                    function plus(id) {
                       all[id].mod = "plus";
                       populate();
                    }

                    function minus(id) {
                       all[id].mod = "minus";
                       populate();
                    }
                    function clear(id) {
                       all[id].mod = null;
                       populate();
                    }

                    function load() {
                        var changes =  JSON.parse(document.getElementById("log").value).changes;
                        var epsilon = 0.001;
                        for (var j = 0; changes.length > j; j++) {
                            var change = changes[j];
                            for (var i = 0; all.length > i; i++) {
                                var peak = all[i];
                                if (peak != null) {
                                    if (epsilon >  Math.abs(change.mass - peak.mass)) {
                                        peak.reverse = change.reverse;
                                        peak.mod = change.mod;
                                        if (peak.mod == "")
                                            peak.mod = null;
                                    }
                                }
                            }
                        }
                        populate();
                    }

                    function separateAll() {
                        for (var i = 0;  peaks.length > i; i++) {
                            var peak = peaks[i];
                            var a = peak.attach;
                            for (var j = 0; a.length > j; j++) {
                                peaks.push(a[j]);
                            }
                            peak.attach = new Array();
                        }
                        populate();
                    }

                    var score;
                    var delta;

                    function process(id) {
                        var peak = peaks[id];
                        separateAll();
                        doProcess(peak);
                        union();
                        showMessage(score + " matches found!");
                    }

                    function processAll() {
                        separateAll();
                        var best = -1;
                        var bd = 0;
                        var peak;
                        var mod;
                        var check = [null, "plus", "minus"];
                        for (var i = 0;  peaks.length > i; i++) {
                            var p = peaks[i];
                            if (p.matched == 1) {
                                for (var j = 0; check.length > j; j++) {
                                    var m = check[j];
                                    p.mod = m;
                                    doProcess(p);
                                    if ((score > best) || (score == best &amp;&amp; bd > delta)) {
                                        peak = p;
                                        mod = m;
                                        best = score;
                                        bd = delta;
                                    }
                                }
                                showMessage( i + " peaks processed out of " + peaks.length);  //TODO: make it work on all browsers
                            }
                        }
                        peak.mod = mod;
                        doProcess(peak);
                        var message = "The best peak is " + peak.mass + " ";
                        if (mod == "plus")  {
                            message += "+1 ";
                        }
                        if (mod == "minus")  {
                            message += "-1 ";
                        }
                        message += "with " + score + " matches and average difference " + (delta/score);

                        for (var i = 0;  peaks.length > i; i++) {
                            if (!peaks[i].selected) {
                                peaks[i].mod = null;
                            }
                        }
                        union();
                        showMessage(message);
                    }

                    function showMessage(message) {
                        document.getElementById("message").innerHTML = message;
                    }

                    function doProcess(peak) {
                        for (var i = 0;  peaks.length > i; i++) {
                            peaks[i].selected = false;
                        }
                        peak.selected = true;
                        var cur = peak.residueId  + 1;
                        m = getMass(peak);
                        score = 0;
                        delta = 0;
                        var t = m;
                        for (var cur = peak.residueId; cur > 0; cur--) {
                            t -= acidMass(protein.charAt(cur));
                            for (var i = 0;  peaks.length > i; i++) {
                                var p = peaks[i];
                                if (p != peak) {
                                    check(p, t);
                                }
                            }
                        }
                        for (var i = 0;  peaks.length > i; i++) {
                            var p = peaks[i];
                            if (p != peak) {
                                check(p, 0);
                            }
                        }
                        t = m;
                        for (var cur = peak.residueId + 1; protein.length > cur; cur++) {
                            t += acidMass(protein.charAt(cur));
                            for (var i = 0;  peaks.length > i; i++) {
                                var p = peaks[i];
                                if (p != peak) {
                                    check(p, t);
                                }
                            }
                        }
                    }

                    function check(p, t) {
                        var modOrig = p.mod;
                        var reverseOrig = p.reverse;
                        if (checkModes(p, t)) return;
                        if (p.matched == 0) {
                            p.reverse = 1 - p.reverse;
                            if (checkModes(p, t)) return;
                        }

                        p.mod = modOrig;
                        p.reverse = reverseOrig
                    }

                    function checkModes(p, t) {
                        if (checkMod(p, t, null)) return true;
                        if (checkMod(p, t, "plus")) return true;
                        if (checkMod(p, t, "minus")) return true;
                        if (p.matched == 0) {
                            if (checkMod(p, t, "ammonia")) return true;
                            if (checkMod(p, t, "water")) return true;
                        }
                        return false;
                    }
                    function checkMod(p, t, mod) {
                        p.mod = mod;
                        var diff =  Math.abs(getMass(p)-t);
                        if (0.1 > diff) {
                            score++;
                            delta += diff;
                            p.selected = true;
                            return true;
                        }
                        return false;
                    }

                    function separate(id) {
                        var peak = peaks[id];
                        var a = peak.attach;
                        for (var j = 0; a.length > j; j++) {
                            peaks.push(a[j]);
                        }
                        peak.attach = new Array();
                        populate();
                    }

                    function union() {
                        peaks.sort(comparePeaks);
                        var t = peaks[0];
                        var i = 1;
                        while ( peaks.length > i) {
                            var s = peaks[i];
                            if (0.1 > Math.abs(getMass(s) - getMass(t))) {
                                doUnite(t, s);
                                peaks.splice(i,1);
                            } else {
                                t = s;
                                i++;
                            }
                        }
                        populate();
                    }

                </script>


                <a href="#" onclick="processAll();">Find the best</a>
                <br/>
                <h2 id="message"></h2>
                <p>
                    * - peak was matched by the original MS-Align+ algorithm, ! - peak was matched during the last process.
                </p>
                <!--<a href="#" onclick="populate();">Populate</a>-->
                <table id="peaks" border="1">
                    <tr><th></th><th width="70">Mass</th><th width="90">Mod</th><th width="70">Intensity</th><th  width="70">N-mass</th>
                        <th width="40">Context</th><th>Context masses</th>
                        <th  width="70">Diff</th>
                        <th width="70">Delta 1</th><th width="70">Delta 2</th><th width="70">Delta 3</th><th width="70">Delta 4</th><th width="70">Delta 5</th><th width="70">Delta 6</th></tr>

                </table>
                Current state:  <br/>
                             <textarea id="log" cols="100" rows="3">
                             </textarea>
                <br/>
                             <a href="#" onclick="load();">Load changes</a>
                <br/>
                <a href="#" onclick="union();">Union</a>
                <br/>
                <a href="#" onclick="separateAll();">Separate</a>

                    <script>
                        populate();
                    </script>
                <xsl:call-template name="navigation"/>
            </body>
        </html>
    </xsl:template>

    <xsl:template match="peak">
        {
            id: <xsl:value-of select="position()"/>,
            mass: <xsl:value-of select="monoisotopic-mass"/>,
            intensity: <xsl:value-of select="intensity"/>,
            matched: <xsl:choose>
                        <xsl:when test="count(breaks/break)>0">1</xsl:when>
                        <xsl:otherwise>0</xsl:otherwise>
                    </xsl:choose>,
            reverse: <xsl:choose>
                        <xsl:when test="count(breaks/break[type=$endIon])>0">1</xsl:when>
                        <xsl:otherwise>0</xsl:otherwise>
                    </xsl:choose>,
            attach:new Array(),
            selected: false
            <xsl:if test="count(breaks/break)>0">
                ,
                residueId :  <xsl:value-of select="min(breaks/break/residue-id)"/>,
                context :
                    '<xsl:apply-templates select="../../align/residue" >
                        <xsl:with-param name="break"><xsl:value-of select="min(breaks/break/residue-id)"/></xsl:with-param>
                    </xsl:apply-templates>',
                contextFull :
                    '<xsl:apply-templates select="../../align/residue" mode="full">
                        <xsl:with-param name="break"><xsl:value-of select="min(breaks/break/residue-id)"/></xsl:with-param>
                    </xsl:apply-templates>'
            </xsl:if>
        },
    </xsl:template>

    <xsl:template match="residue">
        <xsl:param name="break"/>
        <xsl:if test="residue-id = $break or residue-id = $break + 1">
            <xsl:value-of select="acid"/>
        </xsl:if>
    </xsl:template>

    <xsl:template match="residue" mode="full">
        <xsl:param name="break"/>
        <xsl:if test="(residue-id > ($break -3)) and  (($break + 4)>residue-id)">
            <xsl:value-of select="acid"/>
        </xsl:if>
    </xsl:template>


    <xsl:template match="peak" mode="denovo">
        <tr>
                <td>
                    <xsl:if test="count(breaks/break)>0">
                        <xsl:attribute name="class">break</xsl:attribute>
                    </xsl:if>

                    <xsl:value-of select="monoisotopic-mass"/>
                    <xsl:text> </xsl:text>
                    <xsl:if test="position() > 1">
                        <xsl:value-of
                                select="format-number(monoisotopic-mass - preceding::peak[1]/monoisotopic-mass, '0.000')"/>
                    </xsl:if>
                </td>
                <xsl:apply-templates select="preceding::peak[6>position()]"  mode="preceding">
                    <xsl:with-param name="mass"><xsl:value-of select="monoisotopic-mass"/></xsl:with-param>
                </xsl:apply-templates>
        </tr>
    </xsl:template>

    <xsl:template match="peak" mode="preceding">
        <xsl:param name="mass"/>
        <td>
            <xsl:call-template name="delta">
                <xsl:with-param name="value"><xsl:value-of select="$mass - monoisotopic-mass"/></xsl:with-param>
            </xsl:call-template>
        </td>
    </xsl:template>

    <xsl:template name="delta">
        <xsl:param name="value"/>
        <xsl:apply-templates select="document('acid.xml')/amino_acid_list/amino_acid">
            <xsl:with-param name="mass"><xsl:value-of select="$value"/></xsl:with-param>
            <xsl:sort select="abs($value - mono_mass)" order="ascending"/>
        </xsl:apply-templates>

    </xsl:template>

    <xsl:template match="amino_acid">
        <xsl:param name="mass"/>
        <xsl:if test="position() = 1 or 0.1 > abs(mass - mono_mass)">
        <span>
            <xsl:if test="1.1 > abs($mass - mono_mass)">
                <xsl:attribute name="class">match</xsl:attribute>
            </xsl:if>
        <xsl:text> </xsl:text><xsl:value-of select="one_letter"/>
        <xsl:text> </xsl:text>
        <xsl:value-of select="format-number($mass - mono_mass, '0.000')"/>
        </span>
        </xsl:if>
    </xsl:template>

    <xsl:template match="acid"><xsl:value-of select="."/></xsl:template>

    <xsl:template name="navigation">
        <p>
             <a href="../proteins.html">All proteins</a> /
             <a href="../proteins/protein{sequence-id}.html"><xsl:value-of select="sequence-name"/></a> /
             <a href="../peptides/peptide{peptide-id}.html">Peptide #<xsl:value-of select="peptide-id"/></a>
        </p>
    </xsl:template>
</xsl:stylesheet>
