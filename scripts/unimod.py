#!/bin/python
import xml.etree.ElementTree as ET

tree = ET.parse("unimod.xml")
root = tree.getroot()

mod_lst = root.getchildren()[1].getchildren()

root = ET.Element("ptm_list")


for mod in mod_lst:
	mod_name = mod.attrib["full_name"]
	mod_abbr = mod.attrib["title"]
	mod_id = mod.attrib["record_id"]
	tmp_lst = mod.getchildren()
	mod_mass = ""
	for t in tmp_lst:
		if t.tag == "{http://www.unimod.org/xmlns/schema/unimod_2}delta":
			mod_mass = t.attrib["mono_mass"]
	
	doc = ET.SubElement(root, "ptm")
	ET.SubElement(doc, "name").text = mod_name
	ET.SubElement(doc, "abbreviation").text = mod_abbr
	ET.SubElement(doc, "unimod").text = mod_id
	ET.SubElement(doc, "mono_mass").text = mod_mass


tree = ET.ElementTree(root)
tree.write("unimod_ptm.xml")
