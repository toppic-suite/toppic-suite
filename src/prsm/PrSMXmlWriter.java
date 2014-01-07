package edu.ucsd.proteomics.prsm.writer;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;

import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

import edu.ucsd.proteomics.prsm.base.PrSM;

public class PrSMXmlWriter {

	private XMLOutputter outputter;
	PrintWriter xmlWriter;

	public PrSMXmlWriter(File xmlWriterFile) throws Exception {

		xmlWriter = new PrintWriter(new FileOutputStream(xmlWriterFile));
		xmlWriter.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
		xmlWriter.write("<prsm_list>\n");
		outputter = new XMLOutputter(Format.getPrettyFormat());
	}

	public void write(Element element) throws Exception {
		outputter.output(element, xmlWriter);
		xmlWriter.write("\n");
	}

	public void write(PrSM prsm) throws Exception {
		Element element = PrSMXml.toXml(prsm);
		write(element);
	}

	public void write(PrSM prsms[]) throws Exception {
		for (int i = 0; i < prsms.length; i++) {
			if (prsms[i] != null) {
				Element element = PrSMXml.toXml(prsms[i]);
				write(element);
			}
		}
	}

	public void write(PrSM prsms[][]) throws Exception {
		for (int i = 0; i < prsms.length; i++) {
			for (int j = 0; j < prsms[i].length; j++) {
				if (prsms[i][j] != null) {
					Element element = PrSMXml.toXml(prsms[i][j]);
					write(element);
				}
			}
		}
	}

	public void write(PrSM prsms[][][]) throws Exception {
		for (int i = 0; i < prsms.length; i++) {
			for (int j = 0; j < prsms[i].length; j++) {
				for (int k = 0; k < prsms[i][j].length; k++) {
					if (prsms[i][j][k] != null) {
						Element element = PrSMXml.toXml(prsms[i][j][k]);
						write(element);
					}
				}
			}
		}
	}
	
	   public void write(PrSM prsms[][][][]) throws Exception {
	        for (int i = 0; i < prsms.length; i++) {
	            for (int j = 0; j < prsms[i].length; j++) {
	                if (prsms[i][j] != null) {
	                    for (int k = 0; k < prsms[i][j].length; k++) {
	                        for (int l = 0; l < prsms[i][j][k].length; l++) {
	                            if (prsms[i][j][k][l] != null) {
	                                Element element = PrSMXml.toXml(prsms[i][j][k][l]);
	                                write(element);
	                            }
	                        }
	                    }
	                }
	            }
	        }
	    }

	public void close() {
	    xmlWriter.write("</prsm_list>\n");
		xmlWriter.close();
	}
}
