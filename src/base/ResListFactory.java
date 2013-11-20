/**
 * The class factory to generate residue list.
 *
 * @author  Xiaowen Liu
 * @date    2009-8-26
 */

package com.jap.proteomics.base.residue;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.List;

import org.apache.log4j.Logger;

import org.jdom.Document;
import org.jdom.Element;

import com.jap.proteomics.base.util.XmlUtil;


public class ResListFactory {

	private static Logger logger = Logger.getLogger(ResListFactory.class);
	
	public static ResList getDefinedList() throws Exception {
		InputStream stream = ResMng.getDefinedResXmlStream();
		return ResListFactory.getInstance(stream, true);
	}

	/**
	 * Generate a residue list instance from a default xml file.
	 */
	public static ResList getSystemInstance(String fileName) throws Exception {
		InputStream stream = new ResMng().getClass()
				.getResourceAsStream(ResMng.configDir + fileName);
		return ResListFactory.getInstance(stream, false);
	}
	
	/**
	 * Generate a residue list instance from an xml file.
	 */
	public static ResList getFileInstance(File file) throws Exception {
		InputStream stream = new FileInputStream(file);
		return ResListFactory.getInstance(stream, false);
	}
	
	private static ResList getInstance(InputStream stream, boolean isCompleteList) throws Exception {
		Document doc = XmlUtil.getDocument(stream);
		Element root = doc.getRootElement();
		List<?> elementList = root.getChildren();
		int len = elementList.size();
		if (elementList == null || len <= 0) {
			logger.fatal("ReadXml() error: residue no information");
			throw new Exception("residue xml file is empty.");
		}
		ResList resList = new ResList();
		for (int i = 0; i < len; i++) {
			Element element = (Element) elementList.get(i);
			String acidName = element.getChildText("acid");
			String ptmAbbrName = element.getChildText("ptm");
			if (isCompleteList) {
				resList.add(acidName, ptmAbbrName);
			}
			else {
				// new residues are added to the complete list first, then a reference is added to the new resList */
				Res res = ResList.getCompleteList().add(acidName, ptmAbbrName);
				resList.add(res);
			}
		}
        return resList;
	}
}
