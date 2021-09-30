import os
from psims.mzid import MzIdentMLWriter
import csv
import re

PROTON = 1.007276466879
FIXED_PTM = [{'name':'C57', 'mass':57.021464, 'uni_id':'4'},{'name':'C58', 'mass':58.005479,  'uni_id':'6'}]

def mass_to_charge_ratio(mass, z, charge_carrier=PROTON):
    return abs((mass + (charge_carrier * z)) / z)

class IdentiPyMzIdentMLWriter(object):
    def __init__(self, param, data, db_data, total_ptm, user_fixed_ptm, var_ptm):
        self.param = param
        self.data = data
        self.db_data = db_data
        self.total_ptm = total_ptm
        self.user_fixed_ptm = user_fixed_ptm
        self.var_ptm = var_ptm
        self.spec_id = 0
        self.spec_ident_id = 0

    def trim_seq(self, seq):#remove (), " " and dots from a sequence
        seq = seq.replace('(', '')
        seq = seq.replace(')', '')
        seq = seq.replace('"', '')
        
        new_seq = ''
        
        #remove dots from the beginning or end
        for i in range(0, len(seq)):
            if seq[i] == "." and i > 0:
                try:
                    int(seq[i -1])
                    new_seq = new_seq + seq[i]
                except ValueError:
                    continue
            else:
                if seq[i] != ".":
                    new_seq = new_seq + seq[i]
        return new_seq
        
    def transform_peptide(self, prsm):
        #parse mass shift in a sequence and add it to modification data
        #generate a new sequence without mass shift markers
        modifications = []

        seq = self.trim_seq(prsm['Proteoform'])

        parsedseq = '';
        position = 0;
        massShift = '';
        isMassShift = False;

        for i in range(0,  len(seq)):
            if seq[i] == "[":
                isMassShift = True;
            elif seq[i] == "]":
                isMassShift = False;
                if massShift != '':
                    temp_pos = position;
                    try:
                        int(massShift[-1])
                        mass = float(massShift);
                        spec = {
                            'location': temp_pos,
                            'name': 'unknown modification',
                            'monoisotopic_mass_delta': mass,
                            'residues':parsedseq[temp_pos -1]
                        }
                        modifications.append(spec)
                    except ValueError:
                    #if not in self.total_ptm, there is an error somewhere, so skip it.
                        for ptm in self.total_ptm:
                            if ptm['name'] == massShift:
                                spec = {
                                    'location': temp_pos,
                                    'name': massShift,
                                    'monoisotopic_mass_delta': ptm['mass'],
                                    'residues':parsedseq[temp_pos -1],
                                    'uni_id': ptm['uni_id']#unimod id
                                }
                                modifications.append(spec)
                                break
                    massShift = ''
            else:
                if (isMassShift):
                    massShift = massShift + seq[i];
                else:
                    parsedseq = parsedseq + seq[i];
                    position += 1                
        spec = {
            'id': "peptide_" + str(prsm['Prsm ID']),
            'modifications': modifications,
            'peptide_sequence': parsedseq
        }
        return spec
    
    def transform_search_modifications(self):#generate a list of ptms used for database search
        fmods = self.param['Fixed modifications']
        modifications = []

        #fixed ptm information
        if fmods == "None":#fixed ptm was NONE
            modifications.append({#below skips generation of searchParam tag
                'mass_delta': '',
                'fixed': '',
                'residues':'',
                'params':  {'unknown modification': ''}
            })
        elif fmods.split('.')[-1] == 'txt':#fixed ptm given as a file by a user
            for ptm in self.user_fixed_ptm:
                modifications.append({
                    'mass_delta': ptm['mass'],
                    'fixed': True,
                    'residues': ptm['residue'],
                    'name':ptm['uni_id'],
                    'params':  [{'unknown modification': ''}], 'residues': ptm['residue']}
                    )
        else:#fixed ptm is either C57 or C58
            fixed_ptm = fmods.split(':')[0]
            for ptm in FIXED_PTM:
                if fixed_ptm == ptm['name']:
                    modifications.append({
                        'mass_delta': ptm['mass'],
                        'fixed': True,
                        'residues': 'C',
                        'name':ptm['uni_id'],
                        'params':  [{'unknown modification': ptm['uni_id']}], 'residues': 'C'}
                        )
        #variable PTM information (only for TopMG)
        if len(self.var_ptm) > 0:
            for ptm in self.var_ptm:
                modifications.append({
                    'mass_delta': ptm['mass'],
                    'fixed': False,
                    'residues': ptm['residue'],
                    'name':ptm['uni_id'],
                    'params':  [{'unknown modification': ''}], 'residues': ptm['residue']}
                    )
        return modifications
    
    def build_spectrum_identification_protocol(self):
        return {
            "threshold": [
                {"score threshold": ''},
            ],
            "modification_params": self.transform_search_modifications(),
            "fragment_tolerance": '',
            "parent_tolerance": ''
        }
    
    def transform_search_database(self):#add db information
        return {
            'file_format': 'fasta format',
            'name': os.path.basename(self.param['Protein database file']),
            'id': 1,
            'location': 'file:///' + os.path.abspath(self.param['Protein database file']),
            'params': [
                {"decoy DB accession regexp": ''},
                "decoy DB type randomized"
            ]
        }
    
    def transform_spectra_data(self):#add spectrum file information
        file_path = os.path.abspath(self.param['Spectrum file']).replace('_ms2.msalign', '.mzML')
        spec = {
            'id': 1,
            'name': file_path,
            'location':'',
            'spectrum_id_format': 'multiple peak list nativeID format',
        }
        spec['file_format'] = 'mzML format'
        
        return spec
    
    def transform_spectrum_identification_item(self, prsm):#add matched spectrum information
        items = []
        theo_mass = ''
        charge = prsm['Charge']
        mass_error = ''
        theo_mz = mass_to_charge_ratio(float(prsm['Adjusted precursor mass']), int(charge))
        exper_mz = mass_to_charge_ratio(float(prsm['Precursor mass']), int(charge))
        spec = {
            'charge_state':charge,
            'peptide_id': "peptide_" + prsm['Prsm ID'],
            'peptide_evidence_id': "pepevi_" + prsm['Prsm ID'],
            'id': 'SII_' + str(self.spec_ident_id),
            'experimental_mass_to_charge': exper_mz,
            'calculated_mass_to_charge': theo_mz,
        }
        items.append(spec)
        self.spec_ident_id += 1
        return items
    
    def transform_spectrum_identification_result(self, prsm):#add additional details of matched spectrum
        spec = {
            'id': 'SIR_' + str(self.spec_id),
            'spectrum_id': prsm['Scan(s)'],#this should be scan ID not spectra ID
            'params': [{'name': 'scan start time', 'value': prsm['Retention time'], 'unit_name': 'second'}]
        }
        identifications = []
        identifications.extend(self.transform_spectrum_identification_item(prsm))
        spec['identifications'] = identifications
        self.spec_id += 1
        return spec
    
    def write(self, output_file_path):
        with open(output_file_path, 'wb') as fh:
            writer = MzIdentMLWriter(fh, True)
            with writer:
                writer.controlled_vocabularies()
                writer.providence(software={
                    "name": "TopPIC",
                    "uri": "http://proteomics.informatics.iupui.edu/software/toppic/",
                    "version": self.param['Version'],
                })
                writer.register("SpectraData", 1)
                writer.register("SearchDatabase", 1)
                writer.register("SpectrumIdentificationList", 1)
                writer.register("SpectrumIdentificationProtocol", 1)

                with writer.sequence_collection():#write db sequences
                    for seq in self.db_data:
                        prot_desc = seq['desc']
                        prot_acc = seq['acc']
                        writer.write_db_sequence(prot_acc, None, id=prot_acc, params=[
                            {"protein description": prot_desc.replace(prot_acc + ' ', '')}
                        ])    
                    for prsm in self.data:#write proteoform information
                        writer.write_peptide(self.transform_peptide(prsm))
                        writer.write_peptide_evidence(
                            "peptide_" + prsm['Prsm ID'], prsm['Protein accession'],
                            "pepevi_" + prsm['Prsm ID'],
                            prsm['First residue'], prsm['Last residue']
                        )
                with writer.analysis_collection():
                    writer.SpectrumIdentification([1], [1]).write(writer)

                with writer.analysis_protocol_collection():#parameters used for identification (ex: mod)
                    writer.spectrum_identification_protocol(**self.build_spectrum_identification_protocol())
                with writer.data_collection():
                    writer.inputs([], [self.transform_search_database()],[self.transform_spectra_data()])
                    with writer.analysis_data():
                        with writer.spectrum_identification_list(id=1):
                            for prsm in self.data:#each identification result
                                writer.write_spectrum_identification_result(
                                    **self.transform_spectrum_identification_result(prsm))
            return writer

    @classmethod
    def write_mzid(cls, param, data, db_data, total_ptm, user_fixed_ptm, var_ptm):
        writer = cls(param, data, db_data, total_ptm, user_fixed_ptm, var_ptm)
        mzid_name = param['Spectrum file'].split('_ms2.msalign')[0]
        return writer.write(mzid_name + ".mzid")

def parse_user_ptm(path):#parse ptm txt file
    file_path = path
    data = []

    if path == 'None':
        return data
        
    with open(file_path) as ptm_file:
        while True:
            line = ptm_file.readline()
            if not line:
                break
            line = line.replace('\n','')
            line = line.strip()
            if line != '':
                if line[0] != '#':
                    ptm = {}
                    line = line.split(',')
                    ptm['name'] = line[0]
                    ptm['mass'] = line[1]
                    ptm['residue'] = line[2]
                    ptm['uni_id'] = line[4]
                    data.append(ptm)
    return data
    
def parse_ptm():#parse ptm_data.hpp
    dir_name = os.path.split(os.path.abspath(__file__))[0]
    file_path = os.path.join(dir_name, "ptm_data.hpp")
    data = []
    ptm = {}
    with open(file_path) as ptm_file:
        while True:
            line = ptm_file.readline()
            if not line:
                break
            line = line.replace('\n', '')
            line = line.strip()
            if line[1] == 'p':
                data.append(ptm)
                ptm = {}
            elif line[1] == 'a':
                line = line.split('>')[1]
                line = line.split('<')[0]
                ptm['name'] = line
            elif line[1] == 'm':
                line = line.split('>')[1]
                line = line.split('<')[0]
                ptm['mass'] = line
            elif line[1] == 'u':
                line = line.split('>')[1]
                line = line.split('<')[0]
                ptm['uni_id'] = line
    return data[1:]#data[0] is empty

def parse_fasta(path):
    file_path = path
    sequences= []
    current = {}
    seq = ''
    with open(file_path) as db_file:
        while True:
            line = db_file.readline()
            if not line:
                 break
            line = line.replace('\n','')
            if line[0] == '>':
                if seq != '':
                    current['seq'] = seq
                    sequences.append(current)
                    seq = ''
                    current = {}
                line = line[1:] #remove '>'
                splitted = line.split(' ', 1)
                current['acc'] = splitted[0]
                current['desc'] = splitted[1]
            else:
                seq = seq + line
    return sequences
    
def parse_tsv(path):
    ident_file = path
    isParam = True
    isSkip = False
    parameters = {}
    data =[]
    col_name = []
    with open(ident_file) as tsv_file:
        reader = csv.reader(tsv_file, delimiter="\t", quotechar='|')
        next(reader)
        for row in reader:
            if '*' in row[0]:
                isParam = False#end of parameter
                #read column names. If there is ptm information, skip them
                column_row = next(reader)
                if '*' in column_row[0]:#if it is ptm information, skip
                    new_row = next(reader)
                    while '*' not in new_row[0]:
                        new_row = next(reader)
                        #then new row is the closing ********* for ptm information
                        
                    column_row = next(reader)#go to next line
                
                    if '*' in column_row[0]:#if it is ptm information, skip
                        new_row = next(reader)
                        while '*' not in new_row[0]:
                            new_row = next(reader)

                        column_row = next(reader)#go to next line and read column names

                for col in column_row:
                    col_name.append(col)
            else:
                if isParam:
                    try:
                        tmp = row[1]
                        name = row[0].strip()[:-1]
                        parameters[name] = row[1]
                    except:#fixed ptm file name is not written in tsv format in the toppic used in topmsv server
                        name = row[0].split(',')[0]
                        name = name.strip()[:-1]
                        file_path = (row[0].split(',')[1]) 

                        parameters[name] = file_path
                else:
                    if not isSkip:
                        tmp = {}
                        for i in range(0, len(row)):
                            tmp[col_name[i]] = row[i]
                        data.append(tmp)
    return parameters, data

if __name__ == "__main__":
    import sys
    tsv_file = sys.argv[1]
    db_file = sys.argv[2]
    user_fixed_ptm_path = sys.argv[3]#None is entered if not used 
    user_localize_ptm_path = sys.argv[4]#None entered if not used
    user_var_ptm_path = sys.argv[5]#None entered if not used

    param, data = parse_tsv(tsv_file)#param = parameters in top of tsv file, data = each prsm data
    db_data = parse_fasta(db_file)#original and decoy sequences
    ptm_data = parse_ptm()#all ptm in ptm_data.hpp file
    
    user_fixed_ptm = parse_user_ptm(user_fixed_ptm_path)#empty if no file was used
    common_ptm = parse_user_ptm(user_localize_ptm_path)#empty if no file was used
    var_ptm = parse_user_ptm(user_var_ptm_path)#empty if no file was used

    total_ptm = user_fixed_ptm + ptm_data + common_ptm + var_ptm#all ptm used for sequence annotation
    
    IdentiPyMzIdentMLWriter.write_mzid(param, data, db_data, total_ptm, user_fixed_ptm, var_ptm)

