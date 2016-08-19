import csvGenerator
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIXML
from io import StringIO
import os
search_database = "nucleotide"
root_dir = "NRPSRoot"
dat_dir = "NRPSData"
nrps_dir= "NRPS"
ana_dir = "Analysis"
gen_dir = "DataGenBank"
dat_stan_dir = "DataStandard"


# this method compare a sequence imputted by the user against all or some of the root sequences so that it can then be
# implemented into the decision tree for classification
def compare(subjct, records):
    csvGenerator.create_dir(dat_dir, nrps_dir)
    csvGenerator.create_dir(os.path.join(dat_dir, gen_dir), nrps_dir)
    csvGenerator.create_dir(os.path.join(dat_dir, dat_stan_dir), nrps_dir)
    new_handle = Entrez.efetch(db=search_database, id=subjct.get(), rettype="gb", retmode="text")
    out_handle = open(os.path.join(ana_dir, os.path.join(nrps_dir, os.path.join(dat_dir, os.path.join(gen_dir, subjct.get())))), "w")
    out_handle.write(new_handle.read())
    out_handle.close()
    new_handle.close()
    sbjct = SeqIO.read(os.path.join(ana_dir, os.path.join(nrps_dir, os.path.join(dat_dir, os.path.join(gen_dir, subjct.get())))), format="gb")
    for record in records:
        recrd = SeqIO.read(os.path.join(root_dir, record[:len(record) - 4]) + ".gbk", format="gb")
        result_handle = NcbiblastpCommandline(query=recrd.seq.transcribe().translate(), subject=sbjct.seq, outfmt=5)()[0]
        save_file = open(os.path.join(ana_dir, os.path.join(nrps_dir,
                             os.path.join(dat_dir, subjct.name + ".xml"))), "w")
        blast_result_record = NCBIXML.read(StringIO(result_handle))
        # Print some information on the result
        for alignment in blast_result_record.alignments:
            for hsp in alignment.hsps:
                save_file.write('sequence:' + alignment.title)
                save_file.write('length:' + alignment.length)
                save_file.write('e value:' + hsp.expect)
                save_file.write(hsp.query)
                save_file.write(hsp.match)
                save_file.write(hsp.sbjct)
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()