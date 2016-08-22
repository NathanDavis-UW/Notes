import csvGenerator
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from io import StringIO
import BLASTWriter
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
    new_handle = Entrez.efetch(db=search_database, id=subjct.get(), rettype="gbwithparts", retmode="text")
    out_handle = open(os.path.join(ana_dir, os.path.join(nrps_dir, os.path.join(dat_dir, os.path.join(gen_dir,
                     subjct.get() + ".gbk")))), "w")
    out_handle.write(new_handle.read())
    out_handle.close()
    new_handle.close()
    sbjct = SeqRecord(SeqIO.read(os.path.join(ana_dir, os.path.join(nrps_dir, os.path.join(dat_dir, os.path.join(
                gen_dir, subjct.get() + ".gbk")))), format="gb"), id="sub")
    for record in records:
        recrd = SeqIO.read(os.path.join(root_dir, record[:len(record) - 4]) + ".gbk", format="gb")
        rcrd = SeqRecord(recrd, id="rec")
        SeqIO.write(rcrd.seq, "rec.fasta", "fasta")
        SeqIO.write(sbjct.seq, "sub.fasta", "fasta")
        result_handle = NcbiblastnCommandline(query="rec.fasta", subject="sub.fasta", outfmt=5)()[0]
        os.remove("sub.fasta")
        os.remove("rec.fasta")
        blast_result_record = NCBIXML.read(StringIO(result_handle))
        # Print some information on the result
        i = 0
        for alignment in blast_result_record.alignments:
            for hsp in alignment.hsps:
                BLASTWriter.create_blast_files(open(os.path.join(ana_dir, os.path.join(nrps_dir, os.path.join(dat_dir,
                                                   os.path.join(dat_stan_dir, subjct.get())))), 'w'),
                                                       BLASTWriter.write_blast_standard(i, alignment, hsp))
                i += 1
                return open(os.path.join(ana_dir, os.path.join(nrps_dir, os.path.join(dat_dir,
                                                   os.path.join(dat_stan_dir, subjct.get())))), 'r')