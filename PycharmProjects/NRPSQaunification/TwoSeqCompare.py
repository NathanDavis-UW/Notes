import csvGenerator
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from io import StringIO
import BLASTWriter
import os
import BLASTMatrix

search_database = "nucleotide"
ma_dir = "Matrix"
ma_gen_dir = "MatrixGenBank"
ma_stan_dir = "MatrixStandard"
alt_dir = "Alternative Analysis"
root_dir = "NRPSRoot"
dat_dir = "NRPSData"
alt_root_dir = "Alternative NRPSRoot"
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
        if record is not None:
            recrd = SeqIO.read(os.path.join(root_dir, record + ".gbk"), format="gb")
            rcrd = SeqRecord(recrd, id="rec")
            SeqIO.write(rcrd.seq, "rec.fasta", "fasta")
            SeqIO.write(sbjct.seq, "sub.fasta", "fasta")
            result_handle = NcbiblastnCommandline(query="rec.fasta", subject="sub.fasta", outfmt=5)()[0]
            os.remove("sub.fasta")
            os.remove("rec.fasta")
            blast_result_record = NCBIXML.read(StringIO(result_handle))
            # Print some information on the result
            i = 0
            csvGenerator.create_dir(os.path.join(dat_dir, os.path.join(dat_stan_dir, subjct.get())), nrps_dir)
            k = 0

            for alignment in blast_result_record.alignments:
                csvGenerator.create_dir(os.path.join(dat_dir, os.path.join(dat_stan_dir, os.path.join(subjct.get(),
                                            alignment.title + " - " + str(k)))),
                                                nrps_dir)
                j = 0
                for hsp in alignment.hsps:
                    BLASTWriter.create_blast_files(open(os.path.join(ana_dir, os.path.join(nrps_dir, os.path.join(
                                                        dat_dir, os.path.join(dat_stan_dir, os.path.join(subjct.get(),
                                                            os.path.join(alignment.title + " - " + str(k),
                                                                alignment.title + " - " + str(k) +
                                                                    " - hsp" + str(j))))))), 'w'),
                                                                        BLASTWriter.write_blast_standard(i, alignment,
                                                                            hsp))
                j += 1
                i += 1

            k += 1
            return os.path.join(ana_dir, os.path.join(nrps_dir, os.path.join(dat_dir,
                                                                              os.path.join(dat_stan_dir,
                                                                                           subjct.get()))))
def alt_compare(record, subject):
    print(str(record) + ":" + str(subject))
    BLASTMatrix.simple_dir(ma_dir)
    BLASTMatrix.create_dir(ma_gen_dir, ma_dir)
    BLASTMatrix.create_dir(ma_stan_dir, ma_dir)
    sbject = SeqIO.read(os.path.join(alt_root_dir, subject), format="gb")
    sbjct = SeqRecord(sbject, id="sub")
    recrd = SeqIO.read(os.path.join(alt_root_dir, record), format="gb")
    rcrd = SeqRecord(recrd, id="rec")
    SeqIO.write(rcrd.seq, "rec.fasta", "fasta")
    SeqIO.write(sbjct.seq, "sub.fasta", "fasta")
    result_handle = NcbiblastnCommandline(query="rec.fasta", subject="sub.fasta", outfmt=5)()[0]
    blast_result_record = NCBIXML.read(StringIO(result_handle))
    BLASTMatrix.create_dir(os.path.join(ma_stan_dir, record), ma_dir)
    if len(blast_result_record.alignments) > 0:
        alignment = blast_result_record.alignments[0]
        BLASTMatrix.create_dir(os.path.join(ma_stan_dir, os.path.join(record, alignment.title)), ma_dir)
        e_list = []
        for hsp in alignment.hsps:
            e_list.append(hsp.align_length-hsp.score)
        os.remove("sub.fasta")
        os.remove("rec.fasta")
        sum = 0
        for e in e_list:
            sum += e
        return int((sum/len(e_list))/1000)
    else:
        sum = 0
        f = open("sub.fasta")
        sub_list = f.readlines()
        for line in sub_list[1:]:
            sum += len(line)
        return int(sum/1000)