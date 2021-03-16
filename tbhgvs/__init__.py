import sys
import argparse
from uuid import uuid4
import subprocess as sp
import re
import bisect
import os


module_path = os.path.dirname(os.path.realpath(__file__))


class gene_class:
    def __init__(self,name,gene_id,strand,start,end,feature_start,feature_end,length):
        self.name = name
        self.id = gene_id
        self.strand = strand
        self.start = start
        self.end = end
        self.feature_start = feature_start
        self.feature_end = feature_end
        self.length = length
    def __str__(self):
        return str(vars(self))

class reference_db:
    def __init__(self,directory=None,fasta=None,gff=None):
        if not directory and not fasta and not gff:
            self.fasta = "%s/.tbhgvs/genome.fasta" % os.path.expanduser("~")
            self.gff = "%s/.tbhgvs/genome.gff" % os.path.expanduser("~")
        elif directory:
            self.fasta = "%s/.tbhgvs/genome.fasta" % directory
            self.gff = "%s/.tbhgvs/genome.gff" % directory
        else:
            self.fasta = fasta
            self.gff = gff

        self.genes = load_gff(self.gff)
        self.gene_dict = {g.name:g for g in self.genes}
        self.gene_dict.update({g.id:g for g in self.genes})
    def genome2hgvs(self,variant):
        tmp_vcf = "%s.vcf" % uuid4()
        write_vcf(variant,tmp_vcf)
        hgvs =  vcf2hgvs(tmp_vcf,self.fasta,self.gff,self.genes,self.gene_dict)
        os.remove(tmp_vcf)
        return hgvs
    def hgvs2genome(self,variant,gene_name):
        gene = self.gene_dict[gene_name]
        return hgvs2genome(variant,gene,self.fasta)

def download_files(directory=None):
    import urllib.request
    if not directory:
        directory = "%s/.tbhgvs/" % os.path.expanduser("~")
    if not os.path.isdir(directory):
        os.mkdir(directory)

    urllib.request.urlretrieve('https://raw.githubusercontent.com/jodyphelan/tbdb/master/genome.fasta','%s/genome.fasta' % directory )
    urllib.request.urlretrieve('https://raw.githubusercontent.com/jodyphelan/tbdb/master/genome.gff','%s/genome.gff' % directory )

def load_gff(gff):
    genes = []
    for l in open(gff):
        if l[0]=="#": continue
        fields = l.rstrip().split()
        if fields[2]!="gene" and fields[2]!="rRNA_gene": continue
        strand = fields[6]
        p1 = int(fields[3])
        p2 = int(fields[4])
        gene_length = p2-p1+1
        re_obj = re.search(r"Name=([a-zA-Z0-9\.\-\_]+)",l)
        gene_name = re_obj.group(1) if re_obj else None
        re_obj = re.search(r"ID=gene:([a-zA-Z0-9\.\-\_]+)",l)
        gene_id = re_obj.group(1) if re_obj else None
        feature_start,feature_end = p1,p2
        start = p1 if strand=="+" else p2
        end =  p2 if strand=="+" else p1
        tmp = gene_class(gene_name,gene_id,strand,start,end,feature_start,feature_end,gene_length)
        genes.append(tmp)
    return genes


a2aaa = {
'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln',
'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp',
'Y': 'Tyr', 'V': 'Val', '*': '*', '-': '-'
}

aaa2codon = {
    'Phe': ['TTT', 'TTC'],
    'Leu': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'Ile': ['ATT', 'ATC', 'ATA'],
    'Met': ['ATG'],
    'Val': ['GTT', 'GTC', 'GTA', 'GTG'],
    'Ser': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'Pro': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Thr': ['ACT', 'ACC', 'ACA', 'ACG'],
    'Ala': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Tyr': ['TAT', 'TAC'],
    '*': ['TAA', 'TAG', 'TGA'],
    'His': ['CAT', 'CAC'],
    'Gln': ['CAA', 'CAG'],
    'Asn': ['AAT', 'AAC'],
    'Lys': ['AAA', 'AAG'],
    'Asp': ['GAT', 'GAC'],
    'Glu': ['GAA', 'GAG'],
    'Cys': ['TGT', 'TGC'],
    'Trp': ['TGG'],
    'Arg': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'Gly': ['GGT', 'GGC', 'GGA', 'GGG']
}

def write_vcf(variants,filename):
    with open(filename,"w") as O:
        O.write('##fileformat=VCFv4.2\n')
        O.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        O.write('##contig=<ID=Chromosome,length=4411532>\n')
        O.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest\n')
        variants = sorted(variants, key=lambda x:x["pos"])
        for var in variants:
            O.write("Chromosome\t%(pos)s\t.\t%(ref)s\t%(alt)s\t255\t.\t.\tGT\t1\n" % var)

def find_closest_genes(position,genes):
    left_gene = genes[bisect.bisect_left([g.feature_end for g in genes],position)-1]
    i = bisect.bisect_right([g.feature_start for g in genes],position)
    right_gene = genes[i] if len(genes)>i else None
    return (left_gene,right_gene)

def codon_pos2genome_pos(gene,codon):
    gene_positions = [codon*3-2,codon*3-1,codon*3]
    if gene.strand=="+":
        genome_positions = [x + gene.feature_start - 1 for x in gene_positions]
    else:
        genome_positions = [gene.feature_end - x + 1 for x in gene_positions]
    return genome_positions

def genome2gene_pos(pos,gene):
    if gene.strand=="+":
        return pos - gene.feature_start + 1
    else:
        return gene.feature_end - pos + 1

def gene_pos2genome(gene_pos,gene):
    if gene.strand=="+":
        if gene_pos>0:
            return gene_pos + gene.feature_start - 1
        else:
            return gene_pos + gene.feature_start
    else:
        if gene_pos>0:
            return gene.feature_end - gene_pos + 1
        else:
            return gene.feature_end - gene_pos

def revcom(s):
    """Return reverse complement of a sequence"""
    def complement(s):
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        letters = list(s)
        letters = [basecomplement[base] for base in letters]
        return ''.join(letters)
    return complement(s[::-1])


def vcf2hgvs(infile,ref,gff,genes,gene_dict,promoter_offset=50):
    variants = []
    cmd = sp.Popen("bcftools csq -f %s -g %s %s | bcftools query -f '%%POS\t%%REF\t%%ALT\t%%BCSQ\n'" % (ref,gff,infile),shell=True,stdout=sp.PIPE,stderr=sp.PIPE)
    stderr = cmd.stderr.read().decode()
    if "Error" in stderr:
        quit(stderr)
    for l in cmd.stdout:
        row = l.decode().strip().split()
        pos = int(row[0])
        ref = row[1]
        alt = row[2]
        csqs = [x.split("|") for x in row[3].split(",")]
        for csq in csqs:
            if csq[0][0]=="@":
                pass
            elif csq[0]=="missense" or csq[0]=="stop_gained" or csq[0]=="stop_lost":
                # ['missense', 'gyrB', 'Rv0005', 'protein_coding', '+', '301V>301W', '6140G>T+6141T>G']
                r = re.search("(\\d+)([A-Z\\*])>\\d+([A-Z\\*])",csq[5])
                variants.append({
                    "nucleotide_change":csq[6],
                    "type": csq[0],
                    "gene_name": csq[1],
                    "gene_id": csq[2],
                    "hgvs": "p.%s%s%s" % (a2aaa[r.group(2)],r.group(1),a2aaa[r.group(3)])
                })
            elif csq[0]==".":
                left_gene,right_gene = find_closest_genes(pos,genes)
                if left_gene and left_gene.strand=="-" and (left_gene.feature_end+promoter_offset)>pos:
                    offset = left_gene.feature_end - pos
                    variants.append({
                        "nucleotide_change":"%s%s>%s" % (pos,revcom(ref),revcom(alt)),
                        "type": "promoter",
                        "gene_name": left_gene.name,
                        "gene_id": left_gene.id,
                        "hgvs": "c.%s%s>%s" % (offset,revcom(ref),revcom(alt))
                    })
                elif right_gene and right_gene.strand=="+" and (right_gene.feature_start-promoter_offset)<pos:
                    offset = pos-right_gene.feature_start
                    variants.append({
                        "nucleotide_change":"%s%s>%s" % (pos,ref,alt),
                        "type": "promoter",
                        "gene_name": right_gene.name,
                        "gene_id": right_gene.id,
                        "hgvs": "c.%s%s>%s" % (offset,ref,alt)
                    })
                else:
                    variants.append({
                        "nucleotide_change":"%s%s>%s" % (pos,ref,alt),
                        "type": "intergenic",
                        "gene_name": None,
                        "gene_id": None,
                        "hgvs": "g.%s%s>%s" % (pos,ref,alt)
                    })
            elif csq[0]=="frameshift" or csq[0]=="inframe_deletion" or csq[0]=="inframe_insertion":
                gene = gene_dict[csq[1]]
                if len(ref)>1:
                    ### Deletion
                    del_size = len(ref) - len(alt)
                    if gene.strand=="+":
                        del_start = genome2gene_pos(pos + 1,gene )
                        del_end = genome2gene_pos(pos+del_size,gene )
                    else:
                        del_start = genome2gene_pos(pos+ del_size,gene )
                        del_end = genome2gene_pos(pos+ 1,gene )
                    variants.append({
                        "nucleotide_change":csq[6],
                        "type": csq[0],
                        "gene_name": csq[1],
                        "gene_id": csq[2],
                        "hgvs": "c.%s_%sdel" % (del_start,del_end)
                    })
                else:
                    ### Insertion
                    if gene.strand=="+":
                        ins_seq = alt[1:]
                        ins_start = genome2gene_pos(pos,gene)
                    else:
                        ins_seq = revcom(alt[1:])
                        ins_start = genome2gene_pos(pos+1,gene)

                    variants.append({
                        "nucleotide_change":csq[6],
                        "type": csq[0],
                        "gene_name": csq[1],
                        "gene_id": csq[2],
                        "hgvs": "c.%s_%sins%s" % (ins_start,ins_start+1,ins_seq)
                    })
            elif csq[0]=="non_coding" and csq[3]=="rRNA":
                gene = gene_dict[csq[1]]
                gene_pos = genome2gene_pos(pos,gene)
                variants.append({
                    "nucleotide_change":"%s%s>%s" % (pos,ref,alt),
                    "type": csq[0],
                    "gene_name": csq[1],
                    "gene_id": csq[1],
                    "hgvs": "r.%s%s>%s" % (gene_pos,ref.lower(),alt.lower())
                })
            elif csq[0]=="start_lost":
                gene = gene_dict[csq[1]]
                gene_pos = genome2gene_pos(pos,gene)

                variants.append({
                    "nucleotide_change":"%s%s>%s" % (pos,ref,alt),
                    "type": csq[0],
                    "gene_name": csq[1],
                    "gene_id": csq[1],
                    "hgvs": "c.%s%s>%s" % (gene_pos,ref,alt)
                })
            elif csq[0]=="synonymous" or csq[0]=="stop_retained":
                gene = gene_dict[csq[1]]
                gene_pos = genome2gene_pos(pos,gene)
                if gene.strand=="-":
                    ref_nt = revcom(ref)
                    alt_nt = revcom(alt)
                else:
                    ref_nt = ref
                    alt_nt = alt
                variants.append({
                    "nucleotide_change":"%s%s>%s" % (pos,ref,alt),
                    "type": csq[0],
                    "gene_name": csq[1],
                    "gene_id": csq[1],
                    "hgvs": "c.%s%s>%s" % (gene_pos,ref_nt,alt_nt)
                })
            else:
                quit("Error: cannot parse %s\n" % csq)
    return variants

# def genome2hgvs(variants):
#     tmp_vcf = "%s.vcf" % uuid4()
#     write_vcf(variants,tmp_vcf)
#     fasta = "%s/MTB-h37rv_asm19595v2-eg18.fa" % (module_path)
#     gff = "%s/MTB-h37rv_asm19595v2-eg18.gff" % (module_path)
#     hgvs =  vcf2hgvs(tmp_vcf,fasta,gff)
#     os.remove(tmp_vcf)
#     return hgvs

def fetch_seq(genome,start,end):
    if start>end:
        start,end = end,start
    tmp = []
    for l in sp.Popen("samtools faidx %s Chromosome:%s-%s" % (genome,start,end),shell=True,stdout=sp.PIPE).stdout:
        line = l.decode().strip()
        if line[0]==">":
            continue
        else:
            tmp.append(line)
    return "".join(tmp)

def get_possible_codon_changes(codon,alt_aaa,genome_coords):
    possible_changes = []
    for ref_codon,alt_codon in zip([codon for _ in aaa2codon[alt_aaa]],aaa2codon[alt_aaa]):
        changes = []
        for i in range(3):
            if ref_codon[i]==alt_codon[i]: continue
            changes.append({"pos":genome_coords[i],"ref":ref_codon[i],"alt":alt_codon[i]})
        possible_changes.append(changes)
    return possible_changes


def hgvs2genome(var,gene,refgenome):
    if var[0]=="p":
        #p.Thr40Ile
        re_obj = re.match(r"p.([A-Za-z]+)([0-9]+)([A-Za-z\*]+)",var)
        codon_pos = int(re_obj.group(2))
        alt_aaa = re_obj.group(3)
        genome_positions = codon_pos2genome_pos(gene,codon_pos)
        codon = fetch_seq(refgenome,genome_positions[0],genome_positions[2])
        if gene.strand=="-":
            codon = revcom(codon)
        possible_changes = get_possible_codon_changes(codon,alt_aaa,genome_positions)

        return (possible_changes)
    elif var[0]=="c":
        if "ins" in var:
            #c.192_193insG
            re_obj = re.match("c.([0-9]+)_([0-9]+)ins([A-Za-z]+)",var)
            if gene.strand=="+":
                gene_pos = int(re_obj.group(1))
                ins_seq = re_obj.group(3)
            else:
                gene_pos = int(re_obj.group(2))
                ins_seq = revcom(re_obj.group(3))

            genome_pos = gene_pos2genome(gene_pos,gene)
            ref = fetch_seq(refgenome,genome_pos,genome_pos)

            return [{"pos":genome_pos,"ref":ref,"alt":ref+ins_seq}]

        elif "del" in var:
            if "_" in var:
                #c.785_787del
                re_obj = re.match(r"c.(\-*[0-9]+)_(\-*[0-9]+)del",var)
                genome_pos_from = gene_pos2genome(int(re_obj.group(1)),gene)
                genome_pos_to = gene_pos2genome(int(re_obj.group(2)),gene)
            else:
                #c.884del
                re_obj = re.match(r"c.(\-*[0-9]+)del",var)
                genome_pos_from = gene_pos2genome(int(re_obj.group(1)),gene)
                genome_pos_to = genome_pos_from

            if gene.strand=="+":
                vcf_genome_pos = genome_pos_from - 1
                del_seq = fetch_seq(refgenome, vcf_genome_pos,genome_pos_to)
            else:
                vcf_genome_pos = genome_pos_to - 1
                del_seq = fetch_seq(refgenome, vcf_genome_pos,genome_pos_from)

            return([{"pos":vcf_genome_pos,"ref":del_seq,"alt":del_seq[0]}])

        else:
            #c.-16C>T
            re_obj = re.match("c.(-[0-9]+)([A-Z])>([A-Z])", var)
            gene_pos = int(re_obj.group(1))
            genome_pos = gene_pos2genome(gene_pos,gene)
            return [{"pos":genome_pos,"ref":re_obj.group(2),"alt":re_obj.group(3)}]

    elif var[0]=="r":
        #rrl r.2814g>t
        re_obj = re.match("r.([0-9]+)([A-Za-z]+)>([A-Za-z]+)",var)
        gene_pos = int(re_obj.group(1))
        genome_pos = gene_pos2genome(gene_pos,gene)
        return [{"pos":genome_pos,"ref":re_obj.group(2).upper(),"alt":re_obj.group(3).upper()}]
