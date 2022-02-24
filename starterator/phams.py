from database import DB, get_db
from phamgene import new_PhamGene
from Bio import AlignIO
from Bio import SeqIO
from collections import Counter
import utils
import subprocess
import os
from utils import StarteratorError
import json


def get_pham_number(phage_name, gene_number):
    try:
        pham_no = get_db().get("SELECT PhamID FROM gene WHERE PhageID = '%s' "
                               "AND Name = %s" % (phage_name, str(gene_number)))
        return str(pham_no)
    except:
        raise StarteratorError("Gene %s of Phage %s not found in database!" % (gene_number, phage_name))


def get_pham_colors(phams=None):
    results = get_db().query("SELECT PhamID, Color from pham;")
    pham_colors = {}
    if phams:
        for row in results:
            if row[0] in phams:
                pham_colors[str(row[0])] = row[1]
    else:
        for row in results:
            pham_colors[str(row[0])] = row[1]
    return pham_colors


def get_version():
    version = get_db().get("SELECT Version from version;")
    return int(version)


class Pham(object):
    def __init__(self, pham_no, genes=None):
        self.pham_no = pham_no
        self.stats = {}
        self.file = ""
        self.count = 0
        self.genes = self.get_genes()
        self.color = self.get_color()
        if genes:
            for gene in genes:
                self.add(gene)
            whole = "All" if len(genes) > 1 else "One"
            self.file = "%s%s" % (genes[0].phage_name, whole)
        self.aligner = None

    def get_genes(self):
        """
            Get the genes of the Phamily
        """
        results = get_db().query("SELECT GeneID, PhageID, Start, Stop, "
                                 "Orientation, Name FROM gene WHERE "
                                 "PhamID = %s;" % str(self.pham_no))
        genes = {}
        self.count = len(results)
        for geneid, phageid, start, stop, orientation, name in results:
            gene = new_PhamGene(geneid, start, stop, orientation, phageid, name)
            # Data validation:
            # ignore phages with 'N' nucleotides
            # ignore genes lacking valid starts
            seq = get_db().get("SELECT Sequence FROM phage WHERE "
                               "PhageID = '%s';" % phageid)

            if "N" not in seq and gene.has_valid_start():
                genes[gene.gene_id] = gene

        if len(genes) < 1:
            raise StarteratorError("Pham Number %s not found or all genes fail validation!" % self.pham_no)

        return genes

    def get_phage_genes(self):
        pass

    def add(self, gene):
        """
            Add an unphameratored gene to the pham
        """
        self.genes[gene.gene_id] = gene

    def get_color(self):
        """
            Get the color of the phamily from the database
        """
        try:
            color = get_db().get("SELECT Color FROM pham WHERE PhamID = %s;" %
                                 str(self.pham_no))
            return color
        except:
            raise StarteratorError("Pham number %s not found in database!" % self.pham_no)

    def add_alignment(self, alignment):
        """
            Using the alignment, add the alignment to the each gene in the pham
        """ 
        for record in alignment:
            gene = self.genes[record.id]
            gene.alignment = record
            gene.add_alignment_start_site()
            gene.add_alignment_candidate_starts()
            gene.add_gaps_as_features()

    def call_clustal(self, fasta_file):
        # self.aligner = 'ClustalO'
        self.aligner = 'ClustalW'

        if self.aligner == 'ClustalO':
            outfile = fasta_file.replace(".fasta", ".aln")
            subprocess.check_call(['clustalo', '--infile=%s' % fasta_file, '--outfile=%s' % outfile, '--outfmt=clu'])
        else:
            subprocess.check_call(['clustalw', '-infile=%s' % (fasta_file), '-quicktree'])

        aln_file = fasta_file.replace(".fasta", ".aln")
        alignment = AlignIO.read(aln_file, "clustal")
        return alignment

    def make_fasta(self, file_name=None):
        if file_name is None:
            file_name = os.path.join(utils.INTERMEDIATE_DIR, "%sPham%s" % (self.file, self.pham_no))
        genes = [gene.sequence for gene in self.genes.values()]
        count = SeqIO.write(genes, "%s.fasta" % file_name, "fasta")

    def align(self):
        """
            Makes a fasta file of the genes in the Pham
            if the alignment already exists, uses that .aln as the alignment
            Otherwise, calls Clustalw from the command line and creates alignment
        """
        # files?
        file_name = os.path.join(utils.INTERMEDIATE_DIR, "%sPham%s" % (self.file, self.pham_no))
        genes = [gene.sequence for gene in self.genes.values()]
        count = SeqIO.write(genes, "%s.fasta" % file_name, "fasta")
        if len(self.genes) == 1:
            alignment = [gene.sequence]
        else:
            try:
                alignment = AlignIO.read(file_name + ".aln", "clustal")
            except:
                # cline =  ClustalwCommandline("clustalw", infile=("%s.fasta" % file_name))
                # cline()
                alignment = self.call_clustal(file_name + ".fasta")
                # alignment = AlignIO.read(file_name+".aln", "clustal")
        self.add_alignment(alignment)

    def add_total_possible_starts(self):
        """ Returns a list of all the candidate starts from the alignment
        """
        self.total_possible_starts = []
        for gene in self.genes.values():
            for site in gene.alignment_candidate_starts:
                if site not in self.total_possible_starts:
                    self.total_possible_starts.append(site) 
        self.total_possible_starts = sorted(self.total_possible_starts)
        return self.total_possible_starts

    def add_alignment_stats_to_phamgenes(self):
        for gene in self.genes.values():
            gene.add_alignment_start_stats(self)
        return

    def group_similar_genes(self, start_with=None):
        """
            Groups genes that have the same called start site, the same candidate starts
            and the same alignment (gaps are the same) together
        start_with: phage to be first item of first list
        """
        groups = []
        i = 0
        genes = self.genes.values()

        grouped = [False for gene in genes]
        while i < len(self.genes):
            if not grouped[i]:
                # gene is not in a group yet
                gene = genes[i]
                j = i + 1  # genes before index i have been grouped
                group = []
                group.append(gene)  # add gene to this group - first one
                while j < len(self.genes):  # see if other genes are similar
                    if not grouped[j]:   # skip genes that have already been grouped
                        gene_2 = genes[j]
                        if gene.is_equal(gene_2):  # if similar, then add to the group
                            grouped[i] = True 
                            grouped[j] = True
                            group.append(gene_2) 
                    j += 1
                groups.append(group)
            i += 1

        if start_with:
            for i, gene_list in enumerate(groups):
                if len(gene_list) == 1:
                    if start_with.lower() == gene_list[0].phage_name.lower():
                        split_gene_lists_on = i
                else:
                    for j, gene in enumerate(gene_list):
                        if start_with.lower() == gene.phage_name.lower():
                            split_gene_lists_on = i
                            groups[i] = groups[i][j:] + groups[i][:j]
            groups = groups[split_gene_lists_on:] + groups[:split_gene_lists_on]

            if groups[0][0].subcluster != 'Unassigned':
                sort_from = groups[0][0].cluster_hash
            else:
                sort_from = 1

            if len(groups) > 1:
                remaining = groups[1:]
                remaining.sort(key=lambda x: abs(x[0].cluster_hash-sort_from))
                groups[1:] = remaining

        else:
            groups.sort(key=lambda x: x[0].subcluster)

        return groups

    def find_most_common_start(self, ignore_draft=False):
        """
            From the total candidate strats of each gene in the pham and all the start
            called in each gene, finds the start that is most commonly called.
            Returns a dictionary containing:
                "most_called" : a list of genes currently call the "most common start"
                "most_not_called" : a list of genes that have the "most common start" but do not call it
                "no_most_called" : a list of genes that do no have the "most called start"
                "possible" : a list containing lists of genes with the start of the index of self.total_possible_starts
                "called_start: a list containing lists of genes with the called start of the index of self.total_possible_starts

            Also, for each gene in the pham, a suggested start is given, gene.suggested_start["most_commom"]
            For genes that have the most common start called (or not) a tuple containing the index of the
            most common start and the coordinate of the sequence is given.
            For genes that do not have the most common start, a list of all possible starts, containing the index
            (useful when looking at the graphical output), and the coordinate is given.
        """
        # TODO:
        # add functionality for ignoring DRAFT phages?
        # use term Called_start for all genes irrespective of method to determine location of start codon
        # use term Annotated_start for genes in which manual annotation was used to determine start codon
        # use term predicted_start for gene in which computational prediction was used to determine start codon
        all_start_sites = [gene.alignment_start_site for gene in self.genes.values()]
        all_annotated_start_sites = [gene.alignment_start_site for gene in self.genes.values() if not gene.draftStatus]
        all_predicted_start_sites = [gene.alignment_start_site for gene in self.genes.values() if gene.draftStatus]

        all_start_sites_set = set([gene.alignment_start_site for gene in self.genes.values()])
        start_stats = {}
        # creates two lists each containing a list of gene ids
        # for each candidate start of the pham:
        # start_stats["possible"] contains a list of genes with the candidate starts
        # for the index of each start in the pham
        # start_stats["called_starts"] contains of list of the genes that have the site
        #   of the index called as their start
        start_stats["possible"] = {}
        start_stats["called_starts"] = {}
        # start_stats["most_called"] = {}
        self.add_total_possible_starts()
        for i, site in enumerate(self.total_possible_starts):
            start_stats["possible"][i+1] = []
            # start_stats["most_called"][i+1] = []
            start_stats["called_starts"][i+1] = []
            for gene in self.genes.values():
                if site in gene.alignment_candidate_starts:
                    start_stats["possible"][i+1].append(gene.full_name)
                if site == gene.alignment_start_site:
                    start_stats["called_starts"][i+1].append(gene.full_name)

        all_starts_count = Counter(all_start_sites)
        all_annot_count = Counter(all_annotated_start_sites)
        all_predicted_count = Counter(all_predicted_start_sites)

        called_starts_count = all_starts_count.most_common()
        annot_starts_count = all_annot_count.most_common()
        predicted_starts_count = all_predicted_count.most_common()

        most_called_start_index = self.total_possible_starts.index(called_starts_count[0][0])+1
        if len(annot_starts_count) > 0:  # i.e. at least 1 annotated gene
            most_annot_start_index = self.total_possible_starts.index(annot_starts_count[0][0])+1
        else:
            most_annot_start_index = None

        genes_start_most_called = start_stats["called_starts"][most_called_start_index]
        start_stats["most_called_start"] = most_called_start_index
        start_stats["most_annotated_start"] = most_annot_start_index

        if most_annot_start_index is not None:
            genes_start_most_annot = start_stats["called_starts"][most_annot_start_index]
        else:
            genes_start_most_annot = None

        # start_stats["most_called"] = start_stats["called_starts"][most_called_start_index]
        start_stats["most_called"] = []
        start_stats["most_not_called"] = []
        start_stats["no_most_called"] = []
        start_stats["most_annotated"] = []
        start_stats["most_not_annotated"] = []
        start_stats["no_most_annot"] = []
        start_stats["annot_list"] = [g for g in self.genes.values() if not g.draftStatus]
        start_stats["draft_list"] = [g for g in self.genes.values() if g.draftStatus]

        start_stats['called_counts'] = {}
        for k, l in start_stats['called_starts'].items():
            if len(l) > 0:
                start_stats['called_counts'][k] = len(l)

        start_stats['annot_counts'] = {}
        not_drafts = [pg for pg in self.genes.values() if not pg.draftStatus]
        for gene in not_drafts:
            start_number = [key for key,val in start_stats['called_starts'].items() if gene.full_name in val]
            start_number = start_number[0]
            if start_number not in start_stats['annot_counts'].keys():
                start_stats['annot_counts'][start_number] = 0
            start_stats['annot_counts'][start_number] += 1
             #   start_stats['annot_counts'][k] = len(annotated)

        genes_without_most_called = []
        print "phams.find_most_common_start: genes_start_most_called " + str(genes_start_most_called)
        for gene in self.genes.values():
            # check if the gene even has the most called start
            if gene.full_name in start_stats["possible"][most_called_start_index]:
                if gene.full_name in genes_start_most_called:
                    if gene.orientation == 'F':   # only +1 for forward genes
                        # genes where most called start is present and it is called as the start are "most_called"
                        gene.suggested_start["most_called"] = (most_called_start_index, gene.start+1)
                    else:
                        gene.suggested_start["most_called"] = (most_called_start_index, gene.start)
                    start_stats["most_called"].append(gene.full_name)
                else:
                    # genes where most called start is present but it's not the called start are "most_not_called
                    start_stats["most_not_called"].append(gene.full_name)
                    most_called_alignment_index = self.total_possible_starts[most_called_start_index-1]
                    suggested_start = gene.alignment_index_to_coord(most_called_alignment_index)
                    gene.suggested_start["most_called"] = (most_called_start_index, suggested_start)

            else:
                # genes where the most called start is NOT even present are no_most_called
                start_stats["no_most_called"].append(gene.full_name)
                possible_starts_coords = []
                for start in gene.alignment_candidate_starts:
                    index = self.total_possible_starts.index(start) + 1
                    new_start = gene.alignment_index_to_coord(start) + 1
                    possible_starts_coords.append((index, new_start))
                gene.suggested_start["most_called"] = possible_starts_coords

            if most_annot_start_index is not None:
                if gene.full_name in start_stats["possible"][most_annot_start_index]:
                    if gene.full_name in genes_start_most_annot:
                        # code below used for deprecated "suggested starts" list
                        # if gene.orientation == 'F':  # only +1 for forward genes
                        #     # genes where most annotated start is present and it is called as the start are "most_annotated"
                        #     gene.suggested_start["most_annotated"] = (most_annot_start_index, gene.start + 1)
                        # else:
                        #     gene.suggested_start["most_annotated"] = (most_annot_start_index, gene.start)
                        start_stats["most_annotated"].append(gene.full_name)
                    else:
                        # genes where most annotated start is present but it's not the called start are "most_not_annotated"
                        start_stats["most_not_annotated"].append(gene.full_name)
                        # code below used for deprecated "suggested starts" list
                        most_annot_alignment_index = self.total_possible_starts[most_annot_start_index - 1]
                        suggested_start = gene.alignment_index_to_coord(
                            most_annot_alignment_index)  # +1 issue dealt with in function
                        gene.suggested_start["most_called"] = (most_annot_start_index, suggested_start)

                else:
                    # genes where the most annotated start is NOT even present are no_most_annot
                    start_stats["no_most_annot"].append(gene.full_name)
                    # Code below used for deprecated "suggested starts" list
                    # possible_starts_coords = []
                    # for start in gene.alignment_candidate_starts:
                    #     index = self.total_possible_starts.index(start) + 1
                    #     new_start = gene.alignment_index_to_coord(start) + 1
                    #     possible_starts_coords.append((index, new_start))
                    # gene.suggested_start["most_called"] = possible_starts_coords

            # section to add summary of annotations based only on set of starts found in the gene
            # start by looking through all possible starts in this particular gene and see if
            # there are any annotations that call that start

            alignment_start_coord_with_annotations = []
            for start in gene.alignment_candidate_starts:
                if start in all_annotated_start_sites:
                    alignment_start_coord_with_annotations.append(start)

            gene.suggested_start["alignment_start_coord_with_annotations"] = alignment_start_coord_with_annotations

            alignment_start_indices_with_annotations = []
            for start in alignment_start_coord_with_annotations:
                alignment_start_indices_with_annotations.append(self.total_possible_starts.index(start) + 1)

            gene.suggested_start["alignment_start_indices_with_annotations"] = alignment_start_indices_with_annotations

            alignment_start_counts_with_annotations = []
            for annotated_start in alignment_start_coord_with_annotations:
                for start, count in annot_starts_count:
                    if start == annotated_start:
                        alignment_start_counts_with_annotations.append(count)

            gene.suggested_start["alignment_start_counts_with_annotations"] = alignment_start_counts_with_annotations
            try:  # this happens when annotated start of gene is not one of the three typical start codons (ATG,GTG,TTG)
                gene.suggested_start["current_start_number"] = self.total_possible_starts.index(gene.alignment_start_site) + 1
            except:
                gene.suggested_start["current_start_number"] = None

        self.stats["most_common"] = start_stats

        # now update genes based on start analysis
        self.add_alignment_stats_to_phamgenes()
        self.add_cluster_stats(start_stats)
        return start_stats

    def annot_summary(self):

        summary_dict = {}
        summary_dict['Name'] = self.pham_no
        summary_dict['MemberCount'] = self.count
        summary_dict['AnnotCount'] = len(self.stats['most_common']['annot_list'])
        summary_dict['TotalStarts'] = len(self.total_possible_starts)
        summary_dict['DbVersion'] = get_version()
        summary_dict['Aligner'] = self.aligner

        genelist = []
        for gene in self.genes.values():
            gene_dict = {}
            gene_dict['GeneID'] = gene.gene_id
            gene_dict['Start'] = gene.start
            gene_dict['Stop'] = gene.stop
            gene_dict['Orientation'] = gene.orientation
            gene_dict['AvailableStarts'] = gene.alignment_candidate_start_nums
            gene_dict['AvailableCoord'] = [gene.alignment_index_to_coord(s) for s in gene.alignment_candidate_starts]
            if gene.locustag is None:
                gene.get_locustag()
            if gene.locustag is not "":
                gene_dict['locustag'] = gene.locustag
            if gene.annot_author == 1 and gene.status == 'final': # Pitt "owns" it and is in genbank already
                gene_dict['Editable'] = "True"

            genelist.append(gene_dict)

        summary_dict['Genes'] = genelist
        annotlist = {}
        annotlist['Starts'] = self.stats['most_common']['annot_counts'].keys()
        annotlist['Counts'] = self.stats['most_common']['annot_counts'].values()
        summary_dict['Annots'] = annotlist

        conservationdict = {}
        for i in range(0, summary_dict['TotalStarts']):
            conservationdict[i+1] = float(len(self.stats["most_common"]['possible'][i+1]))/self.count
        summary_dict['Conservation'] = conservationdict

        return summary_dict

    def export_json(self, filename):
        blob = self.annot_summary()
        with open(filename, "w") as outfile:
            json.dump(blob, outfile, indent=4, ensure_ascii=False)

    def add_cluster_stats(self, self_stats):
        clusters_present = set()
        genes_by_cluster = {}
        for g, pg in self.genes.items():
            clusters_present.add(pg.subcluster)
            if pg.subcluster in genes_by_cluster.keys():
                genes_by_cluster[pg.subcluster].append(g)
            else:
                genes_by_cluster[pg.subcluster] = [g]

        self_stats['clusters_present'] = [str(t) for t in clusters_present]
        self_stats['genes_by_cluster'] = genes_by_cluster
