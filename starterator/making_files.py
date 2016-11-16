#!/usr/bin/env python

# Copyright (c) 2013, 2014 All Right Reserved, Hatfull Lab, University of Pittsburgh
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
# PARTICULAR PURPOSE.  USE AT YOUR OWN RISK.
#
# Marissa Pacey
# April 4, 2014
# Functions that create the PDF outputsst

import cPickle
import argparse
import time
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors
import reportlab.lib.pagesizes
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
import PyPDF2
from Bio import SeqIO
import math
import StringIO
import utils
import phams
import phamgene
import os
# from phage import



def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--pham_no', default = -1)
    parser.add_argument('-p' , '--phage', default=None)
    parser.add_argument('-d', '--dir', default=None)
    parser.add_argument('-a', '--one_or_all', default='One')
    parser.add_argument('-f', '--pickle_file', default=None)
    parser.add_argument("-l", "--phage_length", type=int, default=-1)
    parser.add_argument('-m', '--make', default=None)

    return parser.parse_args()

def check_file(file_name):
    try: 
        f = open(file_name, "r")
        return True
    except:
        return False

def output_suggested_starts(pham, phage=None, all_genes=False):
    """
        gives a String list of the output of the suggested start site
        If all_genes is true, finds the suggested start site of all genes in a pham
        Otherwise, finds only the suggested start site of the specific phage gene(s)
    """
    output = []
    if all_genes == False:
        for gene in pham.genes.values():
            if phage in gene.gene_id:
                coord_suggestion = gene.suggested_start["most_called"]
                output.append(""+ gene.gene_id + ", " +str(coord_suggestion))
    else:
        for gene in pham.genes.values():
            coord_suggestion = gene.suggested_start["most_called"]
            output.append(gene.gene_id + ", " +str(coord_suggestion))
    return output


def output_start_sites(stats):
        """
            Writes a report of the start sites statistics
            Returns a list of strings with information about
            start sites for this pham  
        """
        most_called_start = stats["most_called_start"]
        most_annotated_start = stats["most_annotated_start"]
        total_genes = stats["phamCount"]
        annotatedCount = stats["annotCount"]
        draftCount = stats["draftCount"]
        calledCount = len(stats["called_starts"][most_called_start])
        output = []
        output.append("")

        output.append("Info on published annotations and draft predictions (start numbers based on diagram):")

        output.append('"Most Annotated" Start is %s, annotated in %s of %s genes.' % (str(stats["most_annotated_start"]), str(annotatedCount), str(total_genes)))
        output.append('"Most Predicted" Start is %s, predicted in %s of %s genes.' % (str(stats["most_predicted_start"]), str(draftCount), str(total_genes)))
        output.append('"Most Called" Start is %s, called in %s of %s genes.' % (str(stats["most_called_start"]), str(calledCount), str(total_genes)))
        # percent_with_most_annotated = (float(len(stats["most_called"]))
        #                             /total_genes *100 )
        #
        # output.append('Percent of genes that begin at the "Most Annotated" start: %10.1f%%'
        #                 % percent_with_most_called )
        output.append('Genes that call the "Most Annotated" start:')
        s = u'\u2022' + ''
        for gene in stats["called_starts"][most_annotated_start]:
            s += gene+ ", "
        output.append(s)
        output.append("")
        output.append('Genes that have the "Most Annotated" start but do not call it:')
        s = u'\u2022' + ''
        for gene in stats["most_not_annotated"]:
            s += gene + ", "
        output.append(s)
        output.append('')
        output.append('Genes that do not have the "Most Annotated" start:')
        s = u'\u2022' + ""
        for gene in stats["no_most_annotated"]:
            s += gene + ", "
        output.append(s + '')
        output.append('')
        output.append("Starts Called:")
        for start, genes in stats["called_starts"].items():
            if len(genes) == 0:
                continue

            s = ''
            for gene in genes:
                s += gene + ", "
            output.append(u'\u2022' + " Start number " + str(start) + ":\t" + s +'')
            percent = float(len(genes)) / total_genes * 100
            output.append("Percent with start %s called: %10.1f%% \n\t" % (str(start), percent))
        return output

def add_pham_no_title(args, pham_no, first_graph_path, i=""):
    # print i, type(i)
    # print first_graph_path
    packet = StringIO.StringIO()
    can = canvas.Canvas(packet, pagesize=reportlab.lib.pagesizes.letter)
    width, height = reportlab.lib.pagesizes.letter
    # print width, height
    can.drawString(280, 750, 'Pham ' + str(pham_no))
    can.save()

    packet.seek(0)
    new_pdf = PyPDF2.PdfFileReader(packet)
    existing_pdf = PyPDF2.PdfFileReader(file(first_graph_path), 'rb')
    output = PyPDF2.PdfFileWriter()
    print first_graph_path
    page = existing_pdf.getPage(0)
    page.mergePage(new_pdf.getPage(0))
    output.addPage(page)
    print utils.INTERMEDIATE_DIR
    print "old graph?", os.path.join(args.dir, "%sPham%sGraph%s.pdf" % (args.phage+ args.one_or_all, pham_no, i))
    outputStream = file(os.path.join(args.dir, "%sPham%sGraph%s.pdf" % (args.phage+ args.one_or_all, pham_no, i)),'wb')
    # print outputStream
    print outputStream
    os.remove(first_graph_path)
    output.write(outputStream)
    outputStream.close()

def combine_graphs(args, phage, pham_no, num_pages):
    merger = PyPDF2.PdfFileMerger()
    for i in xrange(0, num_pages+1):
        print os.path.join(args.dir, "%sPham%sGraph%d.pdf"  % (phage + args.one_or_all, pham_no, i))
        graph = open(os.path.join(args.dir, "%sPham%sGraph%d.pdf"  % (phage + args.one_or_all, pham_no, i)), "rb")
        merger.append(fileobj=graph)
    merger.write(open(os.path.join(args.dir, "%sPham%sGraph.pdf"  % (phage + args.one_or_all, pham_no)), "wb"))

def make_gene_track(gd_diagram, pham, gene_group, num_on_diagram, total):
    """"""
    colors = ['purple', 'red', 'lightblue', 'orange', 'tan', 'brown']
    gene = gene_group[0]

    #change trackname to name of fist gene in list
    track_name = gene.gene_id
    if len(gene_group) > 1:
        track_name += " + "
        track_name += str(len(gene_group)-1)
    gd_gene_track = gd_diagram.new_track(total - num_on_diagram, label=True,
                                         name=track_name, greytrack=1, greytrack_labels=1)
    gd_seq_set = gd_gene_track.new_set()
    gd_feature_set = gd_gene_track.new_set()

    start_site = gene.alignment_start_site
    start_site_feature = SeqFeature(FeatureLocation(start_site, start_site +1), 
                            strand=None)
    for feature in gene.alignment.features:
        if feature.type == 'seq':
            gd_seq_set.add_feature(feature, color='pink')
    for site in gene.alignment_candidate_starts:
        site_color = pham.total_possible_starts.index(site) % len(colors) 
        possible_site = SeqFeature(FeatureLocation(site, site ), strand=None)
        gd_feature_set.add_feature(possible_site, color=colors[site_color], 
            name=str(pham.total_possible_starts.index(site)+1), label=True)
    end_gene_feature = SeqFeature(FeatureLocation(len(gene.alignment), 
                        len(gene.alignment)+1), strand=None)

    # draw blue called start only if non-draft gene in gene group, if all draft use yellow

    allDraftStatus = True
    for gene in gene_group:
        allDraftStatus = allDraftStatus and gene.draftStatus

    if allDraftStatus:
        startcolor="yellow"
    else:
        startcolor = "green"
    gd_feature_set.add_feature(start_site_feature, color=startcolor, label=True)

    gd_feature_set.add_feature(end_gene_feature, color='purple', label=True)

def graph_start_sites(args, pham, file_path):
    """
        graphs the alignment, creates a PDF file called {Phage Name}{One or All}Pham{Pham Number}.pdf
    """

    # genes = sorted(pham.genes_in_pham.values())
    genes = pham.group_similar_genes()
    # for group in genes:
    #     print group, group.id
    if args.phage == None:
        args.phage = ""
    seq_length = len(genes[0][0].sequence.seq)
    if not args.phage:
        graph_path = os.path.join(file_path,"OnePham%sGraph_%%s.pdf" % (pham.pham_no))
    else:
        graph_path = os.path.join(file_path, "%sPham%sGraph_%%s.pdf" % (args.phage+ args.one_or_all, pham.pham_no))
    if len(genes) > 100:
        for i in xrange(0, int(math.ceil(len(genes)/50.0))):
            gd_diagram = GenomeDiagram.Diagram(pham.pham_no)
            final_graph_path = os.path.join(file_path,"OnePham%sGraph%d.pdf" % (pham.pham_no, i)) if not args.phage else  os.path.join(file_path, "%sPham%sGraph%d.pdf" % (args.phage+args.one_or_all, pham.pham_no, i))
            graph_path = os.path.join(file_path, "OnePham%sGraph_%d.pdf" % ( pham.pham_no, i)) if not args.phage else os.path.join(file_path,"%sPham%sGraph_%d.pdf" % ( args.phage+args.one_or_all, pham.pham_no, i))
            if check_file(final_graph_path):
                continue
            graph_path_svg = os.path.join(file_path, "%sPham%sGraph_%s.svg" % (args.phage+ args.one_or_all, pham.pham_no, i))

            for j in xrange(0, 50):
                if i*50 + j >= len(genes):
                    print i * 50, + j, len(genes)
                    gd_gene_track = gd_diagram.new_track(50-j)
                    gd_feature_set = gd_gene_track.new_set()
                    empty_feature = SeqFeature(FeatureLocation(0,1), 
                            strand=None)
                    gd_feature_set.add_feature(empty_feature, color="black", label=True)
                    print 'blank track added'
                else:
                    gene = genes[i*50 + j][0]
                    make_gene_track(gd_diagram, pham, genes[i*50 + j], j, 50)
            print seq_length, i

            gd_diagram.draw(format="linear", orientation="portrait", pagesize=reportlab.lib.pagesizes.letter,
                            fragments=1, start=0, end=seq_length)
            gd_diagram.write(graph_path, "PDF")
            gd_diagram.write(graph_path_svg, "SVG")

            add_pham_no_title(args, args.pham_no, graph_path, str(i))

        combine_graphs(args, args.phage, pham.pham_no, i)
    else:  
        final_graph_path = os.path.join(file_path, "Pham%sGraph_.pdf" % (pham.pham_no)) if not args.phage else os.path.join(file_path, "%sPham%sGraph_.pdf" % (args.phage+args.one_or_all, pham.pham_no))        
        # graph_path_svg = "%sPham%sGraph.svg" % (file_path+args.phage+ args.one_or_all, pham.pham_no)
        graph_path = os.path.join(file_path, "Pham%sGraph_.pdf" % (pham.pham_no)) if not args.phage else  os.path.join(file_path,"%sPham%sGraph_.pdf" % (args.phage, pham.pham_no))
        print "making_files.graph_start_sites: path to graph is " + str(graph_path)
        if check_file(final_graph_path):
            pass
        else:
            gd_diagram = GenomeDiagram.Diagram(pham.pham_no)
            i = 0
            for gene_group in genes:
                print 'making_files.graph_start_sites: adding group ' + str(i)
                make_gene_track(gd_diagram, pham, gene_group, i, len(genes))
                i += 1
            gd_diagram.draw(format="linear", orientation="portrait", pagesize=reportlab.lib.pagesizes.letter,
                            fragments=1, start=0, end=len(gene_group[0].alignment))
            gd_diagram.write(graph_path, "PDF")
        # gd_diagram.write(graph_path_svg, "SVG")
            add_pham_no_title(args, pham.pham_no, graph_path)

        # gd_diagram.write("%s.svg" % (file_path+pham.pham_no), "SVG")
        # gd_diagram.write("%s.eps" % (file_path+pham.pham_no), "EPS")
        # gd_diagram.write("%s.png" % (file_path+pham.pham_no), "PNG")

def make_pham_text(args, pham, pham_no, output_dir, only_pham=False):
    """
        Creates a PDF Report for the specific pham.
        From Start sites statisitics
        phage specific
    """
    if not args.phage:
        name = os.path.join(output_dir, "Pham%sText.pdf" % (pham_no))
    else:
        name = os.path.join(output_dir,"%sPham%sText.pdf" % (args.phage + args.one_or_all, pham_no))
    if check_file(name):
        return
    doc = SimpleDocTemplate(name, pagesize=reportlab.lib.pagesizes.letter)
    story = []
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="paragraph"))
    styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
    text = '<font size=14> Pham %s Report </font>' % pham_no  #item A
    story.append(Paragraph(text, styles['Center']))
    story.append(Spacer(1, 12))
    currentDate = time.strftime("%x")
    rundate = '<font size=12>This analysis was run %s. </font>' % currentDate  #item B
    story.append(Paragraph(rundate, styles["Normal"]))
    story.append(Spacer(1,12))

    phamCount = len(pham.genes.values())
    draftCount = sum(1 for g in pham.genes.values() if g.draftStatus)
    annotCount = phamCount - draftCount
    summaryText = "<font size=12>Pham number %s has %s members, %s are drafts.</font>" % (
    pham_no, len(pham.genes), draftCount)
    story.append(Paragraph(summaryText, styles["Normal"]))  # item E
    story.append(Spacer(1,12))


    story.append(Paragraph('<font size=12>Phages represented in each track:</font>', styles["Normal"])) #item D start

    groups = pham.group_similar_genes()
    tracks_info = []
    for index in range(len(groups)):
        text = "Track %s : " % (index+1)
        text += ", ".join(gene.gene_id for gene in groups[index])
        tracks_info.append("<font size=12> "+ u'\u2022'+" %s</font>" % text)
    for line in tracks_info:
        story.append(Paragraph(line, styles["Normal"]))
    story.append(Spacer(1, 12))  #item D end
    if only_pham:

        start_stats = pham.stats["most_common"]
        start_stats["phamCount"] = phamCount
        start_stats["annotCount"] = annotCount
        start_stats["draftCount"] = draftCount
        output = output_start_sites(start_stats) #this does items F through ??
        for line in output:
            if line == '':
                story.append(Spacer(1, 12))
            text = '<font size=12> %s </font>' % line
            # if 'Genes' not in line or '':
            story.append(Paragraph(text, styles['Normal']))
            story.append(Spacer(1, 12))
            # else:
            #     story.append(Paragraph(text, styles['Normal']))
        # story.append()
        # story.append(Paragraph("<font size=14>Suggested Starts:</font>", styles["Normal"]))
        # suggested_start = output_suggested_starts(pham, all_genes=True)
        # story.append(Spacer(1, 12))

        # for line in suggested_start:
        #     text = '<font size=12>%s</font>' % line
        #     story.append(Paragraph(text, styles["Normal"]))
    else:
        # story.append(Paragraph("<font size=14>Suggested Starts: </font>", styles["Normal"]))
        # starts = pham.stats["most_common"]
        # suggested_start = output_suggested_starts(pham, args.phage)
        # story.append(Spacer(1, 12))
        #
        # for line in suggested_start:
        #     text = '<font size=12>%s</font>' % line
        #     story.append(Paragraph(text, styles["Normal"]))
        story.append(Paragraph("",styles["Normal"]))
        story.append(Paragraph("<font size=12>Gene Information:</font>", styles["Normal"]))
        pham_possible_starts = pham.total_possible_starts
        for gene in pham.genes.values():
            if args.phage in gene.gene_id:
                candidate_starts = []
                for start in gene.alignment_candidate_starts:
                    candidate_starts.append((pham_possible_starts.index(start)+1, gene.alignment_index_to_coord(start)+1))
                if gene.orientation == "R":
                    story.append(Paragraph("<font size = 12> Gene: %s \n Start: %s, Stop: %s </font>" % (gene.gene_id, gene.start, gene.stop+1), styles["Normal"]))
                else:
                    story.append(Paragraph("<font size = 12> Gene: %s \n Start: %s, Stop: %s </font>" % (gene.gene_id, gene.start+1, gene.stop), styles["Normal"]))
                story.append(Paragraph("<font size = 12> Candidate Starts for %s: </font>" % (gene.gene_id), styles["Normal"]))
                story.append(Paragraph("<font size = 12>"+ str(candidate_starts) + "</font>" , styles["Normal"]))


    #note is item C:
    note = '<font size=12>Note: In the above figure, yellow indicates the location of called starts comprised solely of computational predictions, '
    note += 'green indicates the location of called starts with at least 1 manual gene annotation. In the text below, numbers found inside '
    note += 'square brackets (i.e. []) are derived from biopython and are zero-based, add 1 to the coordinate to find the corresponding location in a 1-based coordinate system. </font>'

    story.append(Paragraph(note, styles["Normal"]))
    doc.build(story)


def make_pham_genome(phage_genes, phage_name, length, file_path):
    file_name = os.path.join(file_path,'%sPhamsGraph.pdf' % (phage_name))
    if check_file(file_name):
        return
    pham_colors = phams.get_pham_colors()
    gd_diagram = GenomeDiagram.Diagram(phage_name)
    gd_track = gd_diagram.new_track(1, name=phage_name, greytrack=1)
    gd_pham_set = gd_track.new_set()
    print "making genome page"
    for gene_id in sorted(phage_genes.iterkeys()):
        phage_gene = phage_genes[gene_id]
        pham_no = phage_gene["pham_no"]
        gene = phage_gene["gene"]
        print pham_no, gene.gene_id
        if pham_no == None:
            pham_no = "None"
            pham_color = 'Black'
        else:
            pham_color = colors.HexColor(pham_colors[str(pham_no)])
        if gene.orientation == 'F':
            gene_location = FeatureLocation(gene.start, gene.stop)
            strand = 1
        else:
            strand = -1
            gene_location = FeatureLocation(gene.stop, gene.start)
        gene_feature = SeqFeature(gene_location, strand=strand)
        gene_number = phamgene.get_gene_number(gene.gene_id)
        # label the gene with the gene number
        gd_pham_set.add_feature(gene_feature, name=str(gene_number), label=True, 
            label_size=6, label_angle=75)
        # label gene with pham color and name
        gd_pham_set.add_feature(gene_feature, color=pham_color, name=str(pham_no), label=True, label_position='middle')
    
    print type(length), length
    gd_diagram.draw(format='linear', orientation='portrait', pagesize=reportlab.lib.pagesizes.letter, fragments=8, start=0, end=length)
    gd_diagram.write(file_name, "PDF")

def make_suggested_starts(phage_genes, phage_name, file_path):
    """
        Creates a PDF page of the suggested starts of a phage
        Genes are list in order
        {Gene Name} is a member of Pham {Number}: {Suggested Start Coordinates}
    """
    file_name = os.path.join(file_path, "%sSuggestedStarts.pdf" % (phage_name))
    text_file_name = os.path.join(file_path, "%sSuggestedStarts.txt" % (phage_name))
    if check_file(file_name):
        return
    doc = SimpleDocTemplate(file_name, pagesize=reportlab.lib.pagesizes.letter)
    story = []
    just_text = []
    print "making suggested starts page"
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="paragraph"))
    styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
    text = '<font size=14> Suggested Start Coordinates</font>'
    just_text.append(' Suggested Start Coordinates')
    story.append(Paragraph(text, styles['Center']))
    story.append(Spacer(1, 12))
    for gene_id in sorted(phage_genes.iterkeys()):
        phage_gene = phage_genes[gene_id]
        pham = phage_gene["pham_no"]
        gene = phage_gene["gene"]
        suggested_start = phage_gene["suggested_start"]
        if pham == None:
            text = '<font size=12> %s is not a member of an existing Pham </font>' % (gene.gene_id)
            simple_text = '%s is not a member of an existing Pham ' % (gene.gene_id)
        else:
            text = '<font size=12> %s is a member of Pham %s:  %s </font>' % (gene.gene_id, pham, suggested_start)
            simple_text = '%s is a member of Pham %s:  %s ' % (gene.gene_id, pham, suggested_start)
        story.append(Paragraph(text, styles['Normal']))
        just_text.append(simple_text)
    doc.build(story)
    print "writing text file"
    with open(text_file_name, 'w' ) as outfile:
        outfile.write("\n".join(just_text))


def make_fasta_file(genes, fasta_file):
    count = SeqIO.write(genes, fasta_file, 'fasta')


def main():
    args = parse_arguments()
    print "making_files:main(); args.make is ", args.make
    if 'graph' in args.make:
        print "making_files.main() make 'graph': args.pickle_file " + args.pickle_file
        pham = cPickle.load(open(args.pickle_file.strip('"'), 'rb'))
        graph_start_sites(args, pham, args.dir)

    if 'starts' in args.make:
        phage_genes = cPickle.load(open(args.pickle_file.strip('"'), 'rb'))

        make_suggested_starts(phage_genes, args.phage, args.dir)

    if 'genome' in args.make:
        phage = cPickle.load(open(args.pickle_file.strip('"'), 'rb'))
        make_pham_genome(phage, args.phage, args.phage_length, args.dir)
        make_suggested_starts(phage, args.phage, args.dir)

    if 'text' in args.make:
        print "making_files.main(): Loading pickle file " + str(args.pickle_file)
        pham = cPickle.load(open(args.pickle_file.strip('"'), 'rb'))
        graph_start_sites(args, pham, args.dir)
        print "making_files.main(): 'text' phage is ", args.phage
        if not args.phage:
            print "making_files.main() 'text': no phage"
            make_pham_text(args, pham, args.pham_no, args.dir, only_pham=True)
        else:
            make_pham_text(args, pham, args.pham_no, args.dir)

    if 'fasta' in args.make:
        pass
        # pickle_file = args.file
        # genes = cPickle.load(open(args.pickle_file.strip('"'), 'rb'))
        # make_fasta_file(genes, (args.dir + '.fasta'))

if __name__ == "__main__":
    main()
