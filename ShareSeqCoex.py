#!/usr/bin/env python3.6

import argparse
import gzip
import pandas
from scipy.stats import pearsonr
from scipy.stats import spearmanr

from util.log.Logger import Logger
from util.time.Timer import Timer

#===============================================================================
DESC_COMMENT = "Script to identify co-expressed gene-peak pairs from single cell a atac-seq and gene expression in same cells."
SCRIPT_NAME = "ShareSeqCoex.py"
#===============================================================================

"""
#===============================================================================
@author: Diogo Ribeiro
@date: 15 December 2020 (last update 13 December 2021)
@copyright: Copyright 2020-2021, University of Lausanne
"Script to identify co-expressed gene pairs from single cell a gene expression sparse matrix"
#===============================================================================
"""

class ShareSeqCoex(object):
        
    WINDOW_SIZE = 1000000
    DEFAULT_PARAM_CORR_METHOD = 1
    
    def __init__(self, geneExpressionMatrix, peakExpressionMatrix, geneModelFile,
                 outputFile, corrMethod, verbosityLevel):

        Timer.get_instance().start_chrono()

        self.geneExpressionMatrix = geneExpressionMatrix
        self.peakExpressionMatrix = peakExpressionMatrix
        self.geneModelFile = geneModelFile
        self.outputFile = outputFile
        self.corrMethod = corrMethod

        Logger.get_instance().set_level( verbosityLevel)        


    def read_sparse_gene_expression_matrix(self):
        """
        Read gene expression matrix in sparse format.
        
        Example format (TSV, header):
            gene    cell    value
            ENSG00000238009 R1.02,R2.10,R3.20,P1.51 1
            ENSG00000238009 R1.17,R2.69,R3.50,P1.51 1
            ENSG00000238009 R1.20,R2.74,R3.17,P1.50 1
            
        Note: calculates how many cells in input (i.e. only genes and cells that have non-zero value may be counted)
        """

        Logger.get_instance().info( "read_sparse_expression_matrix: reading file %s" % ( self.geneExpressionMatrix) )

        if self.geneExpressionMatrix.endswith(".gz"):
            fileHandler = gzip.open(self.geneExpressionMatrix, "rt")
        else:
            fileHandler = open(self.geneExpressionMatrix, "r")

        geneExpression = {} # key -> geneID, value -> set of cells where expressed
        setGenes = set()
        setCells = set()
        _ = fileHandler.readline()
        count = 0
        for line in fileHandler:
            count+=1
            if count % 10000000 == 0: 
                Timer.get_instance().step( "read_sparse_expression_matrix: read %s lines.." % (count) )  

            line = line.strip()
            spl = line.split("\t")
            gene = spl[0]
            cell = spl[1]
            
            # Prune last cell barcode tag to map to peak cells
            # Note: this is specific to SHARE-seq cell barcodes
            cell = ",".join(cell.split(",")[0:3])
            
            if gene not in geneExpression:
                geneExpression[gene] = set()
            
            geneExpression[gene].add(cell)
            
            setGenes.add(gene)
            setCells.add(cell)
        
        Logger.get_instance().info( "read_sparse_expression_matrix: unique cells %s" % (len(setCells)) )
        Logger.get_instance().info( "read_sparse_expression_matrix: unique genes %s" % (len(setGenes)) )
                
        self.setCells = setCells
        self.geneExpression = geneExpression


    def read_sparse_peak_expression_matrix(self):
        """
        Read peak expression matrix in sparse format.
        
        Example format (TSV, no header):
            chr21   15352400        15352499        chr21   15352220        15352427        R1.51,R2.57,R3.95,P1.03
            chr21   15352400        15352499        chr21   15352129        15352459        R1.39,R2.22,R3.07,P1.02
            
        Last column must be cell ID, first 3 columns must be BED (coordinates of peak)
        """

        Logger.get_instance().info( "read_sparse_peak_expression_matrix: reading file %s" % ( self.peakExpressionMatrix) )

        if self.peakExpressionMatrix.endswith(".gz"):
            fileHandler = gzip.open(self.peakExpressionMatrix, "rt")
        else:
            fileHandler = open(self.peakExpressionMatrix, "r")

        peakExpression = {} # key -> peak, val -> cell
        chromosomePeaks = {} # key -> chromosome, val -> peak
        setPeaks = set()
        setCells = set()
        count = 0
        for line in fileHandler:
            count+=1
            if count % 10000000 == 0: 
                Timer.get_instance().step( "read_sparse_peak_expression_matrix: read %s lines.." % (count) )  

            line = line.strip()
            spl = line.split("\t")
            chro = spl[0]
            start = spl[1]
            end = spl[2]
            cell = spl[-1]
            
            if chro.startswith("chr"):
                chro = chro.replace("chr","")
            
            # Note: this is specific to SHARE-seq cell barcodes, to match RNA-seq dataset
            # Going from R1.39,R2.22,R3.07,P1.02 to R1.39,R2.22,R3.07
            cell = ",".join(cell.split(",")[0:3])
            
            # Attributing ID to peak based on genomic location
            peak = "_".join([chro,start,end])
            if peak not in peakExpression:
                peakExpression[peak] = set()
            
            peakExpression[peak].add(cell)
            
            setPeaks.add(peak)
            setCells.add(cell)

            if chro not in chromosomePeaks:
                chromosomePeaks[chro] = set()
            chromosomePeaks[chro].add(peak)
        
        Logger.get_instance().info( "read_sparse_peak_expression_matrix: unique cells %s" % (len(setCells)) )
        Logger.get_instance().info( "read_sparse_peak_expression_matrix: unique peaks %s" % (len(setPeaks)) )
        
        self.setPeakCells = setCells
        self.peakExpression = peakExpression
        self.chromosomePeaks = chromosomePeaks
        

    def read_gene_models(self):
        """
        Read file with cells on rows and information about them on columns.

        Example format (TSV, no header):
            1       .  gene    11869   14362   .       +       .       gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; \
                gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; \
                transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";

        """

        Logger.get_instance().info("read_gene_models: reading %s " % self.geneModelFile )

        if ".gz" in self.geneModelFile:  inFile = gzip.open(self.geneModelFile, "rt" )
        else: inFile = open(self.geneModelFile)

        geneModelInfo = {} # key -> gene name, val -> (chr,start,end,strand,info)
        chromosomeGenes = {} # key -> chromosome, val -> list of genes in this chromosome
        chromosomeTSSs = {} # key -> chromosome, val -> list of tss in this chromosome
        for line in inFile:
            if line.startswith("#"):
                continue
            line = line.strip()
            chro,_,typ,start,end,_,strand,_,info = line.split("\t")
            
            if chro.startswith("chr"):
                chro = chro.replace("chr","")
            
            if typ != "gene":              
                continue   # skip entries that are not relative to a gene

            # info example: gene_id "ENSG00000204481.6"; transcript_id "ENSG00000204481.6"; gene_type "protein_coding";
            geneID = ""
            geneType = ""
            geneName = ""
            for inf in info.split(";"):
                if "gene_id" in inf:
                    geneID = inf.split('"')[1]
                if "gene_type" in inf:
                    geneType = inf.split('"')[1]
                elif "gene_biotype" in inf:
                    geneType = inf.split('"')[1]
                if "gene_name" in inf:
                    geneName = inf.split('"')[1]

            if geneID == "":
                raise Exception("read_gene_models: could not identify the gene ID in info: " % info)
            if geneType == "":
                Logger.get_instance().warning("read_gene_models: gene_type not present in GTF file info: %s" % info)
            if geneName == "":
                raise Exception("read_gene_models: could not identify the gene name in info: " % info)
                
            # Compose string like: L=3726;T=lincRNA;R=chr1:89295-133723;N=AL627309.1
            infoText = "L=%s;T=%s;R=chr%s:%s-%s;N=%s" %  ("NA", geneType, chro, start, end, geneName)

            tag = geneID.split(".")[0]

            geneModelInfo[tag] = (chro, start, end, strand, infoText)
            
            if strand == "+":
                tss = int(start)
            else:
                tss = int(end)
            
            if chro not in chromosomeGenes:
                chromosomeGenes[chro] = []
            chromosomeGenes[chro].append(tag)

            if chro not in chromosomeTSSs:
                chromosomeTSSs[chro] = []
            chromosomeTSSs[chro].append(tss)

        Logger.get_instance().info("read_gene_models: read %s gene entries" % len(geneModelInfo) )

        self.geneModelInfo = geneModelInfo
        self.chromosomeGenes = chromosomeGenes
        self.chromosomeTSSs = chromosomeTSSs


    def build_possible_gene_peak_pairs(self):
        """
        Based on gene TSS and a cis window around it, make a list of all gene-peak pairs to be assessed.
        """
        
        possibleGenePeakPairs = set() # set of all possible gene peak pair combinations within cis window
        
        for chro in sorted(self.chromosomeGenes):
            Logger.get_instance().debug("build_possible_gene_pairs: processing chromosome %s" % chro )

            if chro not in self.chromosomePeaks:
                Logger.get_instance().debug("build_possible_gene_pairs: chromosome not in peak file %s (ignored)" % chro )
                continue
                
            p = [] #peak ID
            s = [] # start coordinate
            e = [] # end coordinate
            for peak in self.chromosomePeaks[chro]:
                p.append(peak)
                spl = peak.split("_")
                s.append(int(spl[1]))
                e.append(int(spl[2]))
            
            peakDF = pandas.DataFrame({"peak": p, "start": s, "end": e })

            geneDF = pandas.DataFrame({"gene": self.chromosomeGenes[chro], "tss": self.chromosomeTSSs[chro] })

            # for each gene
            for _,row in geneDF.iterrows():
                tss = row["tss"]
                gene = row["gene"]
                leftBoundary = tss - ShareSeqCoex.WINDOW_SIZE
                rightBoundary = tss + ShareSeqCoex.WINDOW_SIZE
                
                # require that both peak start and end to be inside window
                l1 = peakDF[peakDF["start"] > leftBoundary][peakDF["end"] > leftBoundary]
                possible = l1[l1["start"] < rightBoundary][peakDF["end"] < rightBoundary]
                
                for peak in possible["peak"].values:
                    pair = "|".join([gene, peak])
                    possibleGenePeakPairs.add(pair)                
            
        Logger.get_instance().info("build_possible_gene_peak_pairs: possible gene-peak pairs = %s" % len(possibleGenePeakPairs) )
        
        self.possibleGenePeakPairs = possibleGenePeakPairs


    def loop_gene_peak_pairs(self):
        """
        Go over each gene-peak and calculate their correlation.
        
        The correlation method used depends on the --corrMethod flag.
        Report also number of cases where both gene and peak are 1, or both zero etc.
        
        """

        # TODO: process cell overlap before storing any cell data, to save memory
        # Overlap between cells of each dataset
        cellOverlap = self.setCells & self.setPeakCells        
        
        outFile = open(self.outputFile, "w")
        outFile.write("pairID\toneOne\toneZero\tzeroOne\tzeroZero\tcorr\tcorrSign\tcorrPval\n")
        
        processing = 0
        count = 0
        geneNotFound = set()
        for pair in sorted(self.possibleGenePeakPairs):
            gene,peak = pair.split("|")

            processing+=1
            if processing % 50000 == 0: 
                Timer.get_instance().step( "loop_gene_peak_pairs: processed %s gene pairs.." % (processing) )  
            
            if gene in self.geneExpression and peak in self.peakExpression:
                                                
                geneExpr = self.geneExpression[gene] & cellOverlap
                peakExpr = self.peakExpression[peak] & cellOverlap
                                
                intersect = len(geneExpr & peakExpr)
                g1MinusG2 = len(geneExpr - peakExpr)
                g2MinusG1 = len(peakExpr - geneExpr)
                empty = len(self.setCells) - (intersect + g1MinusG2 + g2MinusG1)
                
                # calculate normal correlation
                l1 = ([1]*intersect + [1]*g1MinusG2 + [0]*g2MinusG1 + [0]*empty)
                l2 = ([1]*intersect + [0]*g1MinusG2 + [1]*g2MinusG1 + [0]*empty)
                corr, corrSign, corrPval = self.perform_correlation(l1, l2)
                                
                outFile.write("%s\t%i\t%i\t%i\t%i\t%.3f\t%s\t%.3g\n" % (pair, intersect, g1MinusG2, g2MinusG1, empty, corr, corrSign, corrPval))
                
                count+=1
            else:
                if gene not in self.geneExpression:
                    geneNotFound.add(gene)
                    Logger.get_instance().debug("loop_gene_peak_pairs: no expression for gene %s" % gene )
                if peak not in self.peakExpression:
                    geneNotFound.add(peak)
                    Logger.get_instance().debug("loop_gene_peak_pairs: no expression for peak %s" % peak )
                    
                continue

        Logger.get_instance().info("loop_gene_peak_pairs: %s genes not found" % len(geneNotFound) )
        Logger.get_instance().info("loop_gene_peak_pairs: wrote %s entries" % count )

        outFile.close()
                

    def perform_correlation(self, x, y):
        """
        Performs correlation between two arrays
 
        @param x and y: arrays to be compared.
         
        @return correlation coefficient, the sign of the correlation (positive/negative) and p-value
        """
 
        assert len(x)==len(y)         
     
        if self.corrMethod:
            corr,pval = pearsonr(x, y)
        else:
            corr,pval = spearmanr(x,y)

        if corr < 0: sign = "-"
        else: sign = "+"
         
        return abs(corr), sign, pval
                

    def run(self):
        """
        Run functions in order
        """

        Timer.get_instance().step( "Read expression matrix.." )        
        self.read_sparse_gene_expression_matrix()

        Timer.get_instance().step( "Read peak matrix.." )        
        self.read_sparse_peak_expression_matrix()
        
        Timer.get_instance().step( "Read gene model info.." )        
        self.read_gene_models()

        Timer.get_instance().step( "Build gene pairs.." )        
        self.build_possible_gene_peak_pairs()
        
        Timer.get_instance().step( "Processing gene pairs.." )        
        self.loop_gene_peak_pairs()


if __name__ == "__main__":

    try:
    
        # Start chrono
        print ("STARTING " + SCRIPT_NAME)

        #===============================================================================
        # Get input arguments
        #===============================================================================
        parser = argparse.ArgumentParser(description= DESC_COMMENT) 
    
        # positional args
        parser.add_argument('geneExpressionMatrix', metavar='geneExpressionMatrix', type=str,
                             help='Sparse single cell gene expression matrix, each row has gene ID, cell ID and value (optional).')
        parser.add_argument('peakExpressionMatrix', metavar='peakExpressionMatrix', type=str,
                             help='Sparse single cell peak expression matrix, each row has gene ID, cell ID and value (optional).')
        parser.add_argument('geneModelFile', metavar='geneModelFile', type=str,
                             help='Information about each gene, GTF file. Only genes in this file will be processed')
        parser.add_argument('outputFile', metavar='outputFile', type=str,
                             help='File where output will be written.')
        parser.add_argument('--corrMethod', metavar='corrMethod', type=int, default = ShareSeqCoex.DEFAULT_PARAM_CORR_METHOD, 
                             help='If 1, use Pearson correlation, otherwise use Spearman correlation. (default = 1)')
        parser.add_argument('--verbosityLevel', metavar='verbosityLevel', type=str, default = "info", 
                             choices = ["debug", "info", "warning", "error", "critical", "fatal"],
                             help='Level of verbosity. Choices: "debug", "info", "warning", "error", "critical", "fatal"')
                
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        # Initialise class    
        run = ShareSeqCoex( args.geneExpressionMatrix, args.peakExpressionMatrix, args.geneModelFile, args.outputFile, args.corrMethod, args.verbosityLevel)

        run.run()

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except Exception as e:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + str(e))

