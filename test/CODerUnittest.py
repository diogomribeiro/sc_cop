
import unittest
import os
import glob
import pandas
import numpy

from cod.cod_identification.CODer import CODer

pandas.options.display.max_rows = 5
pandas.options.display.max_columns = 3
# pandas.options.display.max_rows = 15
# pandas.options.display.max_columns = 30

class CODerUnittest(unittest.TestCase):
        
    # #
    def setUp(self):
        """Runs before each test. Name of this function needs forcedly to be 'setUp'."""

        #===============================================================================
        # Initialise main script
        #===============================================================================

        # Set the (default) options

        self.matrixFile = "test_input/test_expression_matrix_chr21_chr18.bed"
        self.outputFolder = "test_output/"
        self.permutations = CODer.DEFAULT_PARAM_PERMUTATIONS
        self.rankApproach = CODer.DEFAULT_PARAM_RANK_APPROACH
        self.windowSize = CODer.DEFAULT_PARAM_WINDOW_SIZE
        self.fdrCutoff = CODer.DEFAULT_PARAM_FDR_CUTOFF
        self.minimumCorrelation = CODer.DEFAULT_PARAM_MINIMUM_CORRELATION
        self.negativeCorrelation = CODer.DEFAULT_PARAM_NEGATIVE_CORRELATION
        self.consecutiveRank = CODer.DEFAULT_PARAM_CONSECUTIVE_RANK
        self.correctionProcedure = CODer.DEFAULT_PARAM_CORRECTION_PROCEDURE
        self.nullMaxDistance = CODer.DEFAULT_PARAM_NULL_MAX_DISTANCE
        self.determineTSS = CODer.DEFAULT_PARAM_DETERMINE_TSS
        self.corrMethod = CODer.DEFAULT_PARAM_CORR_METHOD
        self.wantedSeed = CODer.DEFAULT_PARAM_WANTED_SEED
        self.verbosityLevel = "debug"
        self.writeLog = CODer.DEFAULT_PARAM_WRITE_LOG
        self.lowMem = CODer.DEFAULT_PARAM_LOW_MEM
        self.rerun = CODer.DEFAULT_PARAM_RERUN
        self.variableMatch = CODer.DEFAULT_PARAM_VARIBLE_MATCH
        self.variableMatchProp = CODer.DEFAULT_PARAM_VARIBLE_MATCH_PROP
        
        # initialise class
        self.run = CODer( self.matrixFile, self.outputFolder, self.permutations, self.rankApproach,
                           self.windowSize, self.fdrCutoff, self.minimumCorrelation, self.negativeCorrelation, 
                           self.consecutiveRank, self.correctionProcedure, self.nullMaxDistance, self.determineTSS, self.corrMethod,
                           self.wantedSeed, self.verbosityLevel, self.writeLog, self.lowMem, self.rerun,
                           self.variableMatch, self.variableMatchProp)

        # unittest specific arguments
        self.expectedFolder = "test_expected/"


    def test_read_matrix_file(self):
        print ("| test_read_matrix_file | ")
  
        self.run.permutations = 0
  
        self.run.read_matrix_file()
  
        matrix = self.run.matrixData
  
        self.assertTrue( len( matrix[21]) == 358, "verifying number of rows")
        self.assertTrue(len(matrix[21].dtypes) == 390, "verifying number of columns")          
        self.assertTrue( len( matrix[18]) == 358, "verifying number of rows")
        self.assertTrue(len(matrix[18].dtypes) == 494, "verifying number of columns")          
 
 
    def test_determine_tss(self):
        print ("| test_read_matrix_file | ")
  
        self.run.permutations = 0
        self.run.determineTSS = 1
  
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
  
        self.assertTrue(self.run.phenoStart["ENSG00000160307.9_2"] == 48025121)
 
 
    def test_read_simple_matrix_file_permute(self):
        print ("| test_read_simple_matrix_file_permute | ")
 
        self.run.matrixFile = "test_input/simple_matrix.bed"
 
        self.run.permutations = 0
        self.run.read_matrix_file()
         
        originalMatrix = self.run.matrixData[21]
         
        # test permuting individuals
        self.run.permutations = 1000
        self.run.read_matrix_file()
                 
        matrix = self.run.matrixData[21]
        rand1 = self.run.randMatrixData[21][0]
        rand2 = self.run.randMatrixData[21][1]
 
        # Confirm that original matrix hasn't changed
        self.assertTrue( originalMatrix.equals(matrix), "confirm that the original dataframe was not modified with the permutations")
        self.assertTrue(list(matrix.loc["ind1"]) == list(originalMatrix.loc["ind1"]) )
        self.assertTrue(list(matrix["gene1"]) == list(originalMatrix["gene1"]) )
 
        # for each phenotype the total expression values should be the same, just shuffled
        self.assertTrue(sorted(list(matrix["gene1"])) == sorted(list(rand1["gene1"]) ))
        # for each individual the total expression values may be different
        #self.assertFalse(sorted(list(matrix.loc["ind1"])) == sorted(list(rand1.loc["ind1"]) ), "this test may fail often (1 every 3 times), rerun test.")
 
 
    def test_read_matrix_file_permutations_permute_individuals(self):
        print ("| test_read_matrix_file_permutations_permute_individuals | ")
 
        self.run.permutations = 2
 
        self.run.read_matrix_file()
 
        matrix = self.run.matrixData[21]
         
        rand1 = self.run.randMatrixData[21][0]
        rand2 = self.run.randMatrixData[21][1]
                 
        self.assertTrue(matrix.shape == rand1.shape)
        self.assertTrue(len(matrix["ENSG00000279579.3_4"] ) == len(rand2["ENSG00000279579.3_4"] ))
 
        # for each phenotype the total expression values should be the same, just shuffled
        self.assertTrue(sorted(list(matrix["ENSG00000279579.3_4"])) == sorted(list(rand1["ENSG00000279579.3_4"]) ))
 
        # check that they are indeed different
        self.assertFalse(list(rand2["ENSG00000279998.1_5"]) == list(rand1["ENSG00000279998.1_5"]) )
        self.assertFalse(list(matrix["ENSG00000279998.1_5"]) == list(rand1["ENSG00000279998.1_5"]) )
        self.assertFalse(list(matrix.loc["HG00096"]) == list(rand1.loc["HG00096"]) )
 
        # confirm that values per individual are different
        self.assertFalse(sorted(list(matrix.loc["HG00096"])) == sorted(list(rand1.loc["HG00096"])))
        self.assertFalse(sorted(list(rand2.loc["HG00096"])) == sorted(list(rand1.loc["HG00096"])))
         
 
    def test_read_real_coordinates(self):
        print ("| test_read_real_coordinates | ")
         
        self.run.permutations = 0
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
          
        self.assertTrue(self.run.realCoordinates["ENSG00000173213.9_5"] == (47390,49557))
  
  
    def test_perform_correlation(self):
        print ("| perform_correlation | ")
 
        self.run.negativeCorrelation = 0
 
        x = [0,0,2]
        y = [2,1,0]
        self.assertTrue(abs( self.run.perform_correlation(x, y)[0] + 0.866) < 0.0001) 
          
        x = [0,1,2]
        y = [0,1,2]
        self.assertTrue(self.run.perform_correlation(x, y) == (1.0, "+") )
  
        x = [0,1,2]
        y = [2,1,0]
        self.assertTrue(self.run.perform_correlation(x, y) == (-1.0, "-" ) )
  
        # test squaring
        self.run.negativeCorrelation = 1
  
        self.assertTrue(self.run.perform_correlation(x, y) == (1.0, "-" ) )
 
        # test absolute
        self.run.negativeCorrelation = 2
  
        x = [0,0,2]
        y = [2,1,0]
        self.assertTrue( abs(self.run.perform_correlation(x, y)[0] - 0.866) < 0.0001)
        self.assertTrue(self.run.perform_correlation(x, y)[1] == "-" )
 
  
    def test_get_cis_phenotypes(self):
        print ("| test_get_cis_phenotypes | ")
 
        self.run.permutations = 0
        self.run.read_matrix_file()
        self.run.windowSize = 1
         
        self.assertTrue(self.run.get_cis_phenotypes( 21, 1) == [], "nothing should be found for this coordinates" ) 
        self.assertTrue(self.run.get_cis_phenotypes( 21, 9647908) == [], "nothing should be found for this coordinates with window size 1"  )
        self.assertTrue(self.run.get_cis_phenotypes( 21, 9647909) == ["ENSG00000279579.3_4"], "gene in this exact position"  )
        self.assertTrue(self.run.get_cis_phenotypes( 21, 9647910) == [], "nothing should be found for this coordinates with window size 1"  )
 
        self.run.windowSize = 2
        self.assertTrue(self.run.get_cis_phenotypes( 21, 9647910) == ["ENSG00000279579.3_4"], "gene one position away"  )
        self.assertTrue(self.run.get_cis_phenotypes( 21, 9647908) == ["ENSG00000279579.3_4"], "gene one position away"  )
 
        self.assertTrue(self.run.get_cis_phenotypes( 21, 78005428) == [], "no gene in this position in this chromosome")
        self.assertTrue(self.run.get_cis_phenotypes( 18, 78005428) == ["ENSG00000178184.15_2"], "no gene in this position in this chromosome")
 
        self.run.windowSize = 20000
         
        self.assertTrue(sorted(self.run.get_cis_phenotypes( 18, 77806899)) == sorted(["ENSG00000101546.12_3","ENSG00000267127.7_6","ENSG00000141759.14_2","ENSG00000261126.7_5"]))
         
  
    def test_cis_window_rank_approach(self):
        print ("| test_cis_window_rank_approach | ")
  
        self.run.windowSize = 10000
        self.run.matrixFile = "test_input/test_expression_matrix_chr21_50genes.bed"
  
        self.run.permutations = 10
        self.run.rankApproach = 1
  
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
  
        self.run.cis_window_rank_approach()
          
        files = glob.glob("test_output/chr21/*")
                 
        self.assertTrue(len(files) == 50)
                  
        self.assertTrue(len(open(files[12],"r").readlines() ) == 13, "should be one line per permutation plus real cases and central pheno")
 
        # ENSG00000188992.11_3 only has 1 phenotype in cis window
        matrix = pandas.read_csv(self.outputFolder + "/chr21/" + "ENSG00000188992.11_3", sep="\t")
        self.assertTrue(abs( matrix["correlation"][10] - 0.2448) < 0.0001, "correlation value validated straight from input file" )
        self.assertTrue(matrix["correlation"][10] != matrix["correlation"][9], "confirm randomisation is different than real value")
        self.assertTrue(matrix["correlation"][8] != matrix["correlation"][9], "confirm two randomisations are different")
  
   
    def test_cop_identification_rank_approach(self):
        print ("| test_cop_identification_rank_approach| ")
   
        self.run.windowSize = 100000
        self.run.matrixFile = "test_input/test_expression_matrix_chr21_50genes.bed"
        self.run.permutations = 100
        self.run.rankApproach = 1
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
        self.run.cis_window_rank_approach()
  
        self.assertTrue( self.run.cop_identification() == 24)
        data = pandas.read_csv(self.outputFolder + "/" + CODer.OUTPUT_COP_RAW_RESULTS, sep="\t", index_col=False)
        data = data[data["rank"] != "central"]
        self.assertTrue(data.shape[0] == 38)
        data = pandas.read_csv(self.outputFolder + "/" + CODer.OUTPUT_COD_IDENTIFICATION_COPS, sep="\t", index_col=False)  
        self.assertTrue(data.shape[0] == 38)
        self.assertTrue(len(set(data["centralPhenotype"])) == 24)
          
  
    def test_cop_identification_first_rank_approach(self):
        print ("| test_cop_identification_first_rank_approach| ")
   
        self.run.windowSize = 100000
        self.run.matrixFile = "test_input/test_expression_matrix_chr21_50genes.bed"
        self.run.permutations = 100
        self.run.rankApproach = 0
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
        self.run.cis_window_rank_approach()
  
        self.assertTrue( self.run.cop_identification() == 24)
          
  
    def test_cop_identification_first_rank_approach_low_mem(self):
        print ("| test_cop_identification_first_rank_approach_low_mem| ")
   
        self.run.windowSize = 100000
        self.run.matrixFile = "test_input/test_expression_matrix_chr21_50genes.bed"
        self.run.permutations = 100
        self.run.rankApproach = 0
          
        self.run.lowMem = 1
          
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
        self.run.cis_window_rank_approach()
  
        self.assertTrue( self.run.cop_identification() == 24)
  
          
    def test_cop_identification_minimum_correlation(self):
        print ("| test_cop_identification_minimum_correlation | ")
   
        self.run.matrixFile = "test_input/test_expression_matrix_chr21_50genes.bed" 
        self.run.windowSize = 100000
        self.run.permutations = 50
        self.run.rankApproach = 1
          
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
  
        self.run.minimumCorrelation = 0
        self.run.cis_window_rank_approach()
        zeroRes1 = self.run.cop_identification()
        self.assertTrue(zeroRes1 > 10 )
  
        self.run.minimumCorrelation = 1
        self.assertTrue(self.run.cop_identification() == 0 )
  
        self.run.minimumCorrelation = 0.5
        self.assertTrue(self.run.cop_identification() < zeroRes1)    
                   
        self.run.minimumCorrelation = -999
        self.run.cis_window_rank_approach()
        self.assertTrue(self.run.cop_identification() == 24 )
  
  
    def test_cop_identification_rank_approach_consecutive(self):
        print ("| test_cop_identification_rank_approach_consecutive| ")
  
        self.run = CODer( "test_input/test_expression_matrix_chr21_50genes.bed", self.outputFolder, 100, 1,
                           1000000, 0.1, self.minimumCorrelation, self.negativeCorrelation, 
                           1, self.correctionProcedure, self.nullMaxDistance, self.determineTSS, self.corrMethod, 
                           54233901, self.verbosityLevel, 
                           self.writeLog, self.lowMem, self.rerun, self.variableMatch, self.variableMatchProp )
          
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
        self.run.cis_window_rank_approach()
        result = self.run.cop_identification()
  
        # some results should be filtered
        self.assertTrue(result < 40)
        #  ENSG00000229306.1_5 should not be filtered 
        self.assertTrue(self.run.outputData[self.run.outputData["centralPhenotype"] == "ENSG00000229306.1_5"].shape[0] == 1)
        # ENSG00000224141.5_4 should be filtered
        self.assertTrue(self.run.outputData[self.run.outputData["centralPhenotype"] == "ENSG00000224141.5_4"].shape[0] < 10)
        # NOTE: this is not using random seed so it may fail from time to time (but not often!)
                  
           
    def test_fdr_filter(self):
        print ("| test_fdr_filter | ") 
          
        x = numpy.array([0.2,0.02,0.3,0.5,0.01,0.9,0.001,1.0])
          
        self.assertTrue(self.run.q_value(x)[0] == 0.4, "confirm the result is the same as in R")
        self.assertTrue(self.run.q_value(x)[6] == 0.008, "confirm the result is the same as in R")
  
        self.run.matrixFile = "test_input/test_expression_matrix_chr21_50genes.bed"
        self.run.permutations = 10
        self.run.rankApproach = 1
  
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
        self.run.cis_window_rank_approach()
          
        CODer.GLOBAL_PVALUE_CORRECTION = 1
          
        self.run.minimumCorrelation = 0.0
           
        noFilterCods = self.run.cop_identification()
       
        self.run.fdrCutoff = 0.5     
        mildFilterCods = self.run.cop_identification()
        self.assertTrue( mildFilterCods < noFilterCods)
           
        self.run.fdrCutoff = 0.001
                           
        self.assertTrue( self.run.cop_identification() < mildFilterCods)
  
          
    def test_fdr_filter_rank0(self):
        print ("| test_fdr_filter_rank0 | ") 
          
        self.run.matrixFile = "test_input/test_expression_matrix_chr21_50genes.bed"
        self.run.permutations = 10
        self.run.rankApproach = 0
  
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
        self.run.cis_window_rank_approach()
          
        noFilterCods = self.run.cop_identification()
          
        self.assertTrue( noFilterCods == 50)
       
        self.run.fdrCutoff = 0.5     
        self.run.cis_window_rank_approach()
        mildFilterCods = self.run.cop_identification()        
        self.assertTrue( mildFilterCods < noFilterCods)
                   
          
    def test_estimate_adjusted_pval_cutoff(self):
        print ("| test_estimate_adjusted_pval_cutoff | ") 
          
        self.run.fdrCutoff = 0.05
          
        x = [0.3,0.2,0.009,0.1,0.09,0.015,0.00099,0.001, 0.025,0.1,0.02,0.5,0.003,1]
        pvals = pandas.DataFrame( {"adjustedPval" : x, "fdr" : self.run.bh_correct(x) } )
  
        self.assertTrue(self.run.estimate_adjusted_pval_cutoff(pvals) == (0.02, 0.04666666666666667) )
  
        self.run.fdrCutoff = 0.01
        self.assertTrue(self.run.estimate_adjusted_pval_cutoff(pvals) == (0.00099, 0.007000) )
          
  
    def test_cod_identification_network(self):
        print ("| test_cod_identification_network | ") 
          
        self.run.matrixFile = "test_input/test_expression_matrix_chr21_50genes.bed"
        self.run.permutations = 1
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
          
        self.run.outputData = pandas.read_csv("test_input/CODer_cod_identification.bed", sep="\t", index_col=False)        
  
        totalEdges,totalNodes,totalComponents = self.run.cod_identification_network()
  
        # values validated using cytoscape
        self.assertTrue( totalEdges == 17866)
        self.assertTrue( totalNodes == 8913)
        self.assertTrue( totalComponents == 2575)
          
          
    def test_cod_identification_network_two(self):
        print ("| test_cod_identification_network_two | ") 
          
        self.run.matrixFile = "test_input/test_expression_matrix_chr21_50genes.bed"
        self.run.permutations = 1
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
        self.run.cis_window_rank_approach()
        self.run.cop_identification()
  
        totalEdges,totalNodes,totalComponents = self.run.cod_identification_network()
          
        self.assertTrue( totalNodes == 50)
        self.assertTrue( totalEdges == 448)
          
        result = open(self.outputFolder + CODer.OUTPUT_COD_IDENTIFICATION_PER_COD, "r").readlines()
          
        self.assertTrue(len( result ) == 3)
        cod1 = result[1].split("\t")
  
        self.assertTrue( cod1[1] == "9647910")        
        self.assertTrue( cod1[2] == "11184046")
        self.assertTrue( cod1[5] == "2.5")
        self.assertTrue( cod1[6] == "7")
        self.assertTrue( cod1[7].strip() == "1536136")
          
  
    def test_build_distance_matched_null(self):
        print ("| test_build_distance_matched_null | ")
  
        self.run.matrixFile = "test_input/test_expression_matrix_chr21_50genes.bed"
        self.run.permutations = 3
        self.run.fdrCutoff = 0.05
        self.run.nullMaxDistance = 1000
           
        self.run.read_matrix_file()
        self.run.read_real_coordinates()
        self.run.cis_window_rank_approach()
        self.run.cop_identification()
        self.run.cod_identification_network()
  
        self.run.build_distance_matched_null()
          
        controlData = pandas.read_csv(self.outputFolder + "/" + CODer.OUTPUT_COP_DISTANCE_CONTROL_NULL, sep="\t", index_col=False)                
          
        positive = controlData[controlData["significant"] == 1]
        negative = controlData[controlData["significant"] == 0]
          
        self.assertTrue(positive.shape[0] == negative.shape[0])
        self.assertTrue(abs(sum(positive["distance"]) - sum(negative["distance"])) <= ((self.run.nullMaxDistance * positive.shape[0] ) / 1.9) )
          
        # confirm that when positive appears once or twice (gene1-gene2, gene2-gene1), the matching null also appears once or twice
        from collections import Counter
        self.assertTrue(sorted(list(Counter(positive["pairID"].values).values())) == sorted(list(Counter(negative["pairID"].values).values())) ) 
  
  
#     def test_build_distance_matched_null_rerun(self):
#         print ("| test_build_distance_matched_null_rerun | ")
#    
#    
#         self.run.outputData = pandas.read_csv("/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/geuvadis/all_chr/final_fdr0.01/test/CODer_cod_identification_cops.bed", sep="\t", index_col=False)
#           
#         self.run.outputFolder = "/scratch/axiom/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/geuvadis/all_chr/final_fdr0.01/test"
#    
#         self.run.build_distance_matched_null()


    def test_build_distance_and_variable_matched_null(self):
        print ("| test_build_distance_and_variable_matched_null | ")
 
        self.rerun = 1
        self.nullMaxDistance = 0.1
        self.variableMatch = "diffExpr"
        self.variableMatchProp = 0.1
        self.run = CODer( self.matrixFile, self.expectedFolder, 1, self.rankApproach, \
                           self.windowSize, self.fdrCutoff, self.minimumCorrelation, self.negativeCorrelation, 
                           self.consecutiveRank, self.correctionProcedure, self.nullMaxDistance, self.determineTSS, self.corrMethod, 
                           self.wantedSeed, self.verbosityLevel, 
                           self.writeLog, self.lowMem, self.rerun, self.variableMatch, self.variableMatchProp )

#         self.run.outputData = pandas.read_csv(self.expectedFolder + "/" + CODer.OUTPUT_COD_IDENTIFICATION_COPS, sep="\t", index_col=False)
          
 
        #self.run.build_distance_and_variable_matched_null()
         
#         controlData = pandas.read_csv(self.outputFolder + "/" + CODer.OUTPUT_COP_DISTANCE_CONTROL_NULL, sep="\t", index_col=False)                
#          
#         positive = controlData[controlData["significant"] == 1]
#         negative = controlData[controlData["significant"] == 0]
#          
#         self.assertTrue(positive.shape[0] == negative.shape[0])
#         self.assertTrue(abs(sum(positive["distance"]) - sum(negative["distance"])) <= ((self.run.nullMaxDistance * positive.shape[0] ) / 1.9) )
#          
#         # confirm that when positive appears once or twice (gene1-gene2, gene2-gene1), the matching null also appears once or twice
#         from collections import Counter
#         self.assertTrue(sorted(list(Counter(positive["pairID"].values).values())) == sorted(list(Counter(negative["pairID"].values).values())) ) 
     
     
    def test_seed(self):
   
        print ("| test_run | ") 
          
        self.run = CODer( "test_input/test_expression_matrix_chr21_50genes.bed", self.outputFolder, 1, self.rankApproach, \
                           self.windowSize, self.fdrCutoff, self.minimumCorrelation, self.negativeCorrelation, 
                           self.consecutiveRank, self.correctionProcedure, self.nullMaxDistance, self.determineTSS,
                           self.corrMethod, 666, self.verbosityLevel, 
                           self.writeLog, self.lowMem, self.rerun, self.variableMatch, self.variableMatchProp  )
          
        self.run.read_matrix_file()
          
        run1 = self.run.randMatrixData
  
        self.run = CODer( "test_input/test_expression_matrix_chr21_50genes.bed", self.outputFolder, 1, self.rankApproach, \
                           self.windowSize, self.fdrCutoff, self.minimumCorrelation, self.negativeCorrelation, 
                           self.consecutiveRank, self.correctionProcedure, self.nullMaxDistance, self.determineTSS,
                           self.corrMethod, 666, self.verbosityLevel, 
                           self.writeLog, self.lowMem, self.rerun, self.variableMatch, self.variableMatchProp  )
          
        self.run.read_matrix_file()
          
        run2 = self.run.randMatrixData
  
        self.assertTrue(list(run1[21][0]["ENSG00000279579.3_4"]) == list(run2[21][0]["ENSG00000279579.3_4"]) )
          
        self.run = CODer( "test_input/test_expression_matrix_chr21_50genes.bed", self.outputFolder, 1, self.rankApproach, \
                           self.windowSize, self.fdrCutoff, self.minimumCorrelation, self.negativeCorrelation, 
                           self.consecutiveRank, self.correctionProcedure, self.nullMaxDistance, self.determineTSS, 
                           self.corrMethod, 667, self.verbosityLevel, 
                           self.writeLog, self.lowMem, self.rerun, self.variableMatch, self.variableMatchProp  )
          
        self.run.read_matrix_file()
        run3 = self.run.randMatrixData
  
        self.assertFalse(list(run1[21][0]["ENSG00000279579.3_4"]) == list(run3[21][0]["ENSG00000279579.3_4"]) )
          
          
    def test_run(self):
   
        print ("| test_run | ") 
        self.run = CODer( "test_input/test_expression_matrix_chr21_50genes.bed", self.outputFolder, 1, self.rankApproach, \
                           self.windowSize, self.fdrCutoff, self.minimumCorrelation, self.negativeCorrelation, 
                           self.consecutiveRank, self.correctionProcedure, self.nullMaxDistance, self.determineTSS,
                           self.corrMethod, self.wantedSeed, self.verbosityLevel, 
                           self.writeLog, self.lowMem, self.rerun, self.variableMatch, self.variableMatchProp )
          
        self.run.run()
  
  
    def test_rerun(self):
   
        print ("| test_run | ") 
        self.run = CODer( "test_input/test_expression_matrix_chr21_50genes.bed", self.outputFolder, 1, self.rankApproach, \
                           self.windowSize, self.fdrCutoff, self.minimumCorrelation, self.negativeCorrelation, 
                           self.consecutiveRank, self.correctionProcedure, self.nullMaxDistance, self.determineTSS, 
                           self.corrMethod, self.wantedSeed, self.verbosityLevel, 
                           self.writeLog, self.lowMem, self.rerun, self.variableMatch, self.variableMatchProp )
          
        self.run.run()
          
        data1 = pandas.read_csv(self.outputFolder + "/" + CODer.OUTPUT_COD_IDENTIFICATION_COPS, sep="\t", index_col=False)
          
        self.run = CODer( "test_input/test_expression_matrix_chr21_50genes.bed", self.outputFolder, 1, self.rankApproach, \
                           self.windowSize, self.fdrCutoff, self.minimumCorrelation, self.negativeCorrelation, 
                           self.consecutiveRank, self.correctionProcedure, self.nullMaxDistance, self.determineTSS, 
                           self.corrMethod, self.wantedSeed, self.verbosityLevel, 
                           self.writeLog, self.lowMem, 1, self.variableMatch, self.variableMatchProp )
        self.run.run()
  
        data2 = pandas.read_csv(self.outputFolder + "/" + CODer.OUTPUT_COD_IDENTIFICATION_COPS, sep="\t", index_col=False)
  
        self.assertTrue(data1.shape == data2.shape)
        
     
    def tearDown(self):
        """Runs after each test"""
                                                     
        # Wipe output folder
        cmd = "rm -r %s/*" % self.outputFolder
        os.system(cmd)


if __name__ == '__main__':
    unittest.main()

    