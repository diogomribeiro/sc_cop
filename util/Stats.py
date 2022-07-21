#!/usr/bin/env python3.6

"""
#===============================================================================
@author: Diogo Ribeiro
@date: 03 March 2019
@copyright: Copyright 2019, University of Lausanne

Collection of python implementations of statistical methods.

#===============================================================================
"""


class Stats( object):
    """
    Class containing utilities regarding statistic methods.
    """

    @staticmethod
    def fisher_exact_test( matrix, alt = "greater"):
        """
        Perform fishers exact test using scipy

        Equivalent of doing in R:
          fisher.test( matrix(c(TP,FN,FP,TN), nrow = 2), alternative = "greater")

          However: The calculated odds ratio is different from the one R uses.
          This scipy implementation returns the (more common) unconditional Maximum Likelihood Estimate
          while R uses the conditional Maximum Likelihood Estimate
        
        @param matrix : list - bi-dimensional list in the format [[TP,FP], [FN,TN]]

        @returns oddsRatio : Fisher's exact test odds ratio
        @returns pvalue : Fisher's exact test p-value

        """
        from scipy import stats
        
        oddsRatio, pvalue = stats.fisher_exact( matrix, alternative = alt)
            
        return oddsRatio, pvalue


    @staticmethod
    def empirical_pvalue( list_random, observed, aboveTail = 1):
        """
        Calculates empirical pvalue by comparing how many times the observed (real) value is above or below the list of values obtained from a random experiment.
   
        Extract from North et al. 2002 "A Note on the Calculation of Empirical P Values from Monte Carlo Procedures", doi:  10.1086/341527
        "Typically, the estimate of the P value is obtained as p = r/n, where n is the number of replicate samples that have been simulated and r is the 
        number of these replicates that produce a test statistic greater than or equal to that calculated for the actual data. 
        However, Davison and Hinkley (1997) give the correct formula for obtaining an empirical P value as (r+1)/(n+1)."
    
        @param list_random : list - list of numerical values obtained from randomisation
        @param observed : 'numerical' - numerical value (e.g. float, int) to be compared to list of random values
        @param aboveTail : boolean - which tail of pvalue to be calculated, if 1 counts times observed is above random, if 0, counts times below.   

        @return pval : float - empirical pvalue
        @return times : int - how many times observed value is above/below random values

        """
        
        # number of replicate samples
        n = len( list_random)

        # number of replicates that produce a test statistic greater than or equal to that calculated for the actual data.
        # can be lesser or equal if aboveTail is not True
        r = 0

        for val in sorted( list_random):
            if aboveTail:
                if val >= observed:
                    r +=1
            else:
                if val <= observed:
                    r +=1

        # empirical p-value based on North et al. 2002
        pval = ( r + 1) / float( n + 1)
        
        # number of times the observed value is "better" than random
        times = n - r
                
        return pval, times


    @staticmethod
    def multiple_test_correction( pvalues, meth="fdr_bh"):
        """
        Run multipletest correction and return pvalues

        @param pvalues : list - list of pvalues to be corrected. All items in the list need to be numeric.
        @param meth : string - name of pvalue correction (adjustment) method. 
        
        @return : list - list of corrected pvalues (in same order as input)
        
        Choice of methods:
            `bonferroni` : one-step correction
            `sidak` : one-step correction
            `holm-sidak` : step down method using Sidak adjustments
            `holm` : step-down method using Bonferroni adjustments
            `simes-hochberg` : step-up method  (independent)
            `hommel` : closed method based on Simes tests (non-negative)
            `fdr_bh` : Benjamini/Hochberg  (non-negative)
            `fdr_by` : Benjamini/Yekutieli (negative)
            `fdr_tsbh` : two stage fdr correction (non-negative)
            `fdr_tsbky` : two stage fdr correction (non-negative)
        
        See more information here: 
        https://www.statsmodels.org/stable/generated/statsmodels.sandbox.stats.multicomp.multipletests.html#statsmodels.sandbox.stats.multicomp.multipletests
        
        """

        from statsmodels.stats.multitest import multipletests

        return multipletests(pvalues, method=meth)[1]


    @staticmethod
    def qvalue(pv, m=None, verbose=False, lowmem=False, pi0=None):
        """
        Copyright of alexisboukouvalas & nfusi (Nicolo Fusi) https://github.com/nfusi/qvalue
        
        29-March-2019

        Function taken from https://github.com/nfusi/qvalue/blob/master/qvalue/qvalue.py and renamed

        =====    
        Estimates q-values from p-values
        Args
        =====
        m: number of tests. If not specified m = pv.size
        verbose: print verbose messages? (default False)
        lowmem: use memory-efficient in-place algorithm
        pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
             For most GWAS this is not necessary, since pi0 is extremely likely to be
             1
        """
        
        from scipy import interpolate
        import scipy as sp
        
        assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"
    
        original_shape = pv.shape
        pv = pv.ravel()  # flattens the array in place, more efficient than flatten()
    
        if m is None:
            m = float(len(pv))
        else:
            # the user has supplied an m
            m *= 1.0
    
        # if the number of hypotheses is small, just set pi0 to 1
        if len(pv) < 100 and pi0 is None:
            pi0 = 1.0
        elif pi0 is not None:
            pi0 = pi0
        else:
            # evaluate pi0 for different lambdas
            pi0 = []
            lam = sp.arange(0, 0.90, 0.01)
            counts = sp.array([(pv > i).sum() for i in sp.arange(0, 0.9, 0.01)])
            for l in range(len(lam)):
                pi0.append(counts[l]/(m*(1-lam[l])))
    
            pi0 = sp.array(pi0)
        
            # fit natural cubic spline
            tck = interpolate.splrep(lam, pi0, k=3)
            pi0 = interpolate.splev(lam[-1], tck)
                        
            if verbose:
                print("qvalues pi0=%.3f, estimated proportion of null features " % pi0)
    
            if pi0 > 1:
                if verbose:
                    print("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
                pi0 = 1.0
    
        assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0
    
        if lowmem:
            # low memory version, only uses 1 pv and 1 qv matrices
            qv = sp.zeros((len(pv),))
            last_pv = pv.argmax()
            qv[last_pv] = (pi0*pv[last_pv]*m)/float(m)
            pv[last_pv] = -sp.inf
            prev_qv = last_pv
            for i in range(int(len(pv))-2, -1, -1):
                cur_max = pv.argmax()
                qv_i = (pi0*m*pv[cur_max]/float(i+1))
                pv[cur_max] = -sp.inf
                qv_i1 = prev_qv
                qv[cur_max] = min(qv_i, qv_i1)
                prev_qv = qv[cur_max]
    
        else:
            p_ordered = sp.argsort(pv)
            pv = pv[p_ordered]
            qv = pi0 * m/len(pv) * pv
            qv[-1] = min(qv[-1], 1.0)
    
            for i in range(len(pv)-2, -1, -1):
                qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])
    
            # reorder qvalues
            qv_temp = qv.copy()
            qv = sp.zeros_like(qv)
            qv[p_ordered] = qv_temp
    
        # reshape qvalues
        qv = qv.reshape(original_shape)
        return qv


    @staticmethod
    def qvalue_lambda(pv, lam_report = 0.5, m=None, verbose=False, lowmem=False, pi0=None):
        """
        Diogo Ribeiro modified version to report pi0 at a given lambda.
        
        Copyright of alexisboukouvalas & nfusi (Nicolo Fusi) https://github.com/nfusi/qvalue
        
        29-March-2019

        Function taken from https://github.com/nfusi/qvalue/blob/master/qvalue/qvalue.py and renamed

        =====    
        Estimates q-values from p-values
        Args
        =====
        m: number of tests. If not specified m = pv.size
        verbose: print verbose messages? (default False)
        lowmem: use memory-efficient in-place algorithm
        pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
             For most GWAS this is not necessary, since pi0 is extremely likely to be
             1
        """
        
        from scipy import interpolate
        import scipy as sp
        
        assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"
    
        original_shape = pv.shape
        pv = pv.ravel()  # flattens the array in place, more efficient than flatten()
    
        if m is None:
            m = float(len(pv))
        else:
            # the user has supplied an m
            m *= 1.0
    
        # if the number of hypotheses is small, just set pi0 to 1
        if len(pv) < 100 and pi0 is None:
            pi0 = 1.0
        elif pi0 is not None:
            pi0 = pi0
        else:
            # evaluate pi0 for different lambdas
            pi0 = []
            lam = sp.arange(0, 0.90, 0.01)
            counts = sp.array([(pv > i).sum() for i in sp.arange(0, 0.9, 0.01)])
            for l in range(len(lam)):
                pi0.append(counts[l]/(m*(1-lam[l])))
    
            pi0 = sp.array(pi0)
            
            pi0_lam = pi0[int(lam_report * 100)]
            if pi0_lam > 1:
                pi0_lam = 1
        
            # fit natural cubic spline
            tck = interpolate.splrep(lam, pi0, k=3)
            pi0 = interpolate.splev(lam[-1], tck)
                        
            if verbose:
                print("qvalues pi0=%.3f, estimated proportion of null features " % pi0)
    
            if pi0 > 1:
                if verbose:
                    print("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
                pi0 = 1.0
    
        assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0
    
        if lowmem:
            # low memory version, only uses 1 pv and 1 qv matrices
            qv = sp.zeros((len(pv),))
            last_pv = pv.argmax()
            qv[last_pv] = (pi0*pv[last_pv]*m)/float(m)
            pv[last_pv] = -sp.inf
            prev_qv = last_pv
            for i in range(int(len(pv))-2, -1, -1):
                cur_max = pv.argmax()
                qv_i = (pi0*m*pv[cur_max]/float(i+1))
                pv[cur_max] = -sp.inf
                qv_i1 = prev_qv
                qv[cur_max] = min(qv_i, qv_i1)
                prev_qv = qv[cur_max]
    
        else:
            p_ordered = sp.argsort(pv)
            pv = pv[p_ordered]
            qv = pi0 * m/len(pv) * pv
            qv[-1] = min(qv[-1], 1.0)
    
            for i in range(len(pv)-2, -1, -1):
                qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])
    
            # reorder qvalues
            qv_temp = qv.copy()
            qv = sp.zeros_like(qv)
            qv[p_ordered] = qv_temp
    
        # reshape qvalues
        qv = qv.reshape(original_shape)
        return qv, pi0_lam


    @staticmethod
    def binomial_test(successes,failures, prob = 0.5):
        """
        Perform binomial using scipy.
        Reports two-way test p-value that provided data is different than 'prob'.
        
        """
        from scipy.stats import binom_test
        return binom_test([successes,failures], p = prob)


