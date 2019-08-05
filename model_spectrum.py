#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 22:11:20 2019

@author: misiak
"""

import scipy.stats as st

class double_norm(st.rv_continuous):
    """ Double Gaussian distribution. """

    def _cdf(self, x, f, loc1, scale1, loc2, scale2):
        cdf1 = (1-f) * st.norm.cdf(x, loc=loc1, scale=scale1)
        cdf2 = f * st.norm.cdf(x, loc=loc2, scale=scale2)
        cdf = cdf1 + cdf2
        return cdf

    def _pdf(self, x, f, loc1, scale1, loc2, scale2):
        pdf1 = (1-f)* st.norm.pdf(x, loc=loc1, scale=scale1)
        pdf2 = f* st.norm.pdf(x, loc=loc2, scale=scale2)
        pdf = pdf1 + pdf2
        return pdf

    def _argcheck(self, *args):
        """Default check for correct values on args and keywords.
        Returns condition array of 1's where arguments are correct and
         0's where they are not.
        """
        cond = 1
#        for arg in args:
#            cond = np.logical_and(cond, (np.asarray(arg) > 0))
        return cond

class fid_mixture(st.rv_continuous):
    """ Double Gaussian distribution plus uniform distribution """

    def _cdf(self, x, f, loc1, scale1, loc2, scale2, fu, loc3, scale3):
        cdf1 = (1-fu) * (1-f) * st.norm.cdf(x, loc=loc1, scale=scale1)
        cdf2 = (1-fu) * f * st.norm.cdf(x, loc=loc2, scale=scale2)
        cdf3 = fu * st.uniform.cdf(x, loc=loc3, scale=scale3)
        cdf = cdf1 + cdf2 + cdf3
        return cdf

    def _pdf(self, x, f, loc1, scale1, loc2, scale2, fu, loc3, scale3):
        pdf1 = (1-fu) * (1-f)* st.norm.pdf(x, loc=loc1, scale=scale1)
        pdf2 = (1-fu) * f* st.norm.pdf(x, loc=loc2, scale=scale2)
        pdf3 = fu * st.uniform.pdf(x, loc=loc3, scale=scale3)
        pdf = pdf1 + pdf2 + pdf3
        return pdf

    def _argcheck(self, *args):
        """Default check for correct values on args and keywords.
        Returns condition array of 1's where arguments are correct and
         0's where they are not.
        """
        cond = 1
#        for arg in args:
#            cond = np.logical_and(cond, (np.asarray(arg) > 0))
        return cond