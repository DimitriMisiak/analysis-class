#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@author: misiak

"""
import numpy as np

from core_classes import Analysis, Artifact
from model_spectrum import fid_mixture, double_norm
#from model_spectrum import double_norm

class Analysis_red80(Analysis):
    
    def define_channels(self):
        """ Meta info concerning channels. """
        run_tree = self.all.run_tree
        
        run_tree.chan_label = np.array([
                'heat_a', 'heat_b', 'ion_a', 'ion_b', 'ion_c', 'ion_d'
        ])
        
        run_tree.chan_valid = np.array([
                0, 2, 3, 4, 5
        ])
    
        run_tree.chan_heat = np.array([0, 1])
        run_tree.chan_ion = np.array([2, 3, 4, 5])
        run_tree.chan_collect = np.array([3, 5])
        run_tree.chan_veto = np.array([2, 4])
        
        run_tree.chan_signal = np.append([0,], run_tree.chan_collect)
        
    def temporal_data_extraction(self):
        """ self explanatory """
        run_tree = self.all.run_tree
        
        # acquisition frequency
        freq_root = np.ravel((run_tree.f_max_heat, run_tree.f_max_ion))
        assert np.all(freq_root == freq_root[0])
        run_tree.freq = freq_root[0]

        # maintenance info
        maint_cycle = run_tree.MaintenanceCycle
        assert np.all(maint_cycle == maint_cycle[0])
        run_tree.maint_cycle = maint_cycle[0] / 3600 # in hours
        
        maint_duration = run_tree.MaintenanceDuration
        assert np.all(maint_duration == maint_duration[0])
        run_tree.maint_duration = maint_duration[0] / 3600 # in hours

    def _time_stamp(self, microstp_root, numpart_root):
        """ time_array in hours """
        run_tree = self.all.run_tree
        secondes_array = microstp_root / run_tree.freq
        hours_array =  secondes_array / 3600 + numpart_root
        return hours_array
  
    def _maintenance_cut(self, time_array):
        """ checking if the considered time is within a maintenance """
        run_tree = self.all.run_tree
        full_cycle = run_tree.maint_duration + run_tree.maint_cycle
        remainder_array = time_array % full_cycle
        # thruth array, True if time not in a maintenance
        return remainder_array > run_tree.maint_duration

    def quality_cut_events(self,
                           thresh_chi2_heat = 300,
                           thresh_chi2_ion = 300,
                           thresh_offset_ion = 14000 ):
        """ Define the quality cut """
        
        trig = self.all.trig
        noise = self.all.noise
        etypes = [trig, noise]
        
        for etype in etypes:
        
            if etype is noise:
                chi2 = etype.filt_decor.chi2_OF_t0        
            
            if etype is trig:
                chi2 = etype.filt_decor.chi2_OF
            
            offset = etype.raw.Off
            etype.time = self._time_stamp(etype.filt_decor.MicroStp,
                                          etype.filt_decor.NumPart)
        
            # CUT Maintenance time
            etype.cut.new_cut('maintenance',
                              self._maintenance_cut(etype.time))
        
            # CUT Chi2 heat
            etype.cut.new_cut('chi2_heat',
                              chi2[:, 0]<thresh_chi2_heat)
            
            # CUT Chi2 Ion
            etype.cut.new_cut('chi2_ion',
                              np.all(chi2[:, 2:]<thresh_chi2_ion, axis=1))    
            
            # CUT Offset Ion
            etype.cut.new_cut('offset_ion',
                              np.all(offset[:, 2:]<thresh_offset_ion, axis=1))      
            
            # CUT Quality (all cuts)
            quality_truth_array = np.all((
                    etype.cut.maintenance,
                    etype.cut.chi2_heat,
                    etype.cut.chi2_ion,
                    etype.cut.offset_ion
            ), axis=0)
            etype.cut.new_cut('quality', quality_truth_array)
            
            # number of events passing quality cuts
            etype.nsamples_quality = np.count_nonzero(etype.cut.quality)
    
    def baseline_resolution(self):
        """ Compute baseline resolution from noise events (in adu). """
        noise = self.all.noise
        run_tree = self.all.run_tree
        noise.sigma0 = Artifact('sigma0')
        
        energy = noise.filt_decor.Energy_OF_t0[noise.cut.quality]
        
        # for all the channels 
        for ind in run_tree.chan_valid:
            chan = run_tree.chan_label[ind]
            sigma = np.std(energy[:, ind])
            setattr(noise.sigma0, chan, sigma)
    
    def fiducial_cut_events(self):
        """ Computes the fiducial cut """
        trig = self.all.trig
        noise = self.all.noise
        run_tree = self.all.run_tree
        
        energy = trig.filt_decor.Energy_OF
        
        cond_veto = []
        
        # fiducial condition
        for ind in run_tree.chan_veto:
            lab = run_tree.chan_label[ind]
            sigma0 = getattr(noise.sigma0, lab)
            
            # consider cut at 2-sigma
            cond = np.abs(energy[:, ind]) < 2*sigma0
        
            cond_veto.append(cond)
        
        trig.cut.fiducial_raw = np.logical_and(*cond_veto)
        trig.cut.fiducial = np.logical_and(trig.cut.fiducial_raw, trig.cut.quality)
        
        # number of events passing quality and fiducial cut
        trig.nsamples_fiducial = np.count_nonzero(trig.cut.fiducial)
      
    def init_spectrum_model(self, model='fid'):
        
        self.model = Artifact('model')
        
        dgaussian_terms = ['double gaussian', 'double-gaussian', 'dgaussian',
                           'double norm', 'double-norm', 'dnorm',
                           'double_gaussian', 'double_norm']
        mixture_terms = ['fid', 'fid mixture', 'fid_mixture', 'fid-mixture']

        if model.lower() in dgaussian_terms:
            self.model.dist = double_norm(name='double_gaussian')     
        elif model.lower() in mixture_terms:
            self.model.dist = fid_mixture(name='fid_mixture') 
        else:
            raise Exception('\"{}\" model is not recognized.'.format(model))
        
    def sensitivity_estimation(self):
        """ Estimate the sesitivity with calibration peak and
        given spectrum model.
        """
        run_tree = self.all.run_tree
        trig = self.all.trig
        noise = self.all.noise
        
        self.calibration_peak = Artifact('calibration_peak')
        self.calibration_peak.energy = 10.37 * 1e3 # Ge 10.37keV
        self.calibration_peak.sigma = Artifact('sigma')
        
        self.sensitivity = Artifact('sensitivity')
        
        self.model.popt = Artifact('popt')
        self.model.pinit = Artifact('pinit')
        
        # check if the run is heat only, and choose a good cut for the data
        if np.all(run_tree.Polar_Ion == 0):
            # heat only run
            self.calibration_peak.cut_type = 'quality'
        else:
            # heat + ion run
            self.calibration_peak.cut_type = 'fiducial'
        
        cut_peak = getattr(trig.cut, self.calibration_peak.cut_type)
        energy = trig.filt_decor.Energy_OF[cut_peak]
        
        for ind in run_tree.chan_signal:
            
            lab = run_tree.chan_label[ind]
            data = energy[:, ind]
        
        # XXX define the _fitstart method in model_spectrum.py
        #        # initialization for the double-gaussian model
        #        mu01, mu02 = np.quantile(data, [0.10, 0.90])
        #        sig01 = sig02 = noise.sigma0.heat_a
        #        p0_light = [0.5, mu01, sig01, mu02, sig02]
        #        p0 = np.append(p0_light, [0, 1])
            
            # intialization for the fiducial mixture model
            mu01, mu02 = np.quantile(data, [0.10, 0.90])
            mu03 = np.mean(data)
            sig01 = sig02 = noise.sigma0.heat_a
            sig03 = abs(mu01 - mu02)
            p0_light = [0.50, mu01, sig01, mu02, sig02, 0.05, mu03, sig03]
            p0 = np.append(p0_light, [0, 1])
            
            popt = self.model.dist.fit(data, *p0_light, floc=0, fscale=1)
        
            setattr(self.model.popt, lab, popt)
            # saving also the initialization for debug purpose
            setattr(self.model.pinit, lab, p0)
        
            f, mu1, sig1, mu2, sig2 = popt[:5]
            
            gauss1 = (mu1, sig1)
            gauss2 = (mu2, sig2)
            if abs(mu1) > abs(mu2):
                mu, sig = gauss1
            else:
                mu, sig = gauss2
            
            setattr(self.calibration_peak, lab, mu)
            setattr(self.calibration_peak.sigma, lab, sig)
            
            sens = getattr(self.calibration_peak, lab) / self.calibration_peak.energy
            setattr(self.sensitivity, lab, sens)
    
    def conversion_adu_to_ev(self):
        """ Reconstructed energy in evfrom sensitivity. """
        trig = self.all.trig
        noise = self.all.noise
        run_tree = self.all.run_tree       
        
        # energy in eV, and energy in eV corrected with the sign
        trig.energy_ev = Artifact('energy_ev')
        noise.energy_ev = Artifact('energy_ev')
        noise.sigma0_ev = Artifact('sigma0_ev')
        self.calibration_peak.sigma_ev = Artifact('sigma_ev')
        
        # recovering energy in adu
        energy_adu = trig.filt_decor.Energy_OF
        noise_energy_adu = noise.filt_decor.Energy_OF_t0
        sigma0 = noise.sigma0
        sigma = self.calibration_peak.sigma
        
        # creating the sign correction used when "adding ion channels"
        proto_sign = np.concatenate((-run_tree.Sign_Chal, run_tree.Polar_Ion), axis=1)
        assert np.all(proto_sign == proto_sign[0])
        proto_sign = proto_sign[0]
        
        # no correction for negative polarisation collecting positive signal,
        # -1 correction for the positive polarisation collecting negative signal
        run_tree.sign_corr = np.where(proto_sign<0, +1, -1)
        
        for ind in run_tree.chan_signal:
            lab = run_tree.chan_label[ind]
            
            e_adu = energy_adu[:, ind]
            noise_e_adu = noise_energy_adu[:, ind]
            sig0_adu = getattr(sigma0, lab)
            sig_adu = getattr(sigma, lab)
            
            sens = getattr(self.sensitivity, lab)
            sign = run_tree.sign_corr[ind]        
        
            e_ev = e_adu / sens
            noise_e_ev = noise_e_adu / sens
            sig0_ev = sign * sig0_adu / sens
            sig_ev = sign * sig_adu / sens
        
            setattr(trig.energy_ev, lab, e_ev)
            setattr(noise.energy_ev, lab, noise_e_ev)
            setattr(noise.sigma0_ev, lab, sig0_ev)
            setattr(self.calibration_peak.sigma_ev, lab, sig_ev)

    def virtual_collect_channel(self):
        """ Virtual collect channel : sum of the collect channels """
        trig = self.all.trig
        noise = self.all.noise
        run_tree = self.all.run_tree   
        
        # labels of the channel with energy conversion in eV
        run_tree.chan_label_virtual = ['heat_a', 'collect']
        
        energy_collect = list()
        noise_collect = list()
        for ind in run_tree.chan_collect:
            lab = run_tree.chan_label[ind]
            run_tree.chan_label_virtual.append(lab)
            
            energy = getattr(trig.energy_ev, lab)
            energy_collect.append(energy)
            
            noise_energy = getattr(noise.energy_ev, lab)
            noise_collect.append(noise_energy)
        
        trig.energy_ev.collect = np.sum(energy_collect, axis=0)
        noise.energy_ev.collect = np.sum(noise_collect, axis=0)
        
        noise.sigma0_ev.collect = np.std(noise.energy_ev.collect)
        
#        ## =============================================================================
#        ## DEFINITON OF THE VIRTUAL ELECTRODE (SUMMING THE DATA)
#        ## =============================================================================
#        #    # XXX not really good to do the sum just now
#        #    # beacause the energy is still in ADU
#        #    # we need to do the sum with the energy in eV
#        #    
#        #    # for the virtual all_ion channel
#        #    energy_ion_all = np.sum(energy[:, run_tree.chan_ion], axis=0)
#        #    noise.sigma0.ion_all = np.std(energy_ion_all)
#        #    
#        #    # for the virtual veto channel
#        #    energy_veto = np.sum(energy_corr[:, run_tree.chan_veto], axis=0)
#        #    noise.sigma0.veto = np.std(energy_veto)
#        #        
    
    def __init__(self, run):
        Analysis.__init__(
                self,
                run,
                detector='RED80',
                run_dir='/home/misiak/Data/data_run57'
        )

        self.define_channels()
        self.temporal_data_extraction()
        self.quality_cut_events()
        self.baseline_resolution()
        self.fiducial_cut_events()
        self.init_spectrum_model()
        self.sensitivity_estimation()
        self.conversion_adu_to_ev()
        self.virtual_collect_channel()

if __name__ == '__main__':
    
    run = 'tg25l019'
    ana = Analysis_red80(run)