#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@author: misiak

"""
import uproot
import numpy as np
import re
import os
import abc


class Artifact(object):
    """ Named object. Used to define namespace and later defining attributes.
    """ 
    def __init__(self, name):
        self.name = name
    

class Copse(object):
    """ Defining the cuts in a kinda clean/handy way.
    """

    def new_cut(self, name, thruth_array):            
        setattr(self, name, thruth_array)
        

class Tree(object):
    """ Tree class, mimicking ROOT's Tree.
    """
    
    def __init__(self, root, tree_name):
        
        self._tree = root[tree_name]
        
        self._keys = list()
        for k in self._tree.keys():
            # sanitize keys and values
            key = k.decode('utf-8')
            (self._keys).append(key)
            
            value = np.array(self._tree[k].array())
            if value.ndim > 1 and value.shape[0]==1:
                value = value[0]
            setattr(self, key, value)


class Wisp(object):
    """ Preparing inheritance for Root and Thicket classes.
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self):

        self.run_tree = Artifact('run_tree')
        
        self._event_types = ('trig', 'noise')
        self._processing_types = ('raw', 'filt', 'filt_decor')   

        for etype in self._event_types:
            setattr(self, etype, Artifact(etype))
            
            for ptype in self._processing_types:
                arti = getattr(self, etype)
                setattr(arti, ptype, Artifact(ptype))
        

class Root(Wisp):
    """ Root class, mimicking ROOT's Root.
    """
    def __init__(self, file_path):
        
        Wisp.__init__(self)

        self._root = uproot.open(file_path)
        
        # general info about the run
        self.run_tree = Tree(self._root, 'RunTree_Normal')
        
        for etype in self._event_types:
            arti = getattr(self, etype)
            
            for ptype in self._processing_types:
                tree = Tree(self._root,
                            'EventTree_{}_Normal_{}'.format(etype, ptype))
                setattr(arti, ptype, tree)


class Guardian(Wisp):
    """ Concanetation of Root partitions.
    """
    def __init__(self, root_list):
        
        Wisp.__init__(self)
        
        self.roots = root_list
        
        # self.run_tree for
        # stacking all the run_tree attributes
        # by listing then numpy.stack-ing
        # and
        # self.data_types.tree_labels for
        # concatenate all the data_types.tree_labels attributes
        # by listing then numpy.concatenate-ing
        
        # initializing these attributes with empty lists
        root0 = self.roots[0]
        self.run_tree._keys = root0.run_tree._keys
        for key in self.run_tree._keys:
            setattr(self.run_tree, key, list()) 

        for etype in self._event_types:
            root_arti = getattr(root0, etype)
            guard_arti = getattr(self, etype)
            
            for ptype in self._processing_types:
                root_tree = getattr(root_arti, ptype)
                guard_tree = getattr(guard_arti, ptype)
                setattr(guard_tree, '_keys', root_tree._keys)

                for key in guard_tree._keys:                
                    setattr(guard_tree, key, list())            
        
        # appending the attributes of all root files to the empty list
        for root in self.roots:
            for key in self.run_tree._keys:
                root_attr = getattr(root.run_tree, key)
                attr_list = getattr(self.run_tree, key)
                attr_list.append(root_attr)

            for etype in self._event_types:
                root_arti = getattr(root, etype)
                guard_arti = getattr(self, etype)
                
                for ptype in self._processing_types:
                    root_tree = getattr(root_arti, ptype)
                    guard_tree = getattr(guard_arti, ptype)
                    
                    for key in guard_tree._keys:
                        root_attr = getattr(root_tree, key)
                        guard_attr = getattr(guard_tree, key)                 
                        guard_attr.append(root_attr)                      
                
        # stacking for self.run_tree
        for key in self.run_tree._keys:
            attr_list = getattr(self.run_tree, key)
            attr_array = np.stack(attr_list, axis=0)
            setattr(self.run_tree, key, attr_array)
   
        # concatenating for self.data_types.tree_labels
        for etype in self._event_types:
            arti = getattr(self, etype)
            
            for ptype in self._processing_types:
                tree = getattr(arti, ptype)
                    
                for key in tree._keys:
                    attr_list = getattr(tree, key)
                    attr_array = np.concatenate(attr_list, axis=0)
                    setattr(tree, key, attr_array)
        
            arti.nsamples = attr_array.shape[0]
                
class Analysis(object):
    """ Analysis class.
    """
    
    def __init__(self, run, detector='RED80',
                 run_dir='/home/misiak/Data/data_run57'):
        
        self.run = run
        self.detector = detector
        
        self.dir = '/'.join((run_dir, self.run, self.detector))
    
        self._files_in_data_dir = os.listdir(self.dir)
        
        self.files = []
        for fname in os.listdir(self.dir):
            if 'ProcessedData' in fname:
                self.files.append(fname)
        
        (self.files).sort()
        
        self._part_num = [re.findall('S([0-9]+)_', f)[0] for f in self.files]
        
        self._part_label = ['part_'+ num for num in self._part_num]
        
        self.roots = list()
        
        for fname, label in zip(self.files, self._part_label):
            fpath = '/'.join((self.dir, fname))
            part_root = Root(fpath)
            (self.roots).append(part_root)

        self.all = Guardian(self.roots)

        self.all.trig.ref = self.all.trig.filt_decor

        # definng quality cuts
        self.all.trig.cut = Copse()
        self.all.noise.cut = Copse()
