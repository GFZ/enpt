# -*- coding: utf-8 -*-

from geomultisens.model.dataset import Dataset
from geomultisens.model.METADATA import METADATA


class _EnMAP_Image(Dataset):

    def __init__(self, pathImage=''):
        """This is the base class for all kinds of EnMAP images."""

        # get all attributes of base class "Dataset"
        super(_EnMAP_Image, self).__init__()


        # add EnMAP specific attributes
        #self.example_attr = 1


        # handle initialization arguments
        if pathImage:
            self.arr = pathImage # run the setter for 'arr' of the base class 'Dataset' which creates an Instance of GeoArray


    @property
    def MetaObj(self):
        if self._MetaObj is None:
            self._MetaObj = METADATA(self.identifier)
            #self._MetaObj.read_meta()

        return self._MetaObj


    @MetaObj.setter
    def MetaObj(self, MetaObj):
        assert isinstance(MetaObj, METADATA), "'MetaObj' can only be set to an instance of METADATA class. " \
                                              "Got %s." %type(MetaObj)
        self._MetaObj = MetaObj


    @MetaObj.deleter
    def MetaObj(self):
        self._MetaObj = None


    def from_disk(self, tuple_path_subset):
        """Fills an already instanced EnMAP_Image object with data from disk. Excludes array attributes in Python mode.

        :param tuple_path_subset:    <tuple> e.g. ('/path/gms_file.gms', ['cube', None])
        """
        pass # TODO




class EnMAP_L1B(_EnMAP_Image):
    """This class represents an EnMAP L1B image including all metadata and associated aux-data (masks, DEM, etc.).
    
    All attributes commonly used among different EnMAP images are inerited from the _EnMAP_Image class. 
    L1B specific modifications are to be implemented here."""

    pass



class EnMAP_L2A(_EnMAP_Image):
    """This class represents an EnMAP L2A image including all metadata and associated aux-data (masks, DEM, etc.).
    
    All attributes commonly used among different EnMAP images are inerited from the _EnMAP_Image class. 
    L2A specific modifications are to be implemented here."""

    pass