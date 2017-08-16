"""Handling of EnMAP Level-1B products."""

from xml.etree import ElementTree


class EnMAPL1BProduct(object):
    """Reader for EnMAP Level-1B products.

    Attributes:
        - this is all preliminary since the official format is not yet fixed -> document this if this
          becomes more stable

    """

    def __init__(self, header: str):
        """Load EnMAP Level-1B products from disk to memory.

        :param header: Parh to EnMAP Level-1B header file.
        """
        prefix = header[:-10]
        self.header = header
        self.detector1 = prefix + 'D1.bsq'
        self.detector2 = prefix + 'D2.bsq'
        self.cloudmask1 = prefix + 'D1_cloudmask.tif'
        self.cloudmask2 = prefix + 'D1_cloudmask.tif'
        self.deadpixelmap1 = prefix + 'D1_deadpixelmap.tif'
        self.deadpixelmap2 = prefix + 'D2_deadpixelmap.tif'
        self.quicklook1 = prefix + 'D1_quicklook.png'
        self.quicklook2 = prefix + 'D2_quicklook.png'
        self.xmlRoot = ElementTree.parse(self.header).getroot()

    def __repr__(self):
        """Nice string representing the EnMAP Level-1B product."""
        return '\n'.join(['{} = {}'.format(k, self.__dict__[k]) for k in sorted(self.__dict__) if k != 'xmlRoot'])

    def getHeaderNode(self, xpath: str):
        """Get node for given path.

        :param xpath: string representing a node, e.g. 'detector1/centre_wavelength'
        """
        nodes = self.xmlRoot.findall(xpath)
        if len(nodes) == 1:
            return nodes[0]
        elif len(nodes) == 0:
            raise NoMatchError(xpath)
        elif len(nodes) > 1:
            raise MultipleMatchError(xpath)


class NoMatchError(Exception):
    """Raised of given path in xml could not be matched."""

    pass


class MultipleMatchError(Exception):
    """Raised if more than one match in xml is found."""

    pass
