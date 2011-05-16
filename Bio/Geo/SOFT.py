# Copyright 2011 Phillip Garland  <pgarland@gmail.com> All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# For more documentation of the GEO SOFT file format, see:
# http://www.ncbi.nlm.nih.gov/geo/info/soft2.html#SOFTformat

# Would be useful to have a function for fetching platforms/series/samples
# referenced in a file's metadata.

import string
import re

from utils import _read_key_value, stringIsType, maybeConvertToNumber

# The SOFT class should never be directly instantiated; it just contains
# variables and methods common to some or all of it's subclasses (GSM, GSE, GPL,
# and GDS)
class SOFT(object):
    def __init__(self):
        # These attributes are used to generate a dictionary to hold the values
        # of the attributes found a file. 'count' gives the possible number of
        # occurances of an attribute in a file. A literal number, e.g. '0' or
        # '1' stands for itself, '|' means 'or' and "+" means more. So, for
        # example, '1|+' means that the attribute occurs on 1 or more lines.
        self._SOFTattributes = {'geo_accession': {'count': '0|+'},
                                'status': {'count': '0|+'},
                                'submission_date': {'count': '0|+'},
                                'last_update_date': {'count': '0|+'},
                                'row_count': {'count': '0|+', 'processFunc': int},
                                'contact_name': {'count': '0|+'},
                                'contact_email': {'count': '0|+'},
                                'contact_phone': {'count': '0|+'},
                                'contact_department': {'count': '0|+'},
                                'contact_institute': {'count': '0|+'},
                                'contact_city': {'count': '0|+'},
                                'contact_phone': {'count': '0|+'},
                                'contact_fax': {'count': '0|+'},
                                'contact_web_link': {'count': '0|+'}}
        # Each subclass has it's own attibutes. _attributes will hold all
        # attributes defined in _SOFTattibutes, as well as the attibutes defined
        # in the subclasses.
        self._attributes = {}
        # Some submitters include attributes not specified by NCBI. Biopython
        # users can check this if they want to access anything nonstandard.
        self.nonStdAttributes = []
        # columns has to be a list rather than a dictionary to preserve the
        # order of the columns. The ordering of columns is important because
        # once a delimiter is specified for a column, that is the delimiter for
        # adjacent columns until another delimiter is specified.
        self.columns = []
        # The mapping from column names to column numbers
        self.columnDict = {}
        # The delimiter within each column in the data table, if any
        self._delimiters = []
        # The entitie's metadata
        self.meta = {}
        # The mapping from ID to table row number
        self.idDict = {}
        # The data table. This will be a 2-D array indexed by column and row.
        self._table = []

    # Prints the main metadata associated with a file
    def __str__(self):
        s = ''
        for attr in self.meta.keys():
            if self.meta[attr]:
                if isinstance(self.meta[attr], list):
                    s += attr + ": "
                    for i in range(min(len(self.meta[attr]), 10)):
                        s += str(self.meta[attr][i]) + "\n"
                    if len(self.meta[attr]) > 10:
                        s += '...' + "\n"
                else:
                    s += attr + ": "  + str(self.meta[attr]) + "\n"
        return s

    # Used by GSM, GPL, and GSE subclasses to set default values for each
    # attribute
    def _setMetaDefaultValues(self):
        for attr in self._attributes.keys():
            # If this attribute can occur on multiple lines, make it's value a list
            if self._attributes[attr]['count'][-1] == '+':
                self.meta[attr] = []
            else:
                self.meta[attr] = None

    # Tests for each of the 4 line types found in GEO SOFT files There are 4
    # line types in GEO SOFT files:
    # * entity indicators
    # * entity attributes
    # * data table head descriptions
    # * data table rows

    # Each line type is distinguished by its first character. Entity indicator,
    # entity attribute, and table header description lines each start with a
    # unique character, while data table rows are any line that does not start
    # with any of the three characters used to mark the other line types.
    @classmethod
    def _isEntityIndicator(cls, line):
        return line[0] == '^'
    @classmethod
    def _isEntityAttribute(cls, line):
        return line[0] == '!'
    @classmethod
    def _isTableHeaderDescription(cls, line):
        return line[0] == '#'

    @classmethod
    def _isDataTableRow(cls, line):
        return (not SOFT._isEntityIndicator(line) and \
                    (not SOFT._isEntityAttribute(line)) and \
                    (not SOFT._isTableHeaderDescription(line)))

    # Each attribute label is prefixed with the '!' character and the entity
    # type and a '_' character - e.g. '!Sample_'. The string following the
    # prefix is the attribute label's name. This function extracts the name from
    # the label.
    @classmethod
    def _getAttributeLabelName(cls, attributeLabel):
        """Return the label's name, with the entity name removed."""
        return attributeLabel.split('_', 1)[1]

    @classmethod
    def _splitAttributeLine(cls, line):
        label, value = _read_key_value(line)
        label = SOFT._getAttributeLabelName(label)
        return label, value

    # Table header description lines can have a more complicated layout that the
    # other 3 line types. They usually contain a description of the values found in
    # the columns, but they can also contain 'linking commands' for generating URLs
    # that point to more information about each entry in the column, or for
    # separating subvalues within in a column.
    @classmethod
    def _parseTableHeaderValue(cls, value):
        """Parse the value portion of a data table header line- a line that
        describes a column in the following data table. It returns a list containing
        the column's description and a dictionary with the information needed for
        creating URLs from the values in the column or using an alternative
        delimiter"""

        cmds = {'LINK_PRE': '', 'LINK_SUF': '', 'DELIMIT': None}
    
        # Capture the description, as well as any remaining portion following the
        # description
        descriptionMatch = re.search('(.+?)((?:LINK_PRE|LINK_SUF|DELIMIT):\".+?\")*$', value)

        # There are some files with table header description lines that do not
        # contain descriptions. In this case, make the description an empty string.

        # In the future I might have a 'stricter' mode that flags mistakes like
        # this, and a 'strictest' mode that throws an exception.
        if descriptionMatch == None:
            description = ''
        else:
            description = descriptionMatch.group(1)
            # group(2) contains anything following the description
            if descriptionMatch.group(2) != None:
                # find all the commands
                linkCmdPairs = re.findall('(LINK_PRE|LINK_SUF|DELIMIT):\"(.+?)\"', descriptionMatch.group(2))
                for pair in linkCmdPairs:
                    cmdType, cmdString = pair
                    cmds[cmdType] = cmdString
        return [description, cmds]

    @classmethod
    def _makeLinkFn(cls, pre, suf):
        """Creates a function that generates a hyperlink from a data table value"""
        return lambda(tableValue): pre + tableValue + suf

    # Read in the entity- the metadata about the platform, series, or sample. The
    # GDS class overrides this method because Datasets contain multiple
    # entities, where as platforms (GPL), series (GSE), and samples (GSM)
    # contain a single entity
    def readEntities(self, handle):
        # In GSM, GSE, and GPL files, there should be a entity indicator only on
        # the first line, so we process that line, then break.
        for line in handle:
            if SOFT._isEntityIndicator(line):
                entityType, self.entityID = _read_key_value(line)
                break

        for line in handle:
            line = line.strip('\n').strip('\r')

            if SOFT._isEntityAttribute(line):
                # Normal label-value line; add the attribute's value to the
                # metadata
                label, value = SOFT._splitAttributeLine(line)
                # If there is a processFunc for further parsing the value, call
                # it
                try:
                    value = self._attributes[label]['processFunc'](value)
                except KeyError:
                    # Don't do any additonal processing on value
                    next
                try:
                    if isinstance(self.meta[label], list):
                        self.meta[label].append(value)
                    else:
                        self.meta[label] = value
                except KeyError: # Catch non standard attributes
                    self.nonStdAttributes.append(label)                    
                    # Since these attributes aren't standardized by NCBI, we
                    # can't know if there are 1 or many, so we shove them in a list

                    # channel attributes in GSM sample files are also caught
                    # here. This isn't a problem since all channel attributes
                    # can potentially occur on more than one line.
                    self.meta[label] = [value]
                if label == 'data_row_count':
                        # Last line before data table header descriptions in GSM
                        # and GPL files. GSE files just end; they don't contain
                        # a data table.
                        return

    # Handle the data table descriptions; this is used by the GPL, GSM, and GDS
    # subclasses. GSE (Series) do not contain a data table, so the GSE class
    # doesn't need it.
    def readTableHeaderDescriptions(self, handle):
        # Table header descriptions are after entity attributes (in GSM,
        # GSE, and GPL files)
        columnNum = 0
        for line in handle:
            if not SOFT._isTableHeaderDescription(line):
                break
            else:
                columnName, value = _read_key_value(line)
                description, linkCmds = SOFT._parseTableHeaderValue(value)

                # I don't think the behaviour is specified, but at least one
                # platform file (GPL97) specifies a delimiter on one line, then
                # continues to use that delimiter on subsequent lines, so we
                # carry the previous delimiter forward here
                if (linkCmds['DELIMIT'] == None) and (columnNum != 0):
                    linkCmds['DELIMIT'] = self._delimiters[columnNum -1]

                self._delimiters.append(linkCmds['DELIMIT'])
                self.columns.append({'name': columnName,
                                     'description': description,
                                     'link': SOFT._makeLinkFn(linkCmds['LINK_PRE'],
                                                         linkCmds['LINK_SUF'])})
                # Add to the mapping between column names and column
                # numbers
                self.columnDict[columnName] = columnNum
                columnNum += 1

    # Handle the data table. This is used by the GPL, GSM, and GDS
    # subclasses. GSE (Series) do not contain a data table, so the GSE class
    # doesn't need it.
    def readTable(self, handle):
        rowNumber = 0

        for line in handle:
            line = line.strip('\n').strip('\r')
            if not line: continue

            if SOFT._isEntityAttribute(line):
                if SOFT._getAttributeLabelName(line) == 'table_end':
                    return # Reached the end of the data table
            elif SOFT._isDataTableRow(line):
                row = line.split("\t")
                for columnNum in range(len(self.columns)):
                    if self._delimiters[columnNum]: # Split the fields with the column
                        row[columnNum] = re.split(self._delimiters[columnNum], row[columnNum])

                self._table.append(row)
                # Add to the mapping between ID and row number
                self.idDict[row[0]] = rowNumber
                rowNumber += 1

    def getTableValue(self, id, column):
        idNum = self.idDict[id]
        columnNum = self.columnDict[column]
        value = self._table[idNum][columnNum]

        if isinstance(value, list):
            return [maybeConvertToNumber(elt) for elt in value]
        else:    
            return maybeConvertToNumber(value)

class GPL(SOFT):
    def __init__(self):
        super(GPL, self).__init__()
        self._GPLattributes= {'title': {'count': '1'},
                             'distribution': {'count': '1'},
                             'technology': {'count': '1'},
                             'organism': {'count': '1|+'},
                             'manufacturer': {'count': '1'},
                             'manufacture_protocol': {'count': '1|+'},
                             'catalog_number': {'count': '0|+'},
                             'web_link': {'count': '0|+'},
                             'support': {'count': '0|1'},
                             'coating': {'count': '0|1'},
                             'description': {'count': '0|+'},
                             'contributor': {'count': '0|+'},
                             'pubmed_id': {'count': '0|+'},
                             'geo_accession': {'count': '0|1'}}
        # Merge the attributes of the superclass and this class
        self._attributes = dict(self._SOFTattributes, **self._GPLattributes)
        self._setMetaDefaultValues()

    def read(self, fileName):
        handle = open(fileName)
        self.readEntities(handle)
        self.readTableHeaderDescriptions(handle)
        # Between the table header descriptions and the data table itself,
        # there's a line that contains the column names. Since the information
        # in this column is redundant with table header descriptions, we skip
        # past it.
        handle.next()
        self.readTable(handle)
        handle.close()

class GSM(SOFT):
    def __init__(self):
        super(GSM, self).__init__()
        self._GSMattributes = {'title': {'count': '1'},
                               'supplementary_file': {'count':'1|+'},
                               'table': {'count': '0|1'},
                               'hyb_protocol': {'count': '1|+'},
                               'scan_protocol': {'count': '1|+'},
                               'data_processing': {'count': '1|+'},
                               'description': {'count': '0|+'},
                               'platform_id': {'count': '1'},
                               'geo_accession': {'count': '1'},
                               'channel_count': {'count': '1', 'processFunc': int},
                               # The following four attributes are used only for
                               # SAGE (Serial Analysis of Gene Expression)
                               # submissions. Maybe ; 'count' should be '0|1'?
                               'anchor': {'count': '1'},
                               'type': {'count': '1'},
                               'tag_count': {'count': '1'},
                               'tag_length': {'count': '1'}}
        self._channelAttributes = {'source_name': {'count': '1|+'},
                                   'organism': {'count': '1|+'},
                                   'characteristics': {'count': '1|+'},
                                   'biomaterial_provider': {'count': '0|+'},
                                   'treatment_protocol': {'count': '0|+'},
                                   'growth_protocol': {'count': '0|+'},
                                   'molecule': {'count': '1|+'},
                                   'extract_protocol': {'count': '1|+'},
                                   'label': {'count': '1|+'},
                                   'label_protocol': {'count': '1|+'}}
        # Merge the attributes of the superclass and this class
        self._attributes = dict(self._SOFTattributes, **self._GSMattributes)
        self._setMetaDefaultValues()

    def read(self, fileName):
        self.readEntities(fileName)
        handle = open(fileName)
        self.readEntities(handle)
        self.readTableHeaderDescriptions(handle)
        # Between the table header descriptions and the data table itself,
        # there's a line that contains the column names. Since the information
        # in this column is redundant with table header descriptions, we skip
        # past it.
        handle.next()        
        self.readTable(handle)
        handle.close()

        # Group the attributes by channel by filling meta.['channels']
        channel_count = int(self.meta['channel_count'])
        # Channels seem to be labeled starting at 1, so we fill self.meta['channels'][0] with None
        self.meta['channels'] = [None]
        for channel in range(1, (channel_count  + 1)):
            channelDict = {}
            for attr in self._channelAttributes.keys():
                try:
                    channelDict[attr] = self.meta[attr + '_ch' + str(channel)]
                except KeyError:
                    channelDict[attr] = None
            self.meta['channels'].append(channelDict)
        # Remove NCBI-defined channel attributes from nonStdAttributes:
        # First generate all possible NCBI-defined channel attributes for this file
        stdChanAttrs = [attr + '_ch' + str(chan) for attr in self._channelAttributes.keys() for chan in range(1, channel_count + 1)]
        # Then remove standard channel attributes from nonStdAttributes
        self.nonStdAttributes = [attr for attr in self.nonStdAttributes if attr not in stdChanAttrs]

class GSE(SOFT):
    def __init__(self):
        super(GSE, self).__init__()
        self._GSEattributes = {'title': {'count': '1'},
                              'summary': {'count': '1|+'},
                              'overall_design': {'count': '1'},
                              'pubmed_id': {'count': '0|+'},
                              'web_link': {'count': '0|+'},
                              'contributor': {'count': '0|+'},
                              'sample_id': {'count': '1|+'},
                              'geo_accession': {'count': '0|1'},
                              'type': {'count': '1'}}
        # Merge the attributes of the superclass and this class
        self._attributes = dict(self._SOFTattributes, **self._GSEattributes)
        self._setMetaDefaultValues()

    def read(self, filename):
        """Read a GSE file"""
        handle = open(filename)
        self.readEntities(handle)
        self.readTableHeaderDescriptions(handle)
        # GSE files do not contain a table, so we just read the metadata
        handle.close()

class GDS(SOFT):
    def __init__(self):
        super(GDS, self).__init__()
        # These attributes are taken from GDS files. 
        self._attributes = {'database': {'entityValue': {},
                                         'name': {},
                                         'institute': {},
                                         'web_link': {},
                                         'email': {},
                                         'ref': {}},
                            # The dataset entry will be deleted after
                            # parsing. Use self.meta['attribute'] to get these
                            # values.
                            'dataset': {'type': {},
                                        'pubmed_id': {},
                                        'platform': {},
                                        'platform_organism': {},
                                        'platform_technology_type': {},
                                        'feature_count': {'processFunc': int},
                                        'sample_organism': {},
                                        'sample_type': {},
                                        'channel_count': {'processFunc': int},
                                        'sample_count': {'processFunc': int},
                                        'value_type': {},
                                        'reference_series': {},
                                        'order': {},
                                        'update_date': {}},
                            # I don't understand why the platform* attributes
                            # exist- is there any reason for them to have
                            # different values than the same attributes in
                            # 'dataset'?
                            'annotation': {'date': {},
                                           'platform': {},
                                           'platform_title': {},
                                           'platform_organism': {}},
                            # The subset entry will be deleted after
                            # parsing. Use self.subsets to get these value.
                            'subset':  {'subset': {},
                                        'dataset_id': {},
                                        'description': {},
                                        'sample_id': {},
                                        'type': {}}}
        
        # Initialize the attributes
        self.entities = {'database': {}, 'dataset': {}, 'subset': {}, 'annotation': {}}
        for attr in self._attributes['database'].keys():
            self.entities['database'][attr] = None
        for attr in self._attributes['dataset'].keys():
            self.entities['dataset'][attr] = None
        for attr in self._attributes['annotation'].keys():
            self.entities['annotation'][attr] = None
        for attr in self._attributes['subset'].keys():
            self.entities['subset'][attr] = None
        # FIXME: Maybe make this a dictionary, indexed by the subset's name
        self.subsets = []

    def readEntities(self, handle):
        """Read in all the entities from a GDS file"""
        entity = ''
        dataSetEntityCount = 0

        line = handle.next()
        while line:
            line = line.strip('\n').strip('\r')
            if not line: continue
           
            if SOFT._isEntityIndicator(line):
                entity, value = _read_key_value(line)
                entity = string.lower(entity)

                if entity == 'dataset':
                    dataSetEntityCount += 1
                    # In GDS files, a second '^DATASET = GDS...' line occurs
                    # immediately before the data table header description
                    # lines.
                    if dataSetEntityCount == 2:
                        # Put the dataset attributes in a more standard place
                        self.meta = self.entities['dataset']
                        del self.entities['dataset']
                        # Get rid of the temporary subset store
                        del self.entities['subset']
                        break

                self.entities[entity]['entityValue'] = value
                line = handle.next()
                while SOFT._isEntityAttribute(line):
                    label, value = SOFT._splitAttributeLine(line)
                    try:
                        value = self._attributes[entity][label]['processFunc'](value)
                    except KeyError:
                        # Don't do any additonal processing on value
                        next
                    self.entities[entity][label] = value        
                    line = handle.next()
                if entity == 'subset':
                    self.subsets.append(self.entities['subset'])

    def readMeta(self, fileName):
        """Read in the entities and data table header descriptions from a GDS file"""
        handle = open(fileName)
        self.readEntities(handle)
        self.readTableHeaderDescriptions(handle)
        handle.close()
                            
    def read(self, fileName):
        """Read in the metadata and data table from a GDS file"""
        handle = open(fileName)
        self.readEntities(handle)
        self.readTableHeaderDescriptions(handle)
        handle.next()
        self.readTable(handle)
        handle.close()
