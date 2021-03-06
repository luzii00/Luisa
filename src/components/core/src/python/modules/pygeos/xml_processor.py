from lxml import etree as ElementTree
from lxml.etree import XMLSyntaxError 
import re
import sys
import os
from numpy import *
from . import UnitManager, DictRegexHandler, symbolicMathRegexHandler, regexConfig


def MergeIncludedXMLFiles(root, fname, includeCount, maxInclude=100):
  # Expand the input path
  pwd = os.getcwd()
  includePath, fname = os.path.split(os.path.abspath(os.path.expanduser(fname)))
  os.chdir(includePath)

  # Check to see if the code has fallen into a loop
  includeCount += 1
  if (includeCount > maxInclude):
    raise Exception('Reached maximum recursive includes...  Is there an include loop?')

  # Check to make sure the file exists
  if (not os.path.isfile(fname)):
    print('Included file does not exist: %s' % (fname))
    raise Exception('Check included file path!')

  # Load target xml
  try:
    parser = ElementTree.XMLParser(remove_comments=True, remove_blank_text=True)
    includeTree = ElementTree.parse(fname, parser)
    includeRoot = includeTree.getroot()
  except XMLSyntaxError as err:
    print('\nCould not load included file: %s' % (fname))
    print err.msg
    raise Exception('\nCheck included file!')

  # Recursively add the includes:
  for includeNode in includeRoot.findall('Included'):
    for f in includeNode.findall('File'):
      MergeIncludedXMLFiles(root, f.get('name'), includeCount)

  # Merge the results into the xml tree
  for topLevelNode in list(includeRoot):
    rootMatchingNodes = root.findall(topLevelNode.tag)
    if (rootMatchingNodes):
      for secondLevelNode in list(topLevelNode):
        rootMatchingNodes[0].insert(-1, secondLevelNode)
    else:
      root.insert(-1, topLevelNode)
  os.chdir(pwd)


def generateRandomName(prefix='', suffix='.xml'):
  from hashlib import md5
  from time import time
  from os import getpid

  return '%s%s%s' % (prefix, md5(str(time())+str(getpid())).hexdigest(), suffix)


def PreprocessGEOSXML(inputFile, schema='/g/g17/sherman/GEOS/geosx/src/components/core/src/schema/gpac_new.xsd', verbose=1):
  
  if (verbose > 0):
    print('\nReading input xml parameters and parsing symbolic math...')

  # Expand the input path
  pwd = os.getcwd()
  rootPath, inputFile = os.path.split(os.path.abspath(os.path.expanduser(inputFile)))
  os.chdir(rootPath)

  # Load the xml files and merge includes
  try:
    parser = ElementTree.XMLParser(remove_comments=True, remove_blank_text=True)
    tree = ElementTree.parse(inputFile, parser=parser)
    root = tree.getroot()
  except XMLSyntaxError as err:
    print('\nCould not load input file: %s' % (inputFile))
    print err.msg
    raise Exception('\nCheck input file!')

  includeCount = 0
  for includeNode in root.findall('Included'):
    for f in includeNode.findall('File'):
      MergeIncludedXMLFiles(root, f.get('name'), includeCount)
  os.chdir(pwd)

  # Build the parameter map, convert function
  Pmap = {}
  for parameters in root.findall('Parameters'):
    for p in parameters.findall('Parameter'):
      Pmap[p.get('name')] = p.get('value')
  tmp_fname_a = generateRandomName(prefix='int_')
  tree.write(tmp_fname_a, pretty_print=True)
  parameterHandler = DictRegexHandler()
  parameterHandler.target = Pmap

  # Parse the raw xml file:  
  tmp_fname_b = generateRandomName(prefix='prep_')
  unitManager = UnitManager()
  with open(tmp_fname_a, 'r') as ifile, open(tmp_fname_b, 'w') as ofile:
    for line in ifile:
      # Fill in any paramters (format:  $Parameter or $:Parameter)
      if ('$' in line):
        line = re.sub(regexConfig.parameters, parameterHandler, line)

      # Parse any units (format: 9.81[m**2/s] or 1.0 [bbl/day])
      if ('[' in line):
        line = re.sub(regexConfig.units, unitManager.regexHandler, line)

      # Evaluate symbolic math (format: {1 + 2.34e5*2 * ...})
      if ('{' in line):
        line = re.sub(regexConfig.symbolic, symbolicMathRegexHandler, line)
      
      ofile.write(line)

  # Check for un-matched special characters
  os.remove(tmp_fname_a)
  with open(tmp_fname_b, 'r') as ofile:
    for line in ofile:
      if any([sc in line for sc in ['$', '[', ']', '{', '}']]):
        raise Exception('Found un-matched special characters in the pre-processed input file on line:\n%s\n Check your input xml for errors!' % (line))

  if (verbose > 0):
    print('Preprocessed xml file stored in %s\n' % (tmp_fname_b))

    # Validate against the schema 
    print('Validating the xml against the schema...')
    try:
      ofile = ElementTree.parse(tmp_fname_b) 
      sfile = ElementTree.XMLSchema(ElementTree.parse(schema))
      sfile.assertValid(ofile)
    except ElementTree.DocumentInvalid as err:
      print('\nWarning: input XML contains potentially invalid input parameters:')
      print('-'*20+'\n')
      print sfile.error_log
      print('\n'+'-'*20)
      print('(Total schema warnings: %i)\n' % (len(sfile.error_log)))

  return tmp_fname_b

if (__name__ == "__main__"):
  ifile = sys.argv[1]
  PreprocessGEOSXML(ifile)