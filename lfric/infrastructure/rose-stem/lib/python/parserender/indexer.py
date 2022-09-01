#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Examine the output directory from a nightly suite run and determine things
about its contents.
'''
from __future__ import print_function

from __future__ import absolute_import
from abc import ABCMeta, abstractmethod
import datetime
import glob
import hashlib
import os
import os.path
import xml.etree.ElementTree as et

##############################################################################
class IndexerException(Exception):
  pass

##############################################################################
class TreeNode(object):
  __MetaClass__ = ABCMeta

  def __init__( self ):
    self._parent = None
    self._children = {}

  @abstractmethod
  def getType( self ):
    pass

  def setParent( self, parent ):
    self._parent = parent

  def getParent( self ):
    return self._parent

  def addChild( self, child ):
    self._children[child.name] = child
    child.setParent( self )

  def getChild( self, name ):
    if name in self._children:
      return self._children[name]
    else:
      return None

  def getChildren( self ):
    return list(self._children.values())

  def containsChild( self, name ):
    for node in self._children.values():
      if node.name == name:
        return True
    return False

  def descend( self ):
    def descendRecursor( node, breadcrumbs ):
      yield node, breadcrumbs
      children = node.getChildren()
      if children is not None:
        breadcrumbs.append(node)
        pages = []
        directories = []
        for child in children:
          if isinstance(child, Page):
            pages.append( child )
          elif isinstance(child, Directory):
            directories.append( child )
          else:
            raise IndexerException( 'Unknown file node type: ' + child.__class__ )
        for page in sorted(pages):
          for x in descendRecursor( page, breadcrumbs ): yield x
        for directory in sorted(directories):
          for x in descendRecursor( directory, breadcrumbs ): yield x
        breadcrumbs.pop()

    breadcrumbs = []
    for x in descendRecursor( self, breadcrumbs ): yield x

##############################################################################
class Root(TreeNode):
  def __init__( self ):
    super(Root,self).__init__()

  def getType( self ):
    return 'root'

##############################################################################
class FileNode(object):
  __MetaClass__ = ABCMeta

  def __init__( self, name, timestamp ):
    super(FileNode, self).__init__()
    self.name      = name
    self._timestamp = timestamp

  @abstractmethod
  def urlPath( self ):
    pass

  def getTimestamp( self ):
    return self._timestamp

##############################################################################
class Directory(FileNode, TreeNode):
  '''
  Holds details about a directory.
  '''
  def __init__( self, name, timestamp ):
    super(Directory, self).__init__( name, timestamp )
    self._indexPage = None
    self._buildPage = None
    self._runPage   = None

  def getType( self ):
    return 'directory'

  def urlPath( self ):
    path = []
    node = self
    while node is not None:
      path.insert( 0, node.name )
      node = node.getParent()
      if isinstance( node, Root ): node = None
    return '/'.join( path )

  def setIndexPage( self, name ):
    self._indexPage = self.getChild( name )
    if self._indexPage is None:
      raise IndexerException( 'Index page does not exist: ' + name )
    self._indexPage.setHidden()

  def getIndexPage( self ):
    return self._indexPage

  def setBuildPage( self, fileObject ):
    self._buildPage = fileObject
    self._buildPage.setHidden()

  def getBuildPage( self ):
      return self._buildPage

  def setRunPage( self, fileObject ):
    self._runPage = fileObject
    self._runPage.setHidden()

  def getRunPage( self ):
    return self._runPage

  def __eq__(self, other):
      return self.name == other.name

  def __lt__(self, other):
      return self.name < other.name

##############################################################################
class Page(FileNode, TreeNode):
  '''
  Holds details about a page.
  '''
  def __init__( self, name, timestamp, digest ):
    super(Page, self).__init__( name, timestamp )
    self._currentDigest  = digest
    self._previousDigest = None
    self._hidden         = False

  def addChild( self, child ):
    raise IndexerException( 'Can''t add children to a page' )

  def getType( self ):
    return 'page'

  def urlPath( self ):
    return self.getParent().urlPath() + '/' + self.name

  def setHidden( self ):
    self._hidden = True

  def isHidden( self ):
    return self._hidden

  def getCurrentDigest( self ):
    return self._currentDigest

  def setPreviousDigest( self, digest ):
    self._previousDigest = digest

  def getPreviousDigest( self ):
    return self._previousDigest

  def hasChanged( self ):
    if self._currentDigest is None or self._previousDigest is None:
      return None
    return self._currentDigest != self._previousDigest

  def __eq__(self, other):
      return self.name == other.name

  def __lt__(self, other):
      return self.name < other.name

##############################################################################
class Indexer(object):
  '''
  Indexers inherit from this class.
  '''
  __MetaClass__ = ABCMeta

##############################################################################
class LFRicIndexer(Indexer):
  '''
  Examine an LFRic site.
  '''
  ############################################################################
  def __init__( self, pathname, ignorePrefixes=[] ):
    self._rootPath         = pathname
    self._ignorePrefixes   = ignorePrefixes
    self._tree             = None

  ############################################################################
  def addPreviousHash( self, filename, digest ):
    node = self._tree
    path, leaf = filename.rsplit( os.path.sep, 1 )
    for directory in path.split( os.path.sep ):
      node = node.getChild( directory )
    node.getChild( leaf ).setPreviousDigest( digest )

  ############################################################################
  def getPreviousHashes( self ):
    def _recurseHashes( node ):
      if isinstance(node, Page):
        digest = node.getPreviousDigest()
        if digest:
          yield node.urlPath(), digest
      else: # node is a branch
        for child in node.getChildren():
          for x in _recurseHashes( child ): yield x

    for x in _recurseHashes( self._tree): yield x

  ############################################################################
  def getCurrentHashes( self ):
    def _recurseHashes( node ):
      if isinstance(node, Page):
        digest = node.getCurrentDigest()
        if digest:
          yield node.urlPath(), digest
      else: # node is a branch
        for child in node.getChildren():
          for x in _recurseHashes( child ): yield x

    for x in _recurseHashes( self._tree): yield x

  ############################################################################
  def getTree( self ):
    return self._tree

  ############################################################################
  def _timestamp( self, filename ):
    timestamp = os.path.getmtime( filename )
    return datetime.datetime.utcfromtimestamp( timestamp )

  ############################################################################
  def examine( self ):
    '''
    Whizz through the directory finding out what's in it and commenting
    on the result.
    '''
    cronLogFilename = os.path.join( self._rootPath, 'cron.out' )
    if os.path.exists( cronLogFilename ):
      self.cronOut = 'cron.out'
      timestamp = os.path.getmtime( cronLogFilename )
      self.cronTimestamp = datetime.datetime.utcfromtimestamp( timestamp )
    else:
      self.cronOut = None
      self.cronTimestamp = None

    self._tree = Root()
    for path, directories, filenames in os.walk( self._rootPath ):
      relativePath = os.path.relpath( path, self._rootPath )
      absolutePath = os.path.join( self._rootPath, path )

      parent = self._tree
      if relativePath != '.':
        for node in relativePath.split( os.sep ):
          parent = parent.getChild( node )

      # In the special case of an index file we descend no further and register
      # the index file with out parent...
      #
      if relativePath != '.' and 'index.html' in filenames:
        directories[:] = []
        parent.addChild( Page( 'index.html', None, None ) )
        parent.setIndexPage( 'index.html' )
      else:
        self._examineFiles( relativePath, absolutePath, filenames, parent )

      # Check for builder pages and notify our parent...
      #
      if 'compile.html' in filenames:
        parent.setBuildPage( parent.getChild( 'compile.html' ) )

      if 'run.html' in filenames:
        parent.setRunPage( parent.getChild( 'run.html' ) )

      for directory in directories:
        relativeFilename = os.path.join( relativePath, directory )
        absoluteFilename = os.path.join( absolutePath, directory )

        timestamp = self._timestamp( absoluteFilename )
        newDirectory = Directory( directory, timestamp )

        # Where a build stage creates a directory, the build output shares
        # its name with that directory.
        #
        if directory + '.html' in filenames:
          newDirectory.setBuildPage( parent.getChild( directory + '.html' ) )

        parent.addChild( newDirectory )

  ############################################################################
  def _examineFiles( self, relativePath, absolutePath, filenames, parent ):
    for filename in filenames:
      relativeFilename = os.path.join( relativePath, filename )
      absoluteFilename = os.path.join( absolutePath, filename )

      # First check to see if we are supposed to be ignoring this file
      #
      ignore = False
      for prefix in self._ignorePrefixes:
        if relativeFilename.startswith( prefix ):
          ignore = True
          break
      if ignore: break

      timestamp = self._timestamp( absoluteFilename )

      digest = None

      if filename.endswith( '.html' ):
        hasher = hashlib.sha1()
        try:
          document = et.parse( absoluteFilename )
          compilerNode = document.getroot().find( './/span[@id=\'compiler\']' )
          if compilerNode is not None:
            compiler = compilerNode.text
            hasher.update( compiler.encode("utf-8") )
          contextNode  = document.getroot().find( './/span[@id=\'context\']' )
          if contextNode is not None:
            context = contextNode.text
            hasher.update( context.encode("utf-8") )
          timestampNode = document.getroot().find( './/span[@id=\'timestamp\']' )
          if timestampNode is not None:
            timestamp = datetime.datetime.strptime( timestampNode.text,
                                                    '%Y-%m-%dT%H:%M:%SZ' )
          eventNodes = document.getroot().find( './/div[@id=\'events\']' )
          if eventNodes is not None:
            for eventNode in eventNodes:
              for node in eventNode:
                if node.text and node.text.strip() != '':
                  hasher.update( node.text.strip().encode("utf-8") )
        except et.ParseError as ex:
          message = 'Unable to parse file "{filename}": {exception}'
          usefulname = os.path.join( relativePath, filename )
          raise IndexerException( message.format( exception=ex,
                                                  filename=usefulname) )
        digest = hasher.hexdigest()

      newPage = Page( filename, timestamp, digest )
      parent.addChild( newPage )
