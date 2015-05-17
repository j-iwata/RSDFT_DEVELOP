#!/usr/bin/env python

import fnmatch
import os
import sys
import numpy as np
import csv
import re
import math
import itertools
from numpy import array

def rmItemsFromList( items, lists ):
  for item in items:
    try:
      lists.remove( item )
    except ValueError:
      pass
def add( name, catagory_list ):
  if name not in catagory_list:
    catagory_list.append( name )
def writeList( f, list ):
  for item in list:
    f.write( ' {0:}\n'.format(item) )
def findKeyword( keyword, data, startline ):
  for idx in xrange(startline,len(data)):
    items = data[idx].split()
    try:
      if items[0].lower() == keyword.lower():
        itemname = items[1].split(",")[0].split("(")[0]
        return [itemname,idx+1]
    except IndexError:
      pass
  return ['',len(data)]
def findUse( source, data ):
  findline = 0
  while findline < len( data ):
    usename = ''
    [usename, findline] = findKeyword( 'use', data, findline )
    if usename is not '':
      source.addUse( usename )
def findSubroutine( source, data ):
  findline = 0
  while findline < len(data):
    subroutinename = ''
    [subroutinename, findline] = findKeyword( 'subroutine', data, findline )
    if subroutinename is not '':
      source.addSubroutine(subroutinename)
def findCall( source, data ):
  findline = 0
  while findline < len(data):
    callname=''
    [callname, findline] = findKeyword( 'call', data, findline )
    if callname is not '':
      source.addCall(callname)

class ModuleFile:
  def __init__( self, file_name, module_name ):
    self.file_name       = file_name
    self.module_name     = module_name
    self.use_list        = []
    self.subroutine_list = []
    self.call_list       = []
    self.contains_mpi    = False
    self.hierarchy_level = 0
  def addUse( self, use_name ):
    add( use_name, self.use_list )
  def addSubroutine( self, subroutine_name ):
    add( subroutine_name, self.subroutine_list )
  def addCall( self, call_name ):
    add( call_name, self.call_list )
  def getMPI( self ):
    tmp = []
    for item in self.call_list:
      if item.lower().startswith( 'mpi_' ):
        tmp.append( item )
        self.contains_mpi = True
    rmItemsFromList( tmp, self.call_list )
  def rmSelfSubroutines( self ):
    rmItemsFromList( self.subroutine_list, self.call_list )
  def isElementaryModule( self ):
    if not self.use_list:
      self.hierarchy_level = 1
      return True
    else:
      self.hierarchy_level = 999
      return False
  def getHierarchyLevel( self ):
    return self.hierarchy_level
  def raiseHierarchyLevel( self, module_name ):
    if module_name.hierarchy_level >= self.hierarchy_level or self.hierarchy_level > 100:
      self.hierarchy_level = module_name.hierarchy_level + 1
  def objectFilename( self ):
    return self.file_name.split( '.' )[0] + '.o'
  def writeDependency( self ):
    with open( make_common_dep, 'a' ) as f:
      basename = self.file_name.split( '.' )[0]
      f.write( self.module_name.lower() + '.mod: ' + basename + '.o' + '\n' )
      f.write( '	@true\n' )
      if self.hierarchy_level == 1:
        f.write( basename + '.o: ' + basename + '.f90' + '\n' )
      else:
        f.write( basename + '.o: ' + basename + '.f90 \\' + '\n' )
        for use_item in self.use_list[:-1]:
          f.write( '                 ' + use_item.lower() + '.mod \\\n' )
        f.write( '                 ' + self.use_list[-1].lower() + '.mod\n' )
      f.write( '	$(FC) $(FFLAGS) -c ' + self.file_name + ' -o $@\n' )
      f.write( '\n' )
  def show( self, filename ):
    with open(filename,'a') as f:
      f.write( '<file_name>       :  {0:}\n'.format(self.file_name) )
      f.write( '<module_name>     :  {0:}\n'.format(self.module_name) )
      f.write( '<hierarchy_level> :  {0:}\n'.format(self.hierarchy_level) )
      f.write( '<this module contains MPI calls> : {0:}\n'.format(self.contains_mpi) )
      f.write( '<this module use the following modules> :\n' )
      if not self.use_list:
        f.write( ' NONE\n' )
      else:
        writeList( f, self.use_list )
      f.write( '<this module contains the following subroutines> :\n' )
      if not self.subroutine_list:
        f.write( ' NONE\n' )
      else:
        writeList( f, self.subroutine_list )
      f.write( '<this module contains the following calls from other modules> :\n' )
      if not self.call_list:
        f.write( ' NONE\n' )
      else:
        writeList( f, self.call_list )
      f.write( '\n' )


DIRNAME     = os.getcwd()
MD_DIRNAME  = DIRNAME+'/mdsource'
FILES       = os.listdir(DIRNAME)
MDFILES     = os.listdir(MD_DIRNAME)
make_common_dep ='makefile.common.dep'
make_common     ='makefile.common'
with open(make_common,'w') as f:
  pass
with open(make_common_dep,'w') as f:
  pass

#----- get src files
filelist = []
for filename in FILES:
  if 'f90' not in filename:
    continue
  filelist.append(filename)
num_files = len(filelist)
print '<total number of source file> : ' + str(num_files)

#----- get relationship
count = 0
# module_file is the raw data
module_file  = []
# leveled_file will be appended according to the hierarchy level of the module
leveled_file = []
# dic_module is used for linking the module name and it's index
dic_module   = {}
program_file = []

print "-----started raw reading-----"
#----- clean previous files
with open('dep.log','w') as f:
  pass
with open('src.log','w') as f:
  pass
for filename in filelist:
  with open(filename,'r') as f:
    data = f.readlines()
  modulename = ''
  lastline   = 0
  [modulename,lastline] = findKeyword( 'module', data, 0 )
  if modulename == '':
#    [modulename, lastline] = findKeyword( 'program', data, 0 )
    if modulename == '':
      continue
  module_file.append( ModuleFile(filename, modulename) )
  dic_module[modulename] = len(module_file) - 1

  findUse( module_file[-1], data )
  findSubroutine( module_file[-1], data )
  findCall( module_file[-1], data )

  module_file[-1].getMPI()
  module_file[-1].rmSelfSubroutines()
  dummy=module_file[-1].isElementaryModule()
  module_file[-1].show('src.log')

for filename in filelist:
  with open( filename, 'r' ) as f:
    data = f.readlines()
  programname = ''
  lastline    = 0
  [programname, lastline] = findKeyword( 'program', data, 0 )
  program_file.append( ModuleFile(filename, programname ) )

  findUse( program_file[-1], data )
  findCall( program_file[-1], data )
  program_file[-1].getMPI()
  program_file[-1].show('src.log')

print "-----finished raw reading-----"
print "-----generating hierarchy-----"

i = 0
while i < 100:
  for classModule in module_file:
    if classModule.isElementaryModule():
      continue
    for usedModule in classModule.use_list:
      if module_file[dic_module[usedModule]].getHierarchyLevel() < 100:
        classModule.raiseHierarchyLevel(module_file[dic_module[usedModule]])
  i += 1

i = 0
while i < 100:
  for classModule in module_file:
    if classModule.getHierarchyLevel() is i:
# append to leveled_file according to hierarchy level from small
      leveled_file.append(classModule)
  i += 1

for classModule in leveled_file:
  classModule.show( 'dep.log' )

print "raw dependency       : 'src.log'"
print "leveled dependency   : 'dep.log'"

print "-----generating makefile.common-----"
print "dependency for make  : " + make_common_dep
print "common for make      : " + make_common
for classModule in leveled_file:
  classModule.writeDependency()
with open( make_common, 'a' ) as f:
  f.write( 'MODS1 = \\\n' )
  for classModule in leveled_file[:]:
    f.write( '       ' + classModule.objectFilename() + '\\\n' )
  f.write( '\n' )
sys.exit()
for program in program_file:
  program.writeDependency()
sys.exit()
for classModule in module_file:
  print classModule.module_name+' '+str(classModule.getHierarchyLevel())
