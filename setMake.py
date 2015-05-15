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
    if use_name not in self.use_list:
      self.use_list.append( use_name )
  def addSubroutine( self, subroutine_name ):
    if subroutine_name not in self.subroutine_list:
      self.subroutine_list.append( subroutine_name )
  def addCall( self, call_name ):
    if call_name not in self.call_list:
      self.call_list.append( call_name )
  def getMPI( self ):
    tmp = []
    for item in self.call_list:
      if item.lower().startswith( 'mpi_' ):
        tmp.append( item )
        self.contains_mpi = True
    for item in tmp:
      self.call_list.remove( item )
  def rmSelfSubroutines( self ):
    for item in self.subroutine_list:
      if item in self.call_list:
         self.call_list.remove( item )
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
    with open( make_common, 'a' ) as f:
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
      f.write( '<file_name>   : ' + self.file_name + '\n' )
      f.write( '<module_name> : ' + self.module_name + '\n' )
      f.write( '<hierarchy_level> : ' + str(self.hierarchy_level) + '\n' )
      f.write( '<this module contains MPI calls> : ' + str(self.contains_mpi) + '\n' )
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

def writeList( f, list ):
  for item in list:
    f.write( ' ' + item + '\n' )
def findKeyword( keyword, data, startline ):
  for idx in xrange(startline,len(data)):
    itemlist = data[idx].split()
    if len(itemlist) == 0:
      continue
    if itemlist[0].lower() == keyword:
      itemname = itemlist[1].split(",")[0].split("(")[0]
      return [itemname,idx+1]
  return ['',len(data)]
def findUse( source_list, idx, data ):
  findline = 0
  while findline < len(data):
    usename = ''
    [usename,findline] = findKeyword( 'use', data, findline )
    if usename is not '':
      source_file[idx].addUse( usename )
def findSubroutine( source_list, idx,data ):
  findline = 0
  while findline < len(data):
    subroutinename = ''
    [subroutinename,findline] = findKeyword( 'subroutine', data, findline )
    if subroutinename is not '':
      source_file[idx].addSubroutine(subroutinename)
def findCall( source_list, idx, data ):
  findline = 0
  while findline < len(data):
    callname=''
    [callname,findline] = findKeyword( 'call', data, findline )
    if callname is not '':
      source_file[idx].addCall(callname)

DIRNAME=os.getcwd()
#MD_DIRNAME=DIRNAME+'/mdsource'
FILES=os.listdir(DIRNAME)
#MDFILES=os.listdir(MD_DIRNAME)
make_common='makefile.dep.common'
make_common_obj='makefile.common'
with open(make_common,'w') as f:
  pass
with open(make_common_obj,'w') as f:
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
# source_file is the raw data
source_file  = []
# leveled_file will be appended according to the hierarchy level of the module
leveled_file = []
# dic_source is used for linking the module name and it's index
dic_source   = {}

print "-----started raw reading-----"
#----- clean previous files
with open('dep.log','w') as f:
  pass
with open('src.log','w') as f:
  pass
for idx,filename in enumerate(filelist):
  with open(filename,'r') as f:
    data=f.readlines()

  modulename = ''
  lastline   = 0
  [modulename,lastline] = findKeyword( 'module', data, 0 )
  if modulename is '':
    [modulename,lastline] = findKeyword( 'program', data, 0 )
  source_file.append(ModuleFile(filename,modulename))
  dic_source[modulename] = idx

  findUse( source_file, idx, data )
  findSubroutine( source_file, idx, data )
  findCall( source_file, idx, data )

  source_file[idx].getMPI()
  source_file[idx].rmSelfSubroutines()
  dummy=source_file[idx].isElementaryModule()
  source_file[idx].show('src.log')

print "-----finished raw reading-----"
print "-----generating hierarchy-----"

i = 0
while i < 100:
  for classModule in source_file:
    if classModule.isElementaryModule():
      continue
    for usedModule in classModule.use_list:
      if source_file[dic_source[usedModule]].getHierarchyLevel() < 100:
        classModule.raiseHierarchyLevel(source_file[dic_source[usedModule]])
  i += 1

i = 0
while i < 100:
  for classModule in source_file:
    if classModule.getHierarchyLevel() is i:
# append to leveled_file according to hierarchy level from small
      leveled_file.append(classModule)
  i += 1

for classModule in leveled_file:
  classModule.show( 'dep.log' )

print "raw dependency       : 'src.log'"
print "leveled dependency   : 'dep.log'"

print "-----generating makefile.common-----"
print "dependency for make  : " + make_common
print "dependency for make  : " + make_common_obj
for classModule in leveled_file:
  classModule.writeDependency()
with open( make_common, 'w' ) as f:
  f.write( 'OBJ_ALL = \\\n' )
  for classModule in leveled_file[:]:
    f.write( '       ' + classModule.objectFilename() + '\\\n' )
  f.write( '\n' )
sys.exit()
for classModule in source_file:
  print classModule.module_name+' '+str(classModule.getHierarchyLevel())
