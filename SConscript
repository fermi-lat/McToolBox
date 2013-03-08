# -*- python -*-
# $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/McToolBox/SConscript,v 1.0 2013/02/21 20:00:04 usher Exp $ 
# Authors: Leon Rochester <lsrea@slac.stanford.edu>, Tracy Usher <usher@slac.stanford.edu>
# Version: McToolBox-00-02-00
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('addLinkDeps', package='McToolBox', toBuild='component')
McToolBox=libEnv.SharedLibrary('McToolBox',
                              listFiles(['src/*.cxx',
                                         'src/*.h']))
progEnv.Tool('McToolBoxLib')

test_McToolBox = progEnv.GaudiProgram('test_McToolBox',
                                     listFiles(['src/test/*.cxx']),
                                     test = 1, package='McToolBox')

progEnv.Tool('registerTargets', package='McToolBox',
             libraryCxts = [[McToolBox,libEnv]],
             testAppCxts = [[test_McToolBox, progEnv]],
             includes = listFiles(['McToolBox/*'], recursive = 1),
             xml=listFiles(files = ['xml/*'], recursive = True),
             jo = ['src/test/jobOptions.txt'])
