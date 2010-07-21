#! /usr/bin/env python

VERSION='0.4'
APPNAME='cpdeflec'

top = '.'
out = 'build'

def set_options(opt):
	opt.tool_options('compiler_cc')

def configure(conf):
	conf.check_tool('compiler_cc')

def build(bld):
	bld(
			features = ['cc', 'cprogram'],
			source = ['cpdeflec.c','commath.c','camera.c','pattern.c',
				'fitting.c'],
			target = 'cpdeflec',
			includes = ['.'],
			defines = ['LINUX=1', 'BIDULE'],
			ccflags = ['-O2', '-std=c99'],
			lib = ['m', 'tiff', 'gsl', 'gslcblas'],
			libpath = ['/usr/lib']
	)
