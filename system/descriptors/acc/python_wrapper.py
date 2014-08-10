# -*- coding: utf8 -*-

import ctypes

def main(args):
	lib = ctypes.cdll.LoadLibrary("lib/libcolordescriptors.so")
	ACC = lib.ACC
	print ACC.restype

if __name__ == '__main__':
	import sys
	main(sys.argv)