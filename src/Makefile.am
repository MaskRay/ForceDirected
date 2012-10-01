bin_PROGRAMS = force
force_SOURCES = Core.hh Core.cc FruchtermanReingold.hh Circle.hh KamadaKawai.hh Main.cc Cmdline.h Cmdline.c
force_CXXFLAGS = -std=c++11

EXTRA_DIST = Cmdline.ggo

Cmdline.c Cmdline.h: Cmdline.ggo
	gengetopt --file-name Cmdline < $<
