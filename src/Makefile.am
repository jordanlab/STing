bin_PROGRAMS = $(top_builddir)/bin/indexer \
	$(top_builddir)/bin/typer \
	$(top_builddir)/bin/detector 

__top_builddir__bin_indexer_SOURCES = Indexer.h \
	Indexer.cpp \
	IndexerOptions.h \
	IndexerOptions.cpp \
	types.hpp \
	custom_utils.hpp \
	dir_utils.hpp \
	data_structures.hpp

__top_builddir__bin_typer_SOURCES = Typer.h \
	Typer.cpp \
	SequenceInfo.h \
	SequenceInfo.cpp \
	KmerMatcher.h \
	KmerMatcher.cpp \
	KmerMatcherOptions.h \
	KmerMatcherOptions.cpp \
	Timer.h \
	Timer.cpp \
	types.hpp \
	custom_utils.hpp \
	data_structures.hpp

__top_builddir__bin_detector_SOURCES = Detector.h \
	Detector.cpp \
	KmerMatcher.h \
	KmerMatcher.cpp \
	KmerMatcherOptions.h \
	KmerMatcherOptions.cpp \
	Timer.h \
	Timer.cpp \
	types.hpp \
	custom_utils.hpp \
	data_structures.hpp

# __top_builddir__bin_classifier_SOURCES = Classifier.h \
# 	Classifier.cpp \
# 	SequenceInfo.h \
# 	SequenceInfo.cpp \
# 	KmerMatcher.h \
# 	KmerMatcher.cpp \
# 	KmerMatcherOptions.h \
# 	KmerMatcherOptions.cpp \
# 	Timer.h \
# 	Timer.cpp \
# 	types.hpp \
# 	custom_utils.hpp \
# 	data_structures.hpp

# AM_CXXFLAGS = -I ../../include -W -Wall -O3 -pedantic -std=c++11 -pthread
# AM_CXXFLAGS = -I $(srcdir)/../include -W -Wall -O3 -pedantic -std=c++11 -pthread
AM_CXXFLAGS = -I $(srcdir)/../include -W -Wall -pedantic -std=c++11 -pthread
AM_LDFLAGS = -static
