R=R
# -> you can do    R=R-devel  make ....

PACKAGE=compResidual
VERSION=$(shell sed -n '/^Version: /s///p' ${PACKAGE}/DESCRIPTION)
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip

CPP_SRC = $(PACKAGE)/src/*.cpp

SUBDIRS := $(wildcard test/*/.)

.PHONY: test testseq testone $(SUBDIRS) all check clean install 

ifeq (testone,$(firstword $(MAKECMDGOALS)))
  ARG := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(ARG):;@:)
endif

all:
	make doc-update
	make build-package
	make install
	make pdf

doc-update: $(PACKAGE)/R/*.R
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\")" | $(R) --slave
	@touch doc-update

namespace-update :: $(PACKAGE)/NAMESPACE
$(PACKAGE)/NAMESPACE: $(PACKAGE)/R/*.R
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\")" | $(R) --slave

build-package: $(TARBALL)
$(TARBALL): $(PACKAGE)/NAMESPACE
	$(R) CMD build --resave-data=no $(PACKAGE)

install: $(TARBALL)
	$(R) CMD INSTALL --preclean $(TARBALL)
	@touch install

debug-install: $(PACKAGE)/src/compResidual.so
	$(R) CMD INSTALL $(PACKAGE)

$(PACKAGE)/src/compResidual.so: $(PACKAGE)/src/compResidual.cpp $(CPP_SRC)
	touch $(PACKAGE)/src/compResidual.cpp
	cd $(PACKAGE)/src; echo "library(TMB); compile('compResidual.cpp','-O0 -g', libinit=FALSE)" | $(R) --slave


unexport TEXINPUTS
pdf: $(PACKAGE).pdf
$(PACKAGE).pdf: $(PACKAGE)/man/*.Rd
	rm -f $(PACKAGE).pdf
	$(R) CMD Rd2pdf --no-preview $(PACKAGE)

check:
	$(R) CMD build --resave-data=no $(PACKAGE)
	$(R) CMD check $(TARBALL)

clean:
	\rm -f install doc-update compResidual_* compResidual.pdf compResidual/src/*.o compResidual/src/*.so compResidual/src/*.dll


NPROCS:=1
OS:=$(shell uname -s)

ifeq ($(OS),Linux)
  NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)
endif

MAKEFLAGS += --silent

test:
	$(MAKE) -j $(NPROCS) testseq

testseq: $(SUBDIRS)

testone:
	@$(MAKE) test/$(ARG)/.

$(SUBDIRS):
	@cp test/Makefile $@
	@$(MAKE) -i -s -C $@
	@rm -f $@/Makefile
