PKG_VERSION=$(shell grep -i ^version ./DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME=$(shell grep -i ^package ./DESCRIPTION | cut -d : -d \  -f 2)

R_FILES := $(wildcard ./R/*.R)
PKG_FILES := ./DESCRIPTION ./NAMESPACE $(R_FILES) README.md 

.PHONY: tarball install check clean roxygen


sourcetar:
	zip -r ../$(PKG_NAME)_$(PKG_VERSION).zip $(PKG_FILES) man test  
	#zip -d $(PKG_NAME)_$(PKG_VERSION).zip .Rproj.user/*

tarball: $(PKG_NAME)_$(PKG_VERSION).tar.gz
$(PKG_NAME)_$(PKG_VERSION).tar.gz: $(PKG_FILES)
	R CMD build .

all: check install

check: roxygen $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD check $(PKG_NAME)_$(PKG_VERSION).tar.gz

install: roxygen $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD INSTALL $(PKG_NAME)_$(PKG_VERSION).tar.gz
			
roxygen:
	Rscript -e "library(roxygen2);roxygenize('.')"

clean:
	-rm -f $(PKG_NAME)_*.tar.gz
	-rm -r -f $(PKG_NAME).Rcheck

