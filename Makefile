.PHONY: install-intel install-intel-local

install-intel:
	./install.sh

install-intel-local:
	PREFIX=$(PWD)/.local-intel ./install.sh
