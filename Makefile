.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

try:
	from urllib import pathname2url
except:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT


define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python3 -c "$$BROWSER_PYSCRIPT"
PY?=python3
VENVDIR?=$(WORKDIR)/.venv

ifdef WORKDIR  # Must be an absolute path
WORKDIR:=$(abspath $(WORKDIR))
else
WORKDIR=$(CURDIR)
endif

VENV=$(VENVDIR)/bin
bold := $(shell tput bold)
sgr0 := $(shell tput sgr0)

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

activate: ## prints the command to run to activate the enironment
	@printf "Run the following: $(bold)source $(VENV)/activate$(sgr0)\n"

clean: clean-build clean-pyc clean-test clean-venv ## remove all build, test, coverage, venv and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .pytest_cache
	rm -f .coverage

.PHONY: clean-venv
clean-venv: ## remove virtual environment artifacts
	[ ! -d $(VENVDIR) ] || rm -rf $(VENVDIR)

init: ## initialise the virtual environment
	mkdir -p $(VENVDIR)
	$(PY) -m venv $(VENVDIR)
	$(VENV)/$(PY) -m pip install --upgrade pip

.PHONY: dev-install
dev-install: ## install dev dependencies
	$(VENV)/$(PY) -m pip install black pytest coverage

lint: dev-install ## reformat with black
	$(VENV)/$(PY) -m black snpit tests

test: dev-install ## run tests
	$(VENV)/$(PY) -m pytest

coverage: dev-install ## check code coverage and open in browser
	$(VENV)/$(PY) -m pip install coverage pytest
	$(VENV)/$(PY) -m coverage run --source snpit -m pytest
	$(VENV)/$(PY) -m coverage report -m
	$(VENV)/$(PY) -m coverage html
	@$(BROWSER) htmlcov/index.html

release: clean lint ## package and upload a release
	$(VENV)/$(PY) setup.py sdist upload
	$(VENV)/$(PY) setup.py bdist_wheel upload

dist: clean lint ## builds source and wheel package
	$(VENV)/$(PY) setup.py sdist
	$(VENV)/$(PY) -m twine upload -r pypi dist/`ls -t dist | head -1`

install: clean init ## install the package into the virtual environment
	$(VENV)/$(PY) -m pip install $(WORKDIR)

.PHONY: show-venv
show-venv: ## show the virtual environment details
	@$(VENV)/$(PY) -c "import sys; print('Python ' + sys.version.replace('\n',''))"
	@$(VENV)/pip --version
	@echo venv: $(VENVDIR)