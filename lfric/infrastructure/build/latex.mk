##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Build all latex documents found in the source directory.
#
# Expects the following environment variables to be set:
#   SOURCE_DIR  : Directory containing documents, these may be organised in
#                 sub-directories.
#   WORKING_DIR : Temporary space in which to park build products.
#
.SECONDEXPANSION:

export TEXINPUTS := $(WORKING_DIR)/figures:$(TEX_STUFF):$(TEXINPUTS)

.PHONY: documents
documents: $(patsubst $(SOURCE_DIR)/%.latex,$(WORKING_DIR)/%.pdf,$(DOCUMENTS)) \
           | $(DOCUMENT_DIR)
	$(Q)echo >/dev/null

$(WORKING_DIR)/%.pdf: $(WORKING_DIR)/%.aux $(WORKING_DIR)/%.bbl
	$(call MESSAGE,Cross-referencing,$@)
	$(Q)TEXINPUTS=$(WORKING_DIR)/figures/$(dir $*):$(TEXINPUTS); \
	while ( grep "Rerun to get" $(WORKING_DIR)/$*.log ); \
	do \
	    echo Rerunning build for cross-references; \
	    pdflatex -interaction errorstopmode -output-directory $(dir $@) $(SOURCE_DIR)/$*.latex ; \
	done
	$(Q)mkdir -p $(DOCUMENT_DIR)
	$(Q)cp $@ $(DOCUMENT_DIR)

.PRECIOUS: $(WORKING_DIR)/%.bbl
$(WORKING_DIR)/%.bbl: $(WORKING_DIR)/%.aux \
                      $(WORKING_DIR)/%.bib $(WORKING_DIR)/%.bst
	$(call MESSAGE,Bibliography,$@)
	$(Q)if [ -e $(WORKING_DIR)/$*.bib ]; then cd $(dir $<); bibtex $(notdir $<); else echo "None"; fi

$(WORKING_DIR)/%.bib: $$(shell find $(SOURCE_DIR)/$$(dir $$*) -name *.bib)
	$(call MESSAGE,Bibliography files,$^)
	$(Q)for file in $^; do echo "Copying $$file"; cp $$file $(dir $@); touch $@; done

$(WORKING_DIR)/%.bst: $$(shell find $(SOURCE_DIR)/$$(dir $$*) -name *.bst)
	$(call MESSAGE,Bibliographic styles,$^)
	$(Q)for file in $^; do echo "Copying $$file"; cp $$file $(dir $@); touch $@; done

.PRECIOUS: $(WORKING_DIR)/%.aux)
$(WORKING_DIR)/%.aux: $(SOURCE_DIR)/%.latex \
                      $(WORKING_DIR)/common/figures \
                      $(WORKING_DIR)/$$(dir $$*)figures
	$(call MESSAGE,Laying out,$@)
	$(Q)mkdir -p $(dir $@)
	$(Q)TEXINPUTS=$(WORKING_DIR)/figures/$(dir $*):$(TEXINPUTS); pdflatex -interaction errorstopmode -output-directory $(dir $@) $<

.PHONY: $(WORKING_DIR)/common/figures
$(WORKING_DIR)/common/figures: $(patsubst $(COMMON_FIGURES)/%,$(WORKING_DIR)/figures/%.pdf,$(basename $(wildcard $(COMMON_FIGURES)/*)))
	$(Q)echo >/dev/null

.PRECIOUS: $(WORKING_DIR)/figures/%.pdf | $(WORKING_DIR)/figures
$(WORKING_DIR)/figures/%.pdf: $(COMMON_FIGURES)/%.svg
	$(call MESSAGE,Transcoding,$<)
	$(Q)mkdir -p $(dir $@)
	$Q$(INKSCAPE) $< --export-pdf=$@

$(WORKING_DIR)/figures/%.pdf: $(COMMON_FIGURES)/%.eps | $(WORKING_DIR)/figures
	$(call MESSAGE,Transcoding,$<)
	$(Q)mkdir -p $(dir $@)
	$(Q)eps2pdf --outfile=$@ $<

.PHONY: $(WORKING_DIR)/%/figures
$(WORKING_DIR)/%/figures: $$(patsubst $$(SOURCE_DIR)/$$*/figures/$$(PERCENT),$(WORKING_DIR)/figures/$$*/$$(PERCENT).pdf,$$(basename $$(wildcard $$(SOURCE_DIR)/$$*/figures/*)))
	$(Q)echo >/dev/null

.PRECIOUS: $(WORKING_DIR)/figures/%.pdf
$(WORKING_DIR)/figures/%.pdf: $(SOURCE_DIR)/$$(dir $$*)/figures/$$(notdir $$*).svg | $(WORKING_DIR)/figures/$$(dir $$*)
	$(call MESSAGE,Transcoding,$<)
	$(Q)mkdir -p $(dir $@)
	$Q$(INKSCAPE) $< --export-pdf=$@

$(WORKING_DIR)/figures/%.pdf: $(SOURCE_DIR)/$$(dir $$*)/figures/$$(notdir $$*).eps | $(WORKING_DIR)/figures/$$(dir $$*)
	$(call MESSAGE,Transcoding,$<)
	$(Q)mkdir -p $(dir $@)
	$(Q)epstopdf --outfile=$@ $<

$(DOCUMENT_DIR) $(WORKING_DIR) $(WORKING_DIR)/figures \
$(patsubst $(SOURCE_DIR)/%,$(WORKING_DIR)/figures/%,$(dir $(DOCUMENTS))): ALWAYS
	$(call Message,Creating,$@)
	$(Q)mkdir -p $@

include $(LFRIC_BUILD)/lfric.mk
