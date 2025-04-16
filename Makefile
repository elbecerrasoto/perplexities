DATA = data
TARGZ = data.tar.gz

.PHONY help:
help:
	less Makefile

.PHONY test:
test: $(DATA)
	Rscript asserts.R

$(DATA): $(TARGZ)
	@if [ -d "$(DATA)" ]; then echo "$(DATA) already exists, skipping extraction."; else tar -xvf $<; fi
	@touch $@

# PHONY Avoids circular dependency
.PHONY targz:
targz:
	tar -czvf $(TARGZ) $(DATA)

.PHONY style:
style:
	Rscript -e 'styler::style_dir(".", recursive = FALSE)'

# .PHONY clean:
# dangerous
# as I sometimes modify the input tables
#	rm -rf data
