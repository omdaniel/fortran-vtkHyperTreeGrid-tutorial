# Build the two tutorial programs with gfortran (any Fortran 2008 compiler works).
FC      = gfortran
FFLAGS  ?= -std=f2018 -Wall -O2

all: single_tree forest export_cells

quadtree_mod.o: quadtree_mod.f90
	$(FC) $(FFLAGS) -c $<

field_mod.o: field_mod.f90 quadtree_mod.o
	$(FC) $(FFLAGS) -c $<

single_tree: single_tree.f90 quadtree_mod.o field_mod.o
	$(FC) $(FFLAGS) $^ -o $@

forest: forest.f90 quadtree_mod.o field_mod.o
	$(FC) $(FFLAGS) $^ -o $@

export_cells: export_cells.f90 quadtree_mod.o field_mod.o
	$(FC) $(FFLAGS) $^ -o $@

# Writes cells_single.txt and cells_forest.txt, then the SVGs in figures/
figures: export_cells
	./export_cells
	python3 make_figures.py

# Builds the styled GitHub Pages version of TUTORIAL.md into docs/ (needs pandoc)
site: TUTORIAL.md make_site.py site/template.html site/style.css
	python3 make_site.py

run: all
	./single_tree
	./forest

clean:
	rm -f *.o *.mod single_tree forest export_cells

.PHONY: all run clean figures site
