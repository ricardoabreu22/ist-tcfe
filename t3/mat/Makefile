all: octave.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	$(eval OCTVTABSM=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	
	$(foreach i, $(OCTVTABSM), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_MEND/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	

octave.log: circuitanalysis.m
	octave $< > $@

clean:
	@rm -f octave.log octave-workspace *.eps *~ *.tex

.PHONY: all clean
