all: octave.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	$(eval OCTVTABSN=$(shell grep _NTAB $< | sed 's/_NTAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval OCTVTABSNN=$(shell grep _MNTAB $< | sed 's/_MNTAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval OCTVTABSN21=$(shell grep _NTAB21 $< | sed 's/_NTAB21//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval OCTVTABSN22=$(shell grep _NTAB22 $< | sed 's/_NTAB22//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval OCTVTABSN41=$(shell grep _TAB41 $< | sed 's/_TAB41//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval OCTVTABSN42=$(shell grep _TAB42 $< | sed 's/_TAB42//g' | sed ':a;N;$!ba;s/\n/ /g'))
	
	
	$(foreach i, $(OCTVTABSN), sed -n '/^$i_NTAB/,/^$i_NEND/{p;/^$i_NEND/q}' $< | grep -v $i_NTAB | grep -v $i_NEND | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(OCTVTABSNN), sed -n '/^$i_MNTAB/,/^$i_MNEND/{p;/^$i_MNEND/q}' $< | grep -v $i_MNTAB | grep -v $i_MNEND | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(OCTVTABSN21), sed -n '/^$i_NTAB21/,/^$i_NEND21/{p;/^$i_NEND21/q}' $< | grep -v $i_NTAB21 | grep -v $i_NEND21 | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(OCTVTABSN22), sed -n '/^$i_NTAB22/,/^$i_NEND22/{p;/^$i_NEND22/q}' $< | grep -v $i_NTAB22 | grep -v $i_NEND22 | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(OCTVTABSN41), sed -n '/^$i_TAB41/,/^$i_END41/{p;/^$i_END41/q}' $< | grep -v $i_TAB41 | grep -v $i_END41 | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(OCTVTABSN42), sed -n '/^$i_TAB42/,/^$i_END42/{p;/^$i_END42/q}' $< | grep -v $i_TAB42 | grep -v $i_END42 | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	
octave.log: circuitanalysis.m
	octave $< > $@

clean:
	@rm -f octave.log octave-workspace *.eps *~ *.tex

.PHONY: all clean
