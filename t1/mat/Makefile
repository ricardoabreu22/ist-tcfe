all: octave.log
	$(eval OCTVTABSM=$(shell grep _MTAB $< | sed 's/_MTAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval OCTVTABSN=$(shell grep _NTAB $< | sed 's/_NTAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval OCTVTABSMM=$(shell grep _MMTAB $< | sed 's/_MMTAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval OCTVTABSNN=$(shell grep _MNTAB $< | sed 's/_MNTAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(OCTVTABSM), sed -n '/^$i_MTAB/,/^$i_MEND/{p;/^$i_MEND/q}' $< | grep -v $i_MTAB | grep -v $i_MEND | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(OCTVTABSN), sed -n '/^$i_NTAB/,/^$i_NEND/{p;/^$i_NEND/q}' $< | grep -v $i_NTAB | grep -v $i_NEND | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(OCTVTABSMM), sed -n '/^$i_MMTAB/,/^$i_MMEND/{p;/^$i_MMEND/q}' $< | grep -v $i_MMTAB | grep -v $i_MMEND | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(OCTVTABSNN), sed -n '/^$i_MNTAB/,/^$i_MNEND/{p;/^$i_MNEND/q}' $< | grep -v $i_MNTAB | grep -v $i_MNEND | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)

octave.log: circuitanalysis.m
	octave $< > $@

clean:
	@rm -f octave.log octave-workspace *.eps *~ *.tex

.PHONY: all clean
