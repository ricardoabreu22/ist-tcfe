all: ngspice.log ng2 ng3 ng4
	#$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	#$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	#$(foreach i, $(NGSFIGPS), ps2pdf $i;)


ngspice.log: ../sim/t0.net
	ngspice -b $< -o $@
	
ng2: ngspice2.log
	#$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	#$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS2=$(shell grep _TAB2 $< | sed 's/_TAB2//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS2), sed -n '/^$i_TAB2/,/^$i_END2/{p;/^$i_END2/q}' $< | grep -v $i_TAB2 | grep -v $i_END2 | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab2.tex;)

ngspice2.log: ../sim/t0_2.net
	ngspice -b $< -o $@
	
ng3: ngspice3.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
ngspice3.log: ../sim/t0_3.net
	ngspice -b $< -o $@
	
ng4:ngspice4.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	

ngspice4.log: ../sim/t0_4.net
	ngspice -b $< -o $@

clean:
	@rm -f ngspice.log *.tex *.pdf *.ps *~ *.cir
	@rm -f ngspice2.log *.tex *.pdf *.ps *~
	@rm -f ngspice3.log *.tex *.pdf *.ps *~
	@rm -f ngspice4.log *.tex *.pdf *.ps *~
	
.PHONY: all clean
