.PHONY: report graph RT clean
REPORT:=report.pdf
REPORT_SRC:=report.typ
GRAPHDIR:=graph
Q2S_PDF = $(GRAPHDIR)/Q2S.svg
Q2U_PDF = $(GRAPHDIR)/Q2U.svg
Q4A_PDF = $(GRAPHDIR)/Q4A.svg
Q4B_PDF = $(GRAPHDIR)/Q4B.svg

all: $(REPORT)
report:$(REPORT)

$(REPORT): $(REPORT_SRC) $(Q2S_PDF) $(Q2U_PDF) $(Q4A_PDF) $(Q4B_PDF)
	typst c $(REPORT_SRC) $(REPORT)

$(Q2U_PDF): Q2.py
	python3 Q2.py

$(Q2S_PDF): Q2.py
	python3 Q2.py

$(Q4A_PDF): Q4A.py
	python3 Q4A.py

$(Q4B_PDF): Q4B.py
	python3 Q4B.py

clean:
	rm -f $(REPORT) $(GRAPHDIR)/*

RT:
	watch -n3 -d make report
