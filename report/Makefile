filename=writeup

pdf:
	pdflatex ${filename}
	# makeindex ${filename}.nlo -s nomencl.ist -o ${filename}.nls
	bibtex ${filename}||true
	pdflatex ${filename}
	pdflatex ${filename}

clean:
	rm -f ${filename}.log ${filename}.aux ${filename}.out ${filename}.bbl
	rm -f ${filename}.blg ${filename}.fls ${filename}.ilg ${filename}.nlo
	rm -f ${filename}.nls ${filename}.fdb_latexmk
	rm -f .${filename}.tex.swp