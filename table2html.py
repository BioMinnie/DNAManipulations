#########################################################################

"""
Script built upon the code provide by Rosetta Code at:

http://rosettacode.org/wiki/CSV_to_HTML_translation

Content is available under GNU Free Documentation License 1.2.

BeautifulHTML Package Index Owner: alice1017

## MODIFIED BY: Melinda Ashcroft 10/09/2014
## MELINDA AFFILIATIONS: Beatson lab | SCMB - UQ
"""

#########################################################################

# Import modules required
import sys
from cgi import escape
from Bhtml import BeautifulHTML

# Open files as variables using arguments passed to the sys module via the command line
with open(sys.argv[1], 'r') as fin:       # Open the infile (typically tab delimited file)
     with open(sys.argv[2], 'w') as fout: # Open the outfile (html)

        # Read in the file line by line and store in a variable	
        rawdata = fin.read()
        # Create the html headers, css elements and body tags
      	html = BeautifulHTML()
      	html1 = BeautifulHTML()
      	html.write("<html>")
      	html.indent("<head>")
      	css = html.setcss()
      	css.element("body",{"background":"#fff","color":"#000"})
      	css.element("table",{"width":"100%","background":"#fff","color":"#000", "border":"1px solid #666"})
      	css.element("tbody:last-child",{"display":"none"})
      	css.element("td",{"border-right":"1px solid #666"})
      	css.element("td:last-child",{"border":"0"})
      	html.writecss(css.release(), indent=True, block=2)
      	html.indent("</head>")
      	html.indent("<body>")
      	html1.indent("</body>")
	
        # Function 1 to read in each row of input file and split on the tab delimiter generating rows and columns
       	def _row2trextra(row, attr=None):
      	    cols = escape(row).split('\t')
      	    attr_tr = attr.get('TR', '')
      	    attr_td = attr.get('TD', '')
      	    return (('<TR%s>' % attr_tr)
      		    + ''.join('<TD%s>%s</TD>' % (attr_td, data) for data in cols)
      		    + '</TR>')
 
        # Function 2 to build the html table based upon the input tab delimited file
      	def csv2htmlextra(txt, header=True, attr=None):
      	    ' attr is a dictionary mapping tags to attributes to add to that tag'
      	 
      	    attr_table = attr.get('TABLE', '')
      	    attr_thead = attr.get('THEAD', '')
      	    attr_tbody = attr.get('TBODY', '')
      	    htmltxt = '<TABLE%s>\n' % attr_table
	          for rownum, row in enumerate(txt.split('\n')):
          		htmlrow = _row2trextra(row, attr)
          		rowclass = ('THEAD%s' % attr_thead) if (header and rownum == 0) else ('TBODY%s' % attr_tbody)
          		htmlrow = '<%s>%s</%s>\n' % (rowclass, htmlrow, rowclass[:0])
          		htmltxt += htmlrow
          	    htmltxt += '</TABLE>\n'
          	    return htmltxt

        # Generate the html variable by building a html table from the dictionary of Function 2, using basic html tags	
      	htmltxt = csv2htmlextra(rawdata, True,
      				dict(TABLE=' cellspacing="0", cellpadding="5", summary="Methylation data table"',
      		                     THEAD=' style="text-align: center; font-weight: bold; border: 1px solid #000; background-color: lightgrey;"',
      		                     TBODY='' 
      		                     )
      		                )
      		                
        # Write the html variables to the output html file
	
      	fout.write(str(html))	
      	fout.write(htmltxt)
      	fout.write(str(html1))
      	fout.write("</html>")
