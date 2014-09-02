#!/usr/bin/python

import os

def get_series_list():
    md = open("/home/mpmendenhall/Documents/aCORN/MetaData_06102014.csv").readlines()[10:]
    md = [l.split(",") for l in md]
    nprev = ""
    slist = []
    for m in md:
        if len(m)==10 and int(m[3])==1 and m[2]!=nprev:
            slist.append(int(m[2]))
            nprev = m[2]
    return slist
            
def makePages(tname):
    pgs = ""
    btext = open(tname,"r").read()
    for s in get_series_list():
        page_settings = {"series":s, "sdir":os.environ["ACORN_WISHBONE"]+"/Series_%i/"%s}
        if not os.path.exists(page_settings["sdir"]):
            print "Missing",page_settings["sdir"]
            continue
        pgs += btext%page_settings
    page_settings = {"series":9999, "sdir":os.environ["ACORN_WISHBONE"]}
    pgs += btext%page_settings
    return pgs
    
if __name__=="__main__":
    
    outPath = os.environ["ACORN_WISHBONE"]+"/Report/"
    os.system("mkdir -p "+outPath)
    
    body = r"\clearpage\section{Detector response}"
    body += makePages("Template_Page_1.tex")
    body += r"\clearpage\section{Wishbone}"
    body += makePages("Template_Page_2.tex")
    
    texFile = outPath+"/aCORN_Summary.tex"
    fOut = open(texFile,"w");
    doc_settings = {"title":r"\aCORN\ analysis report", "author":r"\texttt{aCORN\_MPM/Scripts/report/ReportGen.py}", "body":body}
    fOut.write(open("Template_Document.tex","r").read()%doc_settings)
    fOut.close()
    
    os.system("cd %s; pdflatex aCORN_Summary.tex"%outPath)
    