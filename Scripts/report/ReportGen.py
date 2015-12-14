#!/usr/bin/python

import os

class ReportGenerator:
    
    def __init__(self):
        self.missing = []
        self.slist = [3102, 3109, 3116, 3124, 3126, 3127, 3128, 3129, 3135, 3137, 3138,
                      3143, 3147, 3155, 3186, 3191, 3196, 3197, 3228, 3254, 3290, 3356] 
        #self.slist = self. get_series_list()
        self.outPath = os.environ["ACORN_WISHBONE"]+"/Report/"
        os.system("mkdir -p "+self.outPath)
        
    def get_series_list(self):
        md = open("/home/mpmendenhall/Documents/aCORN/Reference/MetaData_06102014.csv").readlines()[10:]
        md = [l.split(",") for l in md]
        nprev = ""
        slist = []
        self.missing = []
        for m in md:
            if len(m)==10 and int(m[3])==1 and m[2]!=nprev:
                s = int(m[2])
                nprev = m[2]
                sdir = os.environ["ACORN_WISHBONE"]+"/Series_%i/"%s
                if not os.path.exists(sdir):
                    print "Missing", sdir
                    self.missing.append(sdir)
                    continue
                slist.append(s)
        return slist
            
    def makePages(self,tname):
        pgs = ""
        btext = open(tname,"r").read()
        for s in self.slist:
            page_settings = {"series":s, "sdir":os.environ["ACORN_WISHBONE"]+"/Series_%i/"%s}
            if not os.path.exists(page_settings["sdir"]):
                print "Missing",page_settings["sdir"]
                os.system("../../WishboneScanner %i"%s)
                continue
            pgs += btext%page_settings
        page_settings = {"series":9999, "sdir":os.environ["ACORN_WISHBONE"]}
        if os.path.exists(page_settings["sdir"]+"/Series_9999"):
            pgs += btext%page_settings
        return pgs
        
    def build_output(self):
        oneseries = len(self.slist) == 1
        
        self.body = ""
        if not oneseries:
            self.body += "\\tableofcontents\n"
        if self.missing:
            self.body += r"\section{Report errors}"+'\n'
            self.body += "Missing series files:\n"+r"\begin{itemize}"+'\n'
            for m in self.missing:
                self.body += '\t'+r"\item \begin{verbatim}"+ m + r"\end{verbatim}"+'\n';
            self.body += r"\end{itemize}"+'\n\n'
        
        if not oneseries:
            self.body += r"\clearpage\section{Detector response}"
        self.body += self.makePages("Template_Page_1.tex")
        if not oneseries:
            self.body += r"\clearpage\section{Event positions}"
        self.body += self.makePages("Template_Page_3.tex")
        if not oneseries:
            self.body += r"\clearpage\section{Wishbone}"
        self.body += self.makePages("Template_Page_2.tex")
        
        texFile = self.outPath+"/aCORN_Summary.tex"
        fOut = open(texFile,"w");
        doc_settings = {"title":r"\aCORN\ analysis report", "author":r"\texttt{aCORN\_MPM/Scripts/report/ReportGen.py}", "body":self.body}
        if oneseries:
            doc_settings["title"] = r"\aCORN\ series %i report"%self.slist[0]
        fOut.write(open("Template_Document.tex","r").read()%doc_settings)
        fOut.close()
        
        os.system("cd %s; pdflatex aCORN_Summary.tex"%self.outPath)
    
if __name__=="__main__":
    
    R = ReportGenerator()
    R.build_output()
    
    
    
    
