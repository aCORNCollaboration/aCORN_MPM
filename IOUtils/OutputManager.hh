#ifndef OUTPUTMANAGER_HH
#define OUTPUTMANAGER_HH

#include "QFile.hh"
#include <TObject.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

/// collection of saved TObjects
class TObjCollector {
public:
    /// destructor
    virtual ~TObjCollector() { clearItems(); }
    virtual void writeItems() {
        printf("Saving registered objects...");
        fflush(stdout);
        for(vector<TObject*>::iterator i = rootItems.begin(); i != rootItems.end(); i++)
            (*i)->Write();
        printf(" Done.\n");
    }
    void clearItems() {
        for(vector<TObject*>::iterator i = rootItems.begin(); i != rootItems.end(); i++)
            delete(*i);
        rootItems.clear();
        for(vector<TObject*>::iterator i = deleteItems.begin(); i != deleteItems.end(); i++)
            delete(*i);
        deleteItems.clear();
    }
    /// register a root object for output (and eventual deletion)
    virtual TObject* addObject(TObject* o, bool noWrite=false) { if(noWrite) deleteItems.push_back(o); else rootItems.push_back(o); return o; }
    vector<TObject*> rootItems;    ///< objects held until deleted
    vector<TObject*> deleteItems;  ///< other objects never written to file
};

/// indicator for level of problem with analysis
enum WarningLevel {
    BENIGN_WARNING,     ///< just FYI for troubleshooting
    MODERATE_WARNING,   ///< might be a problem to look at
    SEVERE_WARNING,     ///< almost certainly something is very wrong
    FATAL_WARNING       ///< data is corrupted and cannot be analyzed
};

/// manages output directory for grouping related information; manages a canvas, output QFile, output ROOT file, recursive subdirectories
class OutputManager: public TObjCollector {
public:
    /// constructor for top-level
    OutputManager(string nm, string bp);
    /// constructor for nested
    OutputManager(string nm, OutputManager* pnt);
    
    /// destructor
    virtual ~OutputManager() {
        if(writeRootOnDestruct) writeROOT();
        clearItems();
        if(rootOut) rootOut->Close();
        if(defaultCanvas && !parent) delete(defaultCanvas); 
    }
    
    /// generate a TH1F registered with this runs output objects list
    TH1F* registeredTH1F(string hname, string htitle, unsigned int nbins, float x0, float x1);
    /// generate a TH2F registered with this runs output objects list
    TH2F* registeredTH2F(string hname, string htitle, unsigned int nbinsx, float x0, float x1, unsigned int nbinsy, float y0, float y1);
    /// print current canvas
    virtual void printCanvas(string fname, string suffix=".pdf") const;
    
    /// put a data quality warning in the parent output file
    void warn(WarningLevel l, string descrip, Stringmap M = Stringmap());
    
    /// write output QFile
    virtual void write(string outName = "");
    /// open output ROOT file for writing (useful if output tree is being created and needs a home)
    void openOutfile();
    /// set whether to write ROOT output when destructed
    void setWriteRoot(bool w) { writeRootOnDestruct = w; }
    
    QFile qOut;                 ///< QFile output
    TFile* rootOut;             ///< ROOT file output
    TCanvas* defaultCanvas;     ///< canvas for drawing plots
    OutputManager* parent;      ///< parent output manager
    string basePath;            ///< general output path
    string plotPath;            ///< specific output path for plots
    string dataPath;            ///< specific output path for output data
    string rootPath;            ///< specific output path for ROOT files
    string name;                ///< name for this subsystem
    bool writeRootOnDestruct;   ///< whether to write ROOT file when destructed
    
    static bool squelchAllPrinting;     ///< whether to cancel all printCanvas output
    
protected:
    
    /// set name of this detector
    virtual void setName(string nm);
    /// write output ROOT file; WARNING: THIS DELETES ALL REGISTERED ITEMS; do last if you reference these.
    void writeROOT();
};

#endif
