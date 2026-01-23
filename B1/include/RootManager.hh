#ifndef RootManager_h
#define RootManager_h

#include <fstream>

// Singleton class which manages output to text files

class MyROOTManager
{
private:
    MyROOTManager();
    virtual ~MyROOTManager();

    static MyROOTManager* theInstance;
    static int Inversion;

    std::ofstream treeFile;  // file for tree-like data
    std::ofstream histFile;  // file for histogram-like data

public:
    static MyROOTManager* GetPointer();
    void Initialize();
    void Finalize();

    void SetInversion(int);

    void FillHist(double x, double y);
    void FillTree(double x, double y, double e, double k, double fe);
};

#endif
