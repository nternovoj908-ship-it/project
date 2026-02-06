#include "RootManager.hh"
#include <iostream>
#include <fstream>

MyROOTManager* MyROOTManager::theInstance = 0;  // initialize "global" variables
int MyROOTManager::Inversion = 0;

MyROOTManager::MyROOTManager()
{
    Initialize();
}

MyROOTManager::~MyROOTManager()
{
    Finalize();
}

MyROOTManager* MyROOTManager::GetPointer()
{
    if (!theInstance) theInstance = new MyROOTManager();
    return theInstance;
}

void MyROOTManager::Initialize()
{
    // Open files for writing
    treeFile.open("MyTree.txt");
    treeFile << "POSX\tPOSY\tEDEP\tEKIN\tFullEne\n";  // header

    histFile.open("image.txt");
    histFile << "X\tY\n";  // header
}

void MyROOTManager::FillHist(double x, double y)
{
    histFile << x << "\t" << y << "\n";
    histFile.flush();  // <-- Добавлено: сброс буфера
}

void MyROOTManager::FillTree(double x, double y, double e, double k, double fe)
{
    treeFile << x << "\t" << y << "\t" << e << "\t" << k << "\t" << fe << "\n";
    treeFile.flush();  // <-- Добавлено: сброс буфера
}

void MyROOTManager::SetInversion(int x)
{
    Inversion = x;
    if (Inversion) {
        std::cout << "Mask inversion ON\n";
    } else {
        std::cout << "Mask inversion OFF\n";
    }
}

void MyROOTManager::Finalize()
{
    treeFile.close();
    histFile.close();

    if (Inversion == 0) {
        std::cout << "Results saved to MyTree.txt and image.txt\n";
    } else {
        std::cout << "Results saved to MyTree.txt and image.txt (with inversion)\n";
    }
}
