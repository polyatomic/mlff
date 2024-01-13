#pragma once

#include <string>
#include <vector>
#include <map>

using std::string;
using std::vector;
using std::map;

class Atom {
public:
   Atom();
   ~Atom();
   void SetElement(const string&);
   unsigned int GetAtomicNumber() const;
   int GetCoordinationNumber() const;
   int NumberOfNeighbors(int) const;
   Atom *GetNeighbor(int) const;
   int GetBondType(int) const;
   int GetBondStereo(int) const;
   void SetBondStereo(int, int);
   int GetFormalCharge() const;
   bool Contains(Atom *) const;

   static const char *GetAtomicSymbol(unsigned int);
   static unsigned int GetAtomicNumber(const char *);
   static double GetAveragedWeight(unsigned int);

private:
   unsigned int m_elemno;
   int m_FormalCharge;
   vector<Atom *> m_sigmaBonds;
   vector<int> m_bondTypes;
   vector<int> m_bondStereos;

   static map<string, unsigned int> GetElemNos() {
      map<string, unsigned int> ret;
      ret["H"] = 1;
      ret["He"] = 2;
      ret["B"] = 5;
      ret["C"] = 6;
      ret["N"] = 7;
      ret["O"] = 8;
      ret["F"] = 9;
      ret["Ne"] = 10;
      ret["Al"] = 13;
      ret["Si"] = 14;
      ret["P"] = 15;
      ret["S"] = 16;
      ret["Cl"] = 17;
      ret["Ar"] = 18;
      ret["Cu"] = 29;
      ret["Ga"] = 31;
      ret["Ge"] = 32;
      ret["As"] = 33;
      ret["Se"] = 34;
      ret["Br"] = 35;
      ret["Kr"] = 36;
      ret["In"] = 49;
      ret["Sn"] = 50;
      ret["Sb"] = 51;
      ret["Te"] = 52;
      ret["I"] = 53;
      ret["Xe"] = 54;
      return ret;
   }
   static map<string, unsigned int> m_elemnos;
   static map<unsigned int, string> GetElemSymbols() {
      map<unsigned int, string> ret;
      ret[1] = "H";
      ret[2] = "He";
      ret[5] = "B";
      ret[6] = "C";
      ret[7] = "N";
      ret[8] = "O";
      ret[9] = "F";
      ret[10] = "Ne";
      ret[13] = "Al";
      ret[14] = "Si";
      ret[15] = "P";
      ret[16] = "S";
      ret[17] = "Cl";
      ret[18] = "Ar";
      ret[29] = "Cu";
      ret[31] = "Ga";
      ret[32] = "Ge";
      ret[33] = "As";
      ret[34] = "Se";
      ret[35] = "Br";
      ret[36] = "Kr";
      ret[49] = "In";
      ret[50] = "Sn";
      ret[51] = "Sb";
      ret[52] = "Te";
      ret[53] = "I";
      ret[54] = "Xe";
      return ret;
   }
   static map<unsigned int, string> m_elemsymbols;
// http://www.webelements.com
   static map<unsigned int, double> GetAveragedWeights() {
      map<unsigned int, double> ret;
      ret[1] = 1.00794;
      ret[2] = 4.002602;
      ret[5] = 10.811;
      ret[6] = 12.0107;
      ret[7] = 14.0067;
      ret[8] = 15.9994;
      ret[9] = 18.9984032;
      ret[10] = 20.1797;
      ret[13] = 26.9815386;
      ret[14] = 28.0855;
      ret[15] = 30.973761;
      ret[16] = 32.065;
      ret[17] = 35.453;
      ret[18] = 39.948;
      ret[29] = 63.546;
      ret[31] = 69.723;
      ret[32] = 72.64;
      ret[33] = 74.9216;
      ret[34] = 78.96;
      ret[35] = 79.904;
      ret[36] = 83.798;
      ret[49] = 114.818;
      ret[50] = 118.71;
      ret[51] = 121.76;
      ret[52] = 127.6;
      ret[53] = 126.90447;
      ret[54] = 131.293;
      return ret;
   }
   static map<unsigned int, double> m_averagedweights;
   static const string m_dummyName;

   friend class Mol;
};
