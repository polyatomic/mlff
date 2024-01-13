#include "Atom.h"

Atom::Atom():
m_elemno(0),
m_FormalCharge(0) {
}

Atom::~Atom() {
}

void Atom::SetElement(const string& element) {
   map<string, unsigned int>::iterator it;
   it = m_elemnos.find(element);
   if (it != m_elemnos.end()) {
      m_elemno = it->second;
   }
}

unsigned int Atom::GetAtomicNumber() const {
   return m_elemno;
}

int Atom::GetCoordinationNumber() const {
   return m_sigmaBonds.size();
}

int Atom::NumberOfNeighbors(int n) const {
   int i, ret=0;
   for (i=0; i < m_sigmaBonds.size(); i++) {
      if (m_sigmaBonds[i]->GetAtomicNumber() == n) ret++;
   }
   return ret;
}

Atom *Atom::GetNeighbor(int i) const {
   if (i >= 0 && i < m_sigmaBonds.size())
      return m_sigmaBonds[i];
   return NULL;
}

int Atom::GetBondType(int i) const {
   if (i >= 0 && i < m_sigmaBonds.size())
      return m_bondTypes[i];
   return 0;
}

int Atom::GetBondStereo(int i) const {
   if (i >= 0 && i < m_sigmaBonds.size())
      return m_bondStereos[i];
   return 0;
}

void Atom::SetBondStereo(int i, int st) {
   if (i >= 0 && i < m_sigmaBonds.size())
      m_bondStereos[i] = st;
}

int Atom::GetFormalCharge() const {
   return m_FormalCharge;
}

bool Atom::Contains(Atom *a) const {
   for (int i=0; i < m_sigmaBonds.size(); i++) {
      Atom *atm = m_sigmaBonds[i];
      if (a == atm) return true;
   }
   return false;
}

const char *Atom::GetAtomicSymbol(unsigned int atomicNumber) {
   map<unsigned int, string>::iterator it;
   it = m_elemsymbols.find(atomicNumber);
   if (it != m_elemsymbols.end()) {
      return (it->second).c_str();
   }
   return m_dummyName.c_str();
}

unsigned int Atom::GetAtomicNumber(const char *symbol) {
   unsigned int elemno = 0;
   map<string, unsigned int>::iterator it;
   it = m_elemnos.find(symbol);
   if (it != m_elemnos.end()) {
      elemno = it->second;
   }
   return elemno;
}

double Atom::GetAveragedWeight(unsigned int atomicNumber) {
   map<unsigned int, double>::iterator it;
   it = m_averagedweights.find(atomicNumber);
   if (it != m_averagedweights.end()) {
      return it->second;
   }
   return 0.0;
}

map<string, unsigned int> Atom::m_elemnos = GetElemNos();
map<unsigned int, string> Atom::m_elemsymbols = GetElemSymbols();
map<unsigned int, double> Atom::m_averagedweights = GetAveragedWeights();
const string Atom::m_dummyName = "";
