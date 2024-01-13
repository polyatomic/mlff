#include "Mol.h"

Mol::Mol():
m_size(0),
m_AtomCoordinates(0),
m_Atoms(0),
m_nBonds(0) {
}

Mol::~Mol() {
   Clear();
}

void Mol::Init(int nAtoms, int nBonds) {
   Clear();
   m_size = nAtoms;
   m_nBonds = nBonds;
   if (nAtoms > 0) {
      m_AtomCoordinates = new double[3*nAtoms];
      m_Atoms = new Atom[nAtoms];
   }
}

void Mol::SetCoords(const double *coords) {
   int i;
   for (i=0; i < 3*m_size; i++) {
      m_AtomCoordinates[i] = coords[i];
   }
}

void Mol::SetCoords(int n, const double *coords) {
   if (n < 0 || n >= m_size) return;
   m_AtomCoordinates[3*n] = coords[0];
   m_AtomCoordinates[3*n+1] = coords[1];
   m_AtomCoordinates[3*n+2] = coords[2];
}

void Mol::SetAtom(int n, const string& element) {
   if (n < 0 || n >= m_size) return;
   m_Atoms[n].SetElement(element);
}

void Mol::SetComment(const string& c) {
   m_comment = c;
}

string Mol::GetComment() const {
   return m_comment;
}

void Mol::CreateBond(int bnda1, int bnda2, int mlp, int st) {
   if (bnda1 < 0 || bnda1 >= m_size || bnda2 < 0 || bnda2 >= m_size) return;
   m_Atoms[bnda1].m_sigmaBonds.push_back(&m_Atoms[bnda2]);
   m_Atoms[bnda1].m_bondTypes.push_back(mlp);
   m_Atoms[bnda1].m_bondStereos.push_back(st);
   m_Atoms[bnda2].m_sigmaBonds.push_back(&m_Atoms[bnda1]);
   m_Atoms[bnda2].m_bondTypes.push_back(mlp);
   m_Atoms[bnda2].m_bondStereos.push_back(0);
}

void Mol::SetFormalCharge(int n, int fc) {
   if (n < 0 || n >= m_size) return;
   m_Atoms[n].m_FormalCharge = fc;
}

int Mol::GetTotalCharge() const {
   int i;
   int c = 0;
   for (i=0; i < m_size; i++) {
      c += m_Atoms[i].m_FormalCharge;
   }
   return c;
}

unsigned int Mol::NAtoms() const {
   return m_size;
}

unsigned int Mol::FindNBonds() {
   int i;
   m_nBonds = 0;
   for (i=0; i < m_size; i++) {
      m_nBonds += m_Atoms[i].GetCoordinationNumber();
   }
   m_nBonds /= 2;
   return m_nBonds;
}

unsigned int Mol::NBonds() const {
   return m_nBonds;
}

void Mol::GetCoordinates(double *ret) const {
   int i;
   for (i=0; i < 3*m_size; i++) {
      ret[i] = m_AtomCoordinates[i];
   }
}

void Mol::GetCoordinates(int n, double *ret) const {
   int i;
   if (n < 0 || n >= m_size) return;
   for (i=0; i < 3; i++) {
      ret[i] = m_AtomCoordinates[3*n+i];
   }
}

const char *Mol::GetAtomicSymbol(int n) const {
   if (n < 0 || n >= m_size) return NULL;
   return Atom::GetAtomicSymbol(m_Atoms[n].GetAtomicNumber());
}

unsigned int Mol::GetAtomicNumber(int n) const {
   if (n < 0 || n >= m_size) return 0;
   return m_Atoms[n].GetAtomicNumber();
}

Atom *Mol::GetAtom(int n) const {
   if (n < 0 || n >= m_size) return NULL;
   return m_Atoms+n;
}

int Mol::GetAtomPos(Atom *a) const {
   int n;
   n = a - m_Atoms;
   if (n < 0 || n >= m_size) return 0;
   return n;
}

Mol& Mol::operator=(const Mol &rhs) {
   int i, j;
   if (this == &rhs)
      return *this;
   Init(rhs.m_size, rhs.m_nBonds);
   for (i=0; i < 3*m_size; i++)
      m_AtomCoordinates[i] = rhs.m_AtomCoordinates[i];
   for (i=0; i < m_size; i++) {
      m_Atoms[i] = rhs.m_Atoms[i];
      for (j=0; j < m_Atoms[i].m_sigmaBonds.size(); j++)
         m_Atoms[i].m_sigmaBonds[j] = m_Atoms + rhs.GetAtomPos(rhs.m_Atoms[i].m_sigmaBonds[j]);
   }
   m_comment = rhs.m_comment;
   return *this;
}

void Mol::Clear() {
   delete [] m_Atoms;
   delete [] m_AtomCoordinates;
   m_size = 0;
   m_nBonds = 0;
   m_Atoms = 0;
   m_AtomCoordinates = 0;
   m_comment = "";
}
