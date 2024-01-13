#pragma once

#include "Atom.h"

class Mol {
public:
   Mol();
   ~Mol();
   void Init(int, int);
   void SetCoords(const double *);
   void SetCoords(int, const double *);
   void SetAtom(int, const string&);
   void SetComment(const string&);
   string GetComment() const;
   void CreateBond(int, int, int, int);
   void SetFormalCharge(int, int);
   int GetTotalCharge() const;
   unsigned int NAtoms() const;
   unsigned int FindNBonds();
   unsigned int NBonds() const;
   void GetCoordinates(double *) const;
   void GetCoordinates(int, double *) const;
   const char *GetAtomicSymbol(int) const;
   unsigned int GetAtomicNumber(int) const;
   Atom *GetAtom(int) const;
   int GetAtomPos(Atom *) const;
   Mol& operator=(const Mol&);

private:
   void Clear();

   int m_size;
   int m_nBonds;
   double *m_AtomCoordinates;
   Atom *m_Atoms;
   string m_comment;
};
