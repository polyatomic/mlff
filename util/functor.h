#pragma once

class FunctorDaDaI {
public:
   virtual void operator()(double* darg1, double* darg2, int iarg) = 0;
   virtual void Call(double* darg1, double* darg2, int iarg) = 0;
};

template<class FClass>
class TFunctorDaDaI : public FunctorDaDaI {
public:
   TFunctorDaDaI() {
      pt2Object = 0;
      fpt = 0;
   }

   TFunctorDaDaI(FClass* _pt2Object, void(FClass::* _fpt)(double* darg1, double* darg2, int iarg)) {
      pt2Object = _pt2Object;
      fpt = _fpt;
   }

   void init(FClass* _pt2Object, void(FClass::* _fpt)(double* darg1, double* darg2, int iarg)) {
      pt2Object = _pt2Object;
      fpt = _fpt;
   }

   virtual void operator()(double* darg1, double* darg2, int iarg) {
      (*pt2Object.*fpt)(darg1, darg2, iarg);
   }

   virtual void Call(double* darg1, double* darg2, int iarg) {
      (*pt2Object.*fpt)(darg1, darg2, iarg);
   }

private:
   void (FClass::* fpt)(double* darg1, double* darg2, int iarg);
   FClass* pt2Object;
};
