#include "tntsupp.h"
#include "geese.h"
#include "famstr.h"
#include "param.h"
#include "geesubs.h"
#include "utils.h"

#include "inter.h"

//#include "lgtdl.h"
//#include "fgee.h"

/*******************************************************************/

Grad & Grad::operator=(const Grad &G) {
  U1_ = G.U1_; U2_ = G.U2_; U3_ = G.U3_;
  return *this;
}

ostream& operator<<(ostream& s, const Grad &G) {
  s << "U1 = " << G.U1() << "U2 = " << G.U2() << "U3 = " << G.U3(); 
  return s;
}

Hess operator-(Hess &H1, Hess &H2) {
  Hess ans(H1);
  ans.dec(H2);
  return ans;
}


Hess inv(Hess &H, IVector &level) {
  Hess ans(H);
  ans.set_A(inv(H.A()));
  if (level(1) == 1) {
    ans.set_C(inv(H.C()));
    ans.set_B(-1.0 * ans.C() * H.B() * ans.A());
  }
  if (level(2) == 1) {
    ans.set_F(inv(H.F()));
    ans.set_E(-1.0 * ans.F() * H.E() * ans.C());
    ans.set_D(-1.0 * ans.F() * (H.D() * ans.A() + H.E() * ans.B()));
  }
  return ans;
}

Hess operator*(const double &x, const Hess &H) {
  Hess ans(H);
  ans.set_A(x * H.A()); ans.set_B(x * H.B()); ans.set_C(x * H.C());
  ans.set_D(x * H.D()); ans.set_E(x * H.E()); ans.set_F(x * H.F());
  return ans;
}

ostream& operator<<(ostream& s, const Hess &H) {
  s << "A = " << H.A() << "B = " << H.B() << "C = " << H.C() 
    << "D = " << H.D()<< "E = " << H.E() << "F = " << H.F();
  return s;
}

DVector genzi(const DVector &PR) {
  int n = PR.size();
  DVector ans(n * (n - 1)/2);
  int k = 1;
  for (int i = 1; i <= n - 1; i++)
    for (int j = i + 1; j <= n; j++) ans(k++) = PR(i) * PR(j);
  return ans;
}

DVector utri(const DMatrix &R) {
  int n = R.dim(1);
  //assert (n > 1);
  DVector ans(n * (n - 1) / 2);
  int k = 1;
  for (int i = 1; i <= n - 1; i++)
    for (int j = i + 1; j <= n; j++) ans(k++) = R(i,j);
  return ans;
}

DMatrix getZ_Beta(DMatrix &D, DVector &PR, 
		  DVector &V, DVector &V_Mu, DVector &z) {
  //note: this is the version which excludes phi in the formula
  DMatrix ans(z.size(), D.dim(2));
  int k = 1, n = PR.size();
  for (int i = 1; i <= n - 1; i++) {
    DMatrix Di = asMat(MatRow(D,i));
    for (int j = i + 1; j <= n; j++) {
      DMatrix Dj = asMat(MatRow(D,j));
      DMatrix foo = V_Mu(i) * reciproot(V(i)) * Di + 
	V_Mu(j) * reciproot(V(j)) * Dj;
      DMatrix bar = - PR(i) * Di - PR(j) * Dj - 0.5 * PR(i) * PR(j) * foo;
      //cout << "bar = " << bar << "k = " << k;
      MatRow(ans, k) = bar;
      //cout << " ans = " << ans;
      k++;
    }
  }
  return ans;
}
      
DMatrix getZ_Gamma(DMatrix &D, DVector &PR, DVector &Phi, DVector &z) {
  DMatrix ans(z.size(), D.dim(2));
  int k = 1, n = PR.size();
  for (int i = 1; i <= n - 1; i++) {
    DMatrix Di = asMat(MatRow(D,i));
    for (int j = i + 1; j <= n; j++) {
      DMatrix Dj = asMat(MatRow(D,j));
      //MatRow(ans, k) = -0.5 * z(k) * (sqrt(Phi(j)/Phi(i)) * Di +
      //sqrt(Phi(i)/Phi(j)) * Dj);
      //This has caused the scale problem; The first time, scale problem was caused by operator * for Hess, where one component did not get touched, in the old geese (LAPACK);
      MatRow(ans, k) = -0.5 * z(k) * (1.0 / Phi(i) * Di +
				      1.0 / Phi(j) * Dj);
      k++;
    }
  }
  return ans;
}

DMatrix getS_Beta(DMatrix &D, DVector &PR, DVector &V, DVector &V_Mu) {
  DMatrix ans(D);
  for (int i = 1; i <= ans.dim(1); i++) {
    DMatrix Di = asMat(MatRow(D,i));
    double f = -2 * PR(i) / sqrt(V(i)) - PR(i) * PR(i)/V(i) * V_Mu(i);
    MatRow(ans, i) = f * Di;
  }
  return ans;
}

void HiandGi(DVector &PRi, DVector &Phii, DMatrix &Di, DMatrix &R,
	     DVector &Vi, DVector &V_Mui, DMatrix &D2i, DMatrix &E,
	     DVector &Wi, IVector &level,
	     //output
	     Hess &H, Grad &G) {
  int s = PRi.size();
  //beta
  DVector rootPhii = sqrt(Phii);
  DMatrix V1 = diag(rootPhii) * R * diag(rootPhii);
  DVector rootWi = sqrt(Wi);
  DMatrix rootWD = SMult(rootWi, Di); 
  DVector rootWPR = SMult(rootWi, PRi);
  H.set_A(AtBiC(rootWD, V1, rootWD));
  G.set_U1(AtBiC(rootWD, V1, rootWPR));
  //H.set_A(AtBiC(Di, V1, Di));
  //G.set_U1(AtBiC(Di, V1, PRi));
  //gamma
  if (level(1) == 1) {//if (par.ScaleFix() != 1) {
    DVector Si = square(PRi);
    DVector WiV2inv = SMult(Wi, recip(2.0 * Phii));
    H.set_C(Transpose_view(D2i) * SMult(WiV2inv, D2i));
    DMatrix S_Beta = getS_Beta(Di, PRi, Vi, V_Mui);
    H.set_B(Transpose_view(D2i) * SMult(-1.0 * WiV2inv, S_Beta));
    G.set_U2(Transpose_view(D2i) * SMult(WiV2inv, Si - Phii));
    //DMatrix V2 = diag(2.0 * Phii);
    //H.set_C(AtBiC(D2i, V2, D2i));
    //DMatrix S_Beta = getS_Beta(Di, PRi, Vi, V_Mui);
    //H.set_B(AtBiC(D2i, V2, S_Beta));
    //G.set_U2(AtBiC(D2i, V2, S - Phii));
  }
  //alpha
  if (level(2) == 1) {//if (cor.nparam() > 0) {
    if (s == 1) return;
    DVector sPRi = SMult(reciproot(Phii), PRi);
    DVector zi = genzi(sPRi);
    DVector rhoi = utri(R); 
    //DMatrix W = ident(s * (s - 1) / 2);
    DVector Sca = genzi(reciproot(Phii)); 
    DVector WiV3inv = genzi(Wi);
    
    //H.set_F(AtBiC(E, W, E));
    H.set_F(Transpose_view(E) * SMult(WiV3inv, E));
    DMatrix Z_Beta = getZ_Beta(Di, PRi, Vi, V_Mui, zi);
    Z_Beta = SMult(Sca, Z_Beta);
    //H.set_D(AtBiC(E, W, Z_Beta));
    H.set_D(Transpose_view(E) * SMult(-1.0 * WiV3inv, Z_Beta));
    //G.set_U3(AtBiC(E, W, zi - rhoi));
    G.set_U3(Transpose_view(E) * SMult(WiV3inv, zi - rhoi));
    if (level(1) == 1) {//if (par.ScaleFix() != 1) {
      DMatrix Z_Gamma = getZ_Gamma(D2i, PRi, Phii, zi);
      //H.set_E(AtBiC(E, W, Z_Gamma));
      H.set_E(Transpose_view(E) * SMult(-1.0 * WiV3inv, Z_Gamma));
    }
  }
}

void PRandD(DVector &Y, DMatrix &X, DVector &Offset,
	    Index1D &I, IVector &LinkWave, 
	    GeeParam &par, GeeStr &geestr,
	    DVector &PRi, DMatrix &Di) {
  DVector Yi = asVec(VecSubs(Y, I));
  DMatrix Xi = asMat(MatRows(X, I));
  DVector Offseti = asVec(VecSubs(Offset, I));
  IVector Wavei = asVec(VecSubs(LinkWave, I));
  DVector Eta = Xi * par.beta() + Offseti;
  DVector Mu = geestr.MeanLinkinv(Eta, Wavei);
  DVector V = geestr.v(Mu, Wavei);
  DVector Mu_Eta = geestr.MeanMu_eta(Eta, Wavei);

  DVector InvRootV = reciproot(V);
  Di = SMult(InvRootV, SMult(Mu_Eta, Xi));
  PRi = SMult(InvRootV, Yi - Mu);
}

void gee_prep(DVector &Y, DMatrix &X, DVector &Offset,
	      Index1D &I, IVector &LinkWave, 
	      GeeParam &par, GeeStr &geestr,
	      DVector &PRi, DMatrix &Di, DVector &Vi, DVector &V_Mui) {
  DVector Yi = asVec(VecSubs(Y, I));
  DMatrix Xi = asMat(MatRows(X, I));
  DVector Offseti = asVec(VecSubs(Offset, I));
  IVector Wavei = asVec(VecSubs(LinkWave, I));
  DVector Eta = Xi * par.beta() + Offseti;
  DVector Mu = geestr.MeanLinkinv(Eta, Wavei);
  DVector V = geestr.v(Mu, Wavei);
  DVector Mu_Eta = geestr.MeanMu_eta(Eta, Wavei);

  DVector InvRootV = reciproot(V);
  Di = SMult(InvRootV, SMult(Mu_Eta, Xi));
  PRi = SMult(InvRootV, Yi - Mu);
  Vi = geestr.v(Mu, Wavei);
  V_Mui = geestr.v_mu(Mu, Wavei);
}

DMatrix getR(DMatrix &Zmat, Index1D &I, Index1D &J, DVector &CorP,
	     GeeParam &par, GeeStr &geestr, Corr &cor) {
  DVector alp = par.alpha(); 
  DVector corp = asVec(VecSubs(CorP, I)); 
  int s = corp.dim(); // corp should determine meta par for R
  if (s == 1) return ident(1); 
  else if (cor.nparam() == 0) //indenpendence or fixed
    return cor.mat(alp, corp); //if fixed, should have ident link
  else {
    DMatrix Zmati = asMat(MatRows(Zmat, J));
    DVector Eta = Zmati * alp;
    DVector Rho = geestr.CorrLinkinv(Eta);
    return cor.mat(Rho, corp);
  }
}

int RandE(DMatrix &Zmat, Index1D &I, Index1D &J, DVector &CorP,
	  GeeParam &par, GeeStr &geestr, Corr &cor,
	  DMatrix &R, DMatrix &E) {
  DVector alp = par.alpha();
  DVector corp = asVec(VecSubs(CorP, I));
  int s = corp.dim();
  if (s == 1) {
    R = ident(1); 
    return 0;
  }
  else if (cor.nparam() == 0) { //no need for E
    R = cor.mat(alp, corp);
    return 0;
  }
  else {
    DMatrix Zmati = asMat(MatRows(Zmat, J)); 
    DVector Eta = Zmati * alp;
    DVector Rho = geestr.CorrLinkinv(Eta);
    R = cor.mat(Rho, corp);
    DVector Rho_Alp = geestr.CorrMu_eta(Eta);
    DMatrix Cor_Rho = cor.cor_rho(Rho, corp);
    E = Cor_Rho * SMult(Rho_Alp,  Zmati);
    return 0;
  }
}

void gm_prep(DVector &PR, Index1D &I, IVector &LinkWave,
	     DVector &Doffset, DMatrix &Zsca, GeeParam &par, GeeStr &geestr,
	     DVector &Phii, DVector &Si, DMatrix &D2i) {
  DMatrix Zi = asMat(MatRows(Zsca, I));
  DVector Doffseti = asVec(VecSubs(Doffset, I));
  IVector Wavei = asVec(VecSubs(LinkWave, I));
  DVector Zeta = Zi * par.gamma() + Doffseti;
  DVector Phi_Zeta = geestr.ScaleMu_eta(Zeta, Wavei);
  DVector PRi = asVec(VecSubs(PR, I));
  
  Phii = geestr.ScaleLinkinv(Zeta, Wavei);
  Si = square(PRi);
  D2i = Phi_Zeta * Zi;
}

void PhiandD2(Index1D &I, IVector &LinkWave,
	      DVector &Doffset, DMatrix &Zsca, GeeParam &par, GeeStr &geestr,
	      DVector &Phii, DMatrix &D2i) {
  DMatrix Zi = asMat(MatRows(Zsca, I));
  DVector Doffseti = asVec(VecSubs(Doffset, I));
  IVector Wavei = asVec(VecSubs(LinkWave, I));
  DVector Zeta = Zi * par.gamma() + Doffseti;
  Phii = geestr.ScaleLinkinv(Zeta, Wavei);
  if (geestr.ScaleFix() == 1) return;
  DVector Phi_Zeta = geestr.ScaleMu_eta(Zeta, Wavei);
  D2i = Phi_Zeta * Zi;
}

DVector getPR(DVector &Y, DMatrix &X, DVector &Offset, IVector &LinkWave,
	      GeeParam &par, GeeStr &geestr) {
  DVector Eta = X * par.beta() + Offset;
  DVector Mu = geestr.MeanLinkinv(Eta, LinkWave);
  DVector V = geestr.v(Mu, LinkWave);
  DVector InvRootV = reciproot(V);
  return SMult(InvRootV, Y - Mu);
}

DVector getPhi(DVector &Doffset, DMatrix &Zsca, IVector &LinkWave,
	       GeeParam &par, GeeStr &geestr) {
  DVector Zeta = Zsca * par.gamma() + Doffset;
  return geestr.ScaleLinkinv(Zeta, LinkWave);
}

void HnandGis(DVector &Ycur, DMatrix &X, 
	      DVector &Offset, DVector &Doffset, DVector &W, 
	      IVector &LinkWave, IVector &Clusz, IVector &ZcorSize,
	      DMatrix &Zsca, DMatrix &Zcor, DVector &CorP,
	      GeeParam &par, GeeStr &geestr, Corr &cor,
	      IVector &Scur, IVector &level, //Scur is the valid data indicator
	      //output
	      Hess &Hn, Vector<Grad> &Gis) {
  Index1D I(0,0), J(0,0);
  int N = Clusz.size();
  int pb = par.p(), pg = par.r(), pa = par.q();
  DVector V0(pb);
  Hess H(par), Hi(par); Grad Gi(par);
  //cout << "N = " << N;
  for (int i = 1; i <= N; i++) {
    int s1 = Clusz(i), s2 = ZcorSize(i);
    I = Index1D(1, s1) + I.ubound();
    if (s2 > 0) J = Index1D(1, s2) + J.ubound();
    if (Scur(i) == 0) {
      //Gis(i).set_U1(V0);
      continue;
    }
    DVector PRi(s1), Vi(s1), V_Mui(s1); DMatrix Di(s1,pb);
    gee_prep(Ycur, X, Offset, I, LinkWave, par, geestr, PRi, Di, Vi, V_Mui);
    DVector Phii(s1); DMatrix D2i(s1, pg);
    PhiandD2(I, LinkWave, Doffset, Zsca, par, geestr, Phii, D2i);
    DMatrix R(s1, s1), E(s2, pa);
    RandE(Zcor, I, J, CorP, par, geestr, cor, R, E);

    DVector Wi = asVec(VecSubs(W, I));
    HiandGi(PRi, Phii, Di, R, Vi, V_Mui, D2i, E, Wi, level, Hi, Gi);
    H.inc(Hi); Gis(i) = Gi;
  }
  Hn = (1.0/(double) N) * H;
}

void HisandGis(DVector &Ycur, DMatrix &X, 
	      DVector &Offset, DVector &Doffset, DVector &W, 
	      IVector &LinkWave, IVector &Clusz, IVector &ZcorSize,
	      DMatrix &Zsca, DMatrix &Zcor, DVector &CorP,
	      GeeParam &par, GeeStr &geestr, Corr &cor,
	      IVector &Scur, IVector &level,
	      //output
	      Vector<Hess> &His, Vector<Grad> &Gis) {
  Index1D I(0,0), J(0,0);
  int N = Clusz.size();
  int pb = par.p(), pg = par.r(), pa = par.q();
  DVector V0(pb);
  Hess H(par), Hi(par); Grad Gi(par);
  //cout << "N = " << N;
  for (int i = 1; i <= N; i++) {
    int s1 = Clusz(i), s2 = ZcorSize(i);
    I = Index1D(1, s1) + I.ubound();
    if (s2 > 0) J = Index1D(1, s2) + J.ubound();
    if (Scur(i) == 0) {
      //Gis(i).set_U1(V0);
      continue;
    }
    DVector PRi(s1), Vi(s1), V_Mui(s1); DMatrix Di(s1,pb);
    gee_prep(Ycur, X, Offset, I, LinkWave, par, geestr, PRi, Di, Vi, V_Mui);
    DVector Phii(s1); DMatrix D2i(s1, pg);
    PhiandD2(I, LinkWave, Doffset, Zsca, par, geestr, Phii, D2i);
    DMatrix R(s1, s1), E(s2, pa);
    RandE(Zcor, I, J, CorP, par, geestr, cor, R, E);

    DVector Wi = asVec(VecSubs(W, I));
    HiandGi(PRi, Phii, Di, R, Vi, V_Mui, D2i, E, Wi, level, Hi, Gi);
    //H.inc(Hi); Gis(i) = Gi;
    His(i) = Hi; Gis(i) = Gi;
  }
  //Hn = (1.0/(double) N) * H;
}

/*
DVector interpprev(double t, Vector<DVector> &VV, DVector &tis) {
  if (t < tis(1)) return VV(1);
  if (t >= tis(tis.size())) return VV(VV.size());
  for (int i = 1; i <= tis.size() - 1; i++) 
    if (tis(i) <= t && tis(i + 1) > t) return VV(i);
}

double interpprev(double t, DVector &v, DVector &tis) {
  if (t < v(1)) return v(1);
  if (t >= tis(tis.size())) return v(v.size());
  for (int i = 1; i <= tis.size() - 1; i++) 
    if (tis(i) <= t && tis(i + 1) > t) return v(i);
}

//double interplinear(double t, DVector &

DVector getY(double t, DVector &Yall) {
  DVector ans(Yall.size());
  for (int i = 1; i <= ans.size(); i++)
    ans(i) = (Yall(i) <= t) ? 1.0 : 0.0;
  return ans;
}

IVector getS(double t, DVector &S) {
  IVector ans(S.size());
  for (int i = 1; i <= ans.size(); i++) 
    ans(i) = (S(i) <= t) ? 0 : 1;
  return ans;
}

void AandCis(Vector<Lgtdl> &Yall, DMatrix &X, 
	     Vector<DVector> &Offset, Vector<DVector> &Doffset,
	     Vector<DVector> &Weight,
	     IVector &LinkWave, 
	     IVector &Clusz, IVector &ZcorSize,
	     DMatrix &Zsca, DMatrix &Zcor, DVector &CorP,
	     Vector<DVector> &Beta, 
	     Vector<DVector> &Gamma, 
	     Vector<DVector> &Alpha,
	     DVector &W, DVector &S, DVector &Tis, DVector &Tlim,
	     int l, DVector &Ita,
	     GeeStr &geestr, Corr &cor, 
	     int ndivs, int fgconf,
	     //output:
	     DMatrix &A, Vector<DVector> &Cis) {
  int N = Clusz.size();
  double det = (Tlim(2) - Tlim(1))/(double) ndivs;
  Fgee fg(fgconf);
  
  for (double t = Tlim(1) + det/2; t <= Tlim(2); t += det) {
    DVector g = fg.g(t, Ita);
    double  w = interpprev(t, W, Tis);

    DVector beta = interpprev(t, Beta, Tis);
    DVector gamma = interpprev(t, Gamma, Tis);
    DVector alpha = interpprev(t, Alpha, Tis);
    DVector offset = interpprev(t, Offset, Tis);
    DVector doffset = interpprev(t, Doffset, Tis);
    DVector weight = interpprev(t, Weight, Tis);
    beta(l) = fg.f(t, Ita);  //l is the lth component in beta;

    //cout << "beta = " << beta;
    GeeParam par(beta, alpha, gamma);
    int pb = par.p(), pg = par.r(), pa = par.q();
    DMatrix L(1, pb); L(1,l) = 1.0;
    //cout << "L = " << L;
    
    DVector Ycur = getY(t, Yall); //Y is the time to event;
    IVector Scur = getS(t, S); //S is the time to screened out;
    // X, Zcor, Zsca are all time independent;

    Hess Hn(par); Vector<Grad> Gis(N); Grad G0(par); Gis = G0;
    IVector level(2); level(1) = 0; level(2) = 0;
    HnandGis(Ycur, X, offset, doffset, weight, 
	     LinkWave, Clusz, ZcorSize, Zsca, Zcor, CorP,
	     par, geestr, cor, Scur, level, Hn, Gis);
    //cout << "H = " << Hn.A();
    Hess Hinv = inv(Hn, level);
    //cout << "Hinv = " << Hinv.A();
    for (int i = 1; i <= N; i++) {
      if (Scur(i) == 0) continue; //screening
      Cis(i) = Cis(i) + w * asColMat(g) * L * Hinv.A() * Gis(i).U1();
    }
    A = A + w * outerprod(g, g);
  }
  A = inv(det * A);
  for (int i = 1; i <= N; i++) Cis(i) = det * Cis(i);
}
*/
