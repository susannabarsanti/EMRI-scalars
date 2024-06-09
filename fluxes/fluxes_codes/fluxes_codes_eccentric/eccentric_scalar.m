(* ::Package:: *)

(* ::Title:: *)
(*Code for the eccentric scalar fluxes computation*)


(* ::CodeText:: *)
(*If you make use of this code, please cite https://arxiv.org/abs/2203.05003, published in: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.106.044029.*)
(**)
(*Moreover, note that: *)
(*- The convergence criteria have been implemented by following Viktor Skoupy: https://arxiv.org/abs/2201.07044. *)
(*- The functions for the homogeneous solutions have been implemented by Gabriel Piovano, i.e. boundary  conditions for the radial Teukolsky equation in horizon penetrating, hyperboloidal slicing coordinates. *)
(*If you make use of these boundary conditions, please acknowledge https://arxiv.org/pdf/2105.07083.  *)
(**)
(*Finally:	*)
(*This notebook makes use  of the BHPToolkit https://bhptoolkit.org/mathematica-install-dev.html*)
(**)
(*(*IMPORTANT: FindRoot errors in the paclet version of the toolkit, use the GitHub one.*)*)


SetSystemOptions["ParallelOptions"->"MKLThreadNumber"->40];
SetSystemOptions["ParallelOptions"->"ParallelThreadNumber"->40];
LaunchKernels[40];


Exit[];


<<SpinWeightedSpheroidalHarmonics`
<<KerrGeodesics`
<<Teukolsky`
ParallelNeeds["SpinWeightedSpheroidalHarmonics`"];
ParallelNeeds["KerrGeodesics`"];
ParallelNeeds["Teukolsky`"];


(* ::Subsection::Closed:: *)
(*Useful orbital functions *)


(*Useful functions: energy, angular momentum, \[Phi](\[Chi]), t(\[Chi]), Vr(\[Chi]), V\[Phi](\[Chi]), Vt(\[Chi]), J(\[Chi]) *)

en[aa_,p_,e_]:= Block[{en,a},
If[aa<0,
a=-aa;
en = KerrGeoEnergy[a,p,e,-1];
,
a=aa;
en= KerrGeoEnergy[a,p,e,1]];
Return[en];
];

L[aa_,p_,e_]:= Block[{el,a},
If[aa<0,
a=-aa;
(*this line was due to a problem in the toolkit's function, which it is now fixed.
If[e==0, el =-KerrGeoAngularMomentum[a,p,e,-1];,el=KerrGeoAngularMomentum[a,p,e,-1];];*)
el=KerrGeoAngularMomentum[a,p,e,-1];
, 
a=aa;
el=KerrGeoAngularMomentum[a,p,e,1]; 
];
Return[el];
];

\[Phi][a_,p_,e_,\[Chi]_?NumberQ]:=NIntegrate[V\[Phi][a,p,e,chi]/(J[a,p,e,chi] (Vr[a,p,e,chi])^(1/2)),{chi, 0, \[Chi]},WorkingPrecision->ppp];
t[a_,p_,e_,\[Chi]_?NumberQ]:=NIntegrate[Vt[a,p,e,chi]/(J[a,p,e,chi] (Vr[a,p,e,chi])^(1/2)),{chi, 0, \[Chi]},WorkingPrecision->ppp];

x[aa_,p_,e_]:=Block[{ics,a},
If[aa<0, a=-aa;, a=aa;];
ics = L[aa,p,e] - a en[aa,p,e];
Return[ics];
]

Vr[aa_,p_,e_,\[Chi]_]:= Block[{vr,a},
If[aa<0, a=-aa;, a=aa;];
vr=x[aa,p,e]^2+a^2+2 a x[aa,p,e] en[aa,p,e]-(2 M x[aa,p,e]^2)/p (3 + e Cos[\[Chi]]) ;
Return[vr];
]

V\[Phi][aa_,p_,e_,\[Chi]_]:=Block[{vphi,a},
If[aa<0, a=-aa;, a=aa;];
vphi=x[aa,p,e] + a en[aa,p,e] - (2 M x[aa,p,e])/p (1+ e Cos[\[Chi]]);
Return[vphi];
]

Vt[aa_,p_,e_,\[Chi]_]:=Block[{vt,a},
If[aa<0, a=-aa;, a=aa;];
vt = a^2 en[aa,p,e]-(2 a M x[aa,p,e] )/p (1+ e Cos[\[Chi]]) + (en[aa,p,e] p^2)/(1+ e Cos[\[Chi]])^2;
Return[vt];
]

J[aa_,p_,e_,\[Chi]_] := Block[{j,a},
If[aa<0, a=-aa;, a=aa;];
j=1 - (2 M)/p (1+ e Cos[\[Chi]])+ (a/p)^2 (1+ e Cos[\[Chi]])^2;
Return[j];
]

erre[e_,p_,\[Chi]_]:= p/(1 + e Cos[\[Chi]]);


ISCOradius[a_,orbit_?StringQ]:=Module[{M=1,Z1,Z2},
  Z1 = 1+(1-a^2/M^2)^(1/3)((1+a/M)^(1/3)+(1-a/M)^(1/3));
  Z2 = Sqrt[(3a^2/M^2+Z1^2)];

  Which[orbit=="pro",
   M(3+Z2-Sqrt[(3-Z1)(3+Z1+2Z2)]),
   orbit=="ret",
   M(3+Z2+Sqrt[(3-Z1)(3+Z1+2Z2)])
  ]
]


(* ::Subsection::Closed:: *)
(*Homogeneous solutions*)


(*Functions for the homogeneus solutions implemented by Gabriel Piovano, *)

TeukolskyHSCoefficients := 
 Module[{M=1, r, m , a, \[Omega], \[Lambda], \[CapitalDelta],s, f,H,Gtilde, Utilde, p, q, coefficients},
  \[CapitalDelta] = r^2-2M r+a^2;
  f = \[CapitalDelta]/(r^2+a^2); (*dr/drstar*)
  Gtilde = a^2 \[CapitalDelta]+(r^2+a^2)(r s (r-M)-I r((r^2+a^2)\[Omega] H+m a));
  Utilde = 2 I s \[Omega] r^2(r \[CapitalDelta] (1-H) -M (r^2-a^2)(1+H))-2 I a r \[CapitalDelta] (m+ a \[Omega] H)+\[CapitalDelta] (2 a^2-r^2 \[Lambda]-2 M r (s+1))-2 m a \[Omega] r^2 (r^2+a^2)(1+H)+r^2 (r^2+a^2)^2 (\[Omega]^2 (1-H^2)+I \[Omega] f D[H,r]);
  p = (r^2+a^2)/\[CapitalDelta] D[f,r]-1/(\[CapitalDelta](r^2+a^2)) (2Gtilde)/r;
  q =  Utilde/(r^2 \[CapitalDelta]^2);

  coefficients = {Function@@{p}/.Thread[{r,H,s,m,a,\[Omega],\[Lambda]}->Array[Slot,7]], Function@@{q}/.Thread[{r,H,s,m,a,\[Omega],\[Lambda]}->Array[Slot,7]]};
  Remove[r,H,s,m,a,\[Omega],\[Lambda]];
  coefficients
]


bchor[workingprecision_,s_,m_,a_,\[Omega]_,\[Lambda]_]:=Module[{M=1,deltarp,A2n,rp,rm,rin,err,i,chor,p,q,phor,qhor,dphor,dqhor,\[Psi]hor,a2n,cInHor},
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
deltarp=(rp-rm)/50;
rin=rp+deltarp; (*Closest r to the horizon *)
chor=2I ( 2M rp)/(rp-rm) (\[Omega]-(a m)/(2 M rp))+ s;

dqhor[n_]:=Which[
n==0,
0,
n==1,
(2 I a m+2 M (-1+s)-2 I a^2 \[Omega]+rp (2+\[Lambda]-4 I rp s \[Omega]))/((rm-rp) rp),
n>1,
2 (-1+n) (-rp)^-n+(rm-rp)^-n (2+\[Lambda]-4 I rm s \[Omega])+1/rp 2 n (rm-rp)^-n (M (-1+s)+I a (m-a \[Omega])) Hypergeometric2F1[1,1-n,2,rm/rp]
];

dphor[n_]:=Which[
n==0,
1-chor,
n==1,
 1/((rm-rp)^2 rp) (-2 rm^2+a^2 (3+2 s+4 I rp \[Omega])+I rp (-2 a m+2 a^2 \[Omega]+I (rp+2 M s+2 I rp^2 \[Omega]))),
n>1,
2 (-rp)^-n-(rm-rp)^-n+(rm-rp)^(-n-1) (rm (+2  s)+2 I rm^2 \[Omega]+2 I (-a m+I M s+a^2 \[Omega]))
];

a2n[0]=1;
A2n[n_]:= -(1/(n(n-chor)))Sum[(j dphor[n-j]+dqhor[n-j])a2n[j],{j,0,n-1}];

err=1;
i=1;
{p,q}=TeukolskyHSCoefficients;
phor=p[rin,-1,s,m,a,\[Omega],\[Lambda]];
qhor=q[rin,-1,s,m,a,\[Omega],\[Lambda]];

While[err > 10^(-workingprecision),
(*Apparently,it works better when comparing to the MST method by adding more terms instead of decreasing the starting radius. I am not sure why.*)
If[Mod[i,30]==0,
deltarp=deltarp/2;
rin=rp+deltarp;
phor=p[rin,-1,s,m,a,\[Omega],\[Lambda]];
qhor=q[rin,-1,s,m,a,\[Omega],\[Lambda]];
];
        a2n[i]=A2n[i];
\[Psi]hor=Evaluate[1+Sum[a2n[k](#-rp)^k,{k,i}]]&;
err=Abs[\[Psi]hor''[rin]+phor \[Psi]hor'[rin]+qhor \[Psi]hor[rin]];
       i++;
];
cInHor=Table[a2n[k],{k,0,i-1}]; (*Coefficients for ingoing waves near horizon (-)*)
{cInHor,rin}
]


bcinf[workingprecision_,s_,m_,a_,\[Omega]_,\[Lambda]_]:=Module[{M=1,B1n,rp,rm,rout,err,i,p,q,pinf,qinf,dpinf,dqinf,\[Psi]inf,b1n,cOutinf},
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
rout=10\[Pi](1/Abs[\[Omega]]+Abs[\[Omega]]/(1+Abs[\[Omega]]));

dqinf[n_]:=Which[
n==0,
0,
n==1,
0,
n==2,
-(4 a m \[Omega]+4 I M s \[Omega]+\[Lambda]),
n>2,
1/(rm-rp)^3 ((rm-rp)^2 (2 rm^(-2+n) rp-2 rm rp^(-2+n)+2 I a m (-rm^(-2+n)+rp^(-2+n))+2 M (-rm^(-2+n)+rp^(-2+n)) (1+s)+(-rm^(-1+n)+rp^(-1+n)) \[Lambda])+2 I (rm-rp)^2 (-rm^(-1+n) rp+rm rp^(-1+n)) \[Omega]+4 (a m (-rp^n (2M+n rm-n rp)+rm^n (2M-n rm+n rp)+a^2 (rp^(-2+n) (rm-n rm+(-3+n) rp)+rm^(-2+n) (-(-3+n) rm+(-1+n) rp)))+I M (-rp^n (2M+n rm-n rp)+rm^n (2M-n rm+n rp)+a^2 (rp^(-2+n) ((-1+n) rm-(-3+n) rp)+rm^(-2+n) ((-3+n) rm+rp-n rp))) s) \[Omega])
];

dpinf[n_]:=Which[
n==0,
2I \[Omega],
n==1,
-2 s+2 I 2M \[Omega] ,
n>1,
rm^(-1+n)+rp^(-1+n)+(2 rm^(-1+n) ((M-rm) s+I (a m+(a^2+rm^2) \[Omega])))/(rm-rp)-(2 rp^(-1+n) ((M-rp) s+I(a m+(a^2+rp^2) \[Omega])))/(rm-rp)
];

err=1;
i=1;
b1n[0]=1;
B1n[n_]:=(n-1)/(2 I \[Omega] ) b1n[n-1]+1/(2 I \[Omega] n) Sum[(dqinf[j+1]-(n-j)dpinf[j])b1n[n-j],{j,1,n}];

{p,q}=TeukolskyHSCoefficients;
pinf=p[rout,1,s,m,a,\[Omega],\[Lambda]];
qinf=q[rout,1,s,m,a,\[Omega],\[Lambda]];

While[err > 10^(-workingprecision),
(*Apparently, it works better when comparing to the MST method by adding more terms instead of increasing the starting radius. Probably because the longest is the integration interval, the greater are the errors.*)
(*If[Mod[i,50]\[Equal]0,
rout=2*rout;
pinf=p[rout,1,s,m,a,\[Omega],\[Lambda]];
qinf=q[rout,1,s,m,a,\[Omega],\[Lambda]];
];
*)
        b1n[i]=B1n[i];
\[Psi]inf=Evaluate[1+Sum[b1n[k](#)^-k,{k,i}]]&;
err=Abs[\[Psi]inf''[rout]+pinf \[Psi]inf'[rout]+qinf \[Psi]inf[rout]];
        i++;
If[i > 100, Break[]]  (*Asymptotic expansions are not convergent*)    
];
cOutinf=Table[b1n[k],{k,0,i-1}]; (*Coefficients for outgoing waves at \[Infinity] (+)*)
{cOutinf,rout}
]


TeukolskyHS[r0_,s_,l_,m_,a_,\[Omega]_,\[Lambda]_]:=Module[{M=1,prec,rin,rout,rp,rm,\[CapitalDelta],p,q,cInH,cOutinf,\[Psi]hor,\[Psi]inf,eqhor,eqinf,r,rtor,X,Y,\[Psi]in,d\[Psi]in,\[Psi]up,d\[Psi]up,Rin,dRin,Rup,dRup,resfac,nmaxhor,nmaxinf,dfacexp},

rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
\[CapitalDelta]=r0^2+a^2-2 M r0;
rtor=((2M rp)/(rp-rm) Log[(r0-rp)/(2M)]-(2M rm)/(rp-rm) Log[(r0-rm)/(2M)]+r0);

prec=Min[Precision[\[Omega]],precODE];

{p,q}=TeukolskyHSCoefficients;

{cInH,rin}=bchor[prec,s,m,a,\[Omega],\[Lambda]];
{cOutinf,rout}=bcinf[prec,s,m,a,\[Omega],\[Lambda]];
nmaxhor=Length[cInH];
nmaxinf=Length[cOutinf];

\[Psi]hor=Evaluate[Sum[cInH[[i]](#-rp)^(i-1),{i,nmaxhor}]]&;
\[Psi]inf[r_]:=Sum[cOutinf[[i]]r^(-i+1),{i,nmaxinf}];

eqhor={
X'[r]== Y[r],
Y'[r]==-p[r,-1,s,m,a,\[Omega],\[Lambda]]Y[r]-q[r,-1,s,m,a,\[Omega],\[Lambda]] X[r],
X[rin]==\[Psi]hor[rin],Y[rin]== \[Psi]hor'[rin]
};

eqinf={
X'[r]== Y[r],
Y'[r]==-p[r,1,s,m,a,\[Omega],\[Lambda]]Y[r]-q[r,1,s,m,a,\[Omega],\[Lambda]] X[r],
X[rout]==\[Psi]inf[rout],Y[rout]==\[Psi]inf'[rout]
};  

{\[Psi]in,d\[Psi]in}={X[r0],Y[r0]}/.First@NDSolve[eqhor,{X,Y},{r,rin,r0},Method->"StiffnessSwitching",WorkingPrecision->prec];
{\[Psi]up,d\[Psi]up}={X[r0],Y[r0]}/.First@NDSolve[eqinf,{X,Y},{r,r0,rout},Method->"StiffnessSwitching",WorkingPrecision->prec];

dfacexp[H_]:=-(1/r0 +(2 s (r0-M))/\[CapitalDelta])+I/\[CapitalDelta] (H (r0^2+a^2)\[Omega]+a m);

Rin=rp Exp[I a m(1/(2M)+1/rp Log[(rp-rm)/(2M)])] \[Psi]in ;
dRin=rp Exp[I a m(1/(2M)+1/rp Log[(rp-rm)/(2M)])](\[Psi]in dfacexp[-1]+ d\[Psi]in);

Rup= \[Psi]up ;
dRup=\[Psi]up dfacexp[1]+d\[Psi]up;

resfac[H_]:=r0^-1 \[CapitalDelta]^-s Exp[H I \[Omega] rtor]Exp[ I m a /(rp-rm) (Log[(r0-rp)/(r0-rm)])];

(* Remove local variables not garbage collected*)
Remove[X,Y,r];
ClearSystemCache[];
{resfac[-1]{Rin,dRin},resfac[1]{Rup,dRup}}
]

TeukolskyHS2[s_,l_,m_,a_,\[Omega]_,\[Lambda]_,ecc_,semilatus_]:=Module[{M=1,rp,prec,rin,r0in,r0up,rout,p,q,cInH,cOutinf,\[Psi]hor,\[Psi]inf,eqhor,eqinf,r,rtor,X,Y,\[Psi]in,d\[Psi]in,\[Psi]up,d\[Psi]up,nmaxhor,nmaxinf,dfacexp},

prec=Min[Precision[\[Omega]],precODE];
rp=M+Sqrt[M^2-a^2];
{p,q}=TeukolskyHSCoefficients;
r0in=erre[ecc,semilatus,-\[Pi]];
r0up=erre[ecc,semilatus,0];

{cInH,rin}=bchor[prec,s,m,a,\[Omega],\[Lambda]];
{cOutinf,rout}=bcinf[prec,s,m,a,\[Omega],\[Lambda]];
nmaxhor=Length[cInH];
nmaxinf=Length[cOutinf];

\[Psi]hor=Evaluate[Sum[cInH[[i]](#-rp)^(i-1),{i,nmaxhor}]]&;
\[Psi]inf[r_]:=Sum[cOutinf[[i]]r^(-i+1),{i,nmaxinf}];

eqhor={
X''[r]== -p[r,-1,s,m,a,\[Omega],\[Lambda]]X'[r]-q[r,-1,s,m,a,\[Omega],\[Lambda]] X[r],
X[rin]==\[Psi]hor[rin],X'[rin]== \[Psi]hor'[rin]
};

eqinf={
X''[r]==-p[r,1,s,m,a,\[Omega],\[Lambda]]X'[r]-q[r,1,s,m,a,\[Omega],\[Lambda]] X[r],
X[rout]==\[Psi]inf[rout],X'[rout]==\[Psi]inf'[rout]
};  

\[Psi]in=X/.First@NDSolve[eqhor,{X,Y},{r,rin,r0in},Method->"StiffnessSwitching",WorkingPrecision->prec,InterpolationOrder->All];
\[Psi]up=X/.First@NDSolve[eqinf,{X,Y},{r,r0up,rout},Method->"StiffnessSwitching",WorkingPrecision->prec,InterpolationOrder->All];

(* Remove local variables not garbage collected*)
Remove[X,Y,r];
ClearSystemCache[];
{\[Psi]in,\[Psi]up}
]



(* ::Subsection::Closed:: *)
(*Source term*)


(*Source term for equatorial eccentric orbits*)
source[p_,e_,aa_,l_,m_]:= Block[{S,a,sorg,\[Theta],mp=1,d=1,M=1,\[Chi],\[CapitalOmega]r,\[CapitalOmega]\[Phi],\[Omega],r,vr, v\[Phi],vt,j,\[CapitalSigma],dt,\[CapitalDelta]},

If[aa<0, 
a=-aa;
, 
a=aa;
];

\[CapitalSigma] = (r^2+a^2 (Cos[\[Theta]])^2)//.{\[Theta]-> \[Pi]/2};
\[CapitalDelta][r_] = r^2+a^2- 2 M r;
dt =  en[aa,p,e]/\[CapitalSigma] ((r^2+a^2)^2/\[CapitalDelta][r]-a^2 (Sin[\[Theta]])^2) +( a L[aa,p,e])/\[CapitalSigma] (1-(r^2+a^2)/\[CapitalDelta][r]) //.{\[Theta]-> \[Pi]/2} ; 

vr = Vr[aa,p,e,\[Chi]];
v\[Phi] = V\[Phi][aa,p,e,\[Chi]];
vt = Vt[aa,p,e,\[Chi]];
j = J[aa,p,e,\[Chi]];

sorg=- 4 \[Pi] mp  d  \[CapitalOmega]r /(2 \[Pi] (a^2+r^2)^(1/2))  Conjugate[S]/dt  vt/(j (vr)^(1/2))  * (Exp[I (\[Omega] t[aa,p,e,\[Chi]]-m \[Phi][aa,p,e,\[Chi]])]+ Exp[-I (\[Omega] t[aa,p,e,\[Chi]]-m \[Phi][aa,p,e,\[Chi]])]) ;

Return[{sorg,dt}];
]


(* ::Subsection::Closed:: *)
(*Integrator*)


(*Integration over the source term. It returns the fluxes and the amplitudes at the horizon and at infinity*)
Int[p_,e_,l_,m_,k_, aa_]:=Block[{s=0,a,drt,sorg, YH,dYH,dY\[Infinity],Y\[Infinity],mm,\[CapitalOmega]r,ssh,ss\[Infinity],\[CapitalOmega]\[Phi],kk,\[Omega]mk,dfacexp,resfac,rtor,rp,rm,\[CapitalOmega],solTeu,ressource,dEd\[Sigma],Aout,r,M=1,qSh,W,RH,Ain,Bin,R\[Infinity],dRH,dR\[Infinity],Dtran,eqShp,eqShm,eqSh,eqS\[Infinity],eqS\[Infinity]m,eqS\[Infinity]p,f,\[CapitalDelta],df,d2f,\[Chi],ip,\[Omega],pi,y,vr,vt,v\[Phi],j,Aout\[Infinity],Aouth,dEd\[Sigma]\[Infinity],dEd\[Sigma]h,dlm,alpha,eps,ci,im,dLd\[Sigma]\[Infinity],dLd\[Sigma]h,\[Lambda],S,dS,V},

If[aa<0,
a=-aa;
If[e==0, 
\[CapitalOmega] = KerrGeoFrequencies[aa,p,e,1];
\[CapitalOmega]\[Phi] = -\[CapitalOmega][[3]];
,
\[CapitalOmega] = KerrGeoFrequencies[a,p,e,-1];
\[CapitalOmega]\[Phi] = \[CapitalOmega][[3]];
];
, 
a=aa;
\[CapitalOmega] = KerrGeoFrequencies[a,p,e,1];
\[CapitalOmega]\[Phi] = \[CapitalOmega][[3]];
];


\[CapitalOmega]r=\[CapitalOmega][[1]];
\[Omega] = m \[CapitalOmega]\[Phi] + k \[CapitalOmega][[1]]; 
If[a==0 & e==0, \[Omega] = SetPrecision[\[Omega],300]];

pi= \[Omega]- (m a)/(2 M (M+ Sqrt[M^2-a^2]));

rtor = ((2M rp)/(rp-rm) Log[(r-rp)/(2M)]-(2M rm)/(rp-rm) Log[(r-rm)/(2M)]+r);

\[CapitalDelta] = r^2-2M r+a^2;

rp = M+Sqrt[M^2-a^2];
rm = M-Sqrt[M^2-a^2];


\[Lambda] = SpinWeightedSpheroidalEigenvalue[s,l,m,a \[Omega]];
S = SpinWeightedSpheroidalHarmonicS[s,l,m,a \[Omega]][\[Pi]/2,0];

(*een if precision of Teuk is Machine Precision, when you evaluate it in a r`prec it becomes prec (a bit less)*)
dfacexp[H_]:= -(1/r +(2 s (r-M))/\[CapitalDelta])+I/\[CapitalDelta] (H (r^2+a^2)\[Omega]+a m);
resfac[H_]:= r^-1 \[CapitalDelta]^-s Exp[H I \[Omega] rtor]Exp[ I m a /(rp-rm) (Log[(r-rp)/(r-rm)])];

solTeu = TeukolskyHS2[s,l,m,a,\[Omega],\[Lambda],e,p];

(*RH is TeukolskyRadial like the toolkit one*)
RH[r_]= rp Exp[I a m(1/(2M)+1/rp Log[(rp-rm)/(2M)])]resfac[-1]solTeu[[1]][r];


dRH[r_]= rp Exp[I a m(1/(2M)+1/rp Log[(rp-rm)/(2M)])]resfac[-1](solTeu[[1]][r] dfacexp[-1]+solTeu[[1]]'[r]);

R\[Infinity][r_]= resfac[1]solTeu[[2]][r];


dR\[Infinity][r_]= resfac[1](solTeu[[2]][r] dfacexp[1]+solTeu[[2]]'[r]);


YH[r_] = RH[r]*Sqrt[r^2+a^2]/Sqrt[rp^2+a^2];
Y\[Infinity][r_] = R\[Infinity][r]*Sqrt[r^2+a^2];

dYH[r_] = (dRH[r]*Sqrt[r^2+a^2]+RH[r]*D[Sqrt[r^2+a^2],r])/Sqrt[rp^2+a^2];
dY\[Infinity][r_] = (dR\[Infinity][r]*Sqrt[r^2+a^2]+R\[Infinity][r]*D[Sqrt[r^2+a^2],r]);

drt=(r^2+a^2)/\[CapitalDelta];

W=(YH[r]dY\[Infinity][r]-Y\[Infinity][r]dYH[r])/drt//.r->p;

If[e==0, 
If[k!=0, Print["k is not zero"],

sorg =  source[p,e,aa,l,m][[2]];

Aouth = (W^-1 Y\[Infinity][p] drt \[CapitalDelta]/(a^2+r^2)^(3/2)  4 \[Pi]   Conjugate[S]/sorg)//.r->p;

Aout\[Infinity] =  (W^-1 YH[p] drt \[CapitalDelta]/(a^2+r^2)^(3/2)  4 \[Pi]  Conjugate[S]/sorg)//.r->p; ;

dEd\[Sigma]\[Infinity]=1/(4\[Pi]) \[Omega]^2 Abs[Aout\[Infinity]]^2;
dEd\[Sigma]h=1/(4\[Pi]) \[Omega] pi Abs[Aouth]^2;

dLd\[Sigma]\[Infinity] = (m/\[Omega])*dEd\[Sigma]\[Infinity]; 
dLd\[Sigma]h = (m/\[Omega])*dEd\[Sigma]h;

]
,

sorg =  source[p,e,aa,l,m][[1]]  //.r-> erre[e,p,\[Chi]];

eqSh= (Y\[Infinity][r] ) //.r-> erre[e,p,\[Chi]];
ssh[\[Chi]_]= eqSh * sorg;

Aouth = W^-1 NIntegrate[ssh[\[Chi]],{\[Chi],0,\[Pi]},MaxRecursion->25,AccuracyGoal->10,WorkingPrecision->Precision[eqSh[0]]-4];

dEd\[Sigma]h=1/(4\[Pi]) \[Omega] pi  Abs[Aouth]^2;

eqS\[Infinity]=(YH[r] )//.r-> erre[e,p,\[Chi]];
ss\[Infinity][\[Chi]_]=eqS\[Infinity] * sorg ;

Aout\[Infinity] = W^-1 NIntegrate[ss\[Infinity][\[Chi]],{\[Chi],0,\[Pi]},MaxRecursion->25,AccuracyGoal->10,WorkingPrecision->Precision[eqS\[Infinity][0]]-4];


dEd\[Sigma]\[Infinity]=1/(4\[Pi]) \[Omega]^2 Abs[Aout\[Infinity]]^2;

dLd\[Sigma]\[Infinity] = (m/\[Omega])*dEd\[Sigma]\[Infinity]; 
dLd\[Sigma]h = (m/\[Omega])*dEd\[Sigma]h; 

];


Return[{dEd\[Sigma]h,dEd\[Sigma]\[Infinity],dLd\[Sigma]h,dLd\[Sigma]\[Infinity],Aouth,Aout\[Infinity]}]
];


ppp=50;
precODE=50;

Print["get"];
griglia = Import["XXX"];


(* ::Subsection::Closed:: *)
(*Sum over the modes*)


(* ::Subsubsection::Closed:: *)
(*prograde orbits a > 0*)


(*Function to sum over the different modes, for prograde orbits (a>0). Another function is implemented for retrograde orbits. Covergence criteria follows Viktor Skoupy criteria.*)

sum[p_,e_,aa_]:=Block[{
arrayklm00,arraykm0,arrayk0,arrayklm11,arraykm1,arrayk1,arrayklm22,arraykm2,arrayk2,arrayk,
listlm00,listlm0,listm0,listlm11,listm1,liste0m1,listlm22,listm2,liste0m,liste0m2,liste0,list,
sumlm00h,sumlm00,sumlm11h,sumlm11,sumlm22,sumlm22h,suml,
fluxe0lm11,fluxe0lm22,fluxl,
matrixe0,matrix1,matrix,
ind,\[Epsilon],a,\[CapitalOmega],\[CapitalOmega]r,\[CapitalOmega]\[Phi],l,m,k,kfin,flux,i,shift,count,
lmax=120,kmin,kmax=130,mmax=25,mmin=-15,mthr=6},

(*computing the fundamental (l,m,k)=(2,2,0) mode*)
\[Epsilon]=10^-7; 
l=2;
m=2;
i=0;

arrayklm22=ConstantArray[{0,0,0,0,0,0,0},kmax];
arraykm2=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk2=ConstantArray[{0,0,0,0,0,0,0},kmax];

arrayklm00=ConstantArray[{0,0,0,0,0,0,0},kmax];
arraykm0=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk0=ConstantArray[{0,0,0,0,0,0,0},kmax];

arrayklm11=ConstantArray[{0,0,0,0,0,0,0},kmax];
arraykm1=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk1=ConstantArray[{0,0,0,0,0,0,0},kmax];

If[aa<0,
a=-aa;
If[e==0, 
\[CapitalOmega] = KerrGeoFrequencies[aa,p,e,1];
\[CapitalOmega]\[Phi] =  - \[CapitalOmega][[3]];
,
\[CapitalOmega] = KerrGeoFrequencies[a,p,e,-1];
\[CapitalOmega]\[Phi] =  \[CapitalOmega][[3]];
];
, 
a=aa;
\[CapitalOmega] = KerrGeoFrequencies[a,p,e,1];
\[CapitalOmega]\[Phi] = \[CapitalOmega][[3]];
];

\[CapitalOmega]r=\[CapitalOmega][[1]];

(* --- circular case first --- *)
If[e==0, 
Print["e=0"];
\[Epsilon]=10^-15;
k=0;
liste0m1 = ConstantArray[{0,0,0,0,0,0,0,0,0},lmax];
liste0m2 = ConstantArray[{0,0,0,0,0,0,0,0,0},lmax];
fluxe0lm11 = Int[p,e,1,1,k, aa];
fluxe0lm22 = Int[p,e,2,2,k, aa];
liste0m1[[1]]=Join[{1,1,0},fluxe0lm11] ;
liste0m2[[1]]=Join[{2,2,0},fluxe0lm22] ;
Clear[m,l];

count=0;
m=1;
For[l=2, l<=100, l++, 
If[EvenQ[l+Abs[m]],
count=count+1;
Print["(l,m,k)=",{l,m,k}];
flux=Int[p,e,l,m,k, aa];
liste0m1[[count+1]]=Join[{l,1,0},flux]; 
Print["flux=",liste0m1[[count+1]][[5]]+liste0m1[[count+1]][[4]]];
Print["check=",fluxe0lm11[[2]]+fluxe0lm11[[1]]];

If[liste0m1[[count+1]][[5]]+liste0m1[[count+1]][[4]] < \[Epsilon]*(fluxe0lm11[[2]]+fluxe0lm11[[1]]), Break[];];
];
];
liste0m1 = Select[liste0m1, #!={0,0,0,0,0,0,0,0,0}&];

count=0;
m=2;
For[l=3, l<=100, l++, 
If[EvenQ[l+Abs[m]],
count=count+1;
Print["(l,m,k)=",{l,m,k}];
flux=Int[p,e,l,m,k, aa];
liste0m2[[count+1]]=Join[{l,2,0},flux]; 
Print["flux=",liste0m2[[count+1]][[5]]+liste0m2[[count+1]][[4]]];
Print["check=",fluxe0lm22[[2]]+fluxe0lm22[[1]]];

If[liste0m2[[count+1]][[5]]+liste0m2[[count+1]][[4]] < \[Epsilon]*Min[(fluxe0lm22[[2]]+fluxe0lm22[[1]]),(fluxe0lm11[[2]]+fluxe0lm11[[1]])], Break[];];
];
];
liste0m2 = Select[liste0m2, #!={0,0,0,0,0,0,0,0,0}&];
count=0;

liste0m = ConstantArray[{0,0,0,0,0,0,0,0,0},lmax];
matrixe0 = ConstantArray[{0,0,0,0,0,0,0,0,0},lmax];
matrixe0[[1]] = Select[Partition[Flatten[liste0m1],9], #!={0,0,0,0,0,0,0,0,0}&];
matrixe0[[2]] = Select[Partition[Flatten[liste0m2],9], #!={0,0,0,0,0,0,0,0,0}&];

For[m=3,m<=30,m++,
For[l=Max[0,Abs[m]],l<=lmax,l++,
If[EvenQ[l+Abs[m]],
count=count+1;
Print["(l,m,k)=",{l,m,k}];
flux=Int[p,e,l,m,k, aa];
liste0m[[count]]=Join[{l,m,0},flux]; 
If[liste0m[[count]][[5]] < \[Epsilon]*Min[fluxe0lm22[[2]],fluxe0lm11[[2]]] , Break[];];
];
];

matrixe0[[m]] = Select[Partition[Flatten[liste0m],9], #!={0,0,0,0,0,0,0,0,0}&]; 
liste0m = ConstantArray[{0,0,0,0,0,0,0,0,0},lmax];
count=0;

If[m>mthr,
If[(Total[matrixe0[[m]][[All,5]]]+Total[matrixe0[[m]][[All,4]]])/(1-((Total[matrixe0[[m]][[All,5]]]+Total[matrixe0[[m]][[All,4]]])/(Total[matrixe0[[m-1]][[All,5]]]+Total[matrixe0[[m-1]][[All,4]]])))< \[Epsilon]/10 (Total[Partition[Flatten[matrixe0],9][[All,5]]]+Total[Partition[Flatten[matrixe0],9][[All,4]]]), Break[];
];
];
];

matrixe0= Partition[Flatten[matrixe0],9];
matrixe0 = Select[matrixe0, #!={0,0,0,0,0,0,0,0,0}&];

Return[matrixe0];

,

(****** eccentric case : e != 0  ***************************************************)

(******** m=0 **********************************************************************)
l=0;
m=0;
ind=0;
For[k=Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]],k<=kmax,k++,
If[k!=0, 
l=0;
m=0;
ind=ind+1;
flux = Int[p,e,l,m,k,aa];
arrayklm00[[ind]]=Flatten[{k,flux}];
kfin=k;
Print["(l,m,k)=",{l,m,k}];
If[k>17,
If[flux[[2]] <  \[Epsilon]/10 Max[arrayklm00[[All,3]]] , 
i= i+1; 
If[i==3, 
i=0;
Break[];
],
i=0;
]; 
];
];
];

sumlm00= Total[arrayklm00[[All,3]]]; 
sumlm00h= Total[arrayklm00[[All,2]]]; 
Print["sumlm00=",sumlm00];
Clear[l];
listlm00 = ConstantArray[0,kmax*lmax*10];
listm0 = ConstantArray[0,kmax*lmax*10];
list = ConstantArray[0,kmax*lmax*10];
listlm00[[1]]=Table[{0,0,arrayklm00[[i,1]],arrayklm00[[i,2]],arrayklm00[[i,3]],arrayklm00[[i,4]],arrayklm00[[i,5]],arrayklm00[[i,6]],arrayklm00[[i,7]]},{i,1,kfin-Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]]}];
Print["here"];
ind=0;
count=0;

For[l=1,l<=lmax,l++,
If[EvenQ[l+Abs[m]],
kmin= Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]];
count=count+1;
For[k=kmin,k<=kmax,k++,
If[k!=0,
Print["(l,m,k)=",{l,m,k}];
ind=ind+1;
flux=Int[p,e,l,m,k, aa];
arraykm0[[ind]]=Flatten[{k,flux}];
kfin=k;
Print["max=",Max[arrayklm00[[All,3]]+arrayklm00[[All,2]]]];
Print["flux[[2]]+flux[[1]]=",flux[[2]]+flux[[1]]];

If[flux[[2]]+flux[[1]] < \[Epsilon]/10 Max[arrayklm00[[All,3]]+arrayklm00[[All,2]]] &&  ind>1 && (flux[[2]]+flux[[1]]<arraykm0[[k-kmin,3]]+arraykm0[[k-kmin,2]] || flux[[2]]+flux[[1]] < \[Epsilon]/10000000 Max[arrayklm00[[All,3]]+arrayklm00[[All,2]]] ) , 
i= i+1; 
Print["i=",i];
If[i==3, 
i=0;
Break[];
],
i=0;
Print["i=",i];
]; 
];
];
Print["ind=", ind];
Print["kfin=", kfin];
listm0[[count]]=Table[{l,0,arraykm0[[i,1]],arraykm0[[i,2]],arraykm0[[i,3]],arraykm0[[i,4]],arraykm0[[i,5]],arraykm0[[i,6]],arraykm0[[i,7]]},{i,1,ind}]; 
Print["count=", count];
Print["listm0[[count]]=", listm0[[count]]];
ind=0;
arraykm0=ConstantArray[{0,0,0,0,0,0,0},kmax];

If[Total[listm0[[count]][[All,5]]]+Total[listm0[[count]][[All,4]]] < \[Epsilon]/10 (sumlm00+sumlm00h) , Break[];]
(*If[Total[listm2[[l-2]][[All,5]]] < \[Epsilon] (sumlm22) , Break[];];*)
];
];

Clear[m];
matrix= ConstantArray[0,lmax*kmax*20];

matrix[[6]]=Select[Partition[Flatten[Join[listlm00[[1]],listm0]],9], #!={0,0,0,0,0,0,0,0,0}&];
Print["matrix[[6]]=",matrix[[6]]];

ind=0;
count=0;


(** ---------------------------------- m=1 ---------------------------------- **)
l=1;
m=1;
For[k=Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]],k<=kmax,k++,
ind=ind+1;
flux = Int[p,e,l,m,k, aa];
arrayklm11[[ind]]=Flatten[{k,flux}];
kfin=k;
Print["(l,m,k)=",{l,m,k}];
If[k>17,
If[flux[[2]] <  \[Epsilon]/10 Max[arrayklm11[[All,3]]] , 
i= i+1; 
If[i==3, 
i=0;
Break[];
],
i=0;
]; 
];
];

sumlm11= Total[arrayklm11[[All,3]]]; 
sumlm11h= Total[arrayklm11[[All,2]]]; 
Print["sumlm11=",sumlm11];
Clear[l];
listlm11 = ConstantArray[0,kmax*lmax*10];
listm1 = ConstantArray[0,kmax*lmax*10];
list = ConstantArray[0,kmax*lmax*10];
listlm11[[1]]=Table[{1,1,arrayklm11[[i,1]],arrayklm11[[i,2]],arrayklm11[[i,3]],arrayklm11[[i,4]],arrayklm11[[i,5]],arrayklm11[[i,6]],arrayklm11[[i,7]]},{i,1,ind}];

ind=0;
count=0;

For[l=2,l<=lmax,l++,
If[EvenQ[l+Abs[m]],
count=count+1;
kmin= Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]];
For[k=kmin,k<=kmax,k++,
ind=ind+1;
flux = Int[p,e,l,m,k, aa];
arraykm1[[ind]]=Flatten[{k,flux}];
kfin=k;
Print["(l,m,k)=",{l,m,k}];

Print["max=",Max[arrayklm11[[All,3]]+arrayklm11[[All,2]]]];
Print["flux[[2]]+flux[[1]]=",flux[[2]]+flux[[1]]];

If[flux[[2]]+flux[[1]] < \[Epsilon]/10 Max[arrayklm11[[All,3]]+arrayklm11[[All,2]]] &&  ind>1 && (flux[[2]]+flux[[1]]<arraykm1[[k-kmin,3]]+arraykm1[[k-kmin,2]] || flux[[2]]+flux[[1]] < \[Epsilon]*10^-40 Max[arrayklm11[[All,3]]+arrayklm11[[All,2]]] ) , 
i= i+1; 
Print["i=",i];
If[i==3, 
i=0;
Break[];
],
i=0;
Print["i=",i];
]; 
];
listm1[[count]]=Table[{l,1,arraykm1[[i,1]],arraykm1[[i,2]],arraykm1[[i,3]],arraykm1[[i,4]],arraykm1[[i,5]],arraykm1[[i,6]],arraykm1[[i,7]]},{i,1,ind}]; 

ind=0;
arraykm1=ConstantArray[{0,0,0,0,0,0,0},kmax];

If[Total[listm1[[count]][[All,5]]]+Total[listm1[[count]][[All,4]]] < \[Epsilon]/10 (sumlm11+sumlm11h) , Break[];]
(*If[Total[listm2[[l-2]][[All,5]]]<\[Epsilon](sumlm22),Break[];];*)
];
];

Clear[m];
(*matrix= ConstantArray[0,lmax*kmax*20];*)


matrix[[7]]=Select[Partition[Flatten[Join[listlm11[[1]],listm1]],9], #!={0,0,0,0,0,0,0,0,0}&];
Print["matrix[[7]]=",matrix[[7]]];


ind=0;
count=0;


(** ---------------------------------- m=2 ---------------------------------- **)
l=2;
m=2;
For[k=Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]],k<=kmax,k++,
flux = Int[p,e,l,m,k, aa];
arrayklm22[[k-Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]]+1]]=Flatten[{k,flux}];
kfin=k;
Print["(l,m,k)=",{l,m,k}];
If[k>17,
If[flux[[2]] <  \[Epsilon]/10 Max[arrayklm22[[All,3]]] , 
i= i+1; 
If[i==3, 
i=0;
Break[];
],
i=0;
]; 
];
];

sumlm22= Total[arrayklm22[[All,3]]]; 
sumlm22h= Total[arrayklm22[[All,2]]]; 
Print["sumlm22=",sumlm22];
Clear[l];
listlm22 = ConstantArray[0,kmax*lmax*10];
listm2 = ConstantArray[0,kmax*lmax*10];
list = ConstantArray[0,kmax*lmax*10];
listlm22[[1]]=Table[{2,2,arrayklm22[[i,1]],arrayklm22[[i,2]],arrayklm22[[i,3]],arrayklm22[[i,4]],arrayklm22[[i,5]],arrayklm22[[i,6]],arrayklm22[[i,7]]},{i,1,kfin+1-Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]]}];

ind=0;
count=0;

For[l=3,l<=lmax,l++,
If[EvenQ[l+Abs[m]],
count=count+1;
kmin= Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]];
For[k=kmin,k<=kmax,k++,
ind=ind+1;
flux = Int[p,e,l,m,k, aa];
arraykm2[[ind]]=Flatten[{k,flux}];
kfin=k;
Print["(l,m,k)=",{l,m,k}];

Print["max=",Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]]];
Print["flux[[2]]+flux[[1]]=",flux[[2]]+flux[[1]]];

If[flux[[2]]+flux[[1]] < \[Epsilon]/10 Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]] &&  ind>1 && (flux[[2]]+flux[[1]]<arraykm2[[k-kmin,3]]+arraykm2[[k-kmin,2]] || flux[[2]]+flux[[1]] < \[Epsilon]*10^-40 Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]] ) , 
i= i+1; 
Print["i=",i];
If[i==3, 
i=0;
Break[];
],
i=0;
Print["i=",i];
]; 
];
listm2[[count]]=Table[{l,2,arraykm2[[i,1]],arraykm2[[i,2]],arraykm2[[i,3]],arraykm2[[i,4]],arraykm2[[i,5]],arraykm2[[i,6]],arraykm2[[i,7]]},{i,1,ind}]; 

ind=0;
arraykm2=ConstantArray[{0,0,0,0,0,0,0},kmax];

If[Total[listm2[[count]][[All,5]]]+Total[listm2[[count]][[All,4]]] < \[Epsilon]/10 (sumlm22+sumlm22h) , Break[];]
(*If[Total[listm2[[l-2]][[All,5]]] < \[Epsilon] (sumlm22) , Break[];];*)
];
];


Clear[m];
(*matrix= ConstantArray[0,lmax*kmax*20];*)

matrix[[8]]=Select[Partition[Flatten[Join[listlm22[[1]],listm2]],9], #!={0,0,0,0,0,0,0,0,0}&];
Print["matrix[[8]]=",matrix[[8]]];


ind=0;
count=0;

(****** m > - 5 &  m != (0,1,2)  ********************************************************************)

For[m=-5,m<=mmax,m++,
If[m!=2 && m!=1 && m!=0, 

For[l=Max[0,Abs[m]],l<=lmax,l++,
If[EvenQ[l+Abs[m]],
count=count+1;
If[m==0, kmin=1, If[m<0, kmin = Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]], If[aa>=0,kmin = Floor[10 m e^2],kmin = Floor[-10 m e^2] ]]];
(*If[m==0, kmin=1, If[m<0, kmin = Floor[10 Abs[m] e^2], kmin = Ceiling[10 m e^2]]]*)

For[k=kmin,k<=kmax,k++,
ind=ind+1;
Print["(l,m,k)=",{l,m,k}];
flux = Int[p,e,l,m,k, aa];
If[m==0,shift=0,shift=1];
arrayk[[ind]]=Flatten[{k,flux}];
kfin=k;
If[flux[[2]]+flux[[1]] < \[Epsilon]/10 Max[Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]],Max[arrayklm11[[All,3]]+arrayklm11[[All,2]]],Max[arrayklm00[[All,3]]+arrayklm00[[All,2]]]]&&  ind>1 && (flux[[2]]+flux[[1]]<arrayk[[ind-1,3]]+arrayk[[ind-1,2]] || flux[[2]]+flux[[1]] < \[Epsilon]*10^-40 Max[Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]],Max[arrayklm11[[All,3]]+arrayklm11[[All,2]]],Max[arrayklm00[[All,3]]+arrayklm00[[All,2]]] ]) ,  
Print["i=",i];
i= i+1; 
If[i==3, 
i=0;
Break[];
],
i=0;
Print["i=",i];
]; 
];

i=0;
ind=0;

For[k=kmin-1,k>=-100,k--,
If[k<0, If[k<Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]], Break[]], If[k<=Ceiling[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]], Break[]] ];
Print["(l,m,k)=",{l,m,k}];
flux = Int[p,e,l,m,k, aa];
If[m==0,shift=0,shift=1];
ind=ind+1;
arrayk2[[ind]]=Flatten[{k,flux}];

If[flux[[1]]+flux[[2]] < \[Epsilon]/10 Max[Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]], Max[arrayklm11[[All,3]]+arrayklm11[[All,2]]], Max[arrayklm00[[All,3]]+arrayklm00[[All,2]]]] && ind>1 && (flux[[2]]+flux[[1]]<arrayk2[[ind-1,3]]+arrayk2[[ind-1,2]] || flux[[1]]+flux[[2]] < \[Epsilon]*10^-40 Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]]) , 
i= i+1; 
Print["i=",i];
If[i==3, 
i=0;
Break[];
],
i=0;
Print["i=",i];
]; 
];

If[m==0,shift=0,shift=1];
list[[count]]=Join[Table[{l,m,arrayk[[i,1]],arrayk[[i,2]],arrayk[[i,3]],arrayk[[i,4]],arrayk[[i,5]],arrayk[[i,6]],arrayk[[i,7]]},{i,1,kfin-kmin+1}],Table[{l,m,arrayk2[[i,1]],arrayk2[[i,2]],arrayk2[[i,3]],arrayk2[[i,4]],arrayk2[[i,5]],arrayk2[[i,6]],arrayk2[[i,7]]},{i,1,ind}]]; 
Print["list[[count]]=", list[[count]]];

arrayk=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk2=ConstantArray[{0,0,0,0,0,0,0},kmax];

i=0;
ind=0;

If[Total[list[[count]][[All,5]]]+Total[list[[count]][[All,4]]] < \[Epsilon]  (sumlm22+sumlm22h), Break[];];
];
];

count=0;
matrix[[m+6]] = Select[Partition[Flatten[list],9], #!={0,0,0,0,0,0,0,0,0}&]; 


list = ConstantArray[0,kmax*lmax*10];

If[m>mthr,
If[(Total[matrix[[m+6]][[All,5]]]+Total[matrix[[m+6]][[All,4]]])/(1-((Total[matrix[[m+6]][[All,5]]]+Total[matrix[[m+6]][[All,4]]])/(Total[matrix[[m+5]][[All,5]]]+Total[matrix[[m+5]][[All,4]]])))< \[Epsilon]/2 (Total[Partition[Flatten[matrix],9][[All,5]]]+Total[Partition[Flatten[matrix],9][[All,4]]]), Break[];];
];
];
];

matrix= Partition[Flatten[matrix],9];
matrix1 = Select[matrix, #!={0,0,0,0,0,0,0,0,0}&];

Return[matrix1]

];
];


(* ::Subsubsection::Closed:: *)
(*retrograde orbits a < 0*)


sumaneg[p_,e_,aa_]:=Block[{
arrayk1,arraykm1,arrayklm11,arrayk0,arraykm0,arrayklm00,arrayklm22,arrayk2,arraykm2,arrayk,
listlm11,listm1,listlm00,listm0,liste0m1,liste0m,liste0m2,listlm22,listm2,list,
sumlm11,sumlm11h,sumlm00,sumlm00h,sumlm22h,sumlm22,suml,
matrixe0,fluxe0lm11,fluxe0lm22,
matrix,matrix1,
flux,fluxl,
count,ind,\[Epsilon],a,\[CapitalOmega],\[CapitalOmega]r,\[CapitalOmega]\[Phi],l,m,k,kfin,i,shift,
lmax=120,kmin,kmax=120,mmax=25,mmin=-15},

\[Epsilon]=10^-7;
l=2;
m=2;
i=0;

arrayklm22=ConstantArray[{0,0,0,0,0,0,0},kmax];
arraykm2=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk2=ConstantArray[{0,0,0,0,0,0,0},kmax];

arrayklm00=ConstantArray[{0,0,0,0,0,0,0},kmax];
arraykm0=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk0=ConstantArray[{0,0,0,0,0,0,0},kmax];

arrayklm11=ConstantArray[{0,0,0,0,0,0,0},kmax];
arraykm1=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk1=ConstantArray[{0,0,0,0,0,0,0},kmax];

If[aa<0,
a=-aa;
If[e==0, 
\[CapitalOmega] = KerrGeoFrequencies[aa,p,e,1];
\[CapitalOmega]\[Phi] =  - \[CapitalOmega][[3]];
,
\[CapitalOmega] = KerrGeoFrequencies[a,p,e,-1];
\[CapitalOmega]\[Phi] =  \[CapitalOmega][[3]];
];
, 
a=aa;
\[CapitalOmega] = KerrGeoFrequencies[a,p,e,1];
\[CapitalOmega]\[Phi] = \[CapitalOmega][[3]];
];

\[CapitalOmega]r=\[CapitalOmega][[1]];

(****** e=0  ********************************************************************)
If[e==0, 
Print["e=0"];
\[Epsilon]=10^-15;
k=0;
liste0m1 = ConstantArray[{0,0,0,0,0,0,0,0,0},lmax];
liste0m2 = ConstantArray[{0,0,0,0,0,0,0,0,0},lmax];
fluxe0lm11 = Int[p,e,1,1,k, aa];
fluxe0lm22 = Int[p,e,2,2,k, aa];
liste0m1[[1]]=Join[{1,1,0},fluxe0lm11] ;
liste0m2[[1]]=Join[{2,2,0},fluxe0lm22] ;
Clear[m,l];

count=0;
m=1;
For[l=2, l<=100, l++, 
If[EvenQ[l+Abs[m]],
count=count+1;
Print["(l,m,k)=",{l,m,k}];
flux=Int[p,e,l,m,k, aa];
liste0m1[[count+1]]=Join[{l,1,0},flux]; 
Print["flux=",liste0m1[[count+1]][[5]]+liste0m1[[count+1]][[4]]];
Print["check=",fluxe0lm11[[2]]+fluxe0lm11[[1]]];

If[liste0m1[[count+1]][[5]]+liste0m2[[count+1]][[4]] < \[Epsilon]*(fluxe0lm11[[2]]+fluxe0lm11[[1]]), Break[];];
];
];
liste0m1 = Select[liste0m1, #!={0,0,0,0,0,0,0,0,0}&];


count=0;
m=2;
For[l=3, l<=100, l++, 
If[EvenQ[l+Abs[m]],
count=count+1;
Print["(l,m,k)=",{l,m,k}];
flux=Int[p,e,l,m,k, aa];
liste0m2[[count+1]]=Join[{l,2,0},flux]; 
Print["flux=",liste0m2[[count+1]][[5]]+liste0m2[[count+1]][[4]]];
Print["check=",fluxe0lm22[[2]]+fluxe0lm22[[1]]];

If[liste0m2[[count+1]][[5]]+liste0m2[[count+1]][[4]] < \[Epsilon]*(fluxe0lm22[[2]]+fluxe0lm22[[1]]), Break[];];
];
];
liste0m2 = Select[liste0m2, #!={0,0,0,0,0,0,0,0,0}&];
count=0;


liste0m = ConstantArray[{0,0,0,0,0,0,0,0,0},lmax];
matrixe0 = ConstantArray[{0,0,0,0,0,0,0,0,0},lmax];
matrixe0[[1]] = Select[Partition[Flatten[liste0m1],9], #!={0,0,0,0,0,0,0,0,0}&];
matrixe0[[2]] = Select[Partition[Flatten[liste0m2],9], #!={0,0,0,0,0,0,0,0,0}&];

For[m=3,m<=30,m++,
For[l=Max[0,Abs[m]],l<=lmax,l++,
If[EvenQ[l+Abs[m]],
count=count+1;
Print["(l,m,k)=",{l,m,k}];
flux=Int[p,e,l,m,k, aa];
liste0m[[count]]=Join[{l,m,0},flux]; 
If[liste0m[[count]][[5]] < \[Epsilon] fluxe0lm22[[2]] , Break[];];
];
];

matrixe0[[m]] = Select[Partition[Flatten[liste0m],9], #!={0,0,0,0,0,0,0,0,0}&]; 
liste0m = ConstantArray[{0,0,0,0,0,0,0,0,0},lmax];
count=0;

If[m>mthr,
If[(Total[matrixe0[[m]][[All,5]]]+Total[matrixe0[[m]][[All,4]]])/(1-((Total[matrixe0[[m]][[All,5]]]+Total[matrixe0[[m]][[All,4]]])/(Total[matrixe0[[m-1]][[All,5]]]+Total[matrixe0[[m-1]][[All,4]]])))< \[Epsilon]/100 (Total[Partition[Flatten[matrixe0],9][[All,5]]]+Total[Partition[Flatten[matrixe0],9][[All,4]]]), Break[];
];
];
];

matrixe0= Partition[Flatten[matrixe0],9];
matrixe0 = Select[matrixe0, #!={0,0,0,0,0,0,0,0,0}&];

Return[matrixe0];

,
(****** e != 0  ********************************************************************)

(******** m=0 **********************************************************************)
l=0;
m=0;
ind=0;
For[k=Floor[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]],k>=-kmax,k--,
If[k!=0, 
flux = Int[p,e,l,m,k, aa];
ind=ind+1;
arrayklm00[[ind]]=Flatten[{k,flux}];
kfin=k;
Print["(l,m,k)=",{l,m,k}];
If[k<-17,
If[flux[[2]] <  \[Epsilon]/10 Max[arrayklm00[[All,3]]] , 
i= i+1; 
If[i==3, 
i=0;
Break[];
],
i=0;
]; 
];
];
];

sumlm00= Total[arrayklm00[[All,3]]]; 
sumlm00h= Total[arrayklm00[[All,2]]]; 
Print["sumlm00=",sumlm00];
Clear[l];
listlm00 = ConstantArray[0,kmax*lmax*10];
listm0 = ConstantArray[0,kmax*lmax*10];
list = ConstantArray[0,kmax*lmax*10];
listlm00[[1]]=Table[{0,0,arrayklm00[[i,1]],arrayklm00[[i,2]],arrayklm00[[i,3]],arrayklm00[[i,4]],arrayklm00[[i,5]],arrayklm00[[i,6]],arrayklm00[[i,7]]},{i,1,ind}];

ind=0;
count=0;

For[l=1,l<=lmax,l++,
If[EvenQ[l+Abs[m]],
kmin= Floor[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]];
count=count+1;
For[k=kmin,k>=-kmax,k--,
If[k!=0,
Print["(l,m,k)=",{l,m,k}];
ind=ind+1;
flux = Int[p,e,l,m,k, aa];
arraykm0[[ind]]=Flatten[{k,flux}];
kfin=k;
Print["max=",Max[arrayklm00[[All,3]]+arrayklm00[[All,2]]]];
Print["flux[[2]]+flux[[1]]=",flux[[2]]+flux[[1]]];

If[flux[[2]]+flux[[1]] < \[Epsilon]/10 Max[arrayklm00[[All,3]]+arrayklm00[[All,2]]] &&  ind>1 && (flux[[2]]+flux[[1]]<arraykm0[[ind-1,3]]+arraykm0[[ind-1,2]] || flux[[2]]+flux[[1]] < \[Epsilon]*10^-40 Max[arrayklm00[[All,3]]+arrayklm00[[All,2]]]), 
i= i+1; 
Print["i=",i];
If[i==3, 
i=0;
Break[];
],
i=0;
Print["i=",i];
]; 
];
];

listm0[[count]]=Table[{l,0,arraykm0[[i,1]],arraykm0[[i,2]],arraykm0[[i,3]],arraykm0[[i,4]],arraykm0[[i,5]],arraykm0[[i,6]],arraykm0[[i,7]]},{i,1,ind}]; 
Print["count=", count];
Print["listm0[[count]]=", listm0[[count]]];
ind=0;
arraykm0=ConstantArray[{0,0,0,0,0,0,0},kmax];

If[Total[listm0[[count]][[All,5]]]+Total[listm0[[count]][[All,4]]] < \[Epsilon] (sumlm00+sumlm00h) , Break[];];
];
];

Clear[m];
matrix= ConstantArray[0,lmax*kmax*50];

matrix[[6]]=Select[Partition[Flatten[Join[listlm00[[1]],listm0]],9], #!={0,0,0,0,0,0,0,0,0}&];
Print["matrix[[6]]=",matrix[[6]]];

ind=0;
count=0;

(************************************************************ m=1 ******)
l=1;
m=1;
ind=0;
For[k=Floor[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]],k>=-kmax,k--, 
flux = Int[p,e,l,m,k, aa];
ind=ind+1;
arrayklm11[[ind]]=Flatten[{k,flux}];
kfin=k;
Print["(l,m,k)=",{l,m,k}];
If[k<-17,
If[flux[[2]] <  \[Epsilon]/10 Max[arrayklm11[[All,3]]] , 
i= i+1; 
If[i==3, 
i=0;
Break[];
],
i=0;
]; 
];
];

sumlm11= Total[arrayklm11[[All,3]]]; 
sumlm11h= Total[arrayklm11[[All,2]]]; 
Print["sumlm11=",sumlm11];
Clear[l];
listlm11 = ConstantArray[0,kmax*lmax*10];
listm1 = ConstantArray[0,kmax*lmax*10];
list = ConstantArray[0,kmax*lmax*10];
listlm11[[1]]=Table[{1,1,arrayklm11[[i,1]],arrayklm11[[i,2]],arrayklm11[[i,3]],arrayklm11[[i,4]],arrayklm11[[i,5]],arrayklm11[[i,6]],arrayklm11[[i,7]]},{i,1,ind}];

ind=0;
count=0;

For[l=2,l<=lmax,l++,
If[EvenQ[l+Abs[m]],
kmin= Floor[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]];
count=count+1;
For[k=kmin,k>=-kmax,k--,
Print["(l,m,k)=",{l,m,k}];
ind=ind+1;
flux = Int[p,e,l,m,k, aa];
arraykm1[[ind]]=Flatten[{k,flux}];
kfin=k;
Print["max=",Max[arrayklm11[[All,3]]+arrayklm11[[All,2]]]];
Print["flux[[2]]+flux[[1]]=",flux[[2]]+flux[[1]]];

If[flux[[2]]+flux[[1]] < \[Epsilon]/10 Max[arrayklm11[[All,3]]+arrayklm11[[All,2]]] &&  ind>1 && (flux[[2]]+flux[[1]]<arraykm1[[ind-1,3]]+arraykm1[[ind-1,2]] || flux[[2]]+flux[[1]] < \[Epsilon]*10^-40 Max[arrayklm11[[All,3]]+arrayklm11[[All,2]]]), 
i= i+1; 
Print["i=",i];
If[i==3, 
i=0;
Break[];
],
i=0;
Print["i=",i];
]; 
];
listm1[[count]]=Table[{l,1,arraykm1[[i,1]],arraykm1[[i,2]],arraykm1[[i,3]],arraykm1[[i,4]],arraykm1[[i,5]],arraykm1[[i,6]],arraykm1[[i,7]]},{i,1,ind}]; 
Print["count=", count];
Print["listm1[[count]]=", listm1[[count]]];
ind=0;
arraykm1=ConstantArray[{0,0,0,0,0,0,0},kmax];

If[Total[listm1[[count]][[All,5]]]+Total[listm1[[count]][[All,4]]] < \[Epsilon] (sumlm11+sumlm11h) , Break[];];
];
];

Clear[m];

matrix[[7]]=Select[Partition[Flatten[Join[listlm11[[1]],listm1]],9], #!={0,0,0,0,0,0,0,0,0}&];
Print["matrix[[7]]=",matrix[[7]]];

ind=0;
count=0;

(**************** m=2 ***************************************************)
l=2;
m=2;
ind=0;
For[k=Floor[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]],k>=-kmax,k--, 
flux = Int[p,e,l,m,k, aa];
ind=ind+1;
arrayklm22[[ind]]=Flatten[{k,flux}];
kfin=k;
Print["(l,m,k)=",{l,m,k}];
If[k<-17,
If[flux[[2]] <  \[Epsilon]/10 Max[arrayklm22[[All,3]]] , 
i= i+1; 
If[i==3, 
i=0;
Break[];
],
i=0;
]; 
];
];

sumlm22= Total[arrayklm22[[All,3]]]; 
sumlm22h= Total[arrayklm22[[All,2]]]; 
Print["sumlm22=",sumlm22];
Clear[l];
listlm22 = ConstantArray[0,kmax*lmax*10];
listm2 = ConstantArray[0,kmax*lmax*10];
list = ConstantArray[0,kmax*lmax*10];
listlm22[[1]]=Table[{2,2,arrayklm22[[i,1]],arrayklm22[[i,2]],arrayklm22[[i,3]],arrayklm22[[i,4]],arrayklm22[[i,5]],arrayklm22[[i,6]],arrayklm22[[i,7]]},{i,1,ind}];

ind=0;
count=0;

For[l=3,l<=lmax,l++,
If[EvenQ[l+Abs[m]],
kmin= Floor[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]];
count=count+1;
For[k=kmin,k>=-kmax,k--,
Print["(l,m,k)=",{l,m,k}];
ind=ind+1;
flux = Int[p,e,l,m,k, aa];
arraykm2[[ind]]=Flatten[{k,flux}];
kfin=k;
Print["max=",Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]]];
Print["flux[[2]]+flux[[1]]=",flux[[2]]+flux[[1]]];

If[flux[[2]]+flux[[1]] < \[Epsilon]/10 Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]] &&  ind>1 && (flux[[2]]+flux[[1]]<arraykm2[[ind-1,3]]+arraykm2[[ind-1,2]] || flux[[2]]+flux[[1]] < \[Epsilon]*10^-40 Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]]), 
i= i+1; 
Print["i=",i];
If[i==3, 
i=0;
Break[];
],
i=0;
Print["i=",i];
]; 
];
listm2[[count]]=Table[{l,2,arraykm2[[i,1]],arraykm2[[i,2]],arraykm2[[i,3]],arraykm2[[i,4]],arraykm2[[i,5]],arraykm2[[i,6]],arraykm2[[i,7]]},{i,1,ind}]; 
Print["count=", count];
Print["listm2[[count]]=", listm2[[count]]];
ind=0;
arraykm2=ConstantArray[{0,0,0,0,0,0,0},kmax];

If[Total[listm2[[count]][[All,5]]]+Total[listm2[[count]][[All,4]]] < \[Epsilon] (sumlm22+sumlm22h) , Break[];];
];
];

Clear[m];

matrix[[8]]=Select[Partition[Flatten[Join[listlm22[[1]],listm2]],9], #!={0,0,0,0,0,0,0,0,0}&];
Print["matrix[[8]]=",matrix[[8]]];

ind=0;
count=0;


(**************** m>-5 & m != (0,1,2)***************************************************)

For[m=-5,m<=mmax,m++,
If[m!=2 && m!=1 && m!=0, 
For[l=Max[0,Abs[m]],l<=lmax,l++,
If[EvenQ[l+Abs[m]],
count=count+1;
If[m==0, kmin=-1, If[m<0, kmin = Floor[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]], If[aa>=0,kmin = Floor[10 m e^2],kmin = Floor[-10 m e^2] ]]];
(*If[m==0, kmin=1, If[m<0, kmin = Floor[10 Abs[m] e^2], kmin = Ceiling[10 m e^2]]]*)

For[k=kmin,k>=-kmax,k--,
ind=ind+1;
Print["(l,m,k)=",{l,m,k}];
flux = Int[p,e,l,m,k, aa];
If[m==0,shift=0,shift=1];
arrayk[[ind]]=Flatten[{k,flux}];
kfin=k;
If[flux[[2]]+flux[[1]] < \[Epsilon]/10 Max[Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]],Max[arrayklm11[[All,3]]+arrayklm11[[All,2]]],Max[arrayklm00[[All,3]]+arrayklm00[[All,2]]]] &&  ind>1 && (flux[[2]]+flux[[1]]<arrayk[[ind-1,3]]+arrayk[[ind-1,2]] || flux[[2]]+flux[[1]] < \[Epsilon]*10^-40 Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]]),  
Print["i=",i];
i= i+1; 
If[i==3, 
i=0;
Break[];
],
i=0;
Print["i=",i];
]; 
];

i=0;
ind=0;

For[k=kmin+1,k<=kmax,k++,
If[k>0, If[k>Floor[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]], Break[]], If[k>=Floor[- m \[CapitalOmega]\[Phi]/\[CapitalOmega][[1]]], Break[]] ];
Print["(l,m,k)=",{l,m,k}];
flux = Int[p,e,l,m,k, aa];
If[m==0,shift=0,shift=1];
ind=ind+1;
arrayk2[[ind]]=Flatten[{k,flux}];

If[flux[[1]]+flux[[2]] < \[Epsilon]/10 Max[Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]],Max[arrayklm11[[All,3]]+arrayklm11[[All,2]]],Max[arrayklm00[[All,3]]+arrayklm00[[All,2]]]] && ind>1 && (flux[[2]]+flux[[1]]<arrayk2[[ind-1,3]]+arrayk2[[ind-1,2]] || flux[[1]]+flux[[2]] < \[Epsilon]*10^-40 Max[arrayklm22[[All,3]]+arrayklm22[[All,2]]]), 
i= i+1; 
Print["i=",i];
If[i==3, 
i=0;
Break[];
],
i=0;
Print["i=",i];
]; 
];

If[m==0,shift=0,shift=1];
list[[count]]=Join[Table[{l,m,arrayk[[i,1]],arrayk[[i,2]],arrayk[[i,3]],arrayk[[i,4]],arrayk[[i,5]],arrayk[[i,6]],arrayk[[i,7]]},{i,1,-kfin+kmin+1}],Table[{l,m,arrayk2[[i,1]],arrayk2[[i,2]],arrayk2[[i,3]],arrayk2[[i,4]],arrayk2[[i,5]],arrayk2[[i,6]],arrayk2[[i,7]]},{i,1,ind}]]; 

arrayk=ConstantArray[{0,0,0,0,0,0,0},kmax];
arrayk2=ConstantArray[{0,0,0,0,0,0,0},kmax];

i=0;
ind=0;

If[Total[list[[count]][[All,5]]]+Total[list[[count]][[All,4]]] < \[Epsilon]  (sumlm22+sumlm22h), Break[];];
];
];

count=0;
matrix[[m+6]] = Select[Partition[Flatten[list],9], #!={0,0,0,0,0,0,0,0,0}&]; 

list = ConstantArray[0,kmax*lmax*50];

If[m>mthr,
If[(Total[matrix[[m+6]][[All,5]]]+Total[matrix[[m+6]][[All,4]]])/(1-((Total[matrix[[m+6]][[All,5]]]+Total[matrix[[m+6]][[All,4]]])/(Total[matrix[[m+5]][[All,5]]]+Total[matrix[[m+5]][[All,4]]])))< \[Epsilon]/2 (Total[Partition[Flatten[matrix],9][[All,5]]]+Total[Partition[Flatten[matrix],9][[All,4]]]), Break[];];
];
];
];

matrix= Partition[Flatten[matrix],9];
matrix1 = Select[matrix, #!={0,0,0,0,0,0,0,0,0}&];

Return[matrix1]

];
];


(* ::Subsubsection::Closed:: *)
(*Exporting data*)


ppp=40;
precODE=40;

Print["get"];
griglia = Import["/Users/susannabarsanti/Documents/EMRI_ECC/CODICI_FINALI/griglia.wdx"];


Print["Precision=",Precision[griglia]];
Print["Length=",Length[griglia]];
Print["spin=",N[griglia[[1,1]],5]];

check1=1;
check15=15;

Print["griglia_1=",griglia[[check1]]];
Print["griglia_15=",griglia[[check15]]];


multipolestab=ParallelTable[{N[griglia[[i,1]],40],N[griglia[[i,2]],40],N[griglia[[i,3]],40],If[griglia[[i,1]] >= 0 , sum[griglia[[i,2]], griglia[[i,3]], griglia[[i,1]]], sumaneg[griglia[[i,2]], griglia[[i,3]], griglia[[i,1]]]]},{i,{197}}];
totaltab=Table[{multipolestab[[i,1]],multipolestab[[i,2]],multipolestab[[i,3]],Total[multipolestab[[i,4,All,4]]],Total[multipolestab[[i,4,All,5]]],Total[multipolestab[[i,4,All,6]]],Total[multipolestab[[i,4,All,7]]],Total[multipolestab[[i,4,All,8]]],Total[multipolestab[[i,4,All,9]]]},{i,1,multipolestab[[All,1]]//Length}];


Export["scal_multipoles_1_200.wdx", multipolestab];
Export["scal_multipoles_1_200.dat", multipolestab];
Export["scal_multipoles_1_200.mx", multipolestab];
Export["scal_total_1_200.wdx", totaltab];
Export["scal_total_1_200.mx", totaltab];
Export["scal_total_1_200.dat", totaltab];


Exit[];
