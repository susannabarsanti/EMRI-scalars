(* ::Package:: *)

SetSystemOptions["ParallelOptions"->"MKLThreadNumber"->5];
SetSystemOptions["ParallelOptions"->"ParallelThreadNumber"->5];

Needs["EmriInspiral`"];
ParallelNeeds["EmriInspiral`"];
LaunchKernels[5];

months = 12;
montfact = 365/12;

{m1,m2,spin,\[Theta]s,\[Phi]s,\[Theta]L,\[Phi]L,\[Phi]0,charge,dL}={10^6,15,0.9`350,\[Pi]/2,\[Pi]/2,\[Pi]/4,\[Pi]/4,0,1/10,1000};

massphi = 0.018`350;

r0 = findr0[m1,m2,spin,charge,massphi,242/100,365];

Print["r0=",N[r0,10]];;

CC = 299792458 10^-3;
MSUN = 1477 10^-3;
DAY = 3600 24;

rmax = OrbitalEvoT[m1,m2,spin,charge,massphi,r0,\[Phi]0,montfact*months][[2]][montfact*months*DAY]; (* final radius after T *)	
fdopp = dopplerf[m1,m2,spin,charge,massphi,r0,\[Phi]0,\[Theta]s,\[Phi]s,montfact*months,montfact*months*DAY]; (* max fdop @ rmax *)
fmax = -1; (*CC KerrOrbit[spin,rmax][[1]]/(m1 MSUN \[Pi])+fdopp;*) (* f @ rmax + fdop @ rmax *)
fmin = 10^-4;

Print["Precision ",r0,"=",Precision[r0]];

r0 = N[r0,200];

Print["Precision tensor flux = ",Precision[dEGRdr[r0,spin]]];
Print["Precision scalar flux = ",Precision[dEddr[r0,spin,massphi]]];

Print["Final precision r0 = ",Precision[r0]];

\[Delta]={(*{10^-4,10^-4,10^-4,10^-2,10^-4,10^-4,10^-4,10^-4,10^-4,10^-4},
5{10^-5,10^-5,10^-5,10^-3,10^-5,10^-5,10^-5,10^-5,10^-5,10^-5},
{10^-5,10^-5,10^-5,10^-3,10^-5,10^-5,10^-5,10^-5,10^-5,10^-5},
5{10^-6,10^-6,10^-6,10^-4,10^-6,10^-6,10^-6,10^-6,10^-6,10^-6},
{10^-6,10^-6,10^-6,10^-4,10^-6,10^-6,10^-6,10^-6,10^-6,10^-6},
5{10^-7,10^-7,10^-7,10^-5,10^-7,10^-7,10^-7,10^-7,10^-7,10^-7},*)
{10^-7,10^-7,10^-7,10^-5,10^-5,10^-7,10^-7,10^-7,10^-7,10^-7,10^-7},
5{10^-8,10^-8,10^-8,10^-6,10^-6,10^-8,10^-8,10^-8,10^-8,10^-8,10^-8},
{10^-8,10^-8,10^-8,10^-6,10^-6,10^-8,10^-8,10^-8,10^-8,10^-8,10^-8},
5{10^-9,10^-9,10^-9,10^-7,10^-7,10^-9,10^-9,10^-9,10^-9,10^-9,10^-9},
{10^-9,10^-9,10^-9,10^-7,10^-7,10^-9,10^-9,10^-9,10^-9,10^-9,10^-9}};

Print["number of shifts =",Length[\[Delta]]];

Print["final radius = ",OrbitalEvoT[m1,m2,spin,charge,massphi,r0,\[Phi]0,365][[2]][365*24*3600]];

der=ConstantArray[0,Length[\[Delta]]];
\[CapitalGamma]=ConstantArray[0,Length[\[Delta]]];

der=ParallelTable[derwave[\[Delta][[i]],m1,m2,spin,charge,massphi,r0,\[Phi]0,\[Theta]s,\[Phi]s,\[Theta]L,\[Phi]L,dL,months*montfact],{i,1,Length[\[Delta]]}];
Export["derivative_119WHPI_Md6m15a09_mu0018_d01_12months_l3.mx",der];

Table[\[CapitalGamma][[i]]=fisherN[der[[i,1]],"B",12,fmin,fmax]+fisherN[der[[i,2]],"B",12,fmin,fmax],{i,1,Length[\[Delta]]}];
Export["fisher_119WHPIBWD_Md6m15a09_mu0018_d01_12months_l3.mx",\[CapitalGamma]];
