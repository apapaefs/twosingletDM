(* ::Package:: *)

(* ::Chapter:: *)
(*Helper functions*)


(*Import helper functions*)
<< (NotebookDirectory[]<>"CreateModel.m")


(* ::Chapter:: *)
(*Implementing TRSM potential*)


(*Define higgs fields*)
higgsbase = {\[Rho]1, \[Eta]1, \[Psi]1, \[Zeta]1, \[Zeta]X, \[Zeta]S};


(*Assign vevs at T=0*)
higgsvev = {0, 0, 0, v1, 0, vs};


(*Assign vevs at T != 0*)
higgsvevFiniteTemp = {0, 0, 0, w1, wx, ws};


(*Replacement list for the vevs*)
VEVRep = Table[i[[1]]->i[[2]],{i,Transpose[{higgsbase,higgsvev}]}];


(*Replacement list set fields zero*)
RepHiggsZero = Table[i->0,{i,higgsbase}];


(*Define number of Higgses*)
nHiggs = Length[higgsbase];


(*Define parameters of the potential*)
par = {mphiSq, lphi, msSq, ls, mxSq, lx, lphis, lphix, lsx};


(*Define doublets and scalar fields*)
Phi = {
  higgsbase[[1]] + I*higgsbase[[2]],
  higgsbase[[4]] + I*higgsbase[[3]]
}/Sqrt[2];
X = higgsbase[[5]];
S = higgsbase[[6]];
(*Print fields*)
Phi//MatrixForm//TraditionalForm
X//MatrixForm//TraditionalForm
S//MatrixForm//TraditionalForm



(*Write the potential*)
$Assumptions = higgsbase \[Element] Reals;
C11 = Simplify[ConjugateTranspose[Phi] . Phi];
CS = S*S;
CX = X*X;
VHiggs = par[[1]]*C11 + par[[2]] * C11^2 + par[[3]] * CS + par[[4]] * CS^2 + par[[5]] * CX + par[[6]] * CX^2 + par[[7]] * C11 * CS + par[[8]] * C11 * CX + par[[9]] * CX * CS //Simplify


(*Calculate Tadpoles for Minimisation conditions*)
VHiggsGrad=D[VHiggs,{higgsbase}]/.VEVRep//Transpose//Simplify;
{"Derivative with respective to" higgsbase,VHiggsGrad}//Transpose//TableForm


(*Define input parameters of the potential*)
InputParameters = {m1, m2, m3, vs, a12, lx, lphix, lsx};
Mhh = m1^2 Cos[a12]^2 + m2^2 Sin[a12]^2;
Mss = m1^2 Sin[a12]^2 + m2^2 Cos[a12]^2;
Mhs = (m2^2 - m1^2) Sin[a12] Cos[a12];


ListOfDependentParameter = Select[RemoveDuplicates[Join[RemoveDuplicates[InputParameters],RemoveDuplicates[Join[par,higgsvev]]]],MemberQ[Join[par,higgsvev],#]&];
Print["List of dependent parameters that you have to calculate is\n", ListOfDependentParameter//TableForm]


(*Which parameters should be replaced through the minimum conditions?*)
AdditionalConditions = {v1 -> vev0};

MassInputRulesRaw = {
  lphi  -> Mhh/(2 v1^2),
  ls    -> Mss/(8 vs^2),
  lphis -> Mhs/(2 v1 vs),
  mxSq  -> (m3^2 - lphix v1^2 - 2 lsx vs^2)/2
};

TadpoleRulesRaw = Solve[
  VHiggsGrad == Table[0, {VHiggsGrad}],
  {mphiSq, msSq}
][[1]];

MassInputRules =
  MassInputRulesRaw //. AdditionalConditions // Simplify;

TadpoleRules =
  TadpoleRulesRaw //. MassInputRulesRaw //. AdditionalConditions // Simplify;

DepententParameters =
  Flatten[{AdditionalConditions, MassInputRules, TadpoleRules}];

DepententParameters // TableForm
(* verify masses *)
evenFields = {higgsbase[[4]], higgsbase[[6]]};

MassEven =
  Table[
    D[VHiggs, evenFields[[i]], evenFields[[j]]],
    {i, 2}, {j, 2}
  ] /. VEVRep //. DepententParameters // Simplify;

R12 = {{Cos[a12], -Sin[a12]}, {Sin[a12], Cos[a12]}};

Simplify[R12 . MassEven . Transpose[R12]]

(* another check*)
D[VHiggs, higgsbase[[5]], higgsbase[[5]]] /. VEVRep //. DepententParameters // Simplify


(*Sanity check*)
If[(DepententParameters[[All,1]]//Sort) == (ListOfDependentParameter//Sort),Print["Everything looks good!\nAll dependent parameters have been defined!"],Print["You still have to define\n",TableForm[RemoveDuplicates[Join[RemoveDuplicates[DepententParameters[[All,1]]],RemoveDuplicates[ListOfDependentParameter]]]]]]


(* ::Chapter:: *)
(*Code*)
(*In this section the model implementation is generated.*)


(* ::Section:: *)
(*MinimizeOrderVEV*)


(*This section calculates the vector VevOrder which is set in the MinimizeOrderVEV function.*)
VEVList={};
Table[If[Not[PossibleZeroQ[higgsvevFiniteTemp[[i]]]],AppendTo[VEVList,i-1]];{i,higgsvevFiniteTemp[[i]]},{i,Length[higgsvevFiniteTemp]}];
Table["VevOrder[" <> ToString[vev-1] <> "] = " <> ToString[VEVList[[vev]]] <> ";" ,{vev,Length[VEVList]}]//TableForm


(* ::Section:: *)
(*Calculate Higgs Curvature L1*)


CurvatureL1 = GenerateCurvature[1];
If[CurvatureL1=={},"No entries for Curvature_Higgs_L1",CurvatureL1//TableForm]





(* ::Section:: *)
(*Calculate Higgs Curvature L2*)


CurvatureL2 = GenerateCurvature[2];
If[CurvatureL2=={},"No entries for Curvature_Higgs_L2",CurvatureL2//TableForm]


(* ::Section:: *)
(*Calculate Higgs Curvature L3*)


CurvatureL3 = GenerateCurvature[3];
If[CurvatureL3=={},"No entries for Curvature_Higgs_L3",CurvatureL3//TableForm]


(* ::Section:: *)
(*Calculate Higgs Curvature L4*)


CurvatureL4 = GenerateCurvature[4];
If[CurvatureL4=={},"No entries for Curvature_Higgs_L4",CurvatureL4//TableForm]


(* ::Section:: *)
(*Counterterm potential*)


(*Define counterterms*)
parCT = Join[Table[Symbol["d"<>ToString[p]],{p,par}], Table[Symbol["dT"<>ToString[i]],{i,Length[higgsbase]}]]


(*Define the counterterm potential*)
VCT = ((VHiggs/.Table[par[[i]]->parCT[[i]],{i,Length[par]}]) + parCT[[(par//Length)+1;;]] . higgsbase)


(*Define the sytem of equation to determine the CTs*)
nCount=0;
EqTar={};
Table[nCount+=1; GL[nCount]=D[VCT,higgsbase[[i]]]/.VEVRep; If[PossibleZeroQ[GL[nCount]],nCount-=1,AppendTo[EqTar,-NCW[i]]],{i,Length[higgsbase]}];
Table[nCount+=1; GL[nCount]=D[VCT,higgsbase[[i]],higgsbase[[j]]]/.VEVRep; If[PossibleZeroQ[GL[nCount]],nCount-=1,AppendTo[EqTar,-HCW[i,j]]],{i,Length[higgsbase]},{j,i,Length[higgsbase]}];
EqMatrix = Table[D[GL[i],j],{i,nCount},{j,parCT}];
Print["EqMatrix rank is ", EqMatrix//MatrixRank]
EqMatrix//MatrixForm
Print["EqTar is ", EqTar//TableForm]


(*Define the sytem of equation to determine the CTs*)
SysOriginal = Append[EqMatrix//Transpose,EqTar]//Transpose;
Print["SysOriginal rank is ", SysOriginal//MatrixRank]
SysOriginal//MatrixForm


(*Find the CW relations that help solve the system*)
SysOriginal = Append[EqMatrix//Transpose,EqTar]//Transpose;
Stmp=RowReduce[SysOriginal, ZeroTest -> (! FreeQ[#, HCW[_,_]] || !FreeQ[#,NCW[_]] &)];
CWRelations =Solve[Select[Stmp,#[[;;-2]]==Table[0,{i,(EqMatrix//Dimensions)[[2]]}]&][[All,-1]]==0,-EqTar][[1]];
Print["Relations found between the CW 1st and 2nd derivative\n", CWRelations//TableForm]


(*Reduce the system of equation so that Mathematica can solve it*)
SysOriginal = Append[EqMatrix//Transpose,EqTar]//Transpose;
SysOriginal = (SysOriginal/.CWRelations);
NewMatrix={SysOriginal[[1]]};
nCount=1;
While[MatrixRank[NewMatrix]!=MatrixRank[SysOriginal] && nCount < (Length[SysOriginal] + 1),
nCount+=1;
If[MatrixRank[Append[NewMatrix,SysOriginal[[nCount]]]]>MatrixRank[NewMatrix],AppendTo[NewMatrix,SysOriginal[[nCount]]]]]


(*Determine the CTs*)
NullSpaceofCT = NullSpace[NewMatrix[[All,;;-2]]];
NullSpaceofCT = Sum[NullSpaceofCT[[i]]*Subscript[t, i],{i,Length[NullSpaceofCT]}];
CTs = LinearSolve[NewMatrix[[All,;;-2]],NewMatrix[[All,-1]]];
CTs += NullSpaceofCT;
{parCT,Table["\t\[RightArrow]\t",{i,parCT}],CTs}//Transpose//TableForm





Transpose[{Range[Length[parCT]], parCT, CTs}] // TableForm


(* ::Subsection:: *)
(*The Subscript[t, i] still need to be fixed *)
(*We chose the following Subscript[t, i] to match the python and Maple implementation, but any other choice is valid .*)


(*Choose the null-space counterterms dlx, dlphix, dlsx, dT6 to vanish.*)
tiChoice = {
  Subscript[t, 1] -> 0,
  Subscript[t, 2] -> 0,
  Subscript[t, 3] -> 0,
  Subscript[t, 4] -> 0
} /. CWRelations;


(*Fixed CTs*)
tiCTs = CTs /. tiChoice // Simplify;
{parCT, Table["\t\[RightArrow]\t", {i, parCT}], tiCTs} // Transpose // TableForm


(* ::Subsection:: *)
(*Check that the counter terms are correct*)


test = (((((EqMatrix) . parCT)/.(Thread[parCT -> tiCTs]))==(EqTar/.CWRelations))//Simplify);
If[test,Print["Success!"],Print["Something went wrong.\nMake sure the missmatch between the CW derivatives and CT cannot bet solved by a CWRelations substitution rule."];test//TableForm]


(* ::Subsection:: *)
(*Print the CT into C++*)


Table[ToString[parCT[[i]]//CForm] <> " = " <> ToString[CForm[tiCTs[[i]]/.{NCW[x_]->NablaWeinberg[x-1],HCW[x_,y_]->HesseWeinberg[x-1,y-1]}]] <> ";",{i,Length[parCT]}]//TableForm


(* ::Section:: *)
(*CTs Curvatures*)


(* ::Section:: *)
(*Calculate Higgs Curvature CT L1*)


CTCurvatureL1 = GenerateCTsCurvature[1];
If[CTCurvatureL1=={},"No entries for Curvature_Higgs_CT_L1",CTCurvatureL1//TableForm]


(* ::Section:: *)
(*Calculate Higgs Curvature CT L2*)


CTCurvatureL2 = GenerateCTsCurvature[2];
If[CTCurvatureL2=={},"No entries for Curvature_Higgs_CT_L2",CTCurvatureL2//TableForm]


(* ::Section:: *)
(*Calculate Higgs Curvature CT L3*)


CTCurvatureL3 = GenerateCTsCurvature[3];
If[CTCurvatureL3=={},"No entries for Curvature_Higgs_CT_L3",CTCurvatureL3//TableForm]


(* ::Section:: *)
(*Calculate Higgs Curvature CT L4*)


CTCurvatureL4 = GenerateCTsCurvature[4];
If[CTCurvatureL4=={},"No entries for Curvature_Higgs_CT_L4",CTCurvatureL4//TableForm]


(* ::Chapter:: *)
(*Gauge interaction*)
(*This section defines the interaction between gauge bosons and scalars*)
(*First define the basis for the gauge bosons*)


(*Define the gauge fields*)
GaugeBasis={Subscript[W, 1],Subscript[W, 2],Subscript[W, 3],Subscript[B, 0]}


(*Define the covariant derivative*)
Subscript[D, \[Mu]]= -I Subscript[C, g]/2 Sum[Subscript[W, i]PauliMatrix[i],{i,3}] -I Subscript[C, gs]/2 Subscript[B, 0] IdentityMatrix[2];
Subscript[D, \[Mu]]//MatrixForm


(*Define the gauge potential*)
Vgauge =
  (Conjugate[(Subscript[D, \[Mu]] . Phi)] .
   (Subscript[D, \[Mu]] . Phi)) // ComplexExpand // Simplify


(* ::Section:: *)
(*Gauge Curvatures*)


(* ::Section:: *)
(*Calculate Gauge Curvature L4*)


GaugeCurvatureL4 = GenerateGaugeCurvature;
If[GaugeCurvatureL4=={},"No entries for "Curvature_Gauge _G2H2 "",GaugeCurvatureL4//TableForm]


(* ::Chapter:: *)
(*Lepton interaction*)
(*This defines the interactions with Leptons. The N2HDM has several types to avoid FCNC. *)
(*I am showing only Type 1 here as an example, meaning the Leptons and down-type quarks couple to the first doublet while the up-type quark couple to the second doublet.*)


Clear[eL, eR, muL, muR, tauL, tauR, veL, vmuL, vtauL,
      ye, ymu, ytau, NuL, ER, EL, LepBase, PiLep,
      VFLep, MassLep, MLep, RepLepMass];

(*First we define the coupling matrix for the leptons*)
PiLep = {{ye, 0, 0}, {0, ymu, 0}, {0, 0, ytau}};
PiLep//MatrixForm


(*Define the components of the lepton sectors*)
NuL = {veL, vmuL, vtauL};
ER  = {eR, muR, tauR};
EL  = {eL, muL, tauL};

LepBase = {eL, eR, muL, muR, tauL, tauR, veL, vmuL, vtauL};

Phi /. VEVRep // Simplify



(*Defining the Lepton potential (=-L_yuk)*)
VFLep = Sum[
   PiLep[[i, j]] ER[[j]] (NuL[[i]] Phi[[1]] + EL[[i]] Phi[[2]]),
   {i, 3}, {j, 3}
] // Expand;


Phi /. VEVRep


(*Leptonic mass matrix*)
MassLep = D[VFLep, {LepBase, 2}] /. VEVRep // Simplify;
MassLep // MatrixForm


(*Fermionic masses*)
MatrixForm[MLep = Eigenvalues[MassLep]]


(*Calculate the Yukawas*)
RepLepMass = Solve[{MLep[[5]] == CMassElectron, MLep[[7]] == CMassMu, MLep[[9]] == CMassTau}, {ye, ymu, ytau}][[1]]


(* ::Section:: *)
(*Leptonic Curvatures*)


LeptonCurvatureL3 = GenerateLeptonCurvature;
If[LeptonCurvatureL3=={},"No entries for Curvature_Lepton_F2H1",LeptonCurvatureL3//TableForm]


(* ::Chapter:: *)
(*Quark interaction*)


(*First we define the quark doublets*)
UL = {uL, cL, tL};
DL = {dL, sL, bL};
UR = {uR, cR, tR};
DR = {dR, sR, bR};

baseQuarks = {uR, cR, tR, dR, sR, bR, uL, cL, tL, dL, sL, bL};


(*This part is the CKM matrix. If you don't want to include it, set it to unity*)
VCKM = {{V11, V12, V13}, {V21, V22, V23}, {V31, V32, V33}};
VCKM//MatrixForm
CKML = {V11 -> Vud, V12 -> Vus, V13 -> Vub, V21 -> Vcd, V22 -> Vcs, V23 -> Vcb, V31 -> Vtd, V32 -> Vts, V33 -> Vtb};
VCKM/.CKML//MatrixForm


(*Yukawa*)
Dd = {{yd, 0, 0}, {0, ys, 0}, {0, 0, yb}};
DU = {{yu, 0, 0}, {0, yc, 0}, {0, 0, yt}};
Dd//MatrixForm
DU//MatrixForm


VF = (
   (UL . VCKM . Dd . DR) Phi[[1]]
 + (DL . Dd . DR) Phi[[2]]
 + (UL . DU . UR) (higgsbase[[4]] - I*higgsbase[[3]])/Sqrt[2]
 - (DL . ConjugateTranspose[VCKM] . DU . UR)
     (higgsbase[[1]] - I*higgsbase[[2]])/Sqrt[2]
) // Simplify;


MQuark = D[VF,{baseQuarks,2}]/.VEVRep;
MQuark//MatrixForm


mQ = Eigenvalues[MQuark];
mQ//MatrixForm


RepMassQuark = {
  yd -> Sqrt[2] CMassDown/v1,
  ys -> Sqrt[2] CMassStrange/v1,
  yb -> Sqrt[2] CMassBottom/v1,
  yu -> Sqrt[2] CMassUp/v1,
  yc -> Sqrt[2] CMassCharm/v1,
  yt -> Sqrt[2] CMassTop/v1
};
RepMassQuark//MatrixForm


(* ::Section:: *)
(*Quark Curvatures*)


QuarkCurvatureL3 = GenerateQuarkCurvature;
If[QuarkCurvatureL3=={},"No entries for Curvature_Quark_F2H1",QuarkCurvatureL3//TableForm]


(* ::Section:: *)
(*Model implementation*)


ListImplementedModels//TableForm


(*Arguments of CreateModelFile*)
ImplementModel[
"TRSM", (*Model name*)
higgsbase, (*List of all field*)
higgsvev, (*VEVs as T = 0*)
higgsvevFiniteTemp, (*VEVs as T != 0*)
VEVList, (*Derivation VEVList (calculated automatically)*)
par, (*Potential parameters *)
InputParameters, (*Input parameters *)
{125.09, 300, 400, 200, 0.2, 0.1, 0.05, 0.05},(*Specific input parameters for the example file*)
DepententParameters, (*Depentent parameters and how to calculate them*)
CurvatureL1, (*Scalar curvatures L1 (calculated automatically)*)
CurvatureL2, (*Scalar curvatures L2 (calculated automatically)*)
CurvatureL3, (*Scalar curvatures L3 (calculated automatically)*)
CurvatureL4, (*Scalar curvatures L4 (calculated automatically)*)
GaugeCurvatureL4, (*Gauge curvatures L4 (calculated automatically)*)
LeptonCurvatureL3, (*Lepton curvatures L3 (calculated automatically)*)
QuarkCurvatureL3, (*Quark curvatures L3 (calculated automatically)*)
parCT, (*Counter terms*)
CTCurvatureL1, (*Counterterm scalar curvatures L1 (calculated automatically)*)
CTCurvatureL2, (*Counterterm scalar curvatures L2 (calculated automatically)*)
CTCurvatureL3, (*Counterterm scalar curvatures L3 (calculated automatically)*)
CTCurvatureL4, (*Counterterm scalar curvatures L4 (calculated automatically)*)
GaugeBasis, (*Gauge fields*)
LepBase, (*Leptonic fields*)
baseQuarks] (*Quark fields*)




