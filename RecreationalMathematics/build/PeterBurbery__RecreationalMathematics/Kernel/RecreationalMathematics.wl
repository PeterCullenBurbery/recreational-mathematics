(* ::Package:: *)

BeginPackage["PeterBurbery`RecreationalMathematics`"];

(* Declare your packages public symbols here. *)

FivePointConic;
NinePointCubic;
EulerLinePoints;
NinePointQuadric;
BalancedTernary;
CatalanUnrank;
AllBalancedGroupingSymbols;
Derangements;
Begin["`Private`"];

(* Define your public and private symbols here. *)
FivePointConic//ClearAll
FivePointConic[pts_]/;MatrixQ[pts]&&Dimensions[pts]==={5,2}:=FivePointConic[pts,{\[FormalX],\[FormalY]}];FivePointConic[pts_,{xx_,yy_}]/;MatrixQ[pts]&&Dimensions[pts]==={5,2}:=Times@@Apply[Power,DeleteCases[FactorList[Det[PadLeft[Apply[Function[{x,y},{x,y,x^2,x y,y^2}],Append[pts,{xx,yy}],{1}],{6,6},1]]],{_?NumericQ,_}],{1}]

NinePointCubic//ClearAll
NinePointCubic[pts_]/;MatrixQ[pts]&&Dimensions[pts]==={9,2}:=NinePointCubic[pts,{\[FormalX],\[FormalY]}];NinePointCubic[pts_,{xx_,yy_}]/;MatrixQ[pts]&&Dimensions[pts]==={9,2}:=Times@@Apply[Power,DeleteCases[FactorList[Det[PadLeft[Apply[Function[{x,y},{x,y,x^2,x y,y^2,x^3,x^2 y,x y^2,y^3}],Append[pts,{xx,yy}],{1}],{10,10},1]]],{_?NumericQ,_}],{1}]


EulerLinePoints//ClearAll;

EulerLinePoints[vert_]/;MatrixQ[vert,Internal`RealValuedNumericQ]:=Module[{o,g,n,h},
Which[Dimensions[vert]==={3,2},
o=First[Circumsphere[vert]],
Dimensions[vert]==={3,3},
o=ResourceFunction["Circumcircle3D"][vert,"Center"],
True,Return[$Failed,Module]];
g=Mean[vert];
h=3 g -2 o;
n=Mean[{o,h}];
{o,g,n,h}];

EulerLinePoints[Triangle[vert_]]:=If[MatrixQ[vert],EulerLinePoints[vert],EulerLinePoints/@vert]

NinePointQuadric//ClearAll;
NinePointQuadric[pts_]/;MatrixQ[pts]&&Dimensions[pts]==={9,3}:=NinePointQuadric[pts,{\[FormalX],\[FormalY],\[FormalZ]}];
NinePointQuadric[pts_,{xx_,yy_,zz_}]/;MatrixQ[pts]&&Dimensions[pts]==={9,3}:=Times@@Apply[Power,DeleteCases[FactorList[Det[PadLeft[Apply[Function[{x,y,z},{x,x^2,y,x y,y^2,z,x z,y z,z^2}],Append[pts,{xx,yy,zz}],{1}],{10,10},1]]],{_?NumericQ,_}],{1}]



BalancedTernary//ClearAll

BalancedTernary[n:(_Integer|_Row)]:=Module[{list, len},
If[Head[n]===Row,
list=Normal[n]/.\!\(\*UnderscriptBox[\(1\), \(_\)]\)-> -1;
FromDigits[list,3],
len=Ceiling[Log[3,1+2Abs[n]]];
Row[(IntegerDigits[-n- Quotient[3^len-1,2],3,len]-1)/.{{}->{0},-1->\!\(\*UnderscriptBox[\(1\), \(_\)]\)}]]]

CatalanUnrank//ClearAll
CatalanUnrank[n_,rank_]:= 
Module[{ lo=0,y=0, a=Table[0,2 n],m},
Do[
m=Binomial[2 n -x, n-(x+y+1)/2]-Binomial[2 n -x, n-1-(x+y+1)/2];
(*these terms make the Catalan triangle, or ballot numbers*)
If[rank<= lo+m-1,
y=y+1;
a[[x]]=0,
lo=lo+m;
y=y-1;
a[[x]]=1],
{x,1,2 n}];
a]



Derangements//ClearAll
Derangements[n_?IntegerQ]:=With[
{perms=Permutations[Range@n]},
Pick[perms,Length/@PermutationSupport/@perms,n]
]/;n>=0

AllBalancedGroupingSymbols//ClearAll
AllBalancedGroupingSymbols[{openingsymbol_,closingsymbol_},numberofgroups_]:=StringJoin[ResourceFunction["CatalanUnrank"][numberofgroups,#]/.{0->openingsymbol,1->closingsymbol}]&/@Range[0,CatalanNumber[numberofgroups]-1]




End[]; (* End `Private` *)

EndPackage[];
