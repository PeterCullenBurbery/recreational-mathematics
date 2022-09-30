(* ::Package:: *)

BeginPackage["PeterBurbery`RecreationalMathematics`"];

(* Declare your packages public symbols here. *)

FivePointConic;
NinePointCubic;
EulerLinePoints;
NinePointQuadric;
BalancedTernary;
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



End[]; (* End `Private` *)

EndPackage[];
