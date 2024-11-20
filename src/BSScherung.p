
Program Scherung (input, output, f1);

{Version 1.0 Februar 1995}

{Dieses Programm dient zur Berechnung eines Scherungsparameter bei DNS.}
{Als Input File in PDB Format benoetigt man}
{	- File mit den Phosphat atome}
{	- File mit den Basenatomen (eins fuer jede Base bzw. Basenpaar)}

{Das Programm lauft wie folgt:}
{1	Zuerst wird eine BSpline Kurve basierend auf die Phosphatatome}
{    im ersten Input File definiert mit dem Koeffizient k=3}
{    (siehe W.M.Newman, R.F.Sproull, Principles of interactive}
{    computer graphic, Second Edition, Mc Graw Hill, Kapitel 21).}
{    Dabei wird der Abstand zwischen Kurve und Phosphatatome fuer}
{    jedes P berechnet und ausgeschrieben.}

{2	Durch die Basenatome (input File Nr 2) wird eine Beste Ebene}
{    gelegt und der Schnittpunkt zwischen dieser Ebene und der BSpline}
{    Kurve wird bestimmt.}

{3	An diesem Schnittpunkt wird ein Tangenzialvektor (zur BSpline}
{   Kurve) gelegt. Dieser Vektor wird dann auf die Ebene projeziert,}
{    die senkrecht zur Basenebene liegt und der oben berechnete}
{    Schnittpunkt sowie der Schwerpunkt der Basenatome enthaelt.}

{4	Der Scherungswinkel wird dann definiert als der Winkel zwischen}
{    dem Vektor normal zur Baseneben und dem bei Pt. 3 berechnete}
{    Vektor.}
    
{Normalvektor zur besten Ebene durch die Basenatome (B), sowie Tangenzialvektor}
{zeigen von 5' in 3' Richtung.}

{Vorzeichen Definition}
{=====================}

{Scherungswinkel ist > 0 	Skalarprodukt zwischen PB -> S Vektor und Projektion }
{							vom Tangenzialvektor auf die Ebene senkrecht zu }
{							derjeniger, die durch die Basenatome definiert}
{							wird (PaPb) ist > 0.}
{-------------------------------------------}
{INPUT:		File 1 P Atome In PDB Format}
{			File 2 Basenatome In PDB Format}

{OUTPUT:	Scherungswinkel und Torsionswinkel zwischen Tangentialvektor und Basenebennormale}

{Aufrufen des Programmes mit folgender shell}
{-------------------------------------------}
{cp   P.pdb      input1  # Input filename }
{cp   BP_A1.pdb  input2  # Input filename }
{}
{a.out                     # Programmname}
{3                         # BSpline Konstante}
{}
{/bin/rm input1            # Aufrauemen}
{/bin/rm input2            # Aufrauemen}
{-------------------------------------------}
{@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@}
{---Definitionen}

{Mb		Normale zu beta}
{Db			Abstand beta Ursprung}
{SS		Schwerpunkt der basenatome in beta}
{PB		Durchstosspunkt BSpline mit beta}
{B1		PB --> SS Vektor}

{PP,PAx	BSpline Gradient am PB}
{T			Tangenzial Vektor am PB = PB + PP}

{Ma		Normale zur Ebene alfa}
{Db			Abstand alfa Ursprung}

{PA		Projektion von T auf alfa}
{PbPa		Vektor von PB nach PA}
{w			Wnkel zwischen Mb und PbPa}

{cinP		Anzahl P Atome}
{cinB		Anzahl Basenatomen}
{InKP		P Koordinaten: PDB Format}
{InkB		Basenatomenkoordinaten: PDB Format}

{kk			BsPline Konstante}

{@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@}

Const 
  MaxAtom = 500;
  xacc = 1E-13;
  split = 500;{BSpline  und Ebene Schnitt}
  {MAXSTRLEN = 16#50;}
  {MAXSTRLEN = 100000;}

Type 

  {stringArray = packed array [1..MAXSTRLEN] Of char;}

  pa6 = packed array[1..6] Of char;

  Vektor = array[1..3] Of double;
  Matrix = array[1..3, 1..3] Of double;

  Koordinaten = Record
    name: pa6;
    x: Vektor;
  End;

  data = array[0..MaxAtom] Of Koordinaten;

Var 

  InKB, InKP: data;
  Mb, SS, PP, PB, PA, B1, PbPa, Ma, T: Vektor;
  Db, Da: double;
  PBx, PBy, PBz: double;
  PAx, PAy, PAz, u: double;
  kk: integer;
  cinB, cinP: integer;
  w,tau: double;
  i: integer;
  junk,B: Vektor;
  solution: boolean;
  knotK, knotN: integer; {BSpline Stuff...}

{--------------------------------------------------------------------------}

Procedure DataRead (inf: string; Var inp: data; Var cin: integer; start: integer);

Var 
  i, f: integer;
  ch: char;
  f1: TextFile;

Begin
  AssignFile(f1, inf);
  reset(f1);

  cin := start;

  While ((Not eof(f1)) And (cin <= MaxAtom)) Do
    Begin
      For i := 1 To 13 Do
        read(f1, ch);
      For f := 1 To 5 Do
        read(f1, Inp[cin].name[f]);
      For i := 1 To 11 Do
        read(f1, ch);
      With inp[cin] Do
        readln(f1, x[1], x[2], x[3]);
      cin := cin + 1;
    End;

  cin := cin - 1;

  close(f1);

End;

{--------------------------------------------------------------------------}

Procedure SetVektor (Var ab: Vektor; a, b: Vektor);

Var 
  i: integer;

Begin
  For i := 1 To 3 Do
    ab[i] := b[i] - a[i];
End;

{--------------------------------------------------------------------------}

Function atg (s, c: double): double;

Const 
  eps = 1E-13;
  rad = 0.01745329;

Var 
  t: double;

Begin
  If abs(c) < eps Then
    If s < 0 Then
      t := -90
  Else
    t := 90
  Else
    t := arctan(s / c) / rad;
  If (t > 0) And (c < 0) Then
    t := t - 180;
  If (t < 0) And (s > 0) Then
    t := t + 180;
  If (abs(t) < eps) And (c < 0) Then
    t := t + 180;
  atg := t
End;

{--------------------------------------------------------------------------}

Procedure Winkel (Var w: double; a, b: Vektor);

Var 
  m, n, mn, sinw, cosw: double;



Begin

  m := sqrt(sqr(a[1]) + sqr(a[2]) + sqr(a[3]));
  n := sqrt(sqr(b[1]) + sqr(b[2]) + sqr(b[3]));
  mn := (a[1]) * (b[1]) + (a[2]) * (b[2]) + (a[3]) * (b[3]);

  cosw := (mn) / (m * n);
  sinw := (sqrt(1 - sqr(cosw)));
  w := atg(sinw, cosw);

End;

{--------------------------------------------------------------------------}

Procedure VektorProdukt (Var c: Vektor; a, b: Vektor);

Begin
  c[1] := a[2] * b[3] - a[3] * b[2];
  c[2] := a[3] * b[1] - a[1] * b[3];
  c[3] := a[1] * b[2] - a[2] * b[1]
End;

{--------------------------------------------------------------------------}

Function Skalarprodukt (a, b: Vektor): double;

Var 
  i: integer;
  c: double;

Begin
  c := 0;
  For i := 1 To 3 Do
    c := c + (a[i] * b[i]);
  Skalarprodukt := c;
End;

{--------------------------------------------------------------------------}

Procedure SetPlane (Var a: Vektor; Var b: double; c, d, e: vektor);

Var 
  i: integer;
  amod: double;

Begin
  VektorProdukt(a, c, d);
  amod := sqrt(Skalarprodukt(a, a));
  For i := 1 To 3 Do
    a[i] := a[i] / amod;
  b := 0.0;
  For i := 1 To 3 Do
    b := b + a[i] * e[i];
End;

{--------------------------------------------------------------------------}

Procedure Projection (Var p: Vektor; a: Vektor; b: double; c, d: Vektor);

Label 
  10;

Var 
  i: integer;
  ac, ad, t: double;

Begin

  ac := Skalarprodukt(a, c);
  ad := Skalarprodukt(a, d);
  If ad <> 0 Then
    t := (b - ac) / (ad)
  Else
    Begin
      writeln(' P5* - P3* Vektor parallel zur Basenebene !!!!!');
      goto 10;
    End;

  For i := 1 To 3 Do
    p[i] := c[i] + t * d[i];

  10:

End;

{--------------------------------------------------------------------------}
Procedure Eigenvalues (Var z, v: Matrix);

Const 
  kRest = 1E-12;
  kEvSum = 2.99999;  {if sum |v[i,i]| > kEvSum, set v=identity}

Var 
  i, j, k: integer;
  th, ta, t, c, s, zz, sum: double;

Begin

{identitaet}
  For i := 1 To 3 Do
    For j := 1 To 3 Do
      If i = j Then
        v[i, j] := 1.0
      Else
        v[i, j] := 0.0;
  Repeat
    sum := ABS(z[1, 2]) + ABS(z[1, 3]) + ABS(z[2, 3]);
    If sum >= kRest Then
      Begin
        For i := 1 To 2 Do
          For j := (i + 1) To 3 Do
            If z[i, j] <> 0 Then
              Begin
                th := 0.5 * (z[i, i] - z[j, j]) / z[i, j];
                ta := SQRT(1 + SQR(th));
                If th < 0 Then
                  t := 1 / (th - ta)
                Else
                  t := 1 / (th + ta);
                If ABS(t) < kRest Then
                  t := 0;
                c := 1 / sqrt(1 + sqr(t));
                s := c * t;
                zz := z[i, i];
                z[i, i] := zz * sqr(c) + 2 * z[i, j] * s * c + z[j, j] * sqr(s);
                z[j, j] := z[j, j] - z[i, i] + zz;
                z[i, j] := 0;
                z[j, i] := 0;
                For k := 1 To 3 Do
                  Begin
                    zz := v[k, i];
                    v[k, i] := zz * c + v[k, j] * s;
                    v[k, j] := zz * s - v[k, j] * c;
                    If (k <> i) And (k <> j) Then
                      Begin
                        zz := z[k, i];
                        z[i, k] := zz * c + z[k, j] * s;
                        z[k, i] := z[i, k];
                        z[j, k] := zz * s - z[k, j] * c;
                        z[k, j] := z[j, k];
                      End
                  End
              End
      End
  Until sum < kRest;
  For i := 1 To 3 Do
    For j := (i + 1) To 3 Do
      If z[i, i] < z[j, j] Then
        Begin
          t := z[i, i];
          z[i, i] := z[j, j];
          z[j, j] := t;
          For k := 1 To 3 Do
            Begin
              t := v[k, i];
              v[k, i] := v[k, j];
              v[k, j] := t
            End
        End;
  sum := v[1, 1] * v[2, 2] * v[3, 3] + v[1, 2] * v[2, 3] * v[3, 1] + v[1, 3] * v[2, 1] * v[3, 2] - v
         [1, 1] * v[2, 3] * v[3, 2] - v[1, 2] * v[2, 1] * v[3, 3] - v[1, 3] * v[2, 2] * v[3, 1];
  If sum < 0 Then
    For i := 1 To 3 Do
      v[i, 3] := -v[i, 3];
  If (abs(v[1, 1]) + abs(v[2, 2]) + abs(v[3, 3])) >= KEvSum Then
    Begin
      For i := 1 To 3 Do
        For j := 1 To 3 Do
          If i <> j Then
            v[i, j] := 0
          Else
            v[i, j] := 1
    End
End;

{--------------------------------------------------------------------------}
Procedure Best_Plane (NPoint: integer; X: Data; Var M, SS: Vektor; Var D: double);

Label 
  10;

Var 
  i, j, k: integer;
  Y: array[1..MaxAtom] Of Vektor;
  c: Vektor;
  z, v: Matrix;
  eval, Mmod, Sig: double;

Begin
  If NPoint < 3 Then
    Begin
      writeln('Zu wenig Punkte um eine Ebene zu definieren !!!');
      goto 10;
    End;

{Calculate center of Mass c[i]}
  For i := 1 To 3 Do
    c[i] := 0.0;

  For i := 1 To NPoint Do
    With X[i] Do
      For j := 1 To 3 Do
        c[j] := c[j] + x[j];

  For i := 1 To 3 Do
    c[i] := c[i] / NPoint;

{Do coord transformation: Y:= x-c}
  For i := 1 To NPoint Do
    With X[i] Do
      For j := 1 To 3 Do
        Y[i, j] := x[j] - c[j];

{Set up the symmetric Matrix z}
  For i := 1 To 3 Do
    For j := 1 To 3 Do
      Begin
        z[i, j] := 0.0;
      End;

  For i := 1 To NPoint Do
    For j := 1 To 3 Do
      For k := 1 To 3 Do
        z[j, k] := z[j, k] + Y[i, j] * Y[i, k];
  z[2, 1] := z[1, 2];
  z[3, 1] := z[1, 3];
  z[3, 2] := z[2, 3];

{Berechne Eigenvektoren und Eigenwerte}
  Eigenvalues(z, v);

{Gleichung der Ebene: Die kolonne  von v sind die 3 Eigenvektoren}
{Matrix z ist zur Diagonalmatrix umgewandelt worden}
{Kleinste Eigewert is the value of the sum of residuals!!!}
  For i := 1 To 3 Do
    M[i] := v[i, 3];
  eval := z[3, 3];
  Sig := sqrt(eval / (NPoint - 3));

  Mmod := sqrt(Skalarprodukt(M, M));
  For i := 1 To 3 Do
    M[i] := M[i] / Mmod;
  D := 0.0;
  For i := 1 To 3 Do
    Begin
      D := D + M[i] * c[i];
      SS[i] := c[i];
    End;

  10:
End;

{--------------------------------------------------------------------------}
{ BSPLINE....}
{--------------------------------------------------------------------------}
Function Knot (i: integer): integer;

Begin
  If i < knotK Then
    Knot := 0
  Else If i > knotN Then
         Knot := knotN - knotK + 2
  Else
    Knot := i - knotK + 1
End;

{--------------------------------------------------------------------------}
Function NBlend (i, k: integer; u: double): double;

Var 
  t: integer;
  v: double;

Begin
  If k = 1 Then
    Begin
      v := 0;
      If (Knot(i) <= u) And (u < Knot(i + 1)) Then
        v := 1
    End
  Else
    Begin
      v := 0;
      t := Knot(i + k - 1) - Knot(i);
      If t <> 0 Then
        v := (u - Knot(i)) * NBlend(i, k - 1, u) / t;
      t := Knot(i + k) - Knot(i + 1);
      If t <> 0 Then
        v := v + (Knot(i + k) - u) * NBlend(i + 1, k - 1, u) / t
    End;
  NBlend := v
End;

{--------------------------------------------------------------------------}
Procedure BSpline (Var x, y, z: double; u: double; n, k: integer; p: data);

Var 
  i: integer;
  b: double;

Begin

  knotK := k;
  knotN := n;

  x := 0;
  y := 0;
  z := 0;

  For i := 0 To n Do
    Begin
      b := NBlend(i, k, u);
      x := x + p[i].X[1] * b;
      y := y + p[i].X[2] * b;
      z := z + p[i].X[3] * b;
    End;
End;

{--------------------------------------------------------------------------}
{ Bisection....}
{--------------------------------------------------------------------------}
Function AbstandPE (u: double): double;

Var 
  d: double;
  a: Vektor;

Begin
  BSpline(a[1], a[2], a[3], u, cinP, kk, InKP);
  d := Skalarprodukt(a, Mb);
  AbstandPE := d - Db;
End;

{--------------------------------------------------------------------------}
Procedure ZBrac (
{Function fx (z: double): double; AbstandPE}
x1, x2: double;
n: integer;

Var xb1, xb2: double);

Var 
  nbb, i: integer;
  x, fp, fc, dx: double;

Begin
  solution := false;
  nbb := 0;
  dx := (x2 - x1) / n;
  fp := AbstandPE(x1);
  x := x1;
  For i := 1 To (n - 1) Do
    Begin
      x := x + dx;
      fc := AbstandPE(x);
      If (fc * fp < 0.0) Then
        Begin
          nbb := nbb + 1;
          If nbb = 1 Then
            Begin
              xb1 := x - dx;
              xb2 := x;
            End
        End;
      fp := fc;
    End;

  Case nbb Of 
    1:
       solution := true;
    0:
    ;
    otherwise
    solution := true;
  End;

End;
{--------------------------------------------------------------------------}

Function RtBis (
{Function func (z: double): double; AbstandPE}
x1, x2, xacc: double): double;

Label 
  10;

Const 
  JMAX = 280;

Var 
  j: integer;
  dx, f, fmid, xmid, rtb: double;

Begin

  f := AbstandPE(x1);
  fmid := AbstandPE(x2);
  If (f * fmid >= 0.0) Then
    writeln('Error: Root must be bracketed for bisection in rtbis !!!!');
  If f < 0.0 Then
    Begin
      rtb := x1;
      dx := x2 - x1;
    End
  Else
    Begin
      rtb := x2;
      dx := x1 - x2;
    End;
  For j := 1 To JMAX Do
    Begin
      dx := dx * 0.5;
      xmid := rtb + dx;
      fmid := AbstandPE(xmid);
      If (fmid <= 0.0) Then
        rtb := xmid;
      If (abs(dx) <= xacc) Or (fmid = 0.0) Then
        Begin
          RtBis := rtb;
          goto 10;
        End
    End;
  writeln(' Error: too many besection points in rtbis!!!!! rtbis := 0.0');
  rtbis := 0.0;
  10:
End;
{--------------------------------------------------------------------------}
Procedure DurchStoss (Var x, y, z, u: double);

Var 
  xb1, xb2: double;

Begin

  ZBrac(0, (cinP - kk + 2), split, xb1, xb2);
  If solution Then
    Begin
      u := RtBis(xb1, xb2, xacc);
      BSpline(x, y, z, u, cinP, kk, InKP);
    End;
End;

{--------------------------------------------------------------------------}
{ BSpline Ableitung ...}
{--------------------------------------------------------------------------}
Function GNBlend (i, k: integer; u: double): double;

Var 
  t: integer;
  v: double;

Begin
  If k = 1 Then
    Begin
      v := 0;
      If (Knot(i) <= u) And (u < Knot(i + 1)) Then
        v := 1
    End
  Else
    Begin
      v := 0;
      t := Knot(i + k - 1) - Knot(i);
      If t <> 0 Then
        v := NBlend(i, k - 1, u) / t;
      t := Knot(i + k) - Knot(i + 1);
      If t <> 0 Then
        v := v - NBlend(i + 1, k - 1, u) / t
    End;
  GNBlend := (k - 1) * v
End;

{--------------------------------------------------------------------------}
Procedure Grad_BSpline (Var x, y, z: double; u: double; n, k: integer; p: data);

Var 
  i: integer;
  b: double;

Begin

  knotK := k;
  knotN := n;

  x := 0;
  y := 0;
  z := 0;

  For i := 0 To n Do
    Begin
      b := GNBlend(i, k, u);
      x := x + p[i].X[1] * b;
      y := y + p[i].X[2] * b;
      z := z + p[i].X[3] * b;
    End;
End;

{--------------------------------------------------------------------------}

Procedure norm (Var a: Vektor);
{Einheitsvektor wird zugeordnet}

Var 
  i: integer;
  da: double;

Begin
  da := 0;
  For i := 1 To 3 Do
    da := da + sqr(a[i]);
  da := sqrt(da);
  For i := 1 To 3 Do
    a[i] := a[i] / da;
End;

{--------------------------------------------------------------------------}

Procedure torsion (a, b, c, d: Vektor; Var tau: double);
{Angegeben 4 Punkte wird torsionswinkel berechnet}

Var 
  v1, v2, v3, vp1, vp2: Vektor;
  s, si, co, t1, t2: double;

Begin

  SetVektor(v1, a, b);
  SetVektor(v2, b, c);
  SetVektor(v3, c, d);

  norm(v1);
  norm(v2);
  norm(v3);

  VektorProdukt(vp1, v1, v2);
  VektorProdukt(vp2, v2, v3);
  t1 := Skalarprodukt(v1, v2);
  t2 := Skalarprodukt(v2, v3);
  s := sqrt((1 - sqr(t1)) * (1 - sqr(t2)));
  co := Skalarprodukt(vp1, vp2);
  si := Skalarprodukt(v1, vp2);
  tau := atg(si / s, co / s);
End;

{--------------------------------------------------------------------------}
{----------------------- H P --------------------}

Begin

  readln(kk);

  DataRead('input1',InKP, cinP, 0);
  DataRead('input2',InKB, cinB, 1);

  Best_Plane(cinB, InKB, Mb, SS, Db);

  Durchstoss(PBx, PBy, PBz, u);
  If solution Then
    Begin

      PB[1] := PBx;
      PB[2] := PBy;
      PB[3] := PBz;

      Grad_BSpline(PAx, PAy, PAz, u, cinP, kk, InKP);

      PP[1] := PAx;
      PP[2] := PAy;
      PP[3] := PAz;

      For i := 1 To 3 Do
        T[i] := PB[i] + PP[i];

      SetVektor(B1, PB, SS);

      SetPlane(Ma, Da, B1, Mb, PB);

      Projection(PA, Ma, Da, T, Ma);

      SetVektor(PbPa, PB, PA);

      SetVektor(Junk, SS, InKP[1].x);

      If Skalarprodukt(Mb, PbPa) < 0 Then
        Begin
          For i := 1 To 3 Do
            Mb[i] := Mb[i] * (-1);
          Db := -Db;
        End;

      For i := 1 To 3 Do
        B[i] := SS[i] + 5 * Mb[i];

      Winkel(w, Mb, PbPa);

      If Skalarprodukt(B1, PbPa) < 0 Then
        w := -w;

      writeln('Scher  ', w : 10 : 2);

      torsion(T, PB, SS, B, tau);

      writeln('tau  ', tau : 10 : 2);

    End
  Else
    Begin
      writeln('f no_sol');
      writeln('f no_sol');
    End

End.
