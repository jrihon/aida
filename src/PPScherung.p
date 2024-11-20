
Program ScherungPP (input, output, f1);

{Version 1.0 Februar 1995}

{Sherungswinkel ohne BSpline Kurve und Tangentialvektor}

{-------------------------------------------}
{INPUT:		File 1 Mit Basenatome und die 2 P Atome vom Rest A und A+1 am ENDE}
{			des Files}

{OUTPUT:	P-Scherungswinkel und Basennormale_PP_Winkel}

{Aufrufen des Programmes mit folgender shell}
{-------------------------------------------}
{cp   BP.pdb      input1  # Input filename }
{}
{a.out                      # Programmname}
{}
{/bin/rm input1             # Aufrauemen}
{-------------------------------------------}

{---Definitionen}

{Mb		Normale zu beta}
{SS	Schwerpunkt der basenatome in beta}
{PP	P5' --> P3' Vektor}
{PB	Durchstosspunkt von P5' --> P3' Vektor mit beta}
{B1	PB --> SS Vektor}
{Ma	Normale zur Ebene alfa}
{PA	Projektion von P3' auf alfa}

Const 
  MaxAtom = 500;
  {MAXSTRLEN = 1650;}
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

  data = array[1..MaxAtom] Of Koordinaten;

Var 

  InK: data;
  Mb, SS, PP, PB, PA, B1, PbPa, Ma: Vektor;
  Db, Da: double;
  cin: integer;
  w: double;
  i: integer;


 {--------------------------------------------------------------------------}

Procedure DataRead (inf: string; Var inp: data; Var cin: integer);

Var 
  i, f: integer;
  ch : char;
  f1: TextFile;

Begin
  AssignFile(f1, inf);
  reset(f1);

{------- pdb Format -------}

  cin := 1;

  While ((Not eof(f1)) And (cin <= MaxAtom)) Do
    Begin
      For i := 1 To 13 Do
        read(f1, ch);
      For f := 1 To 5 Do
        read(f1, Inp[cin].name[f]);
      For i := 1 To 11 Do
        read(f1, ch);
      With inp[cin] Do
        Begin
          readln(f1, x[1], x[2], x[3]);
        End;
      cin := cin + 1;
    End;

  cin := cin - 1;

  close(f1);

End;

{--------------------------------------------------------------------------}
{ Math... }
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

Procedure SetVektor (Var ab: Vektor; a, b: Vektor);

Var 
  i: integer;

Begin
  For i := 1 To 3 Do
    ab[i] := b[i] - a[i];
End;

{--------------------------------------------------------------------------}

Procedure Winkel (Var w: double; a, b: Vektor);

Var 
  m, n, mn, sinw, cosw: double;

{-------------}

Function atg (s, c: double): double;
     (********************)

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

{-------------}

Begin

  m := sqrt(sqr(a[1]) + sqr(a[2]) + sqr(a[3]));
  n := sqrt(sqr(b[1]) + sqr(b[2]) + sqr(b[3]));
  mn := (a[1]) * (b[1]) + (a[2]) * (b[2]) + (a[3]) * (b[3]);

  cosw := (mn) / (m * n);
  sinw := (sqrt(1 - sqr(cosw)));
  w := atg(sinw, cosw);

End;

{--------------------------------------------------------------------------}

Procedure DurchStoss (Var p: Vektor; a: Vektor; b: double; c, d: Vektor);

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

Procedure VektorProdukt (Var c: Vektor; a, b: Vektor);

Begin
  c[1] := a[2] * b[3] - a[3] * b[2];
  c[2] := a[3] * b[1] - a[1] * b[3];
  c[3] := a[1] * b[2] - a[2] * b[1]
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
{ Best Plane ...}
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
  Mmod: double;

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
{----------------------- H P --------------------}

Begin

  DataRead('input1',InK, cin);

  Best_Plane(cin - 2, InK, Mb, SS, Db);

  SetVektor(PP, InK[cin - 1].x, InK[cin].x);

  If Skalarprodukt(Mb, PP) < 0 Then
    Begin
      For i := 1 To 3 Do
        Mb[i] := Mb[i] * (-1);
      Db := -Db;
    End;

  Durchstoss(PB, Mb, Db, InK[cin - 1].x, PP);

  SetVektor(B1, PB, SS);

  SetPlane(Ma, Da, B1, Mb, PB);

  Durchstoss(PA, Ma, Da, InK[cin].x, Ma);

  SetVektor(PbPa, PB, PA);

 {---- Homo-DNS Fall...---}
  If Skalarprodukt(Mb, PbPa) < 0 Then
    SetVektor(PbPa, PA, PB);

  Winkel(w, Mb, PbPa);

  If Skalarprodukt(B1, PbPa) < 0 Then
    w := -w;

  writeln('PP-Sch   ', w : 10 : 2);

  Winkel(w, Mb, PP);

  writeln('N_PP_w  ', w : 10 : 2);

End.
