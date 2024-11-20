
Program BPSPabstand (input, output, f1);

{Version 1.0 Februar 1995}
{Angelegt wird eine BSpline Kurve durch die input Atome}
{Beachte kk Parameter}
{Anschliessend wird der Abstand von Atomen zur Kurve fuer}
{jedes Atom berechnet}

{INPUT:		File in PDB Format mit Atomen wodurch die BSpline gehen soll}
{			File NAme input1}

{OUTPUT:	Liste der Atome mit Abstand zur BSpline sowie Mittelwert und}
{			Standardabweichung}

{Aufrufen des Programmes mit folgender shell}
{-------------------------------------------}
{cp   BP_A1.pdb  input1  # Input filename }
{}
{a.out                     # Programmname}
{3                         # BSpline Konstante}
{}
{/bin/rm input1            # Aufrauemen}
{-------------------------------------------}
{@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@}

Const 
  MaxAtom = 500;
  {MAXSTRLEN = 16#50;}
  {MAXSTRLEN = 100000;}
  maxP = 500;{Maximale Anzahl P Atome}
  tol = 1.0E-10;

Type 

  {stringArray = packed array [1..MAXSTRLEN] Of char;}

  pa6 = packed array[1..20] Of char;

  Vektor = array[1..3] Of double;

  Koordinaten = Record
    name: pa6;
    x: Vektor;
  End;

  data = array[0..MaxAtom] Of Koordinaten;

Var 

  InKP: data;
  cinP: integer;
  knotK, knotN: integer; {BSpline Stuff...}
  kk : integer;{BSpline dimension}

{--------------------------------------------------------------------------}

Procedure DataRead (inf: string; Var inp: data; Var cin: integer; start: integer);

Var 
  i, f: integer;
  ch: char;
  f1: TextFile;

Begin

  {https://smartpascal.github.io/help/assets/assignfile.htm}
  {AssignFile attaches a FileHandler from the filename `inf`, where -}
  {`f1` as a variable will be handled accordingly}
  AssignFile(f1, inf);
  {https://smartpascal.github.io/help/assets/reset.htm}
  {Opens the actually text file from the handler. Overrides the variable `f1`}
  reset(f1);

  cin := start;

  While ((Not eof(f1)) And (cin <= MaxAtom)) Do
    Begin
      For i := 1 To 8 Do
        read(f1, ch);
      For f := 1 To 20 Do
        read(f1, Inp[cin].name[f]);
      With inp[cin] Do
        readln(f1, x[1], x[2], x[3]);
      cin := cin + 1;
    End;

  cin := cin - 1;

  close(f1);

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

Function AbstandPBS (u: double; b: vektor): double;

Var 
  a: Vektor;

Begin
  BSpline(a[1], a[2], a[3], u, cinP, kk, InKP);

  AbstandPBS := sqrt(sqr(b[1] - a[1]) + sqr(b[2] - a[2]) + sqr(b[3] - a[3]));
End;

{--------------------------------------------------------------------------}

Procedure SHIFT2 (Var a, b, c: double);

Begin
  a := b;
  b := c;
End;

{--------------------------------------------------------------------------}

Procedure SHIFT3 (Var a, b, c, d: double);

Begin
  a := b;
  b := c;
  c := d;
End;

{--------------------------------------------------------------------------}

Procedure mnbrak (Var ax, bx, cx, fa, fb, fc: double;
    {Function func (z: double; v: vektor): double;}
    P: vektor);
{function with specific signature is passed here}
{as a typed argument, without naming the specific function}
Const 
  STEP = 0.2;

Begin
  fa := AbstandPBS(ax, P);
  bx := bx + 0.2;
  fb := AbstandPBS(bx, P);
  cx := bx + step; {First guess for c}
  fc := AbstandPBS(cx, P);
  While (fb > fc) Do  {Keep returning here until we bracket.}
    Begin
      ax := bx;
      fa := AbstandPBS(ax, P);
      bx := cx;
      fb := AbstandPBS(bx, P);
      cx := cx + step;
      fc := AbstandPBS(cx, P);
    End;
End;

{--------------------------------------------------------------------------}

Function golden (ax, bx, cx: double;
    {Function f (z: double; v: vektor): double;}
tol: double;

Var xmin: double;
  P: vektor): double;



{W.H.Press, S.A.Teukolsky, W.T.Vetterling, B.P.Flannery, NUMERICAL RECIPES IN C THE ART OF THE COMPUTING}
{second edition, Cambridge University Press 1992, pp 401-402}
{Routine for golden section search}

{Given a function f, and given a braketing triplet of abscissas ax, bx, cx (such that bx is between ax and cx, and}

{f(bx) is less than both f(ax) and f(cx)), this routine performs a golden section search for the minimum, isolating}

{it to a fractional precision of about tol. The abscissa of the minimum is returned as xmin, and the minimum}
{function value is returned as golden, the returned function value.}

{AbstandPBS}
Const 
  R = 0.61803399;
  C = 0.38196601; {1.0 - R}

Var 
  f1, f2, x0, x1, x2, x3: double;
  dd, ee: double;

Begin
  x0 := ax;
  x3 := cx;
  If (abs(cx - bx) > abs(bx - ax)) Then
    Begin
      x1 := bx;
      x2 := bx + C * (cx - bx);
    End
  Else
    Begin
      x2 := bx;
      x1 := bx - C * (bx - ax);
    End;
  f1 := AbstandPBS(x1, P);
  f2 := AbstandPBS(x2, P);
  While ((abs(x3 - x0)) > (tol * (abs(x1) + abs(x2)))) Do
    Begin
      If (f2 < f1) Then
        Begin
          dd := R * x1 + C * x3;
          SHIFT3(x0, x1, x2, dd);
          ee := AbstandPBS(x2, P);
          SHIFT2(f1, f2, ee);
        End
      Else
        Begin
          dd := R * x2 + C * x0;
          SHIFT3(x3, x2, x1, dd);
          ee := AbstandPBS(x1, P);
          SHIFT2(f2, f1, ee);
        End;
    End;
  If (f1 < f2) Then
    Begin
      xmin := x1;
      golden := f1;
    End
  Else
    Begin
      xmin := x2;
      golden := f2;
    End;
End;

{--------------------------------------------------------------------------}

Procedure FitBSP (Max: integer);

Var 
  a, i,j: integer;
  aa, bb, cc, ff, gg, hh, xmin, mean, sig, varianz: double;
  dd: array[1..maxP] Of double;
  rr: vektor;

Begin
  a := 0;
  aa := 0.1;
  bb := 0.2;
  cc := 0.0;
  ff := 0.0;
  gg := 0.0;
  hh := 0.0;
  For i := 1 To (Max - 1) Do
    Begin
      a := a + 1;
      {mnbrak(aa, bb, cc, ff, gg, hh, AbstandPBS, InkP[i].x);}
      mnbrak(aa, bb, cc, ff, gg, hh, InkP[i].x);
      dd[a] := golden(aa, bb, cc, tol, xmin, InkP[i].x);
      {dd[a] := golden(aa, bb, cc, AbstandPBS, tol, xmin, InkP[i].x);}
      BSpline(rr[1], rr[2], rr[3], xmin, cinP, kk, InKP);

      write('Abstand:  ');
      For j:=1 To 20 Do
        write(InKP[i].name[j]);
      writeln('   ',dd[a] : 10 : 4);

    End;

{--- statistik ---}

  mean := 0.0;
  varianz := 0.0;
  For i := 1 To (Max - 1) Do
    mean := mean + dd[i];
  mean := mean / (Max - 1);
  For i := 1 To (Max - 1) Do
    varianz := varianz + sqr(dd[i] - mean);
  sig := sqrt((1 / (Max - 1)) * varianz);
  writeln('-----------------------------------------------------');
  writeln('Mittlerer Abstand :                  ', mean : 6 : 4);
  writeln('Standardabweichung:                  ', sig : 6 : 4);
  writeln;
End;

{--------------------------------------------------------------------------}
{----------------------- H P --------------------}

Begin
  readln(kk);
  DataRead('input1',InKP, cinP, 0);

  FitBSP(cinP);

End.
