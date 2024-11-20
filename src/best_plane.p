
Program bestplane (input, output,f1);

{PL iii/95}
{Berechnet beste Ebene}
{=========================================}
{}
{INPUT  : File1: Atomkoordinaten in PDB Format}
{OUTPUT : rms Abstand der Atome von der besten Ebene}
{}
{Aufrufen des Programmes mit folgender shell}
{-------------------------------------------}
{cp   BP_A1.pdb  input1  # Input filename 1}
{}
{a.out                     # Programmname}
{}
{/bin/rm input1            # Aufrauemen}
{-------------------------------------------}
{---Definitionen}
{Mb		     Normale zu beta}
{Db			 Abstand beta Ursprung}
{Sb		     Schwerpunkt der basenatome in beta}
{cinB		 Anzahl Basenatomen B}
{InkB		 Basenatomenkoordinaten: PDB Format}
{***********************************************}

Const 
  MaxAtom = 500;
  {MAXSTRLEN = 100000;}
  {MAXSTRLEN = 16#50;}

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
  ch: char;
  InKB: data;
  Mb, Sb: Vektor;
  Db: double;
  cinB: integer;
{--------------------------------------------------------------------------}

Procedure DataRead(inf: string; Var inp: data; Var cin: integer);

Var 
  i, f: integer;
  f1: TextFile;


Begin

  cin := 1;

  AssignFile(f1, inf);
  reset(f1);
  While ((Not eof(f1)) And (cin <= MaxAtom)) Do
    Begin
      For i := 1 To 13 Do
        read(f1, ch);
      For f := 1 To 5 Do
        Begin
          read(f1, inp[cin].name[f]);
        End;
      For i := 1 To 11 Do
        read(f1, ch);
      With inp[cin] Do
        Begin
          readln(f1, x[1], x[2], x[3]);
        End;
      cin  := cin  + 1;
    End;

  cin  := cin  - 1;

  close(f1);
End;

{--------------------------------------------------------------------------}
Procedure Eigenvalues (Var z, v: Matrix);

Const 
  kRest = 1E-12;
  kEvSum = 2.99999;  {if sum |v[i,i]| > kEvSum, set v=identity}

Var 
  i, j, k: integer;
  th, ta, t, c, s, zz, sum:  double;

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
Procedure Best_Plane (NPoint: integer; X: Data; Var M, SS: Vektor; Var D:  double);

Label 
  10;

Var 
  i, j, k: integer;
  Y: array[1..MaxAtom] Of Vektor;
  c: Vektor;
  z, v: Matrix;
  Sig,eval:  double;

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
{Kleinste Eigenvalue is the value of the sum of residuals!!!}
  eval := z[3, 3];
  Sig := sqrt(eval / (NPoint - 3));

  Writeln('rms_D  ', Sig : 10 : 4);
10:
End;

{--------------------------------------------------------------------------}

{----------------------- H P --------------------}

Begin

  DataRead('input1',InkB,cinB);
  Best_Plane(cinB, InKB, Mb, Sb, Db);

End.
