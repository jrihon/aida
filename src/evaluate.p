
Program evaluate (input, output, f1);

{Version 1.0 Maerz 1995}

Const 
  maxin = 500;
  {MAXSTRLEN = 1650;}
  {MAXSTRLEN = 100000;}

Type 

  {stringArray = packed array [1..MAXSTRLEN] Of char;}
  pa7 = packed array[1..7] Of char;
  data24 = array[1..24] Of double;

  Base = Record
    B1, B2: pa7;
    BSP: data24;
  End;
  BBase = array[1..maxin] Of Base;

Var 
  inD: BBase;
  mBSP: data24;
  vBSP: data24;
  f1: TextFile;

{--------------------------------------------------------------------------}

Procedure statistik (Var mean, rms: double; data: BBase; max: integer; ref: integer);

Var 
  i, maxN: integer;
  varianz: double;

Begin
  maxN := 0;
  mean := 0.0;
  varianz := 0.0;
  For i := 1 To max Do
    Begin
      If (data[i].BSP[ref]) <> -1000.0 Then
        Begin
          mean := mean + data[i].BSP[ref];
          maxN := maxN + 1;
        End;
    End;
  mean := mean / (maxN);
  For i := 1 To max Do
    Begin
      If (data[i].BSP[ref] <> -1000.0) Then
        varianz := varianz + sqr(data[i].BSP[ref] - mean);
    End;
  rms := sqrt((1 / (maxN - 1)) * varianz);
End;

{--------------------------------------------------------------------------}

Procedure DataProcess(inf: string);

Var 
  i, a, j, k: integer;
  ch: char;


Begin
  AssignFile(f1, inf);
  reset(f1);
  a := 0;

  Repeat
    a := a + 1;
    With inD[a] Do
      Begin

        For j := 1 To 7 Do
          read(f1, B1[j]);
        readln(f1);
        For j := 1 To 7 Do
          read(f1, B2[j]);
        readln(f1);

        For j := 1 To 24 Do
          Begin
            read(f1, ch);
            If (ch <> 'f') Then
              Begin
                For k := 1 To 7 Do
                  read(f1, ch);
                readln(f1, BSP[j]);
              End
            Else
              Begin
                readln(f1);
                BSP[j] := -1000.0;
              End;
          End;
      End;
  Until eof(f1);
  close(f1);
{statistik ...}
  For i := 1 To 24 Do
    statistik(mBSP[i], vBSP[i], inD, a, i);


{Ausschreiben ...}

  writeln(  'Basenebene');
  writeln(  ' ');
  writeln(  'Basenpaar     rms AB   rms A    rms B   Tau');
  writeln(  '---------------------------------------------');
  For k := 1 To a Do
    With inD[k] Do
      Begin
        For i := 1 To 4 Do
          write(  B1[i]);
        write(  '- ');
        For i := 1 To 4 Do
          write(  B2[i]);
        write(  BSP[1] : 8 : 2, ' ');
        write(  BSP[2] : 8 : 2, ' ');
        write(  BSP[3] : 8 : 2, ' ');
        write(  BSP[4] : 8 : 2, ' ');
        writeln;
      End;
  writeln(  '---------------------------------------------');
  write(  'Mean      ');
  writeln(  mBSP[1] : 8 : 2, ' ', mBSP[2] : 8 : 2, ' ', mBSP[3] : 8 : 2, ' ', mBSP[4] : 8 : 2);
  write(  'Sig       ');
  writeln(  vBSP[1] : 8 : 2, ' ', vBSP[2] : 8 : 2, ' ', vBSP[3] : 8 : 2, ' ', vBSP[4] : 8 : 2);
  writeln;
  writeln;
  writeln;
{@@@@@@@@@@@@@@@@@@@@@@@@@}
  writeln(  'BSpline Scherungswinkeln');
  writeln(  ' ');
  writeln(  'Basenpaar     AB1       A1       B1      AB2       A2       B2     Tau1     Tau2');
  writeln(  '---------------------------------------------------------------------------------');
  For k := 1 To a Do
    With inD[k] Do
      Begin
        For i := 1 To 4 Do
          write(  B1[i]);
        write(  '- ');
        For i := 1 To 4 Do
          write(  B2[i]);
        If (BSP[5] <> -1000.0) Then
          write(  BSP[5] : 8 : 2, ' ')
        Else
          write(  '    -    ');
        If (BSP[7] <> -1000.0) Then
          write(  BSP[7] : 8 : 2, ' ')
        Else
          write(  '    -    ');
        If (BSP[9] <> -1000.0) Then
          write(  BSP[9] : 8 : 2, ' ')
        Else
          write(  '    -    ');
        If (BSP[11] <> -1000.0) Then
          write(  BSP[11] : 8 : 2, ' ')
        Else
          write(  '    -    ');
        If (BSP[13] <> -1000.0) Then
          write(  BSP[13] : 8 : 2, ' ')
        Else
          write(  '    -    ');
        If (BSP[15] <> -1000.0) Then
          write(  BSP[15] : 8 : 2, ' ')
        Else
          write(  '    -    ');
        If (BSP[6] <> -1000.0) Then
          write(  BSP[6] : 8 : 2, ' ')
        Else
          write(  '    -    ');
        If (BSP[12] <> -1000.0) Then
          write(  BSP[12] : 8 : 2, ' ')
        Else
          write(  '    -    ');
        writeln;
      End;
  writeln(  '---------------------------------------------------------------------------------');
  write(  'Mean      ');
  writeln(  mBSP[5] : 8 : 2, ' ', mBSP[7] : 8 : 2, ' ', mBSP[9] : 8 : 2, ' ', mBSP[11] : 8 : 2, ' ',
          mBSP[13] : 8 : 2, ' ', mBSP[15] : 8 : 2, ' ', mBSP[6] : 8 : 2, ' ', mBSP[12] : 8 : 2);
  write(  'Sig       ');
  writeln(  vBSP[5] : 8 : 2, ' ', vBSP[7] : 8 : 2, ' ', vBSP[9] : 8 : 2, ' ', vBSP[11] : 8 : 2, ' ',
          vBSP[13] : 8 : 2, ' ', vBSP[15] : 8 : 2, ' ', vBSP[6] : 8 : 2, ' ', vBSP[12] : 8 : 2);
  writeln;
  writeln;
  writeln;
{@@@@@@@@@@@@@@@@@}
  writeln(  'Backbone Scherungswinkel');
  writeln(  ' ');
  writeln(  'Basenpaar     sA1      sAB1     sB2     sAB2     pp1AB    pp2AB');
  writeln(  '---------------------------------------------------------------');
  For k := 1 To a Do
    With inD[k] Do
      Begin
        For i := 1 To 4 Do
          write(  B1[i]);
        write(  '- ');
        For i := 1 To 4 Do
          write(  B2[i]);
        write(  BSP[17] : 8 : 2, ' ');
        write(  BSP[19] : 8 : 2, ' ');
        write(  BSP[21] : 8 : 2, ' ');
        write(  BSP[23] : 8 : 2, ' ');
        write(  BSP[20] : 8 : 2, ' ');
        write(  BSP[24] : 8 : 2, ' ');
        writeln;
      End;
  writeln(  '---------------------------------------------------------------');
  write(  'Mean      ');
  writeln(  mBSP[17] : 8 : 2, ' ', mBSP[19] : 8 : 2, ' ', mBSP[21] : 8 : 2, ' ', mBSP[23] : 8 : 2,
          ' ', mBSP[20] : 8 : 2, ' ', mBSP[24] : 8 : 2);
  write(  'Sig       ');
  writeln(  vBSP[17] : 8 : 2, ' ', vBSP[19] : 8 : 2, ' ', vBSP[21] : 8 : 2, ' ', vBSP[23] : 8 : 2,
          ' ', vBSP[20] : 8 : 2, ' ', vBSP[24] : 8 : 2);
End;

{--------------------------------------------------------------------------}

Begin
  DataProcess('input1');
End.
