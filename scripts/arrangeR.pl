#!/usr/bin/perl
#

 &transform($ARGV[0],$ARGV[1]);

 exit;

sub transform
{

 local ($inp,$out)=@_;

 open (INPUT,"<$inp") || die ("Can't open input file:  $_\n");
 open (OUTPUT,">$out");
  $n=-1; 
  LINE: while (<INPUT>){
	@fields = split;
      if ($fields[0] eq "ATOM" || $fields[0] eq "HETATM"){
	$n++;
        $type[$n]=$fields[0];
	$type[$n]=~s/\s//g;
 	$atnum[$n]=substr($_,7,6);
	$atnum[$n]=~s/\s//g;
	$atomn[$n]=substr($_,12,5);
	$atomn[$n]=~s/\s//g;
	$res[$n]=substr($_,17,4);
	$res[$n]=~s/\s//g;
	$chain[$n]=substr($_,22,2);
	$chain[$n]=~s/\s//g;
	$resnum[$n]=substr($_,23,7);
	$resnum[$n]=~s/\s//g;
	$coordx[$n]=substr($_,30,8);
	$coordx[$n]=~s/\s//g;
	$coordy[$n]=substr($_,38,8);
	$coordy[$n]=~s/\s//g;
	$coordz[$n]=substr($_,46,8);
	$coordz[$n]=~s/\s//g;
      }else{next LINE;};
  };

  $total = $n;

  for $i (0..$total) {
      if ($i==0) 
      {
       $j=0;
       $check=0;
       $atom = {};
       $atom = {
                type    => $type[0],
                atnum   => $atnum[0],
                atom    => $atomn[0],
                residue => $res[0],
                chain   => $chain[0],
                resnum  => $resnum[0],
                x       => $coordx[0],
                y       => $coordy[0],
                z       => $coordz[0],
                };
       push @atom,$atom;
       if ($atomn[0] eq 'CA') {$check++;};
       next;
      }elsif ($resnum[$i] eq $resnum[$i-1]){
       $atom = {};
       $atom = {
                type    => $type[$i],
                atnum   => $atnum[$i],
                atom    => $atomn[$i],
                residue => $res[$i],
                chain   => $chain[$i],
                resnum  => $resnum[$i],
                x       => $coordx[$i],
                y       => $coordy[$i],
                z       => $coordz[$i],
                };
       push @atom,$atom;
       if ($atomn[$i] eq 'CA') {$check++;};
       next;
      }else{
       $residue[$j]= [@atom];
       undef @atom;
       $j++;
       $atom = {};
       $atom = {
                type    => $type[$i],
                atnum   => $atnum[$i],
                atom    => $atomn[$i],
                residue => $res[$i],
                chain   => $chain[$i],
                resnum  => $resnum[$i],
                x       => $coordx[$i],
                y       => $coordy[$i],
                z       => $coordz[$i],
                };
       push @atom,$atom;
       if ($atomn[$i] eq 'CA') {$check++;};
       next;
      };
  };
      $residue[$j]= [@atom];
      undef @atom;
      $j++;

#  die "Different number of CA ($check) than Residues ($j)" unless ($j==$check);
    
  (@residue)=resorder(@residue);

  $outflush=select(OUTPUT);
  $~ = "PDB_FORMAT";
  select($outflush);
  
  $na=0;
  for $n (0..$#residue)
  {
   for $m (0..$#{$residue[$n]})
   {
    $na++;
    $typew   =  $residue[$n][$m]{type};
    $atnumw  =  $na;
    $atomw   =  $residue[$n][$m]{atom};
    $resw    =  $residue[$n][$m]{residue};
    $chainw  =  " ";
    $resnumw =  $n + 1;
    $coordxw =  $residue[$n][$m]{x};
    $coordyw =  $residue[$n][$m]{y};
    $coordzw =  $residue[$n][$m]{z};
    write OUTPUT;
   };
  };   
 
  close (INPUT);
  close (OUTPUT);

             
format PDB_FORMAT=
@<<<<<<@>>>>  @<<<@<< @@>>>    @>>>>>>> @>>>>>>> @>>>>>>>
$typew,$atnumw,$atomw,$resw,$chainw,$resnumw,$coordxw,$coordyw,$coordzw
.

 
}

sub resorder
{
  local (@residue)=@_;

  $n=1;
  for $j (0..@residue)
  {
   for $i (0..$#{$residue[$j]})
   {
    if ( $residue[$j][$i]{atom} eq 'N' ){
     $atom_new={};
     $atom_new=$residue[$j][$i];
     $atom_new{atnum}=$n;
     $n++;
     push @atom_new,$atom_new;
     last;
    }
   }
   for $i (0..$#{$residue[$j]})
   {
    if ( $residue[$j][$i]{atom} eq 'H' ){
     $atom_new={};
     $atom_new=$residue[$j][$i];
     $atom_new{atnum}=$n;
     $n++;
     push @atom_new,$atom_new;
     last;
    }
   }
   for $i (0..$#{$residue[$j]})
   {
    if ($residue[$j][$i]{atom} eq 'CA' ){
     $atom_new={};
     $atom_new=$residue[$j][$i];
     $atom_new{atnum}=$n;
     $n++;
     push @atom_new,$atom_new;
     last;
    }
   }
   for $i (0..$#{$residue[$j]})
   {
    if (   $residue[$j][$i]{atom} ne 'N' 
        && $residue[$j][$i]{atom} ne 'CA' 
        && $residue[$j][$i]{atom} ne 'C'
        && $residue[$j][$i]{atom} ne 'O'
        && $residue[$j][$i]{atom} ne 'H'){
                                          $atom_new={};
                                          $atom_new=$residue[$j][$i];
                                          $atom_new{atnum}=$n;
                                          $n++;
                                          push @atom_new,$atom_new;
                                          next; };
   };
   for $i (0..$#{$residue[$j]})
   {
    if ($residue[$j][$i]{atom} eq 'C' ){
     $atom_new={};
     $atom_new=$residue[$j][$i];
     $atom_new{atnum}=$n;
     $n++;
     push @atom_new,$atom_new;
     last;
    }
   }
   for $i (0..$#{$residue[$j]})
   {
    if ($residue[$j][$i]{atom} eq 'O' ){
     $atom_new={};
     $atom_new=$residue[$j][$i];
     $atom_new{atnum}=$n;
     $n++;
     push @atom_new,$atom_new;
     last;
    }
   }
   $residue_new[$j]=[@atom_new];
   undef @atom_new;
  }

  return @residue_new;
}
    
       
        
  

    
     
    
