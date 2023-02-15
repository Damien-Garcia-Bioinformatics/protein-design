##########################    In script documentation	##########################
=head

Created by Iyanar VETRIVEL in March 2014.
Modified by Iyanar VETRIVEL in April 2014 to change the input and output method.

Description : This script reads chain seperated dssp files and extracts the torsion angles, calculates the RMSDA and assigns PBs.Output is written to PDB_chainID.dssp.torsion, PDB_chainID.dssp.rmsda and PDB_chainID.dssp.pb files respectively.

Usage : example- perl dssp_to_pb_tor_rmsda.pl path/input_file path/output_directory.

=cut
##################################################################################


#!/usr/bin/perl
use File::Basename;
chomp($ARGV[0]);#input dssp file
chomp($ARGV[1]);#output directory where you want the PDB_chainID.dssp.pb,PDB_chainID.dssp.rmsda and PDB_chainID.dssp.torsion file to be created.

($name,$path,$suffix) = fileparse($ARGV[0],qr/\.[^.]*$/);

unless(open(FIN,"<","$ARGV[0]"))
{
	print"Cannot open \"$name$suffix\"\n";
	exit;
}

unless(open(FOUT1,">","$ARGV[1]/$name$suffix.pb"))
{
	print"Cannot open \"$name$suffix.pb\"\n";
	exit;
}
unless(open(FOUT2,">","$ARGV[1]/$name$suffix.rmsda"))
{
	print"Cannot open \"$name$suffix.rmsda\"\n";
	exit;
}
unless(open(FOUT3,">","$ARGV[1]/$name$suffix.torsion"))
{
	print"Cannot open \"$name$suffix.torsion\"\n";
	exit;
}
@dssp_dump=<FIN>;

foreach $line (@dssp_dump)
{
	unless(($line=~ m/\.$/)||($line=~ m/^(\s*#)/))
	{
		push(@dssp_dump_sans_header,$line);
	}
}


#2D array for the torsions 

for($l=0;$l<@dssp_dump_sans_header;$l++)
{
$torsions[$l][0]=substr($dssp_dump_sans_header[$l],103,6);
$torsions[$l][1]=substr($dssp_dump_sans_header[$l],109,6);
}



print FOUT1"ZZ";
print FOUT2"X\nX\n";
for($i=2;$i<@dssp_dump_sans_header-2;$i++)
{
	@eight_torsions=("$torsions[$i-2][1]","$torsions[$i-1][0]","$torsions[$i-1][1]","$torsions[$i][0]","$torsions[$i][1]","$torsions[$i+1][0]","$torsions[$i+1][1]","$torsions[$i+2][0]");	
	#If '!' is found, add exactly the same number of 'X' as many residues missing	
	if(substr($dssp_dump_sans_header[$i],13,1) eq '!')
	{
		$chain_break_start=substr($dssp_dump_sans_header[$i-1],6,4);
		$chain_break_start=~ s/\s+//g;
		$chain_break_end=substr($dssp_dump_sans_header[$i+1],6,4);
		$chain_break_end=~ s/\s+//g;
		$No_of_residues_missing=$chain_break_end-$chain_break_start-1;
		$appropriate_number_of_X= 'X' x $No_of_residues_missing;
		print FOUT1 $appropriate_number_of_X;
		$appropriate_number_of_X_with_newline=join("\n",split('',$appropriate_number_of_X));
		print FOUT2 $appropriate_number_of_X_with_newline,"\n";	
	}	
	#Elsif either phi or psi are 360.0 then add a single 'X'
	elsif(("$eight_torsions[1]" =='360.0') || ("$eight_torsions[2]" =='360.0')|| ("$eight_torsions[5]" =='360.0')|| ("$eight_torsions[6]" =='360.0'))
	{
		
			print FOUT1 "X";
			print FOUT2 "X\n";
	}
	#Else... no problem, go ahed and calculate rmsda and assign the PB !
	else
	{
		($rmsda,$pb)=&calc_rmsda(\@eight_torsions);
		print FOUT1"$pb";
		print FOUT2"$rmsda\n"
	}
}


#Writing the .torsion file
print FOUT1"ZZ";
print FOUT2"X\nX\n";
for($j=0;$j<@dssp_dump_sans_header;$j++)
{
	$phi_and_psi[$j]=substr($dssp_dump_sans_header[$j],103,12);	
	$phi_and_psi[$j]=~ s/(.{6})/$1,/;
	if ($phi_and_psi[$j]=~ m/ 360.0, 360.0/)
	{
		$chain_break_start2=substr($dssp_dump_sans_header[$j-1],7,3);
		$chain_break_end2=substr($dssp_dump_sans_header[$j+1],7,3);
		$No_of_x2=$chain_break_end2-$chain_break_start2-2;
		$appropriate_number_of_X2= "     X,     X\n" x $No_of_x2;		
		print FOUT3 "$appropriate_number_of_X2";
	}
	$phi_and_psi[$j]=~ s/ 360.0/     X/;
	$phi_and_psi[$j]=~ s/\bX, 360.0\b/X,     X/;
	print FOUT3 "$phi_and_psi[$j]\n";
}

###############################################################################################################################

#Sub routines "calc_rmsda" and "ang_diff" are adapted from the script "read_pdb.pl" created by Swapnil MAHAJAN on 8th Aug. 2011

###############################################################################################################################

# To calculate RMSDa and assign PB using standard torsions for PBs. 
sub calc_rmsda
{
  my %std_dihed=(
	"A"=>[41.14,75.53,13.92,-99.80,131.88,-96.27,122.08,-99.68],
	"B"=>[108.24,-90.12,119.54,-92.21,-18.06,-128.93,147.04,-99.90],
	"C"=>[-11.61,-105.66,94.81,-106.09,133.56,-106.93,135.97,-100.63],
	"D"=>[141.98,-112.79,132.20,-114.79,140.11,-111.05,139.54,-103.16],
	"E"=>[133.25,-112.37,137.64,-108.13,133.00,-87.30,120.54,77.40],
	"F"=>[116.40,-105.53,129.32,-96.68,140.72,-74.19,-26.65,-94.51],
	"G"=>[0.40,-81.83,4.91,	-100.59,85.50,-71.65,130.78,84.98],
	"H"=>[119.14,-102.58,130.83,-67.91,121.55,76.25,-2.95,-90.88],
	"I"=>[130.68,-56.92,119.26,77.85,10.42,-99.43,141.40,-98.01],
	"J"=>[114.32,-121.47,118.14,82.88,-150.05,-83.81,23.35,-85.82],
	"K"=>[117.16,-95.41,140.40,-59.35,-29.23,-72.39,-25.08,-76.16],
	"L"=>[139.20,-55.96,-32.70,-68.51,-26.09,-74.44,-22.60,-71.74],
	"M"=>[-39.62,-64.73,-39.52,-65.54,-38.88,-66.89,-37.76,-70.19],
	"N"=>[-35.34,-65.03,-38.12,-66.34,-29.51,-89.10,-2.91,77.90],
	"O"=>[-45.29,-67.44,-27.72,-87.27,5.13,77.49,30.71,-93.23],
	"P"=>[-27.09,-86.14,0.30,59.85,21.51,-96.30,132.67,-92.91]
	);
	
	my @penta_arr=();
	@dihed=@{$_[0]};
	
	my $min_diff=999999; my $min_pb="";

	foreach my $k(sort(keys(%std_dihed)))
	{
		my $diff=0;
		for(my $n=0;$n<@dihed;$n++)
		{
			my $diff1=0;
			$diff1=&ang_diff($dihed[$n],${$std_dihed{$k}}[$n]);
			$diff+=($diff1**2);
		}
		$diff=sqrt($diff/8);
		if($diff<$min_diff){$min_diff=$diff;$min_pb=$k;}
	}
	
	return($min_diff,$min_pb);
}



# This subroutine calculates the real difference in two angles i.e. (-179)-(+179)=-358 but real difference is 2.
sub ang_diff
{
  $secondAngle=$_[1]; $firstAngle=$_[0];
  my $difference = $secondAngle - $firstAngle;
  while ($difference < -180) {$difference += 360;}
  while ($difference > 180) {$difference -= 360;}
  return $difference;
}

